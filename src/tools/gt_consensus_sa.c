/*
  Copyright (c) 2005-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include "libgtcore/array.h"
#include "libgtcore/fa.h"
#include "libgtcore/fptr.h"
#include "libgtcore/option.h"
#include "libgtcore/range.h"
#include "libgtcore/str.h"
#include "libgtcore/strand.h"
#include "libgtcore/unused.h"
#include "libgtcore/xansi.h"
#include "libgtexercise/sspliced_alignment.h"
#include "libgtext/consensus_sa.h"
#include "tools/gt_consensus_sa.h"

#define DELIMITER         ','
#define FORWARDSTRANDCHAR '+'
#define REVERSESTRANDCHAR '-'

static int parse_input_line(SSplicedAlignment **alignment, const char *line,
                            unsigned long line_length, Error *err)
{
  long leftpos, rightpos;
  unsigned long i = 0;
  Range exon;
  Str *id;
  int had_err = 0;
  error_check(err);

#define CHECKLINELENGTH\
        if (!had_err && i >= line_length) {        \
          error_set(err, "incomplete input line\n" \
                       "line=%s", line);           \
          return had_err = -1;                     \
        }

  /* init */
  id = str_new();
  *alignment = NULL;

  /* parsing id */
  while (!had_err) {
    CHECKLINELENGTH;
    if (line[i] == DELIMITER) {
      /* reference id has been saved, skip this character and break */
      i++;
      CHECKLINELENGTH;
      break;
    }
    else {
      /* save this character of the reference id */
      str_append_char(id, line[i]);
    }

    /* increase counter */
    i++;
  }

  /* parsing orientation */
  if (line[i] == FORWARDSTRANDCHAR)
    *alignment = sspliced_alignment_new(str_get(id), true);
  else if (line[i] == REVERSESTRANDCHAR)
    *alignment = sspliced_alignment_new(str_get(id), false);
  else {
    error_set(err, "wrong formatted input line, orientation must be %c or %c\n"
                   "line=%s", FORWARDSTRANDCHAR, REVERSESTRANDCHAR, line);
    had_err = -1;
  }
  i++;
  CHECKLINELENGTH;

  if (!had_err && line[i] != DELIMITER) {
    error_set(err, "incomplete input line\nline=%s", line);
    had_err = -1;
  }

  while (!had_err) {
    if (line[i] == DELIMITER) {
      i++;
      CHECKLINELENGTH;
      if (!had_err && sscanf(line+i, "%ld-%ld", &leftpos, &rightpos) != 2) {
        error_set(err, "incomplete input line\nline=%s", line);
        had_err = -1;
      }
      if (!had_err) {
        /* save exon */
        exon.start = leftpos;
        exon.end   = rightpos;
        sspliced_alignment_add_exon(*alignment, exon);
      }
    }
    i++;
    if (i >= line_length)
      break;
  }

  if (had_err)
    sspliced_alignment_delete(*alignment);
  else {
    /* alignment contains at least one exon */
    assert(sspliced_alignment_num_of_exons(*alignment));
  }
  str_delete(id);

  return had_err;
}

static int parse_input_file(Array *spliced_alignments,
                            const char *file_name, Error *err)
{
  FILE *input_file;
  SSplicedAlignment *sa;
  int had_err = 0;
  Str *line;
  error_check(err);

  line = str_new();
  input_file = fa_xfopen(file_name, "r");

  while (!had_err && str_read_next_line(line, input_file) != EOF) {
    /* parse input line and save result in spliced alignment */
    had_err = parse_input_line(&sa, str_get(line), str_length(line), err);
    if (!had_err) {
      /* store spliced alignment */
      array_add(spliced_alignments, sa);
      /* reset array */
      str_reset(line);
    }
  }

  fa_xfclose(input_file);
  str_delete(line);
  return had_err;
}

static Range get_genomic_range(const void *sa)
{
  SSplicedAlignment *alignment = *(SSplicedAlignment**) sa;
  assert(alignment);
  return sspliced_alignment_genomic_range(alignment);
}

static Strand get_strand(const void *sa)
{
  SSplicedAlignment *alignment = *(SSplicedAlignment**) sa;
  if (sspliced_alignment_is_forward(alignment))
    return STRAND_FORWARD;
  return STRAND_REVERSE;
}

static void get_exons(Array *exon_ranges, const void *sa)
{
  SSplicedAlignment *alignment = *(SSplicedAlignment**) sa;
  Range exon;
  unsigned long i;
  assert(alignment);
  for (i = 0; i < sspliced_alignment_num_of_exons(alignment); i++) {
    exon = sspliced_alignment_get_exon(alignment, i);
    array_add(exon_ranges, exon);
  }
}

static void process_splice_form(Array *spliced_alignments_in_form,
                                UNUSED const void *set_of_sas,
                                UNUSED unsigned long number_of_sas,
                                UNUSED size_t size_of_sa,
                                UNUSED void *userdata)
{
  unsigned long i;

  printf("contains [");
  for (i = 0; i < array_size(spliced_alignments_in_form); i++) {
    if (i)
      xputchar(',');
    printf("%lu", *((unsigned long*) array_get(spliced_alignments_in_form, i)));
  }
  printf("]\n");
}

static OptionParser* gt_consensus_sa_option_parser_new(UNUSED
                                                       void *tool_arguments)
{
  OptionParser *op;
  op = option_parser_new("spliced_alignment_file", "Read file containing "
                         "spliced alingments, compute consensus spliced "
                         "alignments,\nand print them to stdout.");
  option_parser_set_min_max_args(op, 1, 1);
  return op;
}

static int gt_consensus_sa_runner(UNUSED int argc, const char **argv,
                                  int parsed_args, UNUSED void *tool_arguments,
                                  Error *err)
{
  Array *spliced_alignments;
  SSplicedAlignment *sa;
  unsigned long i;
  int had_err = 0;
  error_check(err);

  /* parse input file and store resuilts in the spliced alignment array */
  spliced_alignments = array_new(sizeof (SSplicedAlignment*));
  had_err = parse_input_file(spliced_alignments, argv[parsed_args], err);

  if (!had_err) {
    /* sort spliced alignments */
    qsort(array_get_space(spliced_alignments), array_size(spliced_alignments),
          array_elem_size(spliced_alignments),
          (Compare) sspliced_alignment_compare_ptr);

    /* compute the consensus spliced alignments */
    consensus_sa(array_get_space(spliced_alignments),
                 array_size(spliced_alignments),
                 array_elem_size(spliced_alignments), get_genomic_range,
                 get_strand, get_exons, process_splice_form, NULL);
  }

  /* free */
  for (i = 0; i < array_size(spliced_alignments); i++) {
    sa = *(SSplicedAlignment**) array_get(spliced_alignments, i);
    sspliced_alignment_delete(sa);
  }
  array_delete(spliced_alignments);

  return had_err;
}

Tool* gt_consensus_sa(void)
{
  return tool_new(NULL,
                  NULL,
                  gt_consensus_sa_option_parser_new,
                  NULL,
                  gt_consensus_sa_runner);
}
