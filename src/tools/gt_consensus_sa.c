/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/option.h"
#include "libgtcore/range.h"
#include "libgtcore/str.h"
#include "libgtcore/strand.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/xansi.h"
#include "libgtext/consensus_sa.h"

#define DELIMITER         ','
#define FORWARDSTRANDCHAR '+'
#define REVERSESTRANDCHAR '-'

typedef struct {
  Str* id;
  bool forward;
  Array *exons; /* the exon ranges */
} SimpleSplicedAlignment;

static void initSimpleSplicedAlignment(SimpleSplicedAlignment *sa)
{
  assert(sa);
  sa->id = str_new();
  sa->forward = true;
  sa->exons = array_new(sizeof (Range));
}

static int parse_input_line(SimpleSplicedAlignment *alignment, const char *line,
                            unsigned long line_length, Error *e)
{
  long leftpos, rightpos;
  unsigned long i = 0;
  Range exon;
  error_check(e);

#define CHECKLINELENGTH\
        if (i >= line_length) {                  \
          error_set(e, "incomplete input line\n" \
                       "line=%s", line);         \
          return -1;                             \
        }

  /* parsing id */
  for (;;) {
    CHECKLINELENGTH;
    if (line[i] == DELIMITER) {
      /* reference id has been saved, skip this character and break */
      i++;
      CHECKLINELENGTH;
      break;
    }
    else {
      /* save this character of the reference id */
      str_append_char(alignment->id, line[i]);
    }

    /* increase counter */
    i++;
  }

  /* parsing orientation */
  if (line[i] == FORWARDSTRANDCHAR)
    alignment->forward = true;
  else if (line[i] == REVERSESTRANDCHAR)
    alignment->forward = false;
  else {
    error_set(e, "wrong formatted input line, orientation must be %c or %c\n"
                 "line=%s", FORWARDSTRANDCHAR, REVERSESTRANDCHAR, line);
    return -1;
  }
  i++;
  CHECKLINELENGTH;

  if (line[i] != DELIMITER) {
    error_set(e, "incomplete input line\nline=%s", line);
    return -1;
  }

  for (;;) {
    if (line[i] == DELIMITER) {
      i++;
      CHECKLINELENGTH;
      if (sscanf(line+i, "%ld-%ld", &leftpos, &rightpos) != 2) {
        error_set(e, "incomplete input line\nline=%s", line);
        return -1;
      }
      exon.start = leftpos;
      exon.end   = rightpos;

      /* save exon */
      array_add(alignment->exons, exon);
    }
    i++;
    if (i >= line_length)
      break;
  }

  /* alignment contains at least one exon */
  assert(array_size(alignment->exons));

  return 0;
}

static int parse_input_file(Array *spliced_alignments,
                             const char *file_name, Error *e)
{
  FILE *input_file;
  SimpleSplicedAlignment sa;
  int had_err = 0;
  Str *line;
  error_check(e);

  line = str_new();
  input_file = fa_xfopen(file_name, "r");

  while (!had_err && str_read_next_line(line, input_file) != EOF) {
    /* init new spliced alignment */
    initSimpleSplicedAlignment(&sa);
    /* parse input line and save result in spliced alignment */
    had_err = parse_input_line(&sa, str_get(line), str_length(line), e);
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
  SimpleSplicedAlignment *alignment = (SimpleSplicedAlignment*) sa;
  Range range;
  assert(alignment);
  range.start = ((Range*) array_get_first(alignment->exons))->start;
  range.end   = ((Range*) array_get_last(alignment->exons))->end;
  return range;
}

static Strand get_strand(const void *sa)
{
  SimpleSplicedAlignment *alignment = (SimpleSplicedAlignment*) sa;
  if (alignment->forward)
    return STRAND_FORWARD;
  return STRAND_REVERSE;
}

static void get_exons(Array *exon_ranges, const void *sa)
{
  SimpleSplicedAlignment *alignment = (SimpleSplicedAlignment*) sa;
  assert(alignment);
  array_add_array(exon_ranges, alignment->exons);
}

static void process_splice_form(Array *spliced_alignments_in_form,
                                /*@unused@*/ const void *set_of_sas,
                                /*@unused@*/ unsigned long number_of_sas,
                                /*@unused@*/ size_t size_of_sa,
                                /*@unused@*/ void *userdata)
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

static int range_compare_long_first(Range range_a, Range range_b)
{
  assert(range_a.start <= range_a.end && range_b.start <= range_b.end);

  if ((range_a.start == range_b.start) && (range_a.end == range_b.end))
    return 0; /* range_a == range_b */

  if ((range_a.start < range_b.start) ||
      ((range_a.start == range_b.start) && (range_a.end > range_b.end)))
    return -1; /* range_a < range_b */

  return 1; /* range_a > range_b */
}

static int compare_spliced_alignment(const void *a, const void *b)
{
  SimpleSplicedAlignment *sa_a = (SimpleSplicedAlignment*) a,
                   *sa_b = (SimpleSplicedAlignment*) b;
  Range range_a, range_b;
  range_a.start = ((Range*) array_get_first(sa_a->exons))->start;
  range_a.end   = ((Range*) array_get_last(sa_a->exons))->end;
  range_b.start = ((Range*) array_get_first(sa_b->exons))->start;
  range_b.end   = ((Range*) array_get_last(sa_b->exons))->end;
  return range_compare_long_first(range_a, range_b);
}

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            Error *err)
{
  OptionParser *op;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("spliced_alignment_file", "Read file containing "
                         "spliced alingments, compute consensus spliced "
                         "alignments,\nand print them to stdout.");
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, 1, err);
  option_parser_delete(op);
  return oprval;
}

int gt_consensus_sa(int argc, const char **argv, Env *env)
{
  Array *spliced_alignments;
  SimpleSplicedAlignment *sa;
  unsigned long i;
  int parsed_args, had_err = 0;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, argc, argv, env_error(env))) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args == 1);

  /* parse input file and store resuilts in the spliced alignment array */
  spliced_alignments = array_new(sizeof (SimpleSplicedAlignment));
  had_err = parse_input_file(spliced_alignments, argv[1], env_error(env));

  if (!had_err) {
    /* sort spliced alignments */
    qsort(array_get_space(spliced_alignments), array_size(spliced_alignments),
          sizeof (SimpleSplicedAlignment), compare_spliced_alignment);

    /* compute the consensus spliced alignments */
    consensus_sa(array_get_space(spliced_alignments),
                 array_size(spliced_alignments),
                 sizeof (SimpleSplicedAlignment), get_genomic_range, get_strand,
                 get_exons, process_splice_form, NULL);
  }

  /* free */
  for (i = 0; i < array_size(spliced_alignments); i++) {
    sa = array_get(spliced_alignments, i);
    str_delete(sa->id);
    array_delete(sa->exons);
  }
  array_delete(spliced_alignments);

  return had_err;
}
