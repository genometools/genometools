/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/bioseq.h"
#include "libgtcore/option.h"
#include "libgtcore/unused.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/xansi.h"
#include "libgtext/align.h"
#include "libgtext/alignment.h"
#include "tools/gt_align.h"

static OPrval parse_options(int *parsed_args, bool *all, int argc,
                            const char **argv, Error *err)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("[option ...] seq_file_1 seq_file_2",
                         "Globally align each sequence in seq_file_1 with each "
                         "sequence in seq_file_2.");
  option = option_new_bool("all", "show all optimal alignments instead of just "
                           "one", all, false);
  option_parser_add_option(op, option);
  option_parser_set_min_max_args(op, 2, 2);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

void show_alignment(const Alignment *a, UNUSED void *data)
{
  assert(a && !data);
  alignment_show(a, stdout);
  xputchar('\n');
}

void show_aligns(unsigned long aligns, UNUSED void *data)
{
  assert(aligns && !data);
  printf("number of optimal alignments: %lu\n\n", aligns);
}

int gt_align(int argc, const char **argv, Error *err)
{
  Bioseq *bioseq_1, *bioseq_2 = NULL;
  unsigned long i, j;
  int parsed_args, had_err = 0;
  Alignment *a;
  bool all;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &all, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args+1 < argc);

  /* init */
  bioseq_1 = bioseq_new(argv[parsed_args], err);
  if (!bioseq_1)
    had_err = -1;
  if (!had_err) {
    bioseq_2 = bioseq_new(argv[parsed_args+1], err);
    if (!bioseq_2)
      had_err = -1;
  }

  /* aligning all sequence combinations */
  if (!had_err) {
    for (i = 0; i < bioseq_number_of_sequences(bioseq_1); i++) {
      for (j = 0; j < bioseq_number_of_sequences(bioseq_2); j++) {
        if (all) {
          align_all(bioseq_get_sequence(bioseq_1, i),
                    bioseq_get_sequence_length(bioseq_1, i),
                    bioseq_get_sequence(bioseq_2, j),
                    bioseq_get_sequence_length(bioseq_2, j),
                    show_alignment, show_aligns, NULL);
        }
        else {
          a = align(bioseq_get_sequence(bioseq_1, i),
                    bioseq_get_sequence_length(bioseq_1, i),
                    bioseq_get_sequence(bioseq_2, j),
                    bioseq_get_sequence_length(bioseq_2, j));
          alignment_show(a, stdout);
          xputchar('\n');
          alignment_delete(a);
        }
      }
    }
  }

  /* free */
  bioseq_delete(bioseq_2);
  bioseq_delete(bioseq_1);

  return had_err;
}
