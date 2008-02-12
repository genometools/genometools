/*
  Copyright (c) 2003-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "libgtcore/bioseq.h"
#include "libgtcore/option.h"
#include "libgtcore/unused.h"
#include "libgtcore/versionfunc.h"
#include "libgtext/linearedist.h"
#include "libgtext/upgma.h"
#include "tools/gt_upgma.h"

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            Error *err)
{
  OptionParser *op;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("sequence_file|example", "Compute and show UPGMA tree "
                         "for the sequences in sequence file (using the unit\n"
                         "cost edit distance as distance function). If "
                         "'example' is given as\nsequence_file, a builtin "
                         "example is used.");
  option_parser_set_min_max_args(op, 1, 1);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

static double distfunc(unsigned long i, unsigned long j, void *data)
{
  Bioseq *bioseq= (Bioseq*) data;
  return linearedist(bioseq_get_sequence(bioseq, i),
                     bioseq_get_sequence_length(bioseq, i),
                     bioseq_get_sequence(bioseq, j),
                     bioseq_get_sequence_length(bioseq, j));
}

static double exampledistfunc(unsigned long i, unsigned long j,
                              UNUSED void *data)
{
  static const double exampledistances[5][5] =
    { {0.0   , 0.1715, 0.2147, 0.3091, 0.2326},
      {0.1715, 0.0   , 0.2991, 0.3399, 0.2058},
      {0.2147, 0.2991, 0.0   , 0.2795, 0.3943},
      {0.3091, 0.3399, 0.2795, 0.0   , 0.4289},
      {0.2326, 0.2058, 0.3943, 0.4289, 0.0   } };
  return exampledistances[i][j];
}

int gt_upgma(int argc, const char **argv, Error *err)
{
  bool use_hard_coded_example = false;
  int parsed_args, had_err = 0;
  Bioseq *bioseq = NULL;
  UPGMA *upgma = NULL;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args == 1);

  if (!strcmp(argv[parsed_args], "example"))
    use_hard_coded_example = true;

  if (use_hard_coded_example)
    upgma = upgma_new(5, NULL, exampledistfunc);
  else {
    bioseq = bioseq_new(argv[parsed_args], err);
    if (!bioseq)
      had_err = -1;
    if (!had_err)
      upgma = upgma_new(bioseq_number_of_sequences(bioseq), bioseq, distfunc);
  }

  if (!had_err)
    upgma_show_tree(upgma, stdout);

  bioseq_delete(bioseq);
  upgma_delete(upgma);

  return had_err;
}
