/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/cstr.h"
#include "libgtcore/fileutils.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/xansi.h"
#include "libgtext/qgramdist.h"
#include "tools/gt_qgramdist.h"

static OPrval parse_options(int *parsed_args, unsigned int *q, int argc,
                            const char **argv, Error *err)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("[option ...] seq_file_1 seq_file_2",
                         "Compute q-gram distance for each sequence "
                         "combination.");
  o = option_new_uint_min("q", "set q", q, 3, 1);
  option_parser_add_option(op, o);
  option_parser_set_min_max_args(op, 2, 2);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

int gt_qgramdist(int argc, const char **argv, Error *err)
{
  Bioseq *bioseq_1 = NULL, *bioseq_2 = NULL;
  unsigned long i, j, dist;
  Seq *seq_1, *seq_2;
  int parsed_args, had_err = 0;
  unsigned int q;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &q, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args+1 < argc);

  /* make sure seq_file_1 exists */
  if (!file_exists(argv[parsed_args])) {
    error_set(err, "seq_file_1 \"%s\" does not exist", argv[parsed_args]);
    had_err = -1;
  }

  /* make sure seq_file_2 exists */
  if (!had_err && !file_exists(argv[parsed_args+1])) {
    error_set(err, "seq_file_2 \"%s\" does not exist", argv[parsed_args+1]);
    had_err = -1;
  }

  /* init */
  if (!had_err) {
    bioseq_1 = bioseq_new(argv[parsed_args], err);
    if (!bioseq_1)
      had_err = -1;
    if (!had_err) {
      bioseq_2 = bioseq_new(argv[parsed_args+1], err);
      if (!bioseq_2)
        had_err = -1;
    }

    /* compute q-gram distance for all sequence combinations */
    for (i = 0; i < bioseq_number_of_sequences(bioseq_1); i++) {
      for (j = 0; j < bioseq_number_of_sequences(bioseq_2); j++) {
        seq_1 = bioseq_get_seq(bioseq_1, i);
        seq_2 = bioseq_get_seq(bioseq_2, j);
        dist = qgramdist(seq_1, seq_2, q);
        printf("qgramdist_%u_(", q);
        cstr_show(seq_get_orig(seq_1), seq_length(seq_1), stdout);
        xputchar(',');
        cstr_show(seq_get_orig(seq_2), seq_length(seq_2), stdout);
        printf(")=%lu\n", dist);
      }
    }
  }

  /* free */
  bioseq_delete(bioseq_2);
  bioseq_delete(bioseq_1);

  return had_err;
}
