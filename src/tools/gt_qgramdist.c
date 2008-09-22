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

#include "core/bioseq.h"
#include "core/cstr.h"
#include "core/fileutils.h"
#include "core/option.h"
#include "core/versionfunc.h"
#include "core/xansi.h"
#include "extended/qgramdist.h"
#include "tools/gt_qgramdist.h"

static OPrval parse_options(int *parsed_args, unsigned int *q, int argc,
                            const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOption *o;
  OPrval oprval;
  gt_error_check(err);
  op = gt_option_parser_new("[option ...] gt_seq_file_1 gt_seq_file_2",
                         "Compute q-gram distance for each sequence "
                         "combination.");
  o = gt_option_new_uint_min("q", "set q", q, 3, 1);
  gt_option_parser_add_option(op, o);
  gt_option_parser_set_min_max_args(op, 2, 2);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, versionfunc,
                                  err);
  gt_option_parser_delete(op);
  return oprval;
}

int gt_qgramdist(int argc, const char **argv, GtError *err)
{
  GtBioseq *gt_bioseq_1 = NULL, *gt_bioseq_2 = NULL;
  unsigned long i, j, dist;
  GtSeq *seq_1, *seq_2;
  int parsed_args, had_err = 0;
  unsigned int q;
  gt_error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, &q, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args+1 < argc);

  /* make sure gt_seq_file_1 exists */
  if (!file_exists(argv[parsed_args])) {
    gt_error_set(err, "seq_file_1 \"%s\" does not exist", argv[parsed_args]);
    had_err = -1;
  }

  /* make sure gt_seq_file_2 exists */
  if (!had_err && !file_exists(argv[parsed_args+1])) {
    gt_error_set(err, "seq_file_2 \"%s\" does not exist", argv[parsed_args+1]);
    had_err = -1;
  }

  /* init */
  if (!had_err) {
    gt_bioseq_1 = gt_bioseq_new(argv[parsed_args], err);
    if (!gt_bioseq_1)
      had_err = -1;
    if (!had_err) {
      gt_bioseq_2 = gt_bioseq_new(argv[parsed_args+1], err);
      if (!gt_bioseq_2)
        had_err = -1;
    }

    /* compute q-gram distance for all sequence combinations */
    for (i = 0; i < gt_bioseq_number_of_sequences(gt_bioseq_1); i++) {
      for (j = 0; j < gt_bioseq_number_of_sequences(gt_bioseq_2); j++) {
        seq_1 = gt_bioseq_get_seq(gt_bioseq_1, i);
        seq_2 = gt_bioseq_get_seq(gt_bioseq_2, j);
        dist = qgramdist(seq_1, seq_2, q);
        printf("qgramdist_%u_(", q);
        gt_cstr_show(gt_seq_get_orig(seq_1), gt_seq_length(seq_1), stdout);
        gt_xputchar(',');
        gt_cstr_show(gt_seq_get_orig(seq_2), gt_seq_length(seq_2), stdout);
        printf(")=%lu\n", dist);
      }
    }
  }

  /* free */
  gt_bioseq_delete(gt_bioseq_2);
  gt_bioseq_delete(gt_bioseq_1);

  return had_err;
}
