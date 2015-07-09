/*
  Copyright (c) 2015 Jörg Winkler <joerg.winkler@studium.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#include <limits.h>
#include "core/encseq_api.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "match/diagbandseed.h"
#include "tools/gt_seed_extend.h"

typedef struct {
  unsigned int dbs_seedlength;
  GtUword dbs_logdiagbandwidth;
  GtUword dbs_mincoverage;
  GtUword dbs_maxfreq;
  GtUword se_alignlength;
  GtUword se_maxalilendiff;
  GtUword se_errorpercentage;
  GtUword se_historysize;
  bool mirror;
  bool verify;
  bool benchmark;
} GtSeedExtendArguments;

static void* gt_seed_extend_arguments_new(void)
{
  GtSeedExtendArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  return arguments;
}

static void gt_seed_extend_arguments_delete(void *tool_arguments)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_free(arguments);
  }
}

static GtOptionParser* gt_seed_extend_option_parser_new(void *tool_arguments)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] encseq_basename [encseq_basename]",
                            "Calculate local alignments using the seed and "
                            "extend algorithm.");

  /* -seedlength */
  option = gt_option_new_uint_min("seedlength", "Minimum length of a seed",
                                  &arguments->dbs_seedlength, 14, 2UL);
  gt_option_parser_add_option(op, option);

  /* -diagbandwidth */
  option = gt_option_new_uword("diagbandwidth", "Logarithm of diagonal band "
                               "width (for filter)",
                               &arguments->dbs_logdiagbandwidth, 6);
  gt_option_parser_add_option(op, option);

  /* -mincoverage */
  option = gt_option_new_uword("mincoverage", "Minimum coverage in two "
                               "neighbouring diagonal bands (for filter)",
                               &arguments->dbs_mincoverage, 35);
  gt_option_parser_add_option(op, option);

  /* -maxfreq */
  option = gt_option_new_uword_min("maxfreq", "Maximum frequency of a k-mer "
                                   "(for filter)",
                                   &arguments->dbs_maxfreq, GT_UWORD_MAX, 1);
  gt_option_parser_add_option(op, option);

  /* -err */
  option = gt_option_new_uword_min_max("err", "Error percentage of matches "
                                       "(for greedy extension)",
                                       &arguments->se_errorpercentage, 15, 0,
                                       100);
  gt_option_parser_add_option(op, option);

  /* -alignlength */
  option = gt_option_new_uword_min("alignlength", "Minimum alignment length",
                                   &arguments->se_alignlength, 1000, 1);
  gt_option_parser_add_option(op, option);

  /* -maxalilendiff */
  option = gt_option_new_uword("maxalilendiff", "Maximum difference of align"
                               "ment length (trimming for greedy extension)",
                               &arguments->se_maxalilendiff, 30);
  gt_option_parser_add_option(op, option);

  /* -historysize */
  option = gt_option_new_uword_min("historysize", "Size of (mis)match history "
                                   "(for greedy extension)",
                                   &arguments->se_historysize, 65, 1);
  gt_option_parser_add_option(op, option);

  /* -mirror */
  option = gt_option_new_bool("mirror", "Add reverse complement reads",
                              &arguments->mirror, false);
  gt_option_parser_add_option(op, option);

  /* -verify */
  option = gt_option_new_bool("verify", "Check that k-mer seeds occur in the "
                              "sequences", &arguments->verify, false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  /* -benchmark */
  option = gt_option_new_bool("benchmark", "Measure time of different steps",
                              &arguments->benchmark, false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  return op;
}

static int gt_seed_extend_arguments_check(GT_UNUSED int rest_argc,
                                          void *tool_arguments, GtError *err)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if (rest_argc == 1 && arguments->dbs_maxfreq == 1) {
    gt_error_set(err, "for 1 input file maxfreq must be >= 2 to find matching "
                 "k-mers");
    had_err = -1;
  }

  if (rest_argc > 2) {
    gt_error_set(err, "too many arguments (-help shows correct usage)");
    had_err = -1;
  } else if (rest_argc < 1) {
    gt_error_set(err, "at least one encseq index name must be specified "
      "(-help shows correct usage)");
    had_err = -1;
  }
  return had_err;
}

static int gt_seed_extend_runner(int argc, const char **argv, int parsed_args,
                                 void *tool_arguments, GtError *err)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  GtEncseqLoader *encseq_loader;
  GtEncseq *aencseq, *bencseq;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);
  gt_assert(argc - parsed_args >= 1);

  /* load encseq A */
  encseq_loader = gt_encseq_loader_new();
  gt_encseq_loader_enable_autosupport(encseq_loader);
  aencseq = gt_encseq_loader_load(encseq_loader, argv[parsed_args], err);
  if (!aencseq)
    had_err = -1;

  /* if there is a 2nd read set: load encseq B */
  if (!had_err) {
    if (argc - parsed_args == 2)
      bencseq = gt_encseq_loader_load(encseq_loader, argv[parsed_args+1], err);
    else
      bencseq = gt_encseq_ref(aencseq);
    if (!bencseq) {
      had_err = -1;
      gt_encseq_delete(aencseq);
    }
  }
  gt_encseq_loader_delete(encseq_loader);

  /* start algorithm */
  if (!had_err) {
    gt_diagbandseed_run(aencseq, bencseq, (GtDiagbandseed *)arguments);
    gt_encseq_delete(aencseq);
    gt_encseq_delete(bencseq);
  }
  return had_err;
}

GtTool* gt_seed_extend(void)
{
  return gt_tool_new(gt_seed_extend_arguments_new,
                     gt_seed_extend_arguments_delete,
                     gt_seed_extend_option_parser_new,
                     gt_seed_extend_arguments_check,
                     gt_seed_extend_runner);
}
