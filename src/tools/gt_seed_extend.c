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
#include "core/error_api.h"
#include "core/ma_api.h"
#include "core/str_api.h"
#include "match/diagbandseed.h"
#include "match/seed-extend.h"
#include "match/xdrop.h"
#include "tools/gt_seed_extend.h"

typedef struct {
  unsigned int dbs_seedlength;
  GtUword dbs_logdiagbandwidth;
  GtUword dbs_mincoverage;
  GtUword dbs_maxfreq;
  GtUword se_minalignlength;
  GtUword se_maxalilendiff;
  GtUword se_errorpercentage;
  GtUword se_historysize;
  GtUword se_perc_match_hist;
  GtStr *se_char_access_mode;
  GtXdropscore se_xdropbelowscore;
  bool extendgreedy;
  bool extendxdrop;
  bool mirror;
  bool overlappingseeds;
  bool verify;
  bool benchmark;
} GtSeedExtendArguments;

static void* gt_seed_extend_arguments_new(void)
{
  GtSeedExtendArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->se_char_access_mode = gt_str_new();
  return arguments;
}

static void gt_seed_extend_arguments_delete(void *tool_arguments)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->se_char_access_mode);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_seed_extend_option_parser_new(void *tool_arguments)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *op_gre, *op_xdr, *op_cam, *op_his, *op_dif, *op_pmh,
    *op_len, *op_err, *op_xbe;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] encseq_basename [encseq_basename]",
                            "Calculate local alignments using the seed and "
                            "extend algorithm.");

  /* -seedlength */
  op_len = gt_option_new_uint_min("seedlength", "Minimum length of a seed",
                                  &arguments->dbs_seedlength, 14, 2UL);
  gt_option_parser_add_option(op, op_len);

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
  op_err = gt_option_new_uword_min_max("err", "Error percentage of matches "
                                       "(for greedy extension)",
                                       &arguments->se_errorpercentage, 15, 0,
                                       100);
  gt_option_parser_add_option(op, op_err);

  /* -alignlength */
  op_len = gt_option_new_uword_min("alignlength", "Minimum alignment length",
                                   &arguments->se_minalignlength, 1000, 1);
  gt_option_parser_add_option(op, op_len);

  /* -maxalilendiff */
  op_dif = gt_option_new_uword("maxalilendiff", "Maximum difference of align"
                               "ment length (trimming for greedy extension)",
                               &arguments->se_maxalilendiff, 30);
  gt_option_parser_add_option(op, op_dif);

  /* -history */
  op_his = gt_option_new_uword_min_max("history", "Size of (mis)match history "
                                       "in range [1..64] (trimming for greedy "
                                       "extension)",
                                       &arguments->se_historysize, 60, 0, 64);
  gt_option_parser_add_option(op, op_his);

  /* -percmathistory */
  op_pmh = gt_option_new_uword_min_max("percmathistory", "percentage of matches"
                                       " required in history",
                                       &arguments->se_perc_match_hist, 55, 1,
                                       100);
  gt_option_parser_add_option(op, op_pmh);


  /* -xdropbelow */
  op_xbe = gt_option_new_word("xdropbelow", "Specify xdrop cutoff score",
                              &arguments->se_xdropbelowscore, 5L);
  gt_option_parser_add_option(op, op_xbe);

  /* -cam */
  op_cam = gt_option_new_string("cam", gt_cam_extendgreedy_comment(),
                                arguments->se_char_access_mode,"");
  gt_option_parser_add_option(op, op_cam);
  gt_option_is_development_option(op_cam);

  /* -extendgreedy */
  op_gre = gt_option_new_bool("extendgreedy", "Extend seed to both sides using "
                              "xdrop algorithm",
                              &arguments->extendgreedy, false);
  gt_option_parser_add_option(op, op_gre);

  /* -extendxdrop */
  op_xdr = gt_option_new_bool("extendxdrop", "Extend seed to both sides using "
                              "greedy algorithm with trimming of waves",
                              &arguments->extendxdrop, false);
  gt_option_parser_add_option(op, op_xdr);

  /* -mirror */
  option = gt_option_new_bool("mirror", "Add reverse complement reads",
                              &arguments->mirror, false);
  gt_option_parser_add_option(op, option);

  /* -overlappingseeds */
  option = gt_option_new_bool("overlappingseeds", "allow SeedPairs, which "
                              "overlap", &arguments->overlappingseeds, false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

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

  gt_option_exclude(op_gre, op_xdr);
  gt_option_imply(op_xbe, op_xdr);
  gt_option_imply(op_cam, op_gre);
  gt_option_imply(op_his, op_gre);
  gt_option_imply(op_dif, op_gre);
  gt_option_imply(op_pmh, op_gre);
  gt_option_imply_either_2(op_len, op_xdr, op_gre);
  gt_option_imply_either_2(op_err, op_xdr, op_gre);

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
  GtGreedyextendmatchinfo *grextinfo = NULL;
  GtXdropmatchinfo *xdropinfo = NULL;
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

  if (!had_err && arguments->extendgreedy) {
    GtExtendCharAccess cam = gt_greedy_extend_char_access(gt_str_get
                             (arguments->se_char_access_mode), err);
    if ((int) cam != -1) {
      grextinfo = gt_greedy_extend_matchinfo_new(arguments->se_errorpercentage,
                                                 arguments->se_maxalilendiff,
                                                 arguments->se_historysize,
                                                 arguments->se_perc_match_hist,
                                                 arguments->se_minalignlength,
                                                 cam);
    } else {
      had_err = -1;
    }
  }

  if (!had_err && arguments->extendxdrop) {
    xdropinfo = gt_xdrop_matchinfo_new(arguments->se_minalignlength,
                                       arguments->se_errorpercentage,
                                       arguments->se_xdropbelowscore,
                                       true);
  }

  /* start algorithm */
  if (!had_err) {
    GtDiagbandseed dbsarguments;
    dbsarguments.seedlength = arguments->dbs_seedlength;
    dbsarguments.logdiagbandwidth = arguments->dbs_logdiagbandwidth;
    dbsarguments.mincoverage = arguments->dbs_mincoverage;
    dbsarguments.maxfreq = arguments->dbs_maxfreq;
    dbsarguments.mirror = arguments->mirror;
    dbsarguments.overlappingseeds = arguments->overlappingseeds;
    dbsarguments.verify = arguments->verify;
    dbsarguments.benchmark = arguments->benchmark;
    dbsarguments.extendgreedyinfo = grextinfo;
    dbsarguments.extendxdropinfo = xdropinfo;

    had_err = gt_diagbandseed_run(aencseq, bencseq, &dbsarguments, err);
    gt_encseq_delete(aencseq);
    gt_encseq_delete(bencseq);
    if (arguments->extendgreedy)
      gt_greedy_extend_matchinfo_delete(grextinfo);
    if (arguments->extendxdrop)
      gt_xdrop_matchinfo_delete(xdropinfo);
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
