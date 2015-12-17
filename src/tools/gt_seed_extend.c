/*
  Copyright (c) 2015 Joerg Winkler <joerg.winkler@studium.uni-hamburg.de>
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
#include "core/alphabet_api.h"
#include "core/encseq.h"
#include "core/encseq_api.h"
#include "core/error_api.h"
#include "core/ma_api.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/showtime.h"
#include "core/str_api.h"
#include "match/diagbandseed.h"
#include "match/seed-extend.h"
#include "match/xdrop.h"
#include "match/ft-polish.h"
#include "match/initbasepower.h"
#include "tools/gt_seed_extend.h"

typedef struct {
  /* diagbandseed options */
  GtStr *dbs_indexname;
  GtStr *dbs_queryname;
  unsigned int dbs_seedlength;
  GtUword dbs_logdiagbandwidth;
  GtUword dbs_mincoverage;
  GtUword dbs_maxfreq;
  GtUword dbs_suppress;
  GtUword dbs_memlimit;
  GtStr *dbs_memlimit_str;
  bool dbs_debug_kmer;
  bool dbs_debug_seedpair;
  bool dbs_verify;
  bool weakends;
  bool onlyseeds;
  /* xdrop extension options */
  GtOption *se_option_xdrop;
  GtUword se_extendxdrop;
  GtXdropscore se_xdropbelowscore;
  /* greedy extension options */
  GtOption *se_option_greedy;
  GtUword se_extendgreedy;
  GtUword se_historysize;
  GtUword se_maxalilendiff;
  GtUword se_perc_match_hist;
  GtStr *se_char_access_mode;
  bool bias_parameters;
  bool relax_polish;
  /* general options */
  GtOption *se_option_withali;
  GtUword se_alignlength;
  GtUword se_minidentity;
  GtUword se_alignmentwidth;
  bool norev;
  bool nofwd;
  bool overlappingseeds;
  bool benchmark;
  bool verbose;
  bool seed_display;
  bool extend_last;
} GtSeedExtendArguments;

static void* gt_seed_extend_arguments_new(void)
{
  GtSeedExtendArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->dbs_indexname = gt_str_new();
  arguments->dbs_queryname = gt_str_new();
  arguments->dbs_memlimit_str = gt_str_new();
  arguments->se_char_access_mode = gt_str_new();
  return arguments;
}

static void gt_seed_extend_arguments_delete(void *tool_arguments)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->dbs_indexname);
    gt_str_delete(arguments->dbs_queryname);
    gt_str_delete(arguments->dbs_memlimit_str);
    gt_str_delete(arguments->se_char_access_mode);
    gt_option_delete(arguments->se_option_greedy);
    gt_option_delete(arguments->se_option_xdrop);
    gt_option_delete(arguments->se_option_withali);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_seed_extend_option_parser_new(void *tool_arguments)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *op_gre, *op_xdr, *op_cam, *op_his, *op_dif, *op_pmh,
    *op_len, *op_err, *op_xbe, *op_sup, *op_frq, *op_mem, *op_ali, *op_bia,
    *op_onl, *op_weakends, *op_seed_display, *op_relax_polish,
    *op_norev, *op_nofwd;
  gt_assert(arguments != NULL);

  /* init */
  op = gt_option_parser_new("[option ...] encseq_basename [encseq_basename]",
                            "Calculate local alignments using the seed and "
                            "extend algorithm.");

  /* DIAGBANDSEED OPTIONS */

  /* -ii */
  option = gt_option_new_string("ii",
                                "Input index for encseq encoded sequences",
                                arguments->dbs_indexname,
                                "");
  gt_option_is_mandatory(option);
  gt_option_hide_default(option);
  gt_option_parser_add_option(op, option);

  /* -qii */
  option = gt_option_new_string("qii",
                                "Query input index (encseq)",
                                arguments->dbs_queryname,
                                "");
  gt_option_hide_default(option);
  gt_option_parser_add_option(op, option);

  /* -seedlength */
  op_len = gt_option_new_uint_min_max("seedlength",
                                      "Minimum length of a seed\n"
                                      "default: logarithm of input length "
                                      "to the basis alphabet size",
                                      &arguments->dbs_seedlength,
                                      UINT_MAX, 1UL, 32UL);
  gt_option_hide_default(op_len);
  gt_option_parser_add_option(op, op_len);

  /* -diagbandwidth */
  option = gt_option_new_uword("diagbandwidth",
                               "Logarithm of diagonal band width (for filter)",
                               &arguments->dbs_logdiagbandwidth,
                               6UL);
  gt_option_parser_add_option(op, option);

  /* -mincoverage */
  option = gt_option_new_uword_min("mincoverage",
                                   "Minimum coverage in two neighbouring "
                                   "diagonal bands (for filter)\n"
                                   "default: 2.5 x seedlength",
                                   &arguments->dbs_mincoverage,
                                   GT_UWORD_MAX, 1UL);
  gt_option_hide_default(option);
  gt_option_parser_add_option(op, option);

  /* -maxfreq */
  op_frq = gt_option_new_uword_min("maxfreq",
                                   "Maximum frequency of a k-mer (for filter)",
                                   &arguments->dbs_maxfreq,
                                   GT_UWORD_MAX, 1UL);
  gt_option_parser_add_option(op, op_frq);

  /* -t */
  op_sup = gt_option_new_uword_min("t",
                                   "Suppress k-mers occurring at least t times "
                                   "(for filter)\n"
                                   "alias for maxfreq - 1",
                                   &arguments->dbs_suppress,
                                   GT_UWORD_MAX, 2UL);
  gt_option_exclude(op_sup, op_frq);
  gt_option_hide_default(op_sup);
  gt_option_is_development_option(op_sup);
  gt_option_parser_add_option(op, op_sup);

  /* -memlimit */
  op_mem = gt_option_new_string("memlimit",
                                "Maximum memory usage to determine the maximum "
                                "frequency of a k-mer (for filter)",
                                arguments->dbs_memlimit_str,
                                "");
  gt_option_parser_add_option(op, op_mem);

  /* -debug-kmer */
  option = gt_option_new_bool("debug-kmer",
                              "Output KmerPos lists",
                              &arguments->dbs_debug_kmer,
                              false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -debug-seedpair */
  option = gt_option_new_bool("debug-seedpair",
                              "Output SeedPair lists",
                              &arguments->dbs_debug_seedpair,
                              false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -verify */
  option = gt_option_new_bool("verify",
                              "Check that k-mer seeds occur in the sequences",
                              &arguments->dbs_verify,
                              false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* SEED EXTENSION OPTIONS */

  /* -extendxdrop */
  op_xdr = gt_option_new_uword_min_max("extendxdrop",
                                       "Extend seed to both sides using xdrop "
                                       "algorithm, /noptional parameter "
                                       "specifies sensitivity",
                                       &arguments->se_extendxdrop,
                                       97UL, 90UL, 100UL);
  gt_option_argument_is_optional(op_xdr);
  gt_option_parser_add_option(op, op_xdr);
  arguments->se_option_xdrop = gt_option_ref(op_xdr);

  /* -xdropbelow */
  op_xbe = gt_option_new_word("xdropbelow",
                              "Specify xdrop cutoff score (0 means "
                              "automatically defined depending on minidentity)",
                              &arguments->se_xdropbelowscore,
                              0L);
  gt_option_imply(op_xbe, op_xdr);
  gt_option_parser_add_option(op, op_xbe);

  /* -extendgreedy */
  op_gre = gt_option_new_uword_min_max("extendgreedy",
                                       "Extend seed to both sides using greedy "
                                       "algorithm, \noptional parameter "
                                       "specifies sensitivity",
                                       &arguments->se_extendgreedy,
                                       97UL, 90UL, 100UL);
  gt_option_argument_is_optional(op_gre);
  gt_option_exclude(op_gre, op_xdr);
  gt_option_parser_add_option(op, op_gre);
  arguments->se_option_greedy = gt_option_ref(op_gre);

  /* -only-seeds */
  op_onl = gt_option_new_bool("only-seeds",
                              "Calculate seeds and do not extend",
                              &arguments->onlyseeds,
                              false);
  gt_option_exclude(op_onl, op_xdr);
  gt_option_exclude(op_onl, op_gre);
  gt_option_is_development_option(op_onl);
  gt_option_parser_add_option(op, op_onl);

  /* -history */
  op_his = gt_option_new_uword_min_max("history",
                                       "Size of (mis)match history in range [1"
                                       "..64]\n(trimming for greedy extension)",
                                       &arguments->se_historysize,
                                       60UL, 1UL, 64UL);
  gt_option_exclude(op_his, op_onl);
  gt_option_exclude(op_his, op_xdr);
  gt_option_is_development_option(op_his);
  gt_option_parser_add_option(op, op_his);

  /* -maxalilendiff */
  op_dif = gt_option_new_uword("maxalilendiff",
                               "Maximum difference of alignment length\n"
                               "(trimming for greedy extension)",
                               &arguments->se_maxalilendiff, 0UL);
  gt_option_exclude(op_dif, op_onl);
  gt_option_exclude(op_dif, op_xdr);
  gt_option_hide_default(op_dif);
  gt_option_is_development_option(op_dif);
  gt_option_parser_add_option(op, op_dif);

  /* -percmathistory */
  op_pmh = gt_option_new_uword_min_max("percmathistory",
                                       "percentage of matches required in "
                                       "history \n(for greedy extension)",
                                       &arguments->se_perc_match_hist,
                                       0UL, 1UL, 100UL);
  gt_option_exclude(op_pmh, op_onl);
  gt_option_exclude(op_pmh, op_xdr);
  gt_option_hide_default(op_pmh);
  gt_option_is_development_option(op_pmh);
  gt_option_parser_add_option(op, op_pmh);

  /* -bias-parameters */
  op_bia = gt_option_new_bool("bias-parameters",
                              "Use -maxalilendiff 30 and let percmathistory "
                              "depend on minidentiy and DNA base distribution",
                              &arguments->bias_parameters,
                              false);
  gt_option_exclude(op_bia, op_onl);
  gt_option_exclude(op_bia, op_xdr);
  gt_option_exclude(op_bia, op_pmh);
  gt_option_exclude(op_bia, op_dif);
  gt_option_is_development_option(op_bia);
  gt_option_parser_add_option(op, op_bia);

  /* -cam */
  op_cam = gt_option_new_string("cam",
                                gt_cam_extendgreedy_comment(),
                                arguments->se_char_access_mode,
                                "");
  gt_option_hide_default(op_cam);
  gt_option_is_development_option(op_cam);
  gt_option_parser_add_option(op, op_cam);

  /* -l */
  op_len = gt_option_new_uword_min("l",
                                   "Minimum alignment length "
                                   "(for seed extension)",
                                   &arguments->se_alignlength,
                                   GT_UWORD_MAX, 1UL);
  gt_option_exclude(op_len, op_onl);
  gt_option_parser_add_option(op, op_len);

  /* -minidentity */
  op_err = gt_option_new_uword_min_max("minidentity",
                                       "Minimum identity of matches "
                                       "(for seed extension)",
                                       &arguments->se_minidentity,
                                       80UL, GT_EXTEND_MIN_IDENTITY_PERCENTAGE,
                                       99UL);
  gt_option_exclude(op_err, op_onl);
  gt_option_parser_add_option(op, op_err);

  /* -a */
  op_ali = gt_option_new_uword_min("a",
                                   "show alignments/sequences (optional "
                                   "argument is number of columns per line)",
                                   &arguments->se_alignmentwidth,
                                   70, 20);
  gt_option_exclude(op_ali, op_onl);
  gt_option_argument_is_optional(op_ali);
  gt_option_parser_add_option(op, op_ali);
  arguments->se_option_withali = gt_option_ref(op_ali);

  /* -relax-polish */
  op_relax_polish = gt_option_new_bool("relax-polish",
                                       "do not force alignments to have "
                                       "polished ends",
                                   &arguments->relax_polish,false);
  gt_option_parser_add_option(op, op_relax_polish);
  gt_option_is_development_option(op_relax_polish);
  gt_option_imply(op_relax_polish, op_ali);

  /* -no-reverse */
  op_norev = gt_option_new_bool("no-reverse",
                                "do not compute matches on reverse "
                                "complemented strand",
                                &arguments->norev,
                                false);
  gt_option_parser_add_option(op, op_norev);

  /* -no-forward */
  op_nofwd = gt_option_new_bool("no-forward",
                                "do not compute matches on forward strand",
                                &arguments->nofwd,
                                false);
  gt_option_exclude(op_nofwd, op_norev);
  gt_option_parser_add_option(op, op_nofwd);

  /* -overlappingseeds */
  option = gt_option_new_bool("overlappingseeds",
                              "Allow overlapping SeedPairs",
                              &arguments->overlappingseeds,
                              false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -benchmark */
  option = gt_option_new_bool("benchmark",
                              "Measure total running time and be silent",
                              &arguments->benchmark,
                              false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -seed-display */
  op_seed_display = gt_option_new_bool("seed-display",
                                       "Display seeds in #-line and by "
                                       "character + (instead of |) in middle "
                                       "row of alignment column",
                                       &arguments->seed_display,
                                       false);
  gt_option_exclude(op_seed_display, op_onl);
  gt_option_is_development_option(op_seed_display);
  gt_option_parser_add_option(op, op_seed_display);

  /* -weakends */
  op_weakends = gt_option_new_bool("weakends",
                                   "reduce minidentity for ends of seeded "
                                   "alignments",
                                   &arguments->weakends,
                                   false);
  gt_option_exclude(op_weakends, op_onl);
  gt_option_is_development_option(op_weakends);
  gt_option_parser_add_option(op, op_weakends);

  /* -extend-last */
  option = gt_option_new_bool("extend-last",
                              "Start extension after all SeedPair lists are "
                              "created",
                              &arguments->extend_last,
                              false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_seed_extend_arguments_check(int rest_argc, void *tool_arguments,
                                          GtError *err)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments != NULL);

  /* -t parameter as alias for maxfreq := t - 1 */
  if (arguments->dbs_suppress < GT_UWORD_MAX) {
    arguments->dbs_maxfreq = arguments->dbs_suppress - 1;
  }

  /* no alignment output */
  if (!gt_option_is_set(arguments->se_option_withali)) {
    arguments->se_alignmentwidth = 0;
  }

  /* parse memlimit argument */
  arguments->dbs_memlimit = GT_UWORD_MAX;
  if (strcmp(gt_str_get(arguments->dbs_memlimit_str), "") != 0) {
    had_err = gt_option_parse_spacespec(&arguments->dbs_memlimit,
                                        "memlimit",
                                        arguments->dbs_memlimit_str,
                                        err);
    if (!had_err && arguments->dbs_memlimit == 0) {
      gt_error_set(err,
                   "argument to option \"-memlimit\" must be at least 1MB");
      had_err = -1;
    }
  }

  /* minimum maxfreq value for 1 input file */
  if (!had_err && arguments->dbs_maxfreq == 1 &&
      strcmp(gt_str_get(arguments->dbs_queryname), "") == 0) {
    if (arguments->dbs_suppress == GT_UWORD_MAX) {
      gt_error_set(err, "argument to option \"-maxfreq\" must be >= 2 to "
                   "find matching k-mers");
    } else {
      gt_error_set(err, "argument to option \"-t\" must be >= 3 to find "
                   "matching k-mers");
    }
    had_err = -1;
  }

  /* no extra arguments */
  if (!had_err && rest_argc > 0) {
    gt_error_set(err, "too many arguments (-help shows correct usage)");
    had_err = -1;
  }
  return had_err;
}

static int gt_seed_extend_runner(GT_UNUSED int argc,
                                 GT_UNUSED const char **argv,
                                 GT_UNUSED int parsed_args,
                                 void *tool_arguments,
                                 GtError *err)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  GtEncseqLoader *encseq_loader = NULL;
  GtEncseq *aencseq = NULL, *bencseq = NULL;
  GtGreedyextendmatchinfo *grextinfo = NULL;
  GtXdropmatchinfo *xdropinfo = NULL;
  GtQuerymatchoutoptions *querymatchoutopt = NULL;
  GtTimer *seedextendtimer = NULL;
  GtExtendCharAccess cam = GT_EXTEND_CHAR_ACCESS_ANY;
  GtUword errorpercentage = 0UL;
  double matchscore_bias = GT_DEFAULT_MATCHSCORE_BIAS;
  bool extendgreedy = true;
  Polishing_info *pol_info = NULL;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments != NULL);
  gt_assert(arguments->se_minidentity >= GT_EXTEND_MIN_IDENTITY_PERCENTAGE &&
            arguments->se_minidentity <= 100UL);

  if (arguments->onlyseeds || gt_option_is_set(arguments->se_option_xdrop)) {
    extendgreedy = false;
  }

  if (arguments->dbs_verify || arguments->verbose) {
    int idx;
    bool minid_out = false, history_out = false;

    printf("# Options:");
    for (idx = 1; idx < argc; idx++) {
      if (strcmp(argv[idx],"-minidentity") == 0)
      {
        minid_out = true;
      }
      if (strcmp(argv[idx],"-history") == 0)
      {
        history_out = true;
      }
      printf(" %s", argv[idx]);
    }
    if (!minid_out)
    {
      printf(" -minidentity " GT_WU,arguments->se_minidentity);
    }
    if (!history_out)
    {
      printf(" -history " GT_WU,arguments->se_historysize);
    }
    printf("\n");
  }

  /* Calculate error percentage from minidentity */
  errorpercentage = 100UL - arguments->se_minidentity;

  /* Measure whole running time */
  if (arguments->benchmark || arguments->verbose) {
    gt_showtime_enable();
  }
  if (gt_showtime_enabled())
  {
    seedextendtimer = gt_timer_new();
    gt_timer_start(seedextendtimer);
  }

  /* Load encseq A */
  encseq_loader = gt_encseq_loader_new();
  gt_encseq_loader_enable_autosupport(encseq_loader);
  aencseq = gt_encseq_loader_load(encseq_loader,
                                  gt_str_get(arguments->dbs_indexname),
                                  err);
  if (aencseq == NULL)
    had_err = -1;

  /* If there is a 2nd read set: Load encseq B */
  if (!had_err) {
    if (strcmp(gt_str_get(arguments->dbs_queryname), "") != 0) {
      bencseq = gt_encseq_loader_load(encseq_loader,
                                      gt_str_get(arguments->dbs_queryname),
                                      err);
    } else {
      bencseq = gt_encseq_ref(aencseq);
    }
    if (bencseq == NULL) {
      had_err = -1;
      gt_encseq_delete(aencseq);
    }
  }
  gt_encseq_loader_delete(encseq_loader);

  if (!had_err && !gt_alphabet_is_dna(gt_encseq_alphabet(bencseq))) {
    if (arguments->nofwd) {
      gt_error_set(err, "option -no-forward is only allowed for DNA sequences");
      had_err = -1;
    } else {
      arguments->norev = true;
    }
  }

  /* set character access method */
  if (!had_err && (!arguments->onlyseeds ||
                   arguments->se_alignmentwidth > 0))
  {
    cam = gt_greedy_extend_char_access(gt_str_get
                                       (arguments->se_char_access_mode),
                                       err);
    if ((int) cam == -1) {
      had_err = -1;
      gt_encseq_delete(aencseq);
      gt_encseq_delete(bencseq);
    }
  }

  /* set seedlength */
  if (!had_err && arguments->dbs_seedlength == UINT_MAX) {
    unsigned int nchars = gt_alphabet_num_of_chars(gt_encseq_alphabet(aencseq));
    double totallength = 0.5 * (gt_encseq_total_length(aencseq) +
                                gt_encseq_total_length(bencseq));
    gt_assert(nchars > 0);
    arguments->dbs_seedlength
      = (unsigned int)gt_round_to_long(gt_log_base(totallength,(double)nchars));
    GT_UPDATE_MAX(arguments->dbs_seedlength, 2);
  }

  /* set mincoverage */
  if (!had_err && arguments->dbs_mincoverage == GT_UWORD_MAX) {
    arguments->dbs_mincoverage = (GtUword) (2.5 * arguments->dbs_seedlength);
  }

  /* set minimum alignment length */
  if (!had_err && arguments->se_alignlength == GT_UWORD_MAX) {
    arguments->se_alignlength = arguments->dbs_mincoverage;
  }

  /* Prepare options for greedy extension */
  if (!had_err && extendgreedy) {
    /* Use bias dependent parameters, adapted from E. Myers' DALIGNER */
    if (!had_err && arguments->bias_parameters)
    {
      matchscore_bias = gt_greedy_dna_sequence_bias_get(aencseq);
      arguments->se_maxalilendiff = 30;
      arguments->se_perc_match_hist
        = (GtUword) (100.0 - errorpercentage * matchscore_bias);
    }
    pol_info = polishing_info_new_with_bias(arguments->weakends
                                              ? MAX(errorpercentage,20)
                                              : errorpercentage,
                                            matchscore_bias,
                                            arguments->se_historysize);
    grextinfo = gt_greedy_extend_matchinfo_new(errorpercentage,
                                               arguments->se_maxalilendiff,
                                               arguments->se_historysize,
                                               arguments->se_perc_match_hist,
                                               arguments->se_alignlength,
                                               cam,
                                               arguments->se_extendgreedy,
                                               pol_info);
    if (arguments->benchmark) {
      gt_greedy_extend_matchinfo_silent_set(grextinfo);
    }
  }

  /* Prepare options for xdrop extension */
  if (!had_err && gt_option_is_set(arguments->se_option_xdrop)) {
    xdropinfo = gt_xdrop_matchinfo_new(arguments->se_alignlength,
                                       errorpercentage,
                                       arguments->se_xdropbelowscore,
                                       arguments->se_extendxdrop);
    if (arguments->benchmark) {
      gt_xdrop_matchinfo_silent_set(xdropinfo);
    }
  }

  /* Prepare output options */
  if (!had_err && (arguments->se_alignmentwidth > 0 ||
                   gt_option_is_set(arguments->se_option_xdrop)))
  {
    querymatchoutopt
      = gt_querymatchoutoptions_new(true,false,arguments->se_alignmentwidth);

    if (!arguments->onlyseeds)
    {
      GtUword sensitivity = 97UL;
      if (gt_option_is_set(arguments->se_option_greedy)) {
        sensitivity = arguments->se_extendgreedy;
      } else if (gt_option_is_set(arguments->se_option_xdrop)) {
        sensitivity = 100UL;
      }

      gt_querymatchoutoptions_extend(querymatchoutopt,
                                     errorpercentage,
                                     arguments->se_maxalilendiff,
                                     arguments->se_historysize,
                                     arguments->se_perc_match_hist,
                                     cam,
                                     arguments->weakends,
                                     sensitivity,
                                     matchscore_bias,
                                     !arguments->relax_polish,
                                     arguments->seed_display);
    }
  }

  /* Start algorithm */
  if (!had_err) {
    GtDiagbandseed dbsarguments;
    unsigned int maxseedlength;

    dbsarguments.errorpercentage = errorpercentage;
    dbsarguments.userdefinedleastlength = arguments->se_alignlength;
    dbsarguments.seedlength = arguments->dbs_seedlength;
    dbsarguments.logdiagbandwidth = arguments->dbs_logdiagbandwidth;
    dbsarguments.mincoverage = arguments->dbs_mincoverage;
    dbsarguments.maxfreq = arguments->dbs_maxfreq;
    dbsarguments.memlimit = arguments->dbs_memlimit;
    dbsarguments.norev = arguments->norev;
    dbsarguments.nofwd = arguments->nofwd;
    dbsarguments.overlappingseeds = arguments->overlappingseeds;
    dbsarguments.verify = arguments->dbs_verify;
    dbsarguments.verbose = arguments->verbose;
    dbsarguments.debug_kmer = arguments->dbs_debug_kmer;
    dbsarguments.debug_seedpair = arguments->dbs_debug_seedpair;
    dbsarguments.seed_display = arguments->seed_display;
    dbsarguments.extend_last = arguments->extend_last;
    dbsarguments.extendgreedyinfo = grextinfo;
    dbsarguments.extendxdropinfo = xdropinfo;
    dbsarguments.querymatchoutopt = querymatchoutopt;

    gt_assert(bencseq != NULL);
    if (gt_encseq_has_twobitencoding(aencseq) &&
        gt_encseq_wildcards(aencseq) == 0 &&
        gt_encseq_has_twobitencoding(bencseq) &&
        gt_encseq_wildcards(bencseq) == 0)
    {
      maxseedlength = 32;
    } else
    {
      unsigned int numofchars_a
        = gt_alphabet_num_of_chars(gt_encseq_alphabet(aencseq));
      unsigned int numofchars_b
        = gt_alphabet_num_of_chars(gt_encseq_alphabet(bencseq));
      if (numofchars_a != numofchars_b)
      {
        gt_error_set(err,"encoded sequences have different alphabet "
                         "sizes %u and %u",numofchars_a,numofchars_b);
        had_err = -1;
      }
      if (!had_err)
      {
        maxseedlength = gt_maxbasepower(numofchars_a) - 1;
        if (dbsarguments.seedlength > maxseedlength)
        {
          gt_error_set(err,"maximum seedlength for alphabet of size %u is %u",
                          numofchars_a,maxseedlength);
          had_err = -1;
        }
      }
    }
    if (!had_err)
    {
      had_err = gt_diagbandseed_run(aencseq, bencseq, &dbsarguments, err);
    }

    /* clean up */
    gt_encseq_delete(aencseq);
    gt_encseq_delete(bencseq);
    if (extendgreedy) {
      gt_greedy_extend_matchinfo_delete(grextinfo);
    }
    if (gt_option_is_set(arguments->se_option_xdrop)) {
      gt_xdrop_matchinfo_delete(xdropinfo);
    }
    if (arguments->se_alignmentwidth > 0 ||
        gt_option_is_set(arguments->se_option_xdrop)) {
      gt_querymatchoutoptions_delete(querymatchoutopt);
    }
  }
  polishing_info_delete(pol_info);

  if (gt_showtime_enabled()) {
    if (!had_err) {
      char *keystring
        = gt_seed_extend_params_keystring(extendgreedy,
                                          gt_option_is_set(arguments->
                                                           se_option_xdrop),
                                          arguments->dbs_seedlength,
                                          arguments->se_alignlength,
                                          arguments->se_minidentity,
                                          arguments->se_maxalilendiff,
                                          arguments->se_perc_match_hist,
                                          arguments->se_extendgreedy,
                                          arguments->se_extendxdrop,
                                          arguments->se_xdropbelowscore);
      printf("# TIME seedextend-%s", keystring);
      gt_free(keystring);
      gt_timer_show_formatted(seedextendtimer,
                              " overall " GT_WD ".%06ld\n",
                              stdout);
    }
    gt_timer_delete(seedextendtimer);
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
