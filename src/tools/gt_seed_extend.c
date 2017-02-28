/*
  Copyright (c) 2015-2016 Joerg Winkler <j.winkler@posteo.de>
  Copyright (c) 2016-2017 Stefan Kurtz  <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2015-2017 Center for Bioinformatics, University of Hamburg

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
#include <float.h>
#include "core/alphabet_api.h"
#include "core/arraydef.h"
#include "core/cstr_api.h"
#include "core/cstr_array.h"
#include "core/encseq.h"
#include "core/encseq_api.h"
#include "core/error_api.h"
#include "core/ma_api.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/parseutils_api.h"
#include "core/range_api.h"
#include "core/showtime.h"
#include "core/str_api.h"
#include "match/diagbandseed.h"
#include "match/seed-extend.h"
#include "match/xdrop.h"
#include "match/ft-polish.h"
#include "match/initbasepower.h"
#include "match/seed_extend_parts.h"
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
  GtUword dbs_parts;
  GtRange seedpairdistance;
  GtStr *dbs_pick_str,
        *chainarguments,
        *dbs_memlimit_str;
  bool dbs_debug_kmer;
  bool dbs_debug_seedpair;
  bool dbs_verify;
  bool weakends;
  bool onlyseeds;
  bool overlappingseeds;
  /* xdrop extension options */
  GtUword se_extendxdrop;
  GtXdropscore se_xdropbelowscore;
  GtOption *se_ref_op_xdr;
  /* greedy extension options */
  GtOption *se_ref_op_gre;
  GtUword se_extendgreedy;
  GtUword se_historysize;
  GtUword se_maxalilendiff;
  GtUword se_perc_match_hist;
  GtStr *char_access_mode, *splt_string, *use_apos_string;
  bool bias_parameters;
  bool relax_polish;
  bool verify_alignment;
  bool only_selected_seqpairs;
  bool cam_generic;
  /* general options */
  GtUword se_alignlength;
  GtUword se_minidentity;
  double se_evalue_threshold;
  GtStrArray *display_args;
  bool norev;
  bool nofwd;
  bool benchmark;
  bool verbose;
  bool histogram;
  bool use_kmerfile;
  bool trimstat_on;
  GtUword maxmat, use_apos;
  GtSeedExtendDisplayFlag *display_flag;
  GtOption *se_ref_op_evalue,
           *se_ref_op_maxmat,
           *se_ref_op_use_apos;
} GtSeedExtendArguments;

static void* gt_seed_extend_arguments_new(void)
{
  GtSeedExtendArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->dbs_indexname = gt_str_new();
  arguments->dbs_queryname = gt_str_new();
  arguments->dbs_pick_str = gt_str_new();
  arguments->chainarguments = gt_str_new();
  arguments->dbs_memlimit_str = gt_str_new();
  arguments->char_access_mode = gt_str_new();
  arguments->splt_string = gt_str_new();
  arguments->use_apos_string = gt_str_new();
  arguments->display_args = gt_str_array_new();
  arguments->display_flag = gt_querymatch_display_flag_new();
  return arguments;
}

static void gt_seed_extend_arguments_delete(void *tool_arguments)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->dbs_indexname);
    gt_str_delete(arguments->dbs_queryname);
    gt_str_delete(arguments->dbs_pick_str);
    gt_str_delete(arguments->chainarguments);
    gt_str_delete(arguments->dbs_memlimit_str);
    gt_str_delete(arguments->char_access_mode);
    gt_str_delete(arguments->splt_string);
    gt_str_delete(arguments->use_apos_string);
    gt_option_delete(arguments->se_ref_op_gre);
    gt_option_delete(arguments->se_ref_op_xdr);
    gt_option_delete(arguments->se_ref_op_evalue);
    gt_option_delete(arguments->se_ref_op_maxmat);
    gt_option_delete(arguments->se_ref_op_use_apos);
    gt_str_array_delete(arguments->display_args);
    gt_querymatch_display_flag_delete(arguments->display_flag);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_seed_extend_option_parser_new(void *tool_arguments)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *op_gre, *op_xdr, *op_cam, *op_splt,
    *op_his, *op_dif, *op_pmh,
    *op_seedlength, *op_minlen, *op_minid, *op_evalue, *op_xbe,
    *op_sup, *op_frq,
    *op_mem, *op_bia, *op_onlyseeds, *op_weakends, *op_relax_polish,
    *op_verify_alignment, *op_only_selected_seqpairs, *op_spdist, *op_display,
    *op_norev, *op_nofwd, *op_part, *op_pick, *op_overl, *op_trimstat,
    *op_cam_generic, *op_diagbandwidth, *op_mincoverage, *op_maxmat,
    *op_use_apos, *op_chain;

  static GtRange seedpairdistance_defaults = {1UL, GT_UWORD_MAX};
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
  op_seedlength = gt_option_new_uint_min_max("seedlength",
                                      "Minimum length of a seed\n"
                                      "default: logarithm of input length "
                                      "to the basis alphabet size",
                                      &arguments->dbs_seedlength,
                                      UINT_MAX, 1UL, 32UL);
  gt_option_hide_default(op_seedlength);
  gt_option_parser_add_option(op, op_seedlength);

  /* -diagbandwidth */
  op_diagbandwidth = gt_option_new_uword("diagbandwidth",
                               "Logarithm of diagonal band width (for filter)",
                               &arguments->dbs_logdiagbandwidth,
                               6UL);
  gt_option_parser_add_option(op, op_diagbandwidth);

  /* -mincoverage */
  op_mincoverage = gt_option_new_uword_min("mincoverage",
                                   "Minimum coverage in two neighbouring "
                                   "diagonal bands (for filter)\n"
                                   "default: 2.5 x seedlength",
                                   &arguments->dbs_mincoverage,
                                   GT_UWORD_MAX, 1UL);
  gt_option_hide_default(op_mincoverage);
  gt_option_parser_add_option(op, op_mincoverage);

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

  /* -extendxdrop */
  op_xdr = gt_option_new_uword_min_max("extendxdrop",
                                       "Extend seed to both sides using xdrop "
                                       "algorithm,\noptional parameter "
                                       "specifies sensitivity",
                                       &arguments->se_extendxdrop,
                                       97UL, 90UL, 100UL);
  gt_option_argument_is_optional(op_xdr);
  gt_option_parser_add_option(op, op_xdr);
  arguments->se_ref_op_xdr = gt_option_ref(op_xdr);

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
                                       "algorithm,\noptional parameter "
                                       "specifies sensitivity",
                                       &arguments->se_extendgreedy,
                                       97UL, 90UL, 100UL);
  gt_option_argument_is_optional(op_gre);
  gt_option_exclude(op_gre, op_xdr);
  gt_option_parser_add_option(op, op_gre);
  arguments->se_ref_op_gre = gt_option_ref(op_gre);

  /* -only-seeds */
  op_onlyseeds = gt_option_new_bool("only-seeds",
                              "Calculate seeds and do not extend",
                              &arguments->onlyseeds,
                              false);
  gt_option_exclude(op_onlyseeds, op_xdr);
  gt_option_exclude(op_onlyseeds, op_gre);
  gt_option_is_development_option(op_onlyseeds);
  gt_option_parser_add_option(op, op_onlyseeds);

  /* -history */
  op_his = gt_option_new_uword_min_max("history",
                                       "Size of (mis)match history in range [1"
                                       "..64]\n(trimming for greedy extension)",
                                       &arguments->se_historysize,
                                       60UL, 1UL, 64UL);
  gt_option_exclude(op_his, op_onlyseeds);
  gt_option_exclude(op_his, op_xdr);
  gt_option_is_development_option(op_his);
  gt_option_parser_add_option(op, op_his);

  /* -maxalilendiff */
  op_dif = gt_option_new_uword("maxalilendiff",
                               "Maximum difference of alignment length\n"
                               "(trimming for greedy extension)",
                               &arguments->se_maxalilendiff, 0UL);
  gt_option_exclude(op_dif, op_onlyseeds);
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
  gt_option_exclude(op_pmh, op_onlyseeds);
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
  gt_option_exclude(op_bia, op_onlyseeds);
  gt_option_exclude(op_bia, op_xdr);
  gt_option_exclude(op_bia, op_pmh);
  gt_option_exclude(op_bia, op_dif);
  gt_option_is_development_option(op_bia);
  gt_option_parser_add_option(op, op_bia);

  /* -cam */
  op_cam = gt_option_new_string("cam",
                                gt_cam_extendgreedy_comment(),
                                arguments->char_access_mode,
                                "bytes,bytes");
  gt_option_hide_default(op_cam);
  gt_option_is_development_option(op_cam);
  gt_option_parser_add_option(op, op_cam);

  op_cam_generic = gt_option_new_bool("cam_generic",
                                      "use generic function to access sequence",
                                      &arguments->cam_generic, false);
  gt_option_is_development_option(op_cam_generic);
  gt_option_parser_add_option(op, op_cam_generic);

  /* -splt */
  op_splt = gt_option_new_string("splt",
                                 gt_diagbandseed_splt_comment(),
                                 arguments->splt_string,
                                 "");
  gt_option_hide_default(op_splt);
  gt_option_is_development_option(op_splt);
  gt_option_parser_add_option(op, op_splt);

  /* -trimstat */
  op_trimstat = gt_option_new_bool("trimstat","show trimming statistics",
                                   &arguments->trimstat_on, false);
  gt_option_is_development_option(op_trimstat);
  gt_option_parser_add_option(op, op_trimstat);
  gt_option_exclude(op_trimstat, op_xdr);
  gt_option_exclude(op_trimstat, op_onlyseeds);

  /* -maxmat */
  op_maxmat = gt_option_new_ulong("maxmat",
                                 "compute maximal matches of minimum length "
                                 "specified by option -l",
                                 &arguments->maxmat,
                                 1);
  arguments->se_ref_op_maxmat = gt_option_ref(op_maxmat);
  gt_option_exclude(op_maxmat, op_diagbandwidth);
  gt_option_exclude(op_maxmat, op_mincoverage);
  gt_option_exclude(op_maxmat, op_xdr);
  gt_option_exclude(op_maxmat, op_xbe);
  gt_option_exclude(op_maxmat, op_gre);
  gt_option_exclude(op_maxmat, op_his);
  gt_option_exclude(op_maxmat, op_dif);
  gt_option_exclude(op_maxmat, op_pmh);
  /* will later be included again: gt_option_exclude(op_maxmat, op_bia); */
  gt_option_exclude(op_maxmat, op_cam);
  gt_option_exclude(op_maxmat, op_trimstat);
  gt_option_argument_is_optional(op_maxmat);
  gt_option_parser_add_option(op, op_maxmat);

  op_chain = gt_option_new_string("chain",
                                  "apply local chaining to maximal matches "
                                  "derived from k-mer seeds",
                                  arguments->chainarguments,"");
  gt_option_argument_is_optional(op_chain);
  gt_option_parser_add_option(op, op_chain);

  /* SEED EXTENSION OPTIONS */

  /* -l */
  op_minlen = gt_option_new_uword_min("l",
                                   "Minimum length of aligned sequences ",
                                   &arguments->se_alignlength,
                                   GT_UWORD_MAX, 1UL);
  gt_option_exclude(op_minlen, op_onlyseeds);
  gt_option_parser_add_option(op, op_minlen);
  gt_option_imply(op_maxmat, op_minlen);

  /* -minidentity */
  op_minid = gt_option_new_uword_min_max("minidentity",
                                       "Minimum identity of matches "
                                       "(for seed extension)",
                                       &arguments->se_minidentity,
                                       80UL, GT_EXTEND_MIN_IDENTITY_PERCENTAGE,
                                       99UL);
  gt_option_exclude(op_minid, op_onlyseeds);
  /* will later be included again: gt_option_exclude(op_maxmat, op_minid); */
  gt_option_parser_add_option(op, op_minid);

  /* -evalue */
  op_evalue = gt_option_new_double("evalue","switch on evalue filtering of "
                                            "matches (optional argument "
                                            "specifies evalue threshold)",
                                   &arguments->se_evalue_threshold,
                                   10.0);
  gt_option_exclude(op_evalue, op_onlyseeds);
  gt_option_parser_add_option(op, op_evalue);
  gt_option_argument_is_optional(op_evalue);
  gt_option_exclude(op_maxmat, op_evalue);
  arguments->se_ref_op_evalue = gt_option_ref(op_evalue);

  /* OUTPUT OPTIONS */

  /* -relax-polish */
  op_relax_polish = gt_option_new_bool("relax-polish",
                                       "do not force alignments to have "
                                       "polished ends",
                                   &arguments->relax_polish,false);
  gt_option_parser_add_option(op, op_relax_polish);
  /*gt_option_exclude(op_maxmat, op_relax_polish);*/
  gt_option_is_development_option(op_relax_polish);

  /* -verify-alignment */
  op_verify_alignment
    = gt_option_new_bool("verify-alignment",
                         "verify alignment directly after its construction "
                         "(without knowning the sequences) and later (after the"
                         "sequence is known), in case the alignment is output",
                                   &arguments->verify_alignment,false);
  gt_option_parser_add_option(op, op_verify_alignment);
  gt_option_exclude(op_maxmat, op_verify_alignment);
  gt_option_is_development_option(op_verify_alignment);

  /* -only-selected-seqpairs */
  op_only_selected_seqpairs
    = gt_option_new_bool("only-selected-seqpairs",
                         "output only sequence pair numbers with "
                         "selected seeds",
                         &arguments->only_selected_seqpairs,false);
  gt_option_parser_add_option(op, op_only_selected_seqpairs);
  gt_option_is_development_option(op_only_selected_seqpairs);

  /* -outfmt */
  op_display = gt_option_new_string_array("outfmt",
                                          gt_querymatch_display_help(),
                                          arguments->display_args);
  gt_option_parser_add_option(op, op_display);

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
  op_overl = gt_option_new_bool("overlappingseeds",
                                "Allow overlapping SeedPairs",
                                &arguments->overlappingseeds,
                                false);
  gt_option_is_development_option(op_overl);
  gt_option_parser_add_option(op, op_overl);

  /* -seedpairdistance */
  op_spdist = gt_option_new_range("seedpairdistance",
                                  "Only use SeedPairs whose positions differ "
                                  "in accordance with the specified range",
                                  &arguments->seedpairdistance,
                                  &seedpairdistance_defaults);
  gt_option_exclude(op_spdist, op_overl);
  gt_option_is_development_option(op_spdist);
  gt_option_hide_default(op_spdist);
  gt_option_parser_add_option(op, op_spdist);

  /* -benchmark */
  option = gt_option_new_bool("benchmark",
                              "Measure total running time and be silent",
                              &arguments->benchmark,
                              false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -weakends */
  op_weakends = gt_option_new_bool("weakends",
                                   "reduce minidentity for ends of seeded "
                                   "alignments",
                                   &arguments->weakends,
                                   false);
  gt_option_exclude(op_weakends, op_onlyseeds);
  gt_option_is_development_option(op_weakends);
  /*gt_option_exclude(op_maxmat, op_weakends);*/
  gt_option_parser_add_option(op, op_weakends);

  /* -use-apos */
  op_use_apos = gt_option_new_string("use-apos",
                "Discard a seed only if both apos and bpos overlap with a "
                "previous alignment\nin the default case all previous "
                "successful alignments are considered (which maximizes "
                "sensitivity but is slower); the optional parameter "
                "\"trackall\" means to consider all previous alignments, not "
                "only those that are successful",
                arguments->use_apos_string,
                "");
  /*gt_option_exclude(op_maxmat, op_use_apos);*/
  gt_option_argument_is_optional(op_use_apos);
  arguments->se_ref_op_use_apos = gt_option_ref(op_use_apos);
  gt_option_parser_add_option(op, op_use_apos);

  /* -parts */
  op_part = gt_option_new_uword_min("parts",
                                    "Divide data into specified number of "
                                    "parts",
                                    &arguments->dbs_parts,
                                    1, 1UL);
  gt_option_parser_add_option(op, op_part);

  /* -pick */
  op_pick = gt_option_new_string("pick",
                                 "Choose parts for 1st/2nd sequence set. "
                                 "Format: i,j",
                                 arguments->dbs_pick_str,
                                 "use all combinations successively");
  gt_option_imply(op_pick, op_part);
  gt_option_is_development_option(op_pick);
  gt_option_parser_add_option(op, op_pick);

  /* -histogram */
  option = gt_option_new_bool("histogram",
                              "Calculate histogram to determine size of mlist",
                              &arguments->histogram,
                              true);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -kmerfile */
  option = gt_option_new_bool("kmerfile",
                              "Use pre-calculated k-mers from file (if exist)",
                              &arguments->use_kmerfile,
                              true);
  gt_option_parser_add_option(op, option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);
  gt_option_exclude(op_cam, op_xbe);
  gt_option_exclude(op_cam, op_xdr);
  gt_option_exclude(op_cam_generic, op_xbe);
  gt_option_exclude(op_cam_generic, op_xdr);
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

  /* parse memlimit argument */
  arguments->dbs_memlimit = GT_UWORD_MAX;
  if (arguments->histogram == true) {
    arguments->dbs_memlimit -= 1;
  }
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

  if (!gt_option_is_set(arguments->se_ref_op_evalue))
  {
    arguments->se_evalue_threshold = DBL_MAX;
  }
  if (!gt_option_is_set(arguments->se_ref_op_maxmat))
  {
    arguments->maxmat = 0;
  }
  arguments->use_apos = 0;
  if (gt_option_is_set(arguments->se_ref_op_use_apos))
  {
    if (gt_str_length(arguments->use_apos_string) == 0)
    {
      arguments->use_apos = 1;
    } else
    {
      if (strcmp(gt_str_get(arguments->use_apos_string),"trackall") == 0)
      {
        arguments->use_apos = 2;
      } else
      {
        gt_error_set(err, "optional argument to option -use_apos must be "
                          "\"trackall\"");
        had_err = -1;
      }
    }
  }
  /* no extra arguments */
  if (!had_err && rest_argc > 0) {
    gt_error_set(err, "too many arguments (-help shows correct usage)");
    had_err = -1;
  }
  return had_err;
}

static int gt_seed_extend_runner(int argc,
                                 const char **argv,
                                 GT_UNUSED int parsed_args,
                                 void *tool_arguments,
                                 GtError *err)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  GtEncseq *aencseq = NULL, *bencseq = NULL;
  GtTimer *seedextendtimer = NULL;
  GtExtendCharAccess cam_a = GT_EXTEND_CHAR_ACCESS_ANY,
                     cam_b = GT_EXTEND_CHAR_ACCESS_ANY;
  GtDiagbandseedPairlisttype splt = GT_DIAGBANDSEED_SPLT_UNDEFINED;
  GtUword errorpercentage = 0UL;
  double matchscore_bias = GT_DEFAULT_MATCHSCORE_BIAS;
  bool extendxdrop, extendgreedy = true;
  unsigned int maxseedlength = 0, nchars = 0;
  GtUwordPair pick = {GT_UWORD_MAX, GT_UWORD_MAX};
  GtUword maxseqlength = 0, a_numofsequences, b_numofsequences;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments != NULL);

  /* Define, whether greedy extension will be performed */
  extendxdrop = gt_option_is_set(arguments->se_ref_op_xdr);
  if (arguments->onlyseeds || extendxdrop) {
    extendgreedy = false;
  }

  /* Print verbose option string */
  if (arguments->verbose) {
    int idx;
    bool minid_out = false, history_out = false;

    printf("# Options:");
    for (idx = 1; idx < argc; idx++) {
      if (strcmp(argv[idx],"-minidentity") == 0) {
        minid_out = true;
      }
      if (strcmp(argv[idx],"-history") == 0) {
        history_out = true;
      }
      printf(" %s", argv[idx]);
    }
    if (arguments->maxmat != 1)
    {
      if (!minid_out)
      {
        printf(" -minidentity " GT_WU,arguments->se_minidentity);
      }
      if (!history_out)
      {
        printf(" -history " GT_WU,arguments->se_historysize);
      }
    }
    printf("\n");
  }

  /* Calculate error percentage from minidentity */
  gt_assert(arguments->se_minidentity >= GT_EXTEND_MIN_IDENTITY_PERCENTAGE &&
            arguments->se_minidentity <= 100UL);
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
  had_err = gt_querymatch_display_flag_args_set(arguments->display_flag,
                                                arguments->display_args,
                                                err);
  /* Set character access method */
  if (!had_err && (!arguments->onlyseeds ||
                   gt_querymatch_display_alignmentwidth(
                        arguments->display_flag) > 0))
  {
    if (gt_greedy_extend_char_access(&cam_a,&cam_b,
                                     gt_str_get(arguments->char_access_mode),
                                     err) != 0)
    {
      had_err = -1;
    }
  }
  gt_querymatch_display_seedpos_relative_set(
             arguments->display_flag,
             cam_a == GT_EXTEND_CHAR_ACCESS_DIRECT ? true : false,
             cam_b == GT_EXTEND_CHAR_ACCESS_DIRECT ? true : false);
  if (!had_err)
  {
    splt = gt_diagbandseed_splt_get(gt_str_get(arguments->splt_string),err);
    if ((int) splt == -1) {
      had_err = -1;
    }
  }

  if (!had_err) {
    GtEncseqLoader *encseq_loader = gt_encseq_loader_new();
    gt_encseq_loader_require_multiseq_support(encseq_loader);
    gt_encseq_loader_require_ssp_tab(encseq_loader);
    if (gt_querymatch_seqdesc_display(arguments->display_flag))
    {
      gt_encseq_loader_require_des_tab(encseq_loader);
      gt_encseq_loader_require_sds_tab(encseq_loader);
    }

    /* Load encseq A */
    aencseq = gt_encseq_loader_load(encseq_loader,
                                    gt_str_get(arguments->dbs_indexname),
                                    err);
    if (aencseq == NULL) {
      had_err = -1;
    } else {
      /* If there is a 2nd read set: Load encseq B */
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
  }

  /* Check alphabet sizes */
  if (!had_err) {
    unsigned int nchars_b;
    nchars = gt_alphabet_num_of_chars(gt_encseq_alphabet(aencseq));
    nchars_b = gt_alphabet_num_of_chars(gt_encseq_alphabet(bencseq));
    if (nchars != nchars_b) {
      gt_error_set(err,"encoded sequences have different alphabet "
                   "sizes %u and %u", nchars, nchars_b);
      gt_encseq_delete(aencseq);
      gt_encseq_delete(bencseq);
      had_err = -1;
    }
  }

  if (had_err) {
    if (gt_showtime_enabled()) {
      gt_timer_delete(seedextendtimer);
    }
    return had_err;
  }

  /* Set seedlength */
  gt_assert(aencseq != NULL && bencseq != NULL);
  a_numofsequences = gt_encseq_num_of_sequences(aencseq);
  b_numofsequences = gt_encseq_num_of_sequences(bencseq);
  if (gt_encseq_has_twobitencoding(aencseq) &&
      gt_encseq_wildcards(aencseq) == 0 &&
      gt_encseq_has_twobitencoding(bencseq) &&
      gt_encseq_wildcards(bencseq) == 0) {
    maxseedlength = 32;
  } else {
    maxseedlength = gt_maxbasepower(nchars) - 1;
  }
  maxseqlength = MIN(gt_encseq_max_seq_length(aencseq),
                     gt_encseq_max_seq_length(bencseq));

  if (arguments->dbs_seedlength == UINT_MAX)
  {
    if (arguments->maxmat == 1)
    {
      arguments->dbs_seedlength = MIN(maxseedlength, arguments->se_alignlength);
    } else
    {
      unsigned int seedlength;
      double totallength = 0.5 * (gt_encseq_total_length(aencseq) +
                                  gt_encseq_total_length(bencseq));
      gt_assert(nchars > 0);
      seedlength = (unsigned int)gt_round_to_long(gt_log_base(totallength,
                                                              (double)nchars));
      seedlength = (unsigned int)MIN3(seedlength, maxseqlength, maxseedlength);
      arguments->dbs_seedlength = MAX(seedlength, 2);
    }
  }
  if (arguments->dbs_seedlength > MIN(maxseedlength, maxseqlength)) {
    if (maxseedlength <= maxseqlength) {
      gt_error_set(err, "maximum seedlength for alphabet of size %u is %u",
                   nchars, maxseedlength);
    } else {
      gt_error_set(err, "argument to option \"-seedlength\" must be an integer "
                   "<= " GT_WU " (length of longest sequence).", maxseqlength);
    }
    had_err = -1;
  }

  /* Set mincoverage */
  if (!had_err)
  {
    if (arguments->dbs_mincoverage == GT_UWORD_MAX)
    {
      arguments->dbs_mincoverage = (GtUword) (2.5 * arguments->dbs_seedlength);
    } else
    {
      if (arguments->dbs_mincoverage < arguments->dbs_seedlength)
      {
        gt_error_set(err, "argument to option \"-mincoverage\" must be an "
                          "integer >= %u (seedlength).",
                          arguments->dbs_seedlength);
        had_err = -1;
      }
    }
  }

  /* Set minimum alignment length */
  if (!had_err && arguments->se_alignlength == GT_UWORD_MAX) {
    arguments->se_alignlength = arguments->dbs_mincoverage;
  }

  /* Check alphabet and direction compatibility */
  if (!had_err && !gt_alphabet_is_dna(gt_encseq_alphabet(bencseq))) {
    if (arguments->nofwd) {
      gt_error_set(err, "option -no-forward is only allowed for DNA sequences");
      had_err = -1;
    } else {
      arguments->norev = true; /* reverse is just for DNA */
    }
  }

  /* Parse pick option */
  if (!had_err && strcmp(gt_str_get(arguments->dbs_pick_str),
                         "use all combinations successively") != 0) {
    GtUword apick, bpick;
    char **items = gt_cstr_split(gt_str_get(arguments->dbs_pick_str), ',');
    if (gt_cstr_array_size((const char **)items) != 2 ||
        gt_parse_uword(&apick, items[0]) != 0 ||
        gt_parse_uword(&bpick, items[1]) != 0) {
      gt_error_set(err, "argument to option -pick must satisfy format i,j");
      had_err = -1;
    } else if (apick > arguments->dbs_parts || bpick > arguments->dbs_parts) {
      gt_error_set(err, "arguments to option -pick must not exceed " GT_WU
                   " (number of parts)", arguments->dbs_parts);
      had_err = -1;
    } else if (apick < 1 || bpick < 1) {
      gt_error_set(err, "arguments to option -pick must be at least 1");
      had_err = -1;
    } else if (apick > a_numofsequences) {
      gt_error_set(err, "first argument to option -pick must not be larger than"
                   " " GT_WU ", which is the number of sequences in the first "
                   "set", a_numofsequences);
      had_err = -1;
    } else if (bpick > b_numofsequences) {
      gt_error_set(err, "second argument to option -pick must not be larger "
                   "than " GT_WU ", which is the number of sequences in the "
                   "second set", b_numofsequences);
      had_err = -1;
    } else if (aencseq == bencseq && apick > bpick) {
      pick.a = bpick - 1;
      pick.b = apick - 1;
    } else {
      pick.a = apick - 1;
      pick.b = bpick - 1;
    }
    gt_cstr_array_delete(items);
  }

  /* Use bias dependent parameters, adapted from E. Myers' DALIGNER */
  if (!had_err && extendgreedy && arguments->bias_parameters) {
    matchscore_bias = gt_greedy_dna_sequence_bias_get(aencseq);
    arguments->se_maxalilendiff = 30;
    arguments->se_perc_match_hist = (GtUword) (100.0 - errorpercentage *
                                               matchscore_bias);
  }

  /* Set SeedPair distance according to overlappingseeds flag */
  if (!had_err && arguments->seedpairdistance.end == GT_UWORD_MAX) {
    arguments->seedpairdistance.end -= gt_encseq_max_seq_length(aencseq);
    if (!arguments->overlappingseeds &&
        arguments->seedpairdistance.start == 1UL) {
      arguments->seedpairdistance.start = (GtUword) arguments->dbs_seedlength;
    }
  }

  /* Fill struct of algorithm arguments */
  if (!had_err) {
    GtDiagbandseedExtendParams *extp = NULL;
    GtDiagbandseedInfo *info = NULL;
    GtUword sensitivity = 0;
    GtSequencePartsInfo *aseqranges, *bseqranges;

    if (extendgreedy) {
      sensitivity = arguments->se_extendgreedy;
    } else if (extendxdrop) {
      sensitivity = arguments->se_extendxdrop;
    }

    /* Get sequence ranges */
    aseqranges = gt_sequence_parts_info_new(aencseq,a_numofsequences,
                                            arguments->dbs_parts);
    if (arguments->verbose && gt_sequence_parts_info_number(aseqranges) > 1)
    {
      gt_sequence_parts_info_variance_show(aseqranges);
    }
    if (aencseq == bencseq)
    {
      bseqranges = aseqranges;
    } else
    {
      bseqranges = gt_sequence_parts_info_new(bencseq,b_numofsequences,
                                              arguments->dbs_parts);
      if (arguments->verbose && gt_sequence_parts_info_number(bseqranges) > 1)
      {
        gt_sequence_parts_info_variance_show(bseqranges);
      }
    }
    gt_assert(pick.a < gt_sequence_parts_info_number(aseqranges) ||
              pick.a == GT_UWORD_MAX);
    gt_assert(pick.b < gt_sequence_parts_info_number(bseqranges) ||
              pick.b == GT_UWORD_MAX);

    extp = gt_diagbandseed_extend_params_new(arguments->se_alignlength,
                                             errorpercentage,
                                             arguments->se_evalue_threshold,
                                             arguments->dbs_logdiagbandwidth,
                                             arguments->dbs_mincoverage,
                                             arguments->display_flag,
                                             arguments->use_apos,
                                             arguments->se_xdropbelowscore,
                                             extendgreedy,
                                             extendxdrop,
                                             arguments->se_maxalilendiff,
                                             arguments->se_historysize,
                                             arguments->se_perc_match_hist,
                                             cam_a,
                                             cam_b,
                                             arguments->cam_generic,
                                             sensitivity,
                                             matchscore_bias,
                                             arguments->weakends,
                                             arguments->benchmark,
                                             !arguments->relax_polish,
                                             arguments->verify_alignment,
                                             arguments->only_selected_seqpairs);

    info = gt_diagbandseed_info_new(aencseq,
                                    bencseq,
                                    arguments->dbs_maxfreq,
                                    arguments->dbs_memlimit,
                                    arguments->dbs_seedlength,
                                    arguments->norev,
                                    arguments->nofwd,
                                    &arguments->seedpairdistance,
                                    splt,
                                    arguments->dbs_verify,
                                    arguments->verbose,
                                    arguments->dbs_debug_kmer,
                                    arguments->dbs_debug_seedpair,
                                    arguments->use_kmerfile,
                                    arguments->trimstat_on,
                                    arguments->maxmat,
                                    arguments->chainarguments,
                                    extp);

    /* Start algorithm */
    had_err = gt_diagbandseed_run(info,
                                  aseqranges,
                                  bseqranges,
                                  &pick,
                                  err);

    /* clean up */
    if (bseqranges != aseqranges)
    {
      gt_sequence_parts_info_delete(bseqranges);
    }
    gt_sequence_parts_info_delete(aseqranges);
    gt_diagbandseed_extend_params_delete(extp);
    gt_diagbandseed_info_delete(info);
  }
  gt_encseq_delete(aencseq);
  gt_encseq_delete(bencseq);

  if (!had_err && gt_showtime_enabled()) {
    char *keystring;
    keystring = gt_seed_extend_params_keystring(extendgreedy,
                                                extendxdrop,
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
                            " overall " GT_WD ".%06ld\n", stdout);
  }
  if (gt_showtime_enabled()) {
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
