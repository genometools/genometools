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
#include "core/arraydef_api.h"
#include "core/cstr_api.h"
#include "core/cstr_array_api.h"
#include "core/encseq.h"
#include "core/encseq_api.h"
#include "core/error_api.h"
#include "core/ma_api.h"
#include "core/mathsupport_api.h"
#include "core/minmax_api.h"
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
#include "match/dbs_spaced_seeds.h"
#include "tools/gt_seed_extend.h"
#ifdef GT_THREADS_ENABLED
#include "core/thread_api.h"
#endif

typedef struct {
  /* diagbandseed options */
  GtStr *dbs_indexname;
  GtStr *dbs_queryname;
  unsigned int dbs_spacedseedweight;
  unsigned int dbs_seedlength;
  GtUword dbs_logdiagbandwidth;
  GtUword dbs_mincoverage;
  GtUword dbs_maxfreq;
  GtUword dbs_suppress;
  GtUword dbs_memlimit;
  GtUword dbs_parts;
  GtRange seedpairdistance;
  GtStr *dbs_pick_str,
        *diagband_statistics_arg,
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
  /* greedy extension options */
  GtUword se_extendgreedy;
  GtUword se_historysize;
  GtUword se_maxalilendiff;
  GtUword se_perc_match_hist;
  GtStr *char_access_mode, *splt_string, *kmplt_string;
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
  bool use_apos, use_apos_track_all, compute_ani;
  GtUword maxmat;
  GtOption *se_ref_op_evalue,
           *se_ref_op_spacedseed,
           *se_ref_op_maxmat,
           *ref_diagband_statistics,
           *se_ref_op_gre,
           *se_ref_op_xdr;
} GtSeedExtendArguments;

static void* gt_seed_extend_arguments_new(void)
{
  GtSeedExtendArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->dbs_indexname = gt_str_new();
  arguments->dbs_queryname = gt_str_new();
  arguments->dbs_pick_str = gt_str_new();
  arguments->chainarguments = gt_str_new();
  arguments->diagband_statistics_arg = gt_str_new();
  arguments->dbs_memlimit_str = gt_str_new();
  arguments->char_access_mode = gt_str_new();
  arguments->splt_string = gt_str_new();
  arguments->kmplt_string = gt_str_new();
  arguments->display_args = gt_str_array_new();
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
    gt_str_delete(arguments->diagband_statistics_arg);
    gt_str_delete(arguments->dbs_memlimit_str);
    gt_str_delete(arguments->char_access_mode);
    gt_str_delete(arguments->splt_string);
    gt_str_delete(arguments->kmplt_string);
    gt_option_delete(arguments->se_ref_op_gre);
    gt_option_delete(arguments->se_ref_op_spacedseed);
    gt_option_delete(arguments->se_ref_op_xdr);
    gt_option_delete(arguments->se_ref_op_evalue);
    gt_option_delete(arguments->se_ref_op_maxmat);
    gt_option_delete(arguments->ref_diagband_statistics);
    gt_str_array_delete(arguments->display_args);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_seed_extend_option_parser_new(void *tool_arguments)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *op_gre, *op_xdr, *op_cam, *op_splt, *op_kmplt,
    *op_his, *op_dif, *op_pmh,
    *op_seedlength, *op_spacedseed, *op_minlen, *op_minid, *op_evalue, *op_xbe,
    *op_sup, *op_frq,
    *op_mem, *op_bia, *op_onlyseeds, *op_weakends, *op_relax_polish,
    *op_verify_alignment, *op_only_selected_seqpairs, *op_spdist, *op_outfmt,
    *op_norev, *op_nofwd, *op_part, *op_pick, *op_overl, *op_trimstat,
    *op_cam_generic, *op_diagbandwidth, *op_mincoverage, *op_maxmat,
    *op_use_apos, *op_use_apos_track_all, *op_chain, *op_diagband_statistics,
    *op_ani, *op_benchmark;

  static GtRange seedpairdistance_defaults = {1UL, GT_UWORD_MAX};
  /* When extending the following array, do not forget to update
     the help message accordingly. */
  static const char *diagband_statistics_choices[] = {"sum", NULL};
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
                                      "with alphabet size as log-base",
                                      &arguments->dbs_seedlength,
                                      UINT_MAX, 1UL, 32UL);
  gt_option_hide_default(op_seedlength);
  gt_option_parser_add_option(op, op_seedlength);

  /* -spacedseed */
  op_spacedseed = gt_option_new_uint_min("spacedseed",
                                          "use spaced seed of length "
                                          "specified by option -seedlength\n"
                                          "(optional argument specifies "
                                          " weight of spaced seed)",
                                          &arguments->dbs_spacedseedweight,
                                          0,
                                          1);
  gt_option_argument_is_optional(op_spacedseed);
  gt_option_parser_add_option(op, op_spacedseed);
  arguments->se_ref_op_spacedseed = gt_option_ref(op_spacedseed);

  /* -diagbandwidth */
  op_diagbandwidth = gt_option_new_uword_min_max("diagbandwidth",
                               "Logarithm of diagonal band width in the "
                               "range\nfrom 0 to 10 (for filter)",
                               &arguments->dbs_logdiagbandwidth,
                               6UL,
                               0,
                               10);
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

  /* -diagband-stat */
  op_diagband_statistics = gt_option_new_choice("diagband-stat",
                                   "Compute statistics from diagonal "
                                   "band scores; parameter specifies "
                                   "kind of statistics, possible choices are\n"
                                   "sum",
                                   arguments->diagband_statistics_arg,
                                   diagband_statistics_choices[0],
                                   diagband_statistics_choices);
  gt_option_parser_add_option(op, op_diagband_statistics);
  arguments->ref_diagband_statistics = gt_option_ref(op_diagband_statistics);

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

  /* -kmplt */
  op_kmplt = gt_option_new_string("kmplt",
                                  gt_diagbandseed_kmplt_comment(),
                                  arguments->kmplt_string,
                                  "");
  gt_option_hide_default(op_kmplt);
  gt_option_is_development_option(op_kmplt);
  gt_option_parser_add_option(op, op_kmplt);

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
  /* will later be included again: gt_option_exclude(op_maxmat, op_bia); */
  /*gt_option_exclude(op_maxmat, op_diagbandwidth);
  gt_option_exclude(op_maxmat, op_mincoverage);
  gt_option_exclude(op_maxmat, op_xdr);
  gt_option_exclude(op_maxmat, op_xbe);
  gt_option_exclude(op_maxmat, op_gre);
  gt_option_exclude(op_maxmat, op_his);
  gt_option_exclude(op_maxmat, op_dif);
  gt_option_exclude(op_maxmat, op_pmh);
  gt_option_exclude(op_maxmat, op_cam);
  gt_option_exclude(op_maxmat, op_trimstat);*/
  gt_option_argument_is_optional(op_maxmat);
  gt_option_parser_add_option(op, op_maxmat);

  op_chain = gt_option_new_string("chain",
                                  "apply local chaining to maximal matches "
                                  "derived from k-mer seeds",
                                  arguments->chainarguments,"");
  gt_option_argument_is_optional(op_chain);
  gt_option_is_development_option(op_chain);
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
  op_outfmt = gt_option_new_string_array("outfmt",
                                          gt_querymatch_display_help(),
                                          arguments->display_args);
  gt_option_parser_add_option(op, op_outfmt);
  gt_option_exclude(op_outfmt,op_onlyseeds);

  /* -ani */
  op_ani = gt_option_new_bool("ani",
                              "output average nucleotide identity determined "
                              "from the computed matches "
                              "(which are not output)",
                              &arguments->compute_ani,
                              false);
  gt_option_parser_add_option(op, op_ani);

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
  op_benchmark = gt_option_new_bool("benchmark",
                                    "Measure total running time and be silent",
                                    &arguments->benchmark,
                                    false);
  gt_option_is_development_option(op_benchmark);
  gt_option_parser_add_option(op, op_benchmark);

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
  op_use_apos = gt_option_new_bool("use-apos",
                "Discard a seed only if both apos and bpos overlap with a "
                "previous successful alignment",
                &arguments->use_apos,false);
  gt_option_parser_add_option(op, op_use_apos);

  /* -use-apos-track-all */
  op_use_apos_track_all = gt_option_new_bool("use-apos-track-all",
                "Discard a seed only if both apos and bpos overlap with a "
                "previous alignment",
                &arguments->use_apos_track_all,false);
  gt_option_parser_add_option(op, op_use_apos_track_all);
  gt_option_is_development_option(op_use_apos_track_all);
  gt_option_exclude(op_use_apos, op_use_apos_track_all);

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

  gt_option_exclude(op_diagband_statistics, op_mincoverage);
  gt_option_exclude(op_diagband_statistics, op_xdr);
  gt_option_exclude(op_diagband_statistics, op_xbe);
  gt_option_exclude(op_diagband_statistics, op_gre);
  gt_option_exclude(op_diagband_statistics, op_his);
  gt_option_exclude(op_diagband_statistics, op_dif);
  gt_option_exclude(op_diagband_statistics, op_pmh);
  gt_option_exclude(op_diagband_statistics, op_cam);
  gt_option_exclude(op_diagband_statistics, op_trimstat);
  gt_option_exclude(op_diagband_statistics, op_maxmat);
  gt_option_exclude(op_diagband_statistics, op_use_apos);
  gt_option_exclude(op_diagband_statistics, op_use_apos_track_all);
  gt_option_exclude(op_diagband_statistics, op_minlen);
  gt_option_exclude(op_diagband_statistics, op_ani);
  gt_option_exclude(op_ani, op_outfmt);
  gt_option_exclude(op_ani, op_onlyseeds);
  gt_option_exclude(op_ani, op_verify_alignment);
  gt_option_exclude(op_ani, op_only_selected_seqpairs);
  gt_option_exclude(op_ani, op_benchmark);
  gt_option_exclude(op_cam, op_xbe);
  gt_option_exclude(op_cam, op_xdr);
  gt_option_exclude(op_cam_generic, op_xbe);
  gt_option_exclude(op_cam_generic, op_xdr);
  gt_option_exclude(op_maxmat, op_spacedseed);

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
#ifdef GT_THREADS_ENABLED
  if (!had_err && arguments->compute_ani && gt_jobs > 1)
  {
    gt_error_set(err,"option -ani does not work with multiple threads");
    had_err = -1;
  }
#endif

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
  if (!gt_option_is_set(arguments->ref_diagband_statistics))
  {
    gt_str_set(arguments->diagband_statistics_arg,"");
  }
  /* no extra arguments */
  if (!had_err && rest_argc > 0) {
    gt_error_set(err, "too many arguments (-help shows correct usage)");
    had_err = -1;
  }
  return had_err;
}

static double gt_seed_extend_ani_evaluate(GtUword sum_of_aligned_len,
                                          GtUword sum_of_distance)
{
  return sum_of_aligned_len > 0
             ? (100.0 * (1.0 - (double)
                               (2 * sum_of_distance)/sum_of_aligned_len))
             : 0.0;
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
  GtDiagbandseedBaseListType splt = GT_DIAGBANDSEED_BASE_LIST_UNDEFINED,
                             kmplt = GT_DIAGBANDSEED_BASE_LIST_UNDEFINED;
  GtUword errorpercentage = 0UL;
  double matchscore_bias = GT_DEFAULT_MATCHSCORE_BIAS;
  bool extendxdrop, extendgreedy = true;
  unsigned int maxseedlength = 0, nchars = 0;
  GtUwordPair pick = {GT_UWORD_MAX, GT_UWORD_MAX};
  GtUword maxseqlength = 0, a_numofsequences, b_numofsequences;
  GtSeedExtendDisplayFlag *out_display_flag = NULL;
  int had_err = 0;
  const GtSeedExtendDisplaySetMode setmode
    = GT_SEED_EXTEND_DISPLAY_SET_STANDARD;
  GtAniAccumulate ani_accumulate[2];

  gt_error_check(err);
  gt_assert(arguments != NULL);
  ani_accumulate[0].sum_of_aligned_len = 0;
  ani_accumulate[0].sum_of_distance = 0;
  ani_accumulate[1].sum_of_aligned_len = 0;
  ani_accumulate[1].sum_of_distance = 0;
  /* Define, whether greedy extension will be performed */
  extendxdrop = gt_option_is_set(arguments->se_ref_op_xdr);
  if (arguments->onlyseeds || extendxdrop) {
    extendgreedy = false;
  }

  /* Calculate error percentage from minidentity */
  gt_assert(arguments->se_minidentity >= GT_EXTEND_MIN_IDENTITY_PERCENTAGE &&
            arguments->se_minidentity <= 100UL);
  errorpercentage = 100UL - arguments->se_minidentity;

  /* Measure whole running time */
  if (arguments->benchmark || arguments->verbose)
  {
    gt_showtime_enable();
  }
  if (gt_showtime_enabled())
  {
    seedextendtimer = gt_timer_new();
    gt_timer_start(seedextendtimer);
  }
  if (!arguments->compute_ani)
  {
    out_display_flag = gt_querymatch_display_flag_new(arguments->display_args,
                                                      setmode,err);
    if (out_display_flag == NULL)
    {
      had_err = -1;
    }
  }

  if (!had_err)
  {
    if (!gt_querymatch_gfa2_display(out_display_flag))
    {
      const bool idhistout
        = (arguments->maxmat != 1 &&
           gt_str_length(arguments->diagband_statistics_arg) == 0)
          ? true : false;
      gt_querymatch_Options_output(stdout,argc,argv,idhistout,
                                   arguments->se_minidentity,
                                   arguments->se_historysize);
      if (!arguments->compute_ani  && !arguments->onlyseeds)
      {
        gt_querymatch_Fields_output(stdout,out_display_flag);
      }
    } else
    {
      printf("H\tVN:Z:2.0");
      if (gt_querymatch_trace_display(out_display_flag))
      {
        printf("\tTS:i:" GT_WU "\n",
               gt_querymatch_trace_delta_display(out_display_flag));
      } else
      {
        fputc('\n',stdout);
      }
    }
  }
  /* Set character access method */
  if (!had_err && !arguments->onlyseeds)
  {
    if (gt_greedy_extend_char_access(&cam_a,&cam_b,
                                     gt_str_get(arguments->char_access_mode),
                                     err) != 0)
    {
      had_err = -1;
    }
  }
  if (!had_err)
  {
    splt = gt_diagbandseed_base_list_get(true,
                                         gt_str_get(arguments->splt_string),
                                         err);
    if ((int) splt == -1) {
      had_err = -1;
    }
  }
  if (!had_err)
  {
    kmplt = gt_diagbandseed_base_list_get(false,
                                          gt_str_get(arguments->kmplt_string),
                                          err);
    if ((int) kmplt == -1) {
      had_err = -1;
    }
  }
  if (!had_err) {
    GtEncseqLoader *encseq_loader = gt_encseq_loader_new();
    gt_encseq_loader_require_multiseq_support(encseq_loader);
    gt_encseq_loader_require_ssp_tab(encseq_loader);
    if (out_display_flag != NULL &&
        gt_querymatch_subjectid_display(out_display_flag))
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
    }
    gt_encseq_loader_delete(encseq_loader);
  }
  if (!had_err)
  {
    /* If there is a 2nd read set: Load encseq B */
    if (gt_str_length(arguments->dbs_queryname) == 0) {
      bencseq = gt_encseq_ref(aencseq);
    } else
    {
      GtEncseqLoader *encseq_loader = gt_encseq_loader_new();
      gt_encseq_loader_require_multiseq_support(encseq_loader);
      gt_encseq_loader_require_ssp_tab(encseq_loader);
      if (out_display_flag != NULL &&
          gt_querymatch_queryid_display(out_display_flag))
      {
        gt_encseq_loader_require_des_tab(encseq_loader);
        gt_encseq_loader_require_sds_tab(encseq_loader);
      }
      bencseq = gt_encseq_loader_load(encseq_loader,
                                      gt_str_get(arguments->dbs_queryname),
                                      err);
      if (bencseq == NULL)
      {
        had_err = -1;
      }
      gt_encseq_loader_delete(encseq_loader);
    }
    if (bencseq == NULL) {
      had_err = -1;
      gt_encseq_delete(aencseq);
    }
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
    gt_querymatch_display_flag_delete(out_display_flag);
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
  maxseqlength = GT_MIN(gt_encseq_max_seq_length(aencseq),
                     gt_encseq_max_seq_length(bencseq));

  if (arguments->dbs_seedlength == UINT_MAX)
  {
    if (arguments->maxmat == 1)
    {
      arguments->dbs_seedlength = GT_MIN(maxseedlength,
                                         arguments->se_alignlength);
    } else
    {
      unsigned int local_seedlength, log_avg_totallength;
      double avg_totallength = 0.5 * (gt_encseq_total_length(aencseq) +
                                      gt_encseq_total_length(bencseq));
      gt_assert(nchars > 0);
      log_avg_totallength
        = (unsigned int) gt_round_to_long(gt_log_base(avg_totallength,
                                                      (double) nchars));
      local_seedlength = (unsigned int) GT_MIN3(log_avg_totallength,
                                             maxseqlength,maxseedlength);
      arguments->dbs_seedlength = GT_MAX(local_seedlength, 2);
    }
    if (gt_option_is_set(arguments->se_ref_op_spacedseed))
    {
      arguments->dbs_seedlength = GT_MIN(maxseedlength,
                                      (arguments->dbs_seedlength * 3)/2);
      arguments->dbs_seedlength = GT_MAX(arguments->dbs_seedlength,
                                      GT_SPACED_SEED_FIRST_SPAN);
    }
  }
  if (!had_err && arguments->dbs_seedlength >
                  GT_MIN(maxseedlength, maxseqlength))
  {
    if (maxseedlength <= maxseqlength) {
      gt_error_set(err, "maximum seedlength for alphabet of size %u is %u "
                        "(if the sequences %scontain wildcards)",
                          nchars, maxseedlength,maxseqlength == 32 ?
                                         "do not " : "");
    } else {
      gt_error_set(err, "argument to option \"-seedlength\" must be an integer "
                   "<= " GT_WU " (length of longest sequence).", maxseqlength);
    }
    had_err = -1;
  }
  if (!had_err)
  {
    if (gt_option_is_set(arguments->se_ref_op_spacedseed))
    {
      if (nchars != 4)
      {
        gt_error_set(err,"spaced seeds only work for sequences over an "
                         "alphabet of size 4");
        had_err = -1;
      } else
      {
        if (arguments->dbs_seedlength > maxseedlength ||
            arguments->dbs_seedlength < GT_SPACED_SEED_FIRST_SPAN)
        {
          gt_error_set(err,"illegal seedlength %u: for this set of sequences "
                           "can only handle spaced seeds of span between %d "
                           "and %u",
                           arguments->dbs_seedlength,
                           GT_SPACED_SEED_FIRST_SPAN,
                           maxseedlength);
          had_err = -1;
        } else
        {
          int min_weight, max_weight;

          gt_spaced_seed_weight_range(&min_weight,&max_weight,
                                      (int) arguments->dbs_seedlength);
          if (arguments->dbs_spacedseedweight == 0)
          {
            /* halfway between min and max */
            arguments->dbs_spacedseedweight = min_weight +
                                              (max_weight - min_weight + 1)/2;
          } else
          {
            if (arguments->dbs_spacedseedweight < (GtUword) min_weight ||
                arguments->dbs_spacedseedweight > (GtUword) max_weight)
            {
              gt_error_set(err,"illegal weight %u: for spaced seeds of span %u "
                               "the weight must be in the range from %d to %d",
                               arguments->dbs_spacedseedweight,
                               arguments->dbs_seedlength,min_weight,max_weight);
              had_err = -1;
            }
          }
        }
      }
    }
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
    GtUword use_apos_local = 0;

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

    if (arguments->use_apos)
    {
      gt_assert(!arguments->use_apos_track_all);
      use_apos_local = 1;
    } else
    {
      if (arguments->use_apos_track_all)
      {
        use_apos_local = 2;
      }
    }
    extp = gt_diagbandseed_extend_params_new(arguments->se_alignlength,
                                             errorpercentage,
                                             arguments->se_evalue_threshold,
                                             arguments->dbs_logdiagbandwidth,
                                             arguments->dbs_mincoverage,
                                             out_display_flag,
                                             use_apos_local,
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
                                             arguments->only_selected_seqpairs,
                                             arguments->compute_ani
                                               ? &ani_accumulate[0]
                                               : NULL);

    info = gt_diagbandseed_info_new(aencseq,
                                    bencseq,
                                    arguments->dbs_maxfreq,
                                    arguments->dbs_memlimit,
                                    arguments->dbs_spacedseedweight,
                                    arguments->dbs_seedlength,
                                    arguments->norev,
                                    arguments->nofwd,
                                    &arguments->seedpairdistance,
                                    splt,
                                    kmplt,
                                    arguments->dbs_verify,
                                    arguments->verbose,
                                    arguments->dbs_debug_kmer,
                                    arguments->dbs_debug_seedpair,
                                    arguments->use_kmerfile,
                                    arguments->trimstat_on,
                                    arguments->maxmat,
                                    arguments->chainarguments,
                                    arguments->diagband_statistics_arg,
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
  gt_querymatch_display_flag_delete(out_display_flag);
  if (arguments->compute_ani)
  {
    int idx;

    printf("ANI %s %s",gt_str_get(arguments->dbs_indexname),
                       gt_str_length(arguments->dbs_queryname) > 0
                         ? gt_str_get(arguments->dbs_queryname)
                         : gt_str_get(arguments->dbs_indexname));
    for (idx = 0; idx < 2; idx++)
    {
      printf(" %.4f",gt_seed_extend_ani_evaluate(
                          ani_accumulate[idx].sum_of_aligned_len,
                          ani_accumulate[idx].sum_of_distance));
    }
    printf("\n");
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
