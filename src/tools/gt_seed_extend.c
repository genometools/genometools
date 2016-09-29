/*
  Copyright (c) 2015-2016 Joerg Winkler <j.winkler@posteo.de>
  Copyright (c) 2015-2016 Center for Bioinformatics, University of Hamburg

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
  GtStr *dbs_pick_str;
  GtStr *dbs_memlimit_str;
  bool dbs_debug_kmer;
  bool dbs_debug_seedpair;
  bool dbs_verify;
  bool weakends;
  bool onlyseeds;
  bool overlappingseeds;
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
  GtStr *char_access_mode;
  bool bias_parameters;
  bool relax_polish;
  bool verify_alignment;
  /* general options */
  GtOption *se_option_withali;
  GtUword se_alignlength;
  GtUword se_minidentity;
  GtUword se_alignmentwidth;
  GtStrArray *display_args;
  bool norev;
  bool nofwd;
  bool benchmark;
  bool verbose;
  bool use_apos;
  bool histogram;
  bool use_kmerfile;
  unsigned int display_flag;
} GtSeedExtendArguments;

static void* gt_seed_extend_arguments_new(void)
{
  GtSeedExtendArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->dbs_indexname = gt_str_new();
  arguments->dbs_queryname = gt_str_new();
  arguments->dbs_pick_str = gt_str_new();
  arguments->dbs_memlimit_str = gt_str_new();
  arguments->char_access_mode = gt_str_new();
  arguments->display_args = gt_str_array_new();
  arguments->display_flag = 0;
  return arguments;
}

static void gt_seed_extend_arguments_delete(void *tool_arguments)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->dbs_indexname);
    gt_str_delete(arguments->dbs_queryname);
    gt_str_delete(arguments->dbs_pick_str);
    gt_str_delete(arguments->dbs_memlimit_str);
    gt_str_delete(arguments->char_access_mode);
    gt_option_delete(arguments->se_option_greedy);
    gt_option_delete(arguments->se_option_xdrop);
    gt_option_delete(arguments->se_option_withali);
    gt_str_array_delete(arguments->display_args);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_seed_extend_option_parser_new(void *tool_arguments)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *op_gre, *op_xdr, *op_cam, *op_his, *op_dif, *op_pmh,
    *op_len, *op_err, *op_xbe, *op_sup, *op_frq, *op_mem, *op_ali, *op_bia,
    *op_onl, *op_weakends, *op_relax_polish,
    *op_verify_alignment, *op_spdist, *op_display,
    *op_norev, *op_nofwd, *op_part, *op_pick, *op_overl;

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
                                arguments->char_access_mode,
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

  /* OUTPUT OPTIONS */

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

  /* -verify-alignment */
  op_verify_alignment
    = gt_option_new_bool("verify-alignment",
                         "verify alignment directly after its construction "
                         "(without knowning the sequences) and later (after the"
                         "sequence is known), in case the alignment is output",
                                   &arguments->verify_alignment,false);
  gt_option_parser_add_option(op, op_verify_alignment);
  gt_option_is_development_option(op_verify_alignment);

  /* -display */
  op_display = gt_option_new_string_array("display",
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
  gt_option_exclude(op_weakends, op_onl);
  gt_option_is_development_option(op_weakends);
  gt_option_parser_add_option(op, op_weakends);

  /* -use-apos */
  option = gt_option_new_bool("use-apos",
                              "Discard a seed only if both apos and bpos "
                              "overlap with previous alignment",
                              &arguments->use_apos,
                              false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

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

  /* no extra arguments */
  if (!had_err && rest_argc > 0) {
    gt_error_set(err, "too many arguments (-help shows correct usage)");
    had_err = -1;
  }
  return had_err;
}

static void gt_seed_extend_parts_variance_show(const GtRange *seqranges,
                                               GtUword numranges,
                                               const GtEncseq *encseq)
{
  GtUword variance_sum = 0, idx, avgpartsize;

  gt_assert(encseq != NULL);
  avgpartsize = gt_encseq_total_length(encseq)/numranges;
  for (idx = 0; idx < numranges; idx++)
  {
    GtUword segmentstart, segmentend, width;

    segmentstart = gt_encseq_seqstartpos(encseq,seqranges[idx].start);
    segmentend = gt_encseq_seqstartpos(encseq,seqranges[idx].end) +
                 gt_encseq_seqlength(encseq,seqranges[idx].end) - 1;
    gt_assert(segmentstart < segmentend);
    width = segmentend - segmentstart + 1;
    if (width > avgpartsize)
    {
      GtUword diff = width - avgpartsize;
      variance_sum += diff * diff;
    } else
    {
      GtUword diff = avgpartsize - width;
      variance_sum += diff * diff;
    }
    printf("# Part " GT_WU ": sequence " GT_WU "..." GT_WU " of total length "
           GT_WU "\n",idx+1,seqranges[idx].start,seqranges[idx].end,
                      segmentend - segmentstart + 1);
  }
  printf("# Variance of parts is %.2e\n",(double) variance_sum/numranges);
}

static GtUword gt_encseq_next_larger_width(const GtEncseq *encseq,
                                           GtUword startseqnum,
                                           GtUword width)
{
  GtUword left, right, found = GT_UWORD_MAX,
          numofsequences = gt_encseq_num_of_sequences(encseq);
  GtUword start_segment = gt_encseq_seqstartpos(encseq,startseqnum);

  left = startseqnum;
  gt_assert(numofsequences > 0);
  right = numofsequences - 1;
  while (left <= right)
  {
    GtUword mid = left + (right - left + 1)/2, mid_end, this_width;

    gt_assert(mid < numofsequences);
    mid_end = gt_encseq_seqstartpos(encseq,mid) +
              gt_encseq_seqlength(encseq,mid) - 1;
    gt_assert(mid_end > start_segment);
    this_width = mid_end - start_segment;
    if (this_width > width)
    {
      found = mid;
      if (right == 0)
      {
        break;
      }
      right = mid - 1;
    } else
    {
      if (left == numofsequences - 1)
      {
        break;
      }
      left = mid + 1;
    }
  }
  return found;
}

static GtUword gt_seed_extend_even_parts(GtRange *seqranges,
                                         const GtEncseq *encseq,
                                         GtUword maxseqnum,
                                         GtUword numparts)
{
  GtUword seqnum = 0, idx;
  const GtUword partwidth = gt_encseq_total_length(encseq)/numparts;

  for (idx = 0; idx < numparts && seqnum <= maxseqnum; idx++)
  {
    GtUword seqnum_next_width;

    seqnum_next_width = gt_encseq_next_larger_width(encseq,seqnum,partwidth);
    seqranges[idx].start = seqnum;
    if (seqnum_next_width == GT_UWORD_MAX)
    {
      seqranges[idx].end = maxseqnum;
      idx++;
      break;
    }
    seqranges[idx].end = seqnum_next_width;
    seqnum = seqnum_next_width + 1;
  }
  gt_assert(idx > 0 && seqranges[idx-1].end == maxseqnum);
  return idx;
}

/* Compute sequence ranges for specified number of parts. */
static GtUword gt_seed_extend_compute_parts(GtRange *seqranges,
                                            GtUword numparts,
                                            const GtEncseq *encseq)
{
  const GtUword maxseqnum = gt_encseq_num_of_sequences(encseq) - 1;
  gt_assert(seqranges != NULL);
  if (numparts > maxseqnum) { /* assign one seq for each part */
    GtUword idx;
    for (idx = 0; idx <= maxseqnum; ++idx) {
      seqranges[idx].start = idx;
      seqranges[idx].end = idx;
    }
    return maxseqnum + 1;
  } else {
    return gt_seed_extend_even_parts(seqranges, encseq, maxseqnum, numparts);
  }
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
  GtExtendCharAccess cam = GT_EXTEND_CHAR_ACCESS_ANY;
  GtUword errorpercentage = 0UL;
  double matchscore_bias = GT_DEFAULT_MATCHSCORE_BIAS;
  bool extendxdrop, extendgreedy = true;
  unsigned int maxseedlength = 0, nchars = 0;
  GtUwordPair pick = {GT_UWORD_MAX, GT_UWORD_MAX};
  GtUword maxseqlength = 0;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments != NULL);
  gt_assert(arguments->se_minidentity >= GT_EXTEND_MIN_IDENTITY_PERCENTAGE &&
            arguments->se_minidentity <= 100UL);

  /* Define, whether greedy extension will be performed */
  extendxdrop = gt_option_is_set(arguments->se_option_xdrop);
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
    if (!minid_out) {
      printf(" -minidentity " GT_WU,arguments->se_minidentity);
    }
    if (!history_out) {
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

  had_err = gt_querymatch_eval_display_args(&arguments->display_flag,
                                            arguments->display_args,
                                            err);

  /* Set character access method */
  if (!had_err && (!arguments->onlyseeds || arguments->se_alignmentwidth > 0)) {
    cam = gt_greedy_extend_char_access(gt_str_get(arguments->char_access_mode),
                                       err);
    if ((int) cam == -1) {
      had_err = -1;
    }
  }

  if (!had_err) {
    GtEncseqLoader *encseq_loader = gt_encseq_loader_new();
    gt_encseq_loader_require_multiseq_support(encseq_loader);

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

  if (arguments->dbs_seedlength == UINT_MAX) {
    unsigned int seedlength;
    double totallength = 0.5 * (gt_encseq_total_length(aencseq) +
                                gt_encseq_total_length(bencseq));
    gt_assert(nchars > 0);
    seedlength = (unsigned int)gt_round_to_long(gt_log_base(totallength,
                                                            (double)nchars));
    seedlength = (unsigned int)MIN3(seedlength, maxseqlength, maxseedlength);
    arguments->dbs_seedlength = MAX(seedlength, 2);
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
  if (!had_err && arguments->dbs_mincoverage == GT_UWORD_MAX) {
    arguments->dbs_mincoverage = (GtUword) (2.5 * arguments->dbs_seedlength);
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
    } else if (apick > gt_encseq_num_of_sequences(aencseq)) {
      gt_error_set(err, "first argument to option -pick must not be larger than"
                   " " GT_WU ", which is the number of sequences in the first "
                   "set", gt_encseq_num_of_sequences(aencseq));
      had_err = -1;
    } else if (bpick > gt_encseq_num_of_sequences(bencseq)) {
      gt_error_set(err, "second argument to option -pick must not be larger "
                   "than " GT_WU ", which is the number of sequences in the "
                   "second set", gt_encseq_num_of_sequences(bencseq));
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
      arguments->seedpairdistance.start = (GtUword)arguments->dbs_seedlength;
    }
  }

  /* Fill struct of algorithm arguments */
  if (!had_err) {
    GtDiagbandseedExtendParams *extp = NULL;
    GtDiagbandseedInfo *info = NULL;
    GtUword sensitivity = 0;
    GtUwordPair numparts = {arguments->dbs_parts, arguments->dbs_parts};
    GtRange *aseqranges = (GtRange *)gt_malloc(numparts.a * sizeof *aseqranges);
    GtRange *bseqranges = (GtRange *)gt_malloc(numparts.b * sizeof *bseqranges);

    if (extendgreedy) {
      sensitivity = arguments->se_extendgreedy;
    } else if (extendxdrop) {
      sensitivity = arguments->se_extendxdrop;
    }

    gt_assert(gt_encseq_num_of_sequences(aencseq) > 0);
    gt_assert(gt_encseq_num_of_sequences(bencseq) > 0);

    /* Get sequence ranges */
    numparts.a = gt_seed_extend_compute_parts(aseqranges,
                                              arguments->dbs_parts,
                                              aencseq);
    if (numparts.a > 1 && arguments->verbose)
    {
      gt_seed_extend_parts_variance_show(aseqranges, numparts.a, aencseq);
    }
    if (aencseq == bencseq && pick.a == pick.b)
    {
      gt_assert(numparts.b >= numparts.a);
      memcpy(bseqranges, aseqranges, numparts.a * sizeof *aseqranges);
      numparts.b = numparts.a;
    } else
    {
      numparts.b = gt_seed_extend_compute_parts(bseqranges,
                                                arguments->dbs_parts,
                                                bencseq);
      if (numparts.b > 1 && arguments->verbose)
      {
        gt_seed_extend_parts_variance_show(bseqranges, numparts.b, bencseq);
      }
    }
    gt_assert(pick.a < numparts.a || pick.a == GT_UWORD_MAX);
    gt_assert(pick.b < numparts.b || pick.b == GT_UWORD_MAX);

    extp = gt_diagbandseed_extend_params_new(errorpercentage,
                                             arguments->se_alignlength,
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
                                             cam,
                                             sensitivity,
                                             matchscore_bias,
                                             arguments->weakends,
                                             arguments->benchmark,
                                             arguments->se_alignmentwidth,
                                             !arguments->relax_polish,
                                             arguments->verify_alignment);

    info = gt_diagbandseed_info_new(aencseq,
                                    bencseq,
                                    arguments->dbs_maxfreq,
                                    arguments->dbs_memlimit,
                                    arguments->dbs_seedlength,
                                    arguments->norev,
                                    arguments->nofwd,
                                    &arguments->seedpairdistance,
                                    arguments->dbs_verify,
                                    arguments->verbose,
                                    arguments->dbs_debug_kmer,
                                    arguments->dbs_debug_seedpair,
                                    arguments->use_kmerfile,
                                    extp,
                                    numparts.a,
                                    numparts.b);

    /* Start algorithm */
    had_err = gt_diagbandseed_run(info,
                                  aseqranges,
                                  bseqranges,
                                  &pick,
                                  err);

    /* clean up */
    gt_free(aseqranges);
    gt_free(bseqranges);
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
