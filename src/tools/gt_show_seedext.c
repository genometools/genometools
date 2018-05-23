/*
  Copyright (c) 2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2015 Joerg Winkler <j.winkler@posteo.de>
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "core/ma_api.h"
#include "core/minmax.h"
#include "core/str_api.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "core/encseq.h"
#include "core/showtime.h"
#include "core/timer_api.h"
#include "match/seed-extend.h"
#include "match/seq_or_encseq.h"
#include "match/seed-extend-iter.h"
#include "tools/gt_show_seedext.h"

typedef struct
{
  bool relax_polish,
       sortmatches,
       verify_alignment,
       affine_alignment;
  GtStr *matchfilename;
  GtStrArray *display_args;
} GtShowSeedextArguments;

static void* gt_show_seedext_arguments_new(void)
{
  GtShowSeedextArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->matchfilename = gt_str_new();
  arguments->display_args = gt_str_array_new();
  return arguments;
}

static void gt_show_seedext_arguments_delete(void *tool_arguments)
{
  GtShowSeedextArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_array_delete(arguments->display_args);
    gt_str_delete(arguments->matchfilename);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_show_seedext_option_parser_new(void *tool_arguments)
{
  GtShowSeedextArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option_filename, *op_relax_polish, *op_sortmatches, *op_display,
           *op_verify_alignment, *op_affine_alignment;

  gt_assert(arguments);
  /* init */
  op = gt_option_parser_new("[options] -f <matchfilename>",
                            "Parse output of a seed extension and show/verify "
                            "the alignment.");

  /* -outfmt */
  op_display = gt_option_new_string_array("outfmt",
                                          gt_querymatch_display_help(),
                                          arguments->display_args);
  gt_option_parser_add_option(op, op_display);

  /* -relax-polish */
  op_relax_polish = gt_option_new_bool("relax-polish",
                                       "do not force alignments to have "
                                       "polished ends",
                                       &arguments->relax_polish,false);
  gt_option_parser_add_option(op, op_relax_polish);
  gt_option_is_development_option(op_relax_polish);

  /* -sort */
  op_sortmatches = gt_option_new_bool("sort","sort matches in ascending order "
                                             "of their end position on the "
                                             "query",
                                      &arguments->sortmatches,false);
  gt_option_parser_add_option(op, op_sortmatches);

  /* -verify-alignment */
  op_verify_alignment = gt_option_new_bool("verify-alignment",
                                           "verify correctned of alignment",
                                           &arguments->verify_alignment,false);
  gt_option_parser_add_option(op, op_verify_alignment);
  gt_option_is_development_option(op_verify_alignment);

  /* -affine */
  op_affine_alignment = gt_option_new_bool("affine",
                                           "compute optimal affine alignment "
                                           "for substrings in given "
                                           "coordinates",
                                           &arguments->affine_alignment,false);
  gt_option_parser_add_option(op, op_affine_alignment);
  gt_option_is_development_option(op_affine_alignment);

  /* -f */
  option_filename = gt_option_new_filename("f",
                                          "path to file with match coordinates",
                                          arguments->matchfilename);
  gt_option_is_mandatory(option_filename);
  gt_option_parser_add_option(op, option_filename);

  return op;
}

static int gt_show_seedext_arguments_check(GT_UNUSED int rest_argc,
                                           void *tool_arguments,
                                           GtError *err)
{
  GtShowSeedextArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments != NULL);
  if (arguments->matchfilename == NULL ||
      gt_str_length(arguments->matchfilename) == 0)
  {
    gt_error_set(err,"option -f requires a file name");
    had_err = -1;
  }
  return had_err;
}

static int gt_show_seedext_runner(GT_UNUSED int argc,
                                  GT_UNUSED const char **argv,
                                  GT_UNUSED int parsed_args,
                                  void *tool_arguments,
                                  GtError *err)
{
  int had_err = 0;
  GtShowSeedextArguments *arguments = tool_arguments;
  GtSeedextendMatchIterator *semi = NULL;
  const GtEncseq *aencseq = NULL, *bencseq = NULL;
  GtSeedExtendDisplayFlag *out_display_flag = NULL;
  GtGreedyextendmatchinfo *greedyextendmatchinfo = NULL;
  const GtExtendCharAccess a_extend_char_access = GT_EXTEND_CHAR_ACCESS_ANY;
  const GtExtendCharAccess b_extend_char_access = GT_EXTEND_CHAR_ACCESS_ANY;
  GtTimer *timer = NULL;
  const GtSeedExtendDisplaySetMode setmode
    = GT_SEED_EXTEND_DISPLAY_SET_STANDARD;

  if (gt_showtime_enabled())
  {
    timer = gt_timer_new();
    gt_timer_start(timer);
  }
  gt_error_check(err);
  gt_assert(arguments != NULL);
  out_display_flag
    = gt_querymatch_display_flag_new(arguments->display_args,setmode,err);
  if (out_display_flag == NULL)
  {
    had_err = true;
  }
  if (!had_err)
  {
    semi = gt_seedextend_match_iterator_new(arguments->matchfilename,err);
    if (semi == NULL)
    {
      had_err = -1;
    }
  }
  if (!had_err)
  {
    if (gt_querymatch_subjectid_display(out_display_flag) &&
        !gt_seedextend_match_iterator_has_subjectid(semi))
    {
      gt_error_set(err,"output of subject id requires subject id in input");
      had_err = -1;
    }
  }
  if (!had_err)
  {
    if (gt_querymatch_queryid_display(out_display_flag) &&
        !gt_seedextend_match_iterator_has_queryid(semi))
    {
      gt_error_set(err,"output of query id requires query id in input");
      had_err = -1;
    }
  }
  if (!had_err)
  {
    printf("%s\n",gt_seedextend_match_iterator_Options_line(semi));
    aencseq = gt_seedextend_match_iterator_aencseq(semi);
    bencseq = gt_seedextend_match_iterator_bencseq(semi);
    /* the following are used if seed_extend is set */
    gt_querymatch_Fields_output(stdout,out_display_flag);
    if (arguments->verify_alignment)
    {
      gt_seedextend_match_iterator_verify_alignment_set(semi);
    }
    if (arguments->verify_alignment ||
        gt_querymatch_run_aligner(out_display_flag))
    {
      if (gt_seedextend_match_iterator_querymatchoutoptions_set(
                            semi,
                            !arguments->relax_polish,
                            a_extend_char_access,
                            b_extend_char_access,
                            out_display_flag,
                            err) != 0)
      {
        had_err = -1;
      }
    } else
    {
      if (arguments->affine_alignment)
      {
        gt_error_set(err,"Option -affine requires option -outfmt without any "
                         " of the following arguments: alignment, trace, "
                         "dtrace, cigar, or cigarX");
        had_err = -1;
      }
    }
  }
  if (!gt_seedextend_match_iterator_has_seqnums(semi))
  {
    gt_error_set(err,"match file without sequence numbers cannot be processed");
    had_err = -1;
  }
  if (!had_err)
  {
    GtKarlinAltschulStat *karlin_altschul_stat = NULL;
    const bool match_has_cigar = gt_seedextend_match_iterator_has_cigar(semi),
               match_has_seed = gt_seedextend_match_iterator_has_seed(semi),
               dtrace = gt_seedextend_match_iterator_dtrace(semi);
    const GtUword trace_delta = gt_seedextend_match_iterator_trace_delta(semi);
    GtSequencepairbuffer *seqpairbuf = NULL;
    GtAffineDPreservoir *adpr = NULL;

    if (gt_querymatch_evalue_display(out_display_flag) ||
        gt_querymatch_bitscore_display(out_display_flag))
    {
      karlin_altschul_stat = gt_karlin_altschul_stat_new_gapped(
                                       gt_encseq_total_length(aencseq),
                                       gt_encseq_num_of_sequences(aencseq),
                                       bencseq);
    }
    if (arguments->sortmatches)
    {
      (void) gt_seedextend_match_iterator_all_sorted(semi,true);
    }
    if (arguments->affine_alignment)
    {
      const bool opt_memory = false, keepcolumns = true;

      seqpairbuf = gt_calloc(1,sizeof *seqpairbuf);
      adpr = gt_affine_diagonalband_new(opt_memory, keepcolumns, 0, 0);
    }
    while (true)
    {
      GtQuerymatch *querymatchptr = gt_seedextend_match_iterator_next(semi);
      const double evalue = gt_seedextend_match_iterator_evalue(semi),
                   bitscore = gt_seedextend_match_iterator_bitscore(semi);

      if (querymatchptr == NULL)
      {
        break;
      }
      gt_querymatch_recompute_alignment(querymatchptr,
                                        out_display_flag,
                                        match_has_cigar,
                                        dtrace,
                                        trace_delta,
                                        match_has_seed,
                                        seqpairbuf,
                                        adpr,
                                        aencseq,
                                        bencseq,
                                        karlin_altschul_stat,
                                        evalue,
                                        bitscore);
    }
    gt_greedy_extend_matchinfo_delete(greedyextendmatchinfo);
    if (seqpairbuf != NULL)
    {
      gt_free(seqpairbuf->a_sequence);
      gt_free(seqpairbuf->b_sequence);
      gt_free(seqpairbuf);
      gt_affine_diagonalband_delete(adpr);
    }
    gt_karlin_altschul_stat_delete(karlin_altschul_stat);
  }
  gt_seedextend_match_iterator_delete(semi);
  if (!had_err && gt_showtime_enabled())
  {
    printf("# TIME show_seedext %s alignment ",
           gt_querymatch_alignment_display(out_display_flag)
           ? "with" : "without");
    gt_timer_show_formatted(timer,GT_WD ".%06ld\n",stdout);
  }
  gt_querymatch_display_flag_delete(out_display_flag);
  if (gt_showtime_enabled())
  {
    gt_timer_delete(timer);
  }
  return had_err;
}

GtTool* gt_show_seedext(void)
{
  return gt_tool_new(gt_show_seedext_arguments_new,
                     gt_show_seedext_arguments_delete,
                     gt_show_seedext_option_parser_new,
                     gt_show_seedext_arguments_check,
                     gt_show_seedext_runner);
}
