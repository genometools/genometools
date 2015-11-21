/*
  Copyright (c) 2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "core/ma_api.h"
#include "core/minmax.h"
#include "core/str_api.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "core/encseq.h"
#include "match/revcompl.h"
#include "match/ft-polish.h"
#include "match/seed-extend.h"
#include "match/seq_or_encseq.h"
#include "match/seed-extend-iter.h"
#include "extended/linearalign.h"
#include "tools/gt_show_seedext.h"

typedef struct
{
  bool show_alignment,
       seed_display,
       relax_polish,
       seed_extend;
  GtStr *matchfilename;
} GtShowSeedextArguments;

static void* gt_show_seedext_arguments_new(void)
{
  GtShowSeedextArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->matchfilename = gt_str_new();
  return arguments;
}

static void gt_show_seedext_arguments_delete(void *tool_arguments)
{
  GtShowSeedextArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->matchfilename);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_show_seedext_option_parser_new(void *tool_arguments)
{
  GtShowSeedextArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *op_ali, *option_filename, *op_relax_polish,
           *op_seed_extend;

  gt_assert(arguments);
  /* init */
  op = gt_option_parser_new("[options] -f <matchfilename>",
                            "Parse output of a seed extension and show/verify "
                            "the alignment.");

  /* -a */
  op_ali = gt_option_new_bool("a",
                              "show alignment",
                              &arguments->show_alignment,
                              false);
  gt_option_parser_add_option(op, op_ali);

  /* -seed-display */
  option = gt_option_new_bool("seed-display",
                              "Display seeds in #-line and by "
                              "character + (instead of |) in middle "
                              "row of alignment column",
                              &arguments->seed_display,
                              false);
  gt_option_parser_add_option(op, option);

  /* -seed-extend */
  op_seed_extend = gt_option_new_bool("seed-extend",
                              "read the seeds from the # seed: -lines and "
                              "extend them; match lines are ignored; "
                              "match coordindates are displayed",
                              &arguments->seed_extend,
                              false);
  gt_option_parser_add_option(op, op_seed_extend);

  /* -relax-polish */
  op_relax_polish = gt_option_new_bool("relax-polish",
                                       "do not force alignments to have "
                                       "polished ends",
                                       &arguments->relax_polish,false);
  gt_option_parser_add_option(op, op_relax_polish);
  gt_option_is_development_option(op_relax_polish);
  gt_option_imply(op_relax_polish, op_ali);

  /* -f */
  option_filename = gt_option_new_filename("f",
                                          "path to file with match coordinates",
                                          arguments->matchfilename);
  gt_option_is_mandatory(option_filename);
  gt_option_parser_add_option(op, option_filename);

  gt_option_exclude(op_relax_polish, op_seed_extend);
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

static void gt_show_seed_extend_plain(LinspaceManagement *linspace_spacemanager,
                                      GtScoreHandler *linspace_scorehandler,
                                      GtAlignment *alignment,
                                      GtUchar *alignment_show_buffer,
                                      GtUword alignmentwidth,
                                      const GtUchar *characters,
                                      GtUchar wildcardshow,
                                      const GtEncseq *aencseq,
                                      GtUchar *asequence,
                                      const GtEncseq *bencseq,
                                      GtUchar *bsequence,
                                      const GtQuerymatch *querymatchptr)
{
  GtUword edist;
  GtReadmode query_readmode = gt_querymatch_query_readmode(querymatchptr);
  const GtUword distance = gt_querymatch_distance(querymatchptr),
                dblen = gt_querymatch_dblen(querymatchptr),
                queryseqnum = gt_querymatch_queryseqnum(querymatchptr),
                querystart = gt_querymatch_querystart(querymatchptr),
                querylen = gt_querymatch_querylen(querymatchptr);

  const GtUword apos_ab = gt_querymatch_dbstart(querymatchptr);
  const GtUword bpos_ab = gt_encseq_seqstartpos(bencseq, queryseqnum) +
                          querystart;

  gt_encseq_extract_encoded(aencseq, asequence, apos_ab, apos_ab + dblen - 1);
  gt_encseq_extract_encoded(bencseq, bsequence, bpos_ab,
                            bpos_ab + querylen - 1);
  if (query_readmode != GT_READMODE_FORWARD)
  {
    gt_assert(query_readmode == GT_READMODE_REVCOMPL);
    gt_inplace_reverse_complement(bsequence,querylen);
  }
  edist = gt_computelinearspace_generic(linspace_spacemanager,
                                        linspace_scorehandler,
                                        alignment,
                                        asequence,
                                        0,
                                        dblen,
                                        bsequence,
                                        0,
                                        querylen);
  if (edist < distance)
  {
    printf("# edist=" GT_WU " (smaller by " GT_WU ")\n",edist,distance - edist);
  }
  gt_assert(edist <= distance);
  gt_alignment_show_generic(alignment_show_buffer,
                            false,
                            alignment,
                            stdout,
                            alignmentwidth,
                            characters,
                            wildcardshow);
  gt_alignment_reset(alignment);
}

static void gt_show_seed_extend_encseq(GtQuerymatch *querymatchptr,
                                       const GtEncseq *aencseq,
                                       const GtEncseq *bencseq,
                                       GtUword query_totallength,
                                       GtUword seedlen,
                                       GtUword seedpos1,
                                       GtUword seedpos2)
{
  GtSeqorEncseq bseqorencseq;

  bseqorencseq.seq = NULL;
  bseqorencseq.encseq = bencseq;
  if (gt_querymatch_process(querymatchptr,
                            seedpos1,
                            seedpos2,
                            seedlen,
                            aencseq,
                            &bseqorencseq,
                            query_totallength,
                            false) != 0)
  {
    gt_querymatch_prettyprint(querymatchptr);
  }
}

static int gt_show_seedext_runner(GT_UNUSED int argc,
                                  GT_UNUSED const char **argv,
                                  GT_UNUSED int parsed_args,
                                  void *tool_arguments,
                                  GtError *err)
{
  int had_err = 0;
  GtUword alignmentwidth;
  GtShowSeedextArguments *arguments = tool_arguments;
  GtSeedextendMatchIterator *semi;

  gt_error_check(err);
  gt_assert(arguments != NULL);
  /* Parse option string in first line of file specified by filename. */
  alignmentwidth = arguments->show_alignment ? 70 : 0;
  semi = gt_seedextend_match_iterator_new(arguments->matchfilename,err);
  if (semi == NULL)
  {
    had_err = -1;
  }
  /* Parse seed extensions. */
  if (!had_err)
  {
    const GtEncseq *aencseq = gt_seedextend_match_iterator_aencseq(semi),
                   *bencseq = gt_seedextend_match_iterator_bencseq(semi);
    GtGreedyextendmatchinfo *greedyextendmatchinfo = NULL;
    Polishing_info *pol_info = NULL;
    GtQuerymatchoutoptions *querymatchoutoptions = NULL;
    const GtUchar *characters = gt_encseq_alphabetcharacters(aencseq);
    const GtUchar wildcardshow = gt_encseq_alphabetwildcardshow(aencseq);
    GtUchar *alignment_show_buffer = gt_alignment_buffer_new(alignmentwidth);
    LinspaceManagement *linspace_spacemanager = gt_linspaceManagement_new();
    GtScoreHandler *linspace_scorehandler = gt_scorehandler_new(0,1,0,1);;
    GtAlignment *alignment = gt_alignment_new();
    const GtUword query_totallength = gt_encseq_total_length(bencseq);
    GtUchar *a_sequence = NULL, *b_sequence = NULL;
    GtUword a_allocated = 0, b_allocated = 0;
    double matchscore_bias = GT_DEFAULT_MATCHSCORE_BIAS;
    GtProcessinfo_and_querymatchspaceptr processinfo_and_querymatchspaceptr;

    if (!arguments->relax_polish)
    {
      if (gt_seedextend_match_iterator_bias_parameters(semi))
      {
        matchscore_bias = gt_greedy_dna_sequence_bias_get(aencseq);
      }
      pol_info = polishing_info_new_with_bias(
                          gt_seedextend_match_iterator_errorpercentage(semi),
                          matchscore_bias,
                          gt_seedextend_match_iterator_history_size(semi));
    }
    if (arguments->seed_display)
    {
      gt_seedextend_match_iterator_seed_display_set(semi);
    }
    querymatchoutoptions = gt_querymatchoutoptions_new(alignmentwidth);
    gt_querymatchoutoptions_for_align_only(querymatchoutoptions,
                          gt_seedextend_match_iterator_errorpercentage(semi),
                          matchscore_bias,
                          gt_seedextend_match_iterator_history_size(semi),
                          !arguments->relax_polish,
                          arguments->seed_display);
    gt_seedextend_match_iterator_outoptions_set(semi,querymatchoutoptions);
    if (arguments->seed_extend)
    {
      greedyextendmatchinfo
        = gt_greedy_extend_matchinfo_new(70,
                              GT_MAX_ALI_LEN_DIFF,
                              gt_seedextend_match_iterator_history_size(semi),
                              GT_MIN_PERC_MAT_HISTORY,
                              0, /* userdefinedleastlength */
                              GT_EXTEND_CHAR_ACCESS_ANY,
                              100,
                              pol_info);
    }
    if (pol_info != NULL)
    {
      gt_alignment_polished_ends(alignment,pol_info,false);
    }
    processinfo_and_querymatchspaceptr.processinfo = greedyextendmatchinfo;
    while (true)
    {
      GtQuerymatch *querymatchptr = gt_seedextend_match_iterator_next(semi);
      if (querymatchptr == NULL)
      {
        break;
      }
      if (gt_seedextend_match_iterator_has_seedline(semi))
      {
        const GtUword seedlen = gt_seedextend_match_iterator_seedlen(semi),
                      seedpos1 = gt_seedextend_match_iterator_seedpos1(semi),
                      seedpos2 = gt_seedextend_match_iterator_seedpos2(semi);

        if (arguments->seed_extend)
        {
          if (aencseq == bencseq)
          {
            processinfo_and_querymatchspaceptr.querymatchspaceptr
              = querymatchptr;
            had_err = gt_greedy_extend_selfmatch_with_output(
                                  &processinfo_and_querymatchspaceptr,
                                  aencseq,
                                  seedlen,
                                  seedpos1,
                                  seedpos2,
                                  err);
            if (had_err)
            {
              break;
            }
          } else
          {
            gt_assert(false);
          }
        } else
        {
          gt_show_seed_extend_encseq(querymatchptr,
                                     aencseq,
                                     bencseq,
                                     query_totallength,
                                     seedlen,
                                     seedpos1,
                                     seedpos2);
        }
      } else
      {
        const GtUword dblen = gt_querymatch_dblen(querymatchptr),
                      querylen = gt_querymatch_querylen(querymatchptr);

        if (dblen >= a_allocated)
        {
          a_sequence = gt_realloc(a_sequence, sizeof *a_sequence * dblen);
          a_allocated = dblen;
        }
        if (querylen >= b_allocated)
        {
          b_sequence = gt_realloc(b_sequence,sizeof *b_sequence * querylen);
          b_allocated = querylen;
        }
        gt_show_seed_extend_plain(linspace_spacemanager,
                                  linspace_scorehandler,
                                  alignment,
                                  alignment_show_buffer,
                                  alignmentwidth,
                                  characters,
                                  wildcardshow,
                                  aencseq,
                                  a_sequence,
                                  bencseq,
                                  b_sequence,
                                  querymatchptr);
      }
    }
    polishing_info_delete(pol_info);
    gt_greedy_extend_matchinfo_delete(greedyextendmatchinfo);
    gt_querymatchoutoptions_delete(querymatchoutoptions);
    gt_free(alignment_show_buffer);
    gt_scorehandler_delete(linspace_scorehandler);
    gt_linspaceManagement_delete(linspace_spacemanager);
    gt_free(a_sequence);
    gt_free(b_sequence);
    gt_alignment_delete(alignment);
  }
  gt_seedextend_match_iterator_delete(semi);
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
