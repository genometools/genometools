/*
  Copyright (c) 2007-2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2015 Center for Bioinformatics, University of Hamburg

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

#include "core/ma_api.h"
#include "core/types_api.h"
#include "core/minmax.h"
#include "revcompl.h"
#include "seed-extend.h"
#include "ft-polish.h"
#include "querymatch-align.h"

struct GtQuerymatchoutoptions
{
  GtUword useqbuffer_size,
          vseqbuffer_size;
  GtUchar *useqbuffer, *vseqbuffer;
  GtEoplist *eoplist;
  GtEoplistReader *eoplist_reader, *eoplist_reader_verify;
  GtFronttrace *front_trace;
  GtGreedyextendmatchinfo *ggemi;
  const GtUchar *characters;
  GtEncseqReader *esr_for_align_show;
  GtUchar wildcardshow;
  GtSeqpaircoordinates correction_info;
  bool always_polished_ends,
       generate_eoplist,
       show_eoplist;
  Polishing_info *pol_info;
};

GtQuerymatchoutoptions *gt_querymatchoutoptions_new(bool generate_eoplist,
                                                    bool show_eoplist,
                                                    GtUword alignmentwidth)
{
  GtQuerymatchoutoptions *querymatchoutoptions
    = gt_malloc(sizeof *querymatchoutoptions);

  querymatchoutoptions->generate_eoplist
    = generate_eoplist || alignmentwidth > 0 || show_eoplist;
  querymatchoutoptions->show_eoplist = show_eoplist;
  querymatchoutoptions->front_trace = NULL;
  querymatchoutoptions->ggemi = NULL;
  querymatchoutoptions->useqbuffer = NULL;
  querymatchoutoptions->useqbuffer_size = 0;
  querymatchoutoptions->vseqbuffer = NULL;
  querymatchoutoptions->vseqbuffer_size = 0;
  querymatchoutoptions->eoplist = gt_eoplist_new();
  querymatchoutoptions->eoplist_reader_verify = NULL;
  if (alignmentwidth > 0)
  {
    querymatchoutoptions->eoplist_reader
      = gt_eoplist_reader_new(querymatchoutoptions->eoplist);
    gt_eoplist_reader_reset_width(querymatchoutoptions->eoplist_reader,
                                  alignmentwidth);
  } else
  {
    querymatchoutoptions->eoplist_reader = NULL;
  }
  querymatchoutoptions->esr_for_align_show = NULL;
  querymatchoutoptions->characters = NULL;
  querymatchoutoptions->always_polished_ends = true;
  return querymatchoutoptions;
}

void gt_querymatchoutoptions_extend(
                  GtQuerymatchoutoptions *querymatchoutoptions,
                  GtUword errorpercentage,
                  GtUword maxalignedlendifference,
                  GtUword history_size,
                  GtUword perc_mat_history,
                  GtExtendCharAccess extend_char_access,
                  bool weakends,
                  GtUword sensitivity,
                  double matchscore_bias,
                  bool always_polished_ends,
                  unsigned int display_flag)
{
  if (errorpercentage > 0)
  {
    gt_assert(querymatchoutoptions != NULL);
    if (querymatchoutoptions->generate_eoplist)
    {
      querymatchoutoptions->front_trace = front_trace_new();
    }
    querymatchoutoptions->pol_info
      = polishing_info_new_with_bias(weakends ? MAX(errorpercentage,20)
                                              : errorpercentage,
                                     matchscore_bias,
                                     history_size);
    querymatchoutoptions->ggemi
      = gt_greedy_extend_matchinfo_new(errorpercentage,
                                       maxalignedlendifference,
                                       history_size, /* default value */
                                       perc_mat_history,
                                       0,/* userdefinedleastlength not used */
                                       extend_char_access,
                                       sensitivity,
                                       querymatchoutoptions->pol_info);
    if (querymatchoutoptions->eoplist != NULL)
    {
      if (always_polished_ends)
      {
        gt_eoplist_polished_ends(querymatchoutoptions->eoplist,
                                 querymatchoutoptions->pol_info,true);
      }
      if (gt_querymatch_seed_display(display_flag))
      {
        gt_eoplist_seed_display_set(querymatchoutoptions->eoplist);
      }
    }
    querymatchoutoptions->always_polished_ends = always_polished_ends;
  }
}

void gt_querymatchoutoptions_for_align_only(
                  GtQuerymatchoutoptions *querymatchoutoptions,
                  GtUword errorpercentage,
                  double matchscore_bias,
                  GtUword history_size,
                  bool always_polished_ends,
                  unsigned int display_flag)
{
  gt_querymatchoutoptions_extend(querymatchoutoptions,
                                 errorpercentage,
                                 GT_MAX_ALI_LEN_DIFF,
                                 history_size,
                                 GT_MIN_PERC_MAT_HISTORY,
                                 GT_EXTEND_CHAR_ACCESS_ANY,
                                 false,
                                 100,
                                 matchscore_bias,
                                 always_polished_ends,
                                 display_flag);
}

void gt_querymatchoutoptions_delete(
        GtQuerymatchoutoptions *querymatchoutoptions)
{
  if (querymatchoutoptions != NULL)
  {
    front_trace_delete(querymatchoutoptions->front_trace);
    gt_greedy_extend_matchinfo_delete(querymatchoutoptions->ggemi);
    gt_free(querymatchoutoptions->useqbuffer);
    if (querymatchoutoptions->vseqbuffer_size > 0)
    {
      gt_free(querymatchoutoptions->vseqbuffer);
    }
    gt_eoplist_delete(querymatchoutoptions->eoplist);
    gt_eoplist_reader_delete(querymatchoutoptions->eoplist_reader);
    gt_eoplist_reader_delete(querymatchoutoptions->eoplist_reader_verify);
    gt_encseq_reader_delete(querymatchoutoptions->esr_for_align_show);
    polishing_info_delete(querymatchoutoptions->pol_info);
    gt_free(querymatchoutoptions);
  }
}

static void seededmatch2eoplist(GtQuerymatchoutoptions *querymatchoutoptions,
                                const GtEncseq *encseq,
                         const GtSeqorEncseq *query,
                         GtReadmode query_readmode,
                         GtUword query_seqstartpos,
                         GtUword query_totallength,
                         GtUword dbstart,
                         GtUword dblen,
                         GtUword abs_querystart,
                         GtUword querylen,
                         GtUword seedpos1,
                         GtUword seedpos2,
                         GtUword seedlen,
                         bool verify_alignment,
                         bool greedyextension)
{
  GtUword ulen, vlen, ustart, vstart;
  Polished_point right_best_polished_point = {0,0,0};
  Polished_point left_best_polished_point = {0,0,0};
  GtUword pol_size;
  GtSeqpaircoordinates *coords;
  GtUword leftcolumn, rightcolumn;

  gt_assert(querymatchoutoptions != NULL &&
            querymatchoutoptions->pol_info != NULL);
  pol_size = GT_MULT2(querymatchoutoptions->pol_info->cut_depth);
  gt_eoplist_reset(querymatchoutoptions->eoplist);
  ustart = seedpos1 + seedlen;
  vstart = seedpos2 + seedlen;
  gt_assert(dbstart + dblen >= ustart);
  ulen = dbstart + dblen - ustart;
  gt_assert(abs_querystart + querylen >= vstart);
  vlen = abs_querystart + querylen - vstart;
  if (ulen > 0 && vlen > 0)
  {
    gt_align_front_prune_edist(true,
                               &right_best_polished_point,
                               querymatchoutoptions->front_trace,
                               encseq,
                               query,
                               query_readmode,
                               query_seqstartpos,
                               query_totallength,
                               querymatchoutoptions->ggemi,
                               greedyextension,
                               seedlen,
                               ustart,
                               ulen,
                               vstart,
                               vlen);
    if (querymatchoutoptions->front_trace != NULL)
    {
      front_trace2eoplist(querymatchoutoptions->always_polished_ends,
                          querymatchoutoptions->eoplist,
                          querymatchoutoptions->front_trace,
                          &right_best_polished_point,
                          pol_size,
                          querymatchoutoptions->pol_info->match_score,
                          querymatchoutoptions->pol_info->difference_score,
                          NULL,
                          ulen,
                          NULL,
                          vlen);
      front_trace_reset(querymatchoutoptions->front_trace,ulen+vlen);
    }
  }
  if (querymatchoutoptions->generate_eoplist)
  {
    gt_eoplist_match_add(querymatchoutoptions->eoplist,seedlen);
  }
  if (seedpos1 > dbstart && seedpos2 > abs_querystart)
  {
    ulen = seedpos1 - dbstart;
    vlen = seedpos2 - abs_querystart;
    gt_align_front_prune_edist(false,
                               &left_best_polished_point,
                               querymatchoutoptions->front_trace,
                               encseq,
                               query,
                               query_readmode,
                               query_seqstartpos,
                               query_totallength,
                               querymatchoutoptions->ggemi,
                               greedyextension,
                               seedlen,
                               dbstart,
                               ulen,
                               abs_querystart,
                               vlen);
    if (querymatchoutoptions->front_trace != NULL)
    {
      GtUword previous_eoplistlen
        = gt_eoplist_length(querymatchoutoptions->eoplist);
      front_trace2eoplist(querymatchoutoptions->always_polished_ends,
                          querymatchoutoptions->eoplist,
                          querymatchoutoptions->front_trace,
                          &left_best_polished_point,
                          pol_size,
                          querymatchoutoptions->pol_info->match_score,
                          querymatchoutoptions->pol_info->difference_score,
                          NULL,
                          ulen,
                          NULL,
                          vlen);
      gt_eoplist_reverse_end(querymatchoutoptions->eoplist,previous_eoplistlen);
      front_trace_reset(querymatchoutoptions->front_trace,ulen+vlen);
    }
  }
  coords = &querymatchoutoptions->correction_info;
  gt_assert(dbstart <= seedpos1 - left_best_polished_point.row);
  coords->uoffset = seedpos1 - left_best_polished_point.row - dbstart;
  coords->ulen = seedlen + left_best_polished_point.row +
                 right_best_polished_point.row;
  leftcolumn = left_best_polished_point.alignedlen -
               left_best_polished_point.row;
  rightcolumn = right_best_polished_point.alignedlen -
                right_best_polished_point.row;
  gt_assert(seedpos2 >= leftcolumn &&
            abs_querystart <= seedpos2 - leftcolumn);
  coords->voffset = seedpos2 - leftcolumn - abs_querystart;
  coords->vlen = seedlen + leftcolumn + rightcolumn;
  coords->sumdist = left_best_polished_point.distance +
                    right_best_polished_point.distance;
  coords->sum_max_mismatches = left_best_polished_point.max_mismatches +
                               right_best_polished_point.max_mismatches;
  gt_eoplist_reverse_end(querymatchoutoptions->eoplist,0);
  if (verify_alignment)
  {
    if (querymatchoutoptions->eoplist_reader_verify == NULL)
    {
      querymatchoutoptions->eoplist_reader_verify
        = gt_eoplist_reader_new(querymatchoutoptions->eoplist);
    }
    gt_eoplist_set_sequences(querymatchoutoptions->eoplist,NULL,
                             coords->ulen,NULL,coords->vlen);
    gt_eoplist_verify(querymatchoutoptions->eoplist,
                      querymatchoutoptions->eoplist_reader_verify,
                      coords->sumdist,
                      true);
  }
}

void gt_frontprune2eoplist(GtQuerymatchoutoptions *querymatchoutoptions,
                           const GtEncseq *encseq,
                           const GtSeqorEncseq *query,
                           GtReadmode query_readmode,
                           GtUword query_seqstartpos,
                           GtUword query_totallength,
                           GtUword dbstart,
                           GtUword abs_querystart,
                           GtUword ustart,
                           GtUword ulen,
                           GtUword vstart,
                           GtUword vlen)
{
  Polished_point right_best_polished_point = {0,0,0};
  GtUword pol_size;
  GtSeqpaircoordinates *coords;
  const bool greedyextension = true, rightextension = true;

  gt_assert(querymatchoutoptions != NULL &&
            querymatchoutoptions->pol_info != NULL);
  pol_size = GT_MULT2(querymatchoutoptions->pol_info->cut_depth);
  gt_eoplist_reset(querymatchoutoptions->eoplist);
  gt_assert(ulen > 0 && vlen > 0);
  gt_align_front_prune_edist(rightextension,
                             &right_best_polished_point,
                             querymatchoutoptions->front_trace,
                             encseq,
                             query,
                             query_readmode,
                             query_seqstartpos,
                             query_totallength,
                             querymatchoutoptions->ggemi,
                             greedyextension,
                             0,
                             ustart,
                             ulen,
                             vstart,
                             vlen);
  if (querymatchoutoptions->front_trace != NULL)
  {
    front_trace2eoplist(querymatchoutoptions->always_polished_ends,
                        querymatchoutoptions->eoplist,
                        querymatchoutoptions->front_trace,
                        &right_best_polished_point,
                        pol_size,
                        querymatchoutoptions->pol_info->match_score,
                        querymatchoutoptions->pol_info->difference_score,
                        NULL,
                        ulen,
                        NULL,
                        vlen);
    front_trace_reset(querymatchoutoptions->front_trace,ulen+vlen);
  }
  coords = &querymatchoutoptions->correction_info;
  coords->uoffset = ustart - dbstart;
  coords->ulen = right_best_polished_point.row;
  coords->voffset = vstart - abs_querystart;
  coords->vlen = right_best_polished_point.alignedlen -
                 right_best_polished_point.row;
  coords->sumdist = right_best_polished_point.distance;
}

bool gt_querymatchoutoptions_alignment_prepare(GtQuerymatchoutoptions
                                                *querymatchoutoptions,
                                               const GtEncseq *encseq,
                                               const GtSeqorEncseq *query,
                                               GtReadmode query_readmode,
                                               GtUword query_seqstartpos,
                                               GtUword query_totallength,
                                               GtUword dbstart,
                                               GtUword dblen,
                                               GtUword abs_querystart,
                                               GtUword abs_querystart_fwdstrand,
                                               GtUword querylen,
                                               GtUword edist,
                                               GtUword seedpos1,
                                               GtUword seedpos2,
                                               GtUword seedlen,
                                               bool verify_alignment,
                                               bool greedyextension)
{
  bool seededalignment = false;

  gt_assert(querymatchoutoptions != NULL);
  if (querymatchoutoptions->eoplist_reader != NULL &&
      querymatchoutoptions->characters == NULL)
  {
    querymatchoutoptions->characters
      = gt_encseq_alphabetcharacters(encseq);
    querymatchoutoptions->wildcardshow
      = gt_alphabet_wildcard_show(gt_encseq_alphabet(encseq));
  }
  if (querymatchoutoptions->esr_for_align_show == NULL)
  {
    querymatchoutoptions->esr_for_align_show
      = gt_encseq_create_reader_with_readmode(encseq,
                                              GT_READMODE_FORWARD,
                                              0);
  }
  if (dblen > querymatchoutoptions->useqbuffer_size)
  {
    querymatchoutoptions->useqbuffer
      = gt_realloc(querymatchoutoptions->useqbuffer,
                   sizeof *querymatchoutoptions->useqbuffer * dblen);
    querymatchoutoptions->useqbuffer_size = dblen;
  }
  gt_encseq_extract_encoded_with_reader(
                            querymatchoutoptions->esr_for_align_show,
                            encseq,
                            querymatchoutoptions->useqbuffer,
                            dbstart,
                            dbstart + dblen - 1);
  if ((query == NULL || query->seq == NULL ||
       query_readmode != GT_READMODE_FORWARD) &&
      querylen > querymatchoutoptions->vseqbuffer_size)
  {
    querymatchoutoptions->vseqbuffer
      = gt_realloc(querymatchoutoptions->vseqbuffer,
                   sizeof *querymatchoutoptions->vseqbuffer * querylen);
    querymatchoutoptions->vseqbuffer_size = querylen;
  }
  if (query == NULL || query->seq == NULL)
  {
    gt_assert(query == NULL || query->encseq != NULL);
    gt_encseq_extract_encoded_with_reader(
                            querymatchoutoptions->esr_for_align_show,
                            query == NULL ? encseq : query->encseq,
                            querymatchoutoptions->vseqbuffer,
                            abs_querystart_fwdstrand,
                            abs_querystart_fwdstrand + querylen - 1);
  } else
  {
    if (query_readmode == GT_READMODE_FORWARD)
    {
      querymatchoutoptions->vseqbuffer
        = (GtUchar *) (query->seq + abs_querystart_fwdstrand);
    } else
    {
      memcpy(querymatchoutoptions->vseqbuffer,
             query->seq + abs_querystart_fwdstrand,
             querylen * sizeof *querymatchoutoptions->vseqbuffer);
    }
  }
  if (query_readmode == GT_READMODE_REVERSE)
  {
    gt_inplace_reverse(querymatchoutoptions->vseqbuffer,querylen);
  } else
  {
    if (query_readmode == GT_READMODE_REVCOMPL)
    {
      gt_inplace_reverse_complement(querymatchoutoptions->vseqbuffer,
                                    querylen);
    } else
    {
      if (query_readmode == GT_READMODE_COMPL)
      {
        gt_inplace_complement(querymatchoutoptions->vseqbuffer,querylen);
      }
    }
  }
  if (edist > 0)
  {
    seededmatch2eoplist(querymatchoutoptions,
                        encseq,
                        query,
                        query_readmode,
                        query_seqstartpos,
                        query_totallength,
                        dbstart,
                        dblen,
                        abs_querystart,
                        querylen,
                        seedpos1,
                        seedpos2,
                        seedlen,
                        verify_alignment,
                        greedyextension);
    if (querymatchoutoptions->eoplist_reader != NULL)
    {
      gt_assert(dbstart <= seedpos1);
      gt_eoplist_set_seedoffset(querymatchoutoptions->eoplist,
                                seedpos1 - dbstart,
                                seedlen);
      gt_eoplist_set_sequences(querymatchoutoptions->eoplist,
                               querymatchoutoptions->useqbuffer +
                                 querymatchoutoptions->correction_info.uoffset,
                               querymatchoutoptions->correction_info.ulen,
                               querymatchoutoptions->vseqbuffer +
                                 querymatchoutoptions->correction_info.voffset,
                               querymatchoutoptions->correction_info.vlen);
    }
    seededalignment = true;
  } else
  {
    if (querymatchoutoptions->eoplist_reader != NULL)
    {
      gt_eoplist_set_sequences(querymatchoutoptions->eoplist,
                               querymatchoutoptions->useqbuffer,
                               dblen,
                               querymatchoutoptions->vseqbuffer,
                               querylen);
    }
  }
  return seededalignment;
}

void gt_querymatchoutoptions_alignment_show(const GtQuerymatchoutoptions
                                              *querymatchoutoptions,
                                            GtUword distance,
                                            bool verify_alignment,
                                            FILE *fp)
{
  if (querymatchoutoptions != NULL)
  {
    if (querymatchoutoptions->show_eoplist)
    {
      if (distance > 0)
      {
        gt_eoplist_show_plain(querymatchoutoptions->eoplist,fp);
      } else
      {
        fprintf(fp, "[]\n");
      }
    }
    if (querymatchoutoptions->eoplist_reader != NULL)
    {
      if (distance > 0)
      {
        if (verify_alignment)
        {
          gt_assert(querymatchoutoptions->eoplist_reader_verify != NULL);
          gt_eoplist_verify(querymatchoutoptions->eoplist,
                            querymatchoutoptions->eoplist_reader_verify,
                            distance,true);
        }
        gt_eoplist_format_generic(fp,
                                  querymatchoutoptions->eoplist,
                                  querymatchoutoptions->eoplist_reader,
                                  true,
                                  querymatchoutoptions->characters,
                                  querymatchoutoptions->wildcardshow);
      } else
      {
        gt_eoplist_format_exact(fp,
                                querymatchoutoptions->eoplist,
                                querymatchoutoptions->eoplist_reader,
                                querymatchoutoptions->characters);
      }
    }
    if (querymatchoutoptions->eoplist_reader != NULL ||
        querymatchoutoptions->show_eoplist)
    {
      gt_eoplist_reset(querymatchoutoptions->eoplist);
    }
  }
}

const GtSeqpaircoordinates *gt_querymatchoutoptions_correction_get(
              const GtQuerymatchoutoptions *querymatchoutoptions)
{
  gt_assert(querymatchoutoptions != NULL);
  return &querymatchoutoptions->correction_info;
}
