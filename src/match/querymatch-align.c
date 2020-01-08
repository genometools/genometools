/*
  Copyright (c) 2007-2016 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2016 Center for Bioinformatics, University of Hamburg

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

#include <float.h>
#include "core/ma_api.h"
#include "core/types_api.h"
#include "core/minmax_api.h"
#include "core/encseq_metadata.h"
#include "revcompl.h"
#include "seed-extend.h"
#include "ft-polish.h"
#include "querymatch-align.h"

struct GtQuerymatchoutoptions
{
  GtUword useqbuffer_size,
          vseqbuffer_size,
          trace_delta;
  GtUchar *useqbuffer, *vseqbuffer;
  GtEoplist *eoplist;
  GtEoplistReader *eoplist_reader, *eoplist_reader_verify;
  GtFrontTrace *front_trace;
  GtGreedyextendmatchinfo *ggemi;
  const GtUchar *characters;
  GtEncseqReader *db_esr_for_align_show,
                 *query_esr_for_align_show;
  GtEncseqMetadata *emd;
  GtUchar wildcardshow;
  GtSeqpaircoordinates correction_info;
  bool always_polished_ends;
  GtFtPolishing_info *pol_info;
};

GtQuerymatchoutoptions *gt_querymatchoutoptions_new(const
                                                     GtSeedExtendDisplayFlag
                                                      *out_display_flag,
                                                    const char *indexname,
                                                    GtError *err)
{
  GtQuerymatchoutoptions *querymatchoutoptions;

  if (indexname != NULL)
  {
    GtAlphabet *alphabet;
    GtEncseqMetadata *emd = gt_encseq_metadata_new(indexname,err);

    if (emd == NULL)
    {
      return NULL;
    }
    alphabet = gt_encseq_metadata_alphabet(emd);
    querymatchoutoptions = gt_malloc(sizeof *querymatchoutoptions);
    querymatchoutoptions->emd = emd;
    querymatchoutoptions->characters = gt_alphabet_characters(alphabet);
    querymatchoutoptions->wildcardshow = gt_alphabet_wildcard_show(alphabet);
  } else
  {
    querymatchoutoptions = gt_malloc(sizeof *querymatchoutoptions);
    querymatchoutoptions->emd = NULL;
    querymatchoutoptions->characters = NULL;
    querymatchoutoptions->wildcardshow = 0;
  }
  querymatchoutoptions->front_trace = NULL;
  querymatchoutoptions->ggemi = NULL;
  querymatchoutoptions->useqbuffer = NULL;
  querymatchoutoptions->useqbuffer_size = 0;
  querymatchoutoptions->vseqbuffer = NULL;
  querymatchoutoptions->vseqbuffer_size = 0;
  querymatchoutoptions->eoplist_reader_verify = NULL;
  querymatchoutoptions->pol_info = NULL;
  querymatchoutoptions->trace_delta
    = gt_querymatch_trace_delta_display(out_display_flag);
  querymatchoutoptions->eoplist_reader = gt_eoplist_reader_new();
  querymatchoutoptions->eoplist = gt_eoplist_new();
  if (gt_querymatch_alignment_display(out_display_flag))
  {
    gt_eoplist_reader_reset_width(querymatchoutoptions->eoplist_reader,
                                  gt_querymatch_display_alignmentwidth(
                                       out_display_flag));
  }
  querymatchoutoptions->db_esr_for_align_show = NULL;
  querymatchoutoptions->query_esr_for_align_show = NULL;
  querymatchoutoptions->always_polished_ends = true;
  return querymatchoutoptions;
}

void gt_querymatchoutoptions_reset(GtQuerymatchoutoptions *querymatchoutoptions)
{
  if (querymatchoutoptions != NULL)
  {
    if (querymatchoutoptions->vseqbuffer_size == 0)
    {
      querymatchoutoptions->vseqbuffer = NULL;
    }
    if (querymatchoutoptions->useqbuffer_size == 0)
    {
      querymatchoutoptions->useqbuffer = NULL;
    }
  }
}

void gt_querymatchoutoptions_extend(
                  GtQuerymatchoutoptions *querymatchoutoptions,
                  GtUword errorpercentage,
                  double evalue_threshold,
                  GtUword maxalignedlendifference,
                  GtUword history_size,
                  GtUword perc_mat_history,
                  GtExtendCharAccess a_extend_char_access,
                  GtExtendCharAccess b_extend_char_access,
                  bool cam_generic,
                  bool weakends,
                  GtUword sensitivity,
                  double matchscore_bias,
                  bool always_polished_ends,
                  const GtSeedExtendDisplayFlag *out_display_flag)
{
  if (errorpercentage > 0)
  {
    gt_assert(querymatchoutoptions != NULL);
    querymatchoutoptions->front_trace = front_trace_new();
    querymatchoutoptions->pol_info
      = polishing_info_new_with_bias(weakends ? GT_MAX(errorpercentage,20)
                                              : errorpercentage,
                                     matchscore_bias,
                                     history_size);
    querymatchoutoptions->ggemi
      = gt_greedy_extend_matchinfo_new(maxalignedlendifference,
                                       history_size, /* default value */
                                       perc_mat_history,
                                       0,/* userdefinedleastlength not used */
                                       errorpercentage,
                                       evalue_threshold,
                                       a_extend_char_access,
                                       b_extend_char_access,
                                       cam_generic,
                                       sensitivity,
                                       querymatchoutoptions->pol_info);
    if (always_polished_ends)
    {
      gt_eoplist_polished_ends(querymatchoutoptions->eoplist,
                               querymatchoutoptions->pol_info,true,
                               gt_querymatch_polinfo_display(out_display_flag));
    }
    if (gt_querymatch_seed_in_algn_display(out_display_flag))
    {
      gt_eoplist_display_seed_in_alignment_set(querymatchoutoptions->eoplist);
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
                  GtExtendCharAccess a_extend_char_access,
                  GtExtendCharAccess b_extend_char_access,
                  const GtSeedExtendDisplayFlag *out_display_flag)
{
  const bool weakends = false;
  const bool cam_generic = false;
  const GtUword sensitivity = 100;
  const double evalue_threshold = DBL_MAX;

  gt_querymatchoutoptions_extend(querymatchoutoptions,
                                 errorpercentage,
                                 evalue_threshold,
                                 GT_MAX_ALI_LEN_DIFF,
                                 history_size,
                                 GT_MIN_PERC_MAT_HISTORY,
                                 a_extend_char_access,
                                 b_extend_char_access,
                                 cam_generic,
                                 weakends,
                                 sensitivity,
                                 matchscore_bias,
                                 always_polished_ends,
                                 out_display_flag);
}

void gt_querymatchoutoptions_delete(
        GtQuerymatchoutoptions *querymatchoutoptions)
{
  if (querymatchoutoptions != NULL)
  {
    front_trace_delete(querymatchoutoptions->front_trace);
    gt_greedy_extend_matchinfo_delete(querymatchoutoptions->ggemi);
    if (querymatchoutoptions->useqbuffer_size > 0)
    {
      gt_free(querymatchoutoptions->useqbuffer);
    }
    if (querymatchoutoptions->vseqbuffer_size > 0)
    {
      gt_free(querymatchoutoptions->vseqbuffer);
    }
    gt_eoplist_delete(querymatchoutoptions->eoplist);
    gt_eoplist_reader_delete(querymatchoutoptions->eoplist_reader);
    gt_eoplist_reader_delete(querymatchoutoptions->eoplist_reader_verify);
    gt_encseq_reader_delete(querymatchoutoptions->db_esr_for_align_show);
    gt_encseq_reader_delete(querymatchoutoptions->query_esr_for_align_show);
    polishing_info_delete(querymatchoutoptions->pol_info);
    gt_encseq_metadata_delete(querymatchoutoptions->emd);
    gt_free(querymatchoutoptions);
  }
}

static void gt_querymtch_alignment_verification(
            GtQuerymatchoutoptions *querymatchoutoptions,
            GtUword dbstart,
            GtUword dblen,
            GtUword abs_querystart,
            GtUword querylen,
            GtUword sumdist)
{
  if (querymatchoutoptions->eoplist_reader_verify == NULL)
  {
    querymatchoutoptions->eoplist_reader_verify = gt_eoplist_reader_new();
  }
  gt_eoplist_set_sequences(querymatchoutoptions->eoplist,NULL,
                           dbstart,
                           dblen,
                           NULL,
                           abs_querystart,
                           querylen);
  gt_eoplist_verify(querymatchoutoptions->eoplist,
                    querymatchoutoptions->eoplist_reader_verify,
                    sumdist);
}

void gt_querymatchoutoptions_seededmatch2eoplist(
                                GtQuerymatchoutoptions *querymatchoutoptions,
                                const GtSeqorEncseq *dbes,
                                GtUword dbstart_relative,
                                GtUword db_seqstart,
                                GtUword dblen,
                                GtReadmode query_readmode,
                                const GtSeqorEncseq *queryes,
                                GtUword query_seqstart,
                                GtUword query_seqlen,
                                GtUword querystart_rel,
                                GtUword querylen,
                                GtUword db_seedpos_rel,
                                GtUword query_seedpos_rel,
                                GtUword seedlen,
                                bool verify_alignment,
                                bool greedyextension)
{
  GtUword ulen, vlen, ustart, vstart;
  GtFtPolished_point right_best_polished_point = {0,0,0,0,0},
                     left_best_polished_point = {0,0,0,0,0};
  GtUword pol_size;
  GtSeqpaircoordinates *coords;
  GtUword leftcolumn, rightcolumn;

  gt_assert(querymatchoutoptions != NULL &&
            querymatchoutoptions->pol_info != NULL);
  pol_size = GT_MULT2(querymatchoutoptions->pol_info->cut_depth);
  gt_eoplist_reset(querymatchoutoptions->eoplist);
  ustart = db_seedpos_rel + seedlen;
  vstart = query_seedpos_rel + seedlen;
  gt_assert(dbstart_relative + dblen >= ustart);
  ulen = dbstart_relative + dblen - ustart;
  gt_assert(querystart_rel + querylen >= vstart);
  vlen = querystart_rel + querylen - vstart;
  if (ulen > 0 && vlen > 0)
  {
    gt_align_front_prune_edist(true,
                               &right_best_polished_point,
                               querymatchoutoptions->front_trace,
                               dbes,
                               queryes,
                               query_readmode,
                               query_seqstart,
                               query_seqlen,
                               querymatchoutoptions->ggemi,
                               greedyextension,
                               seedlen,
                               db_seqstart + ustart,
                               ulen,
                               query_seqstart + vstart,
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
  gt_eoplist_match_add(querymatchoutoptions->eoplist,seedlen);
  if (db_seedpos_rel > dbstart_relative && query_seedpos_rel > querystart_rel)
  {
    ulen = db_seedpos_rel - dbstart_relative;
    vlen = query_seedpos_rel - querystart_rel;
    gt_align_front_prune_edist(false,
                               &left_best_polished_point,
                               querymatchoutoptions->front_trace,
                               dbes,
                               queryes,
                               query_readmode,
                               query_seqstart,
                               query_seqlen,
                               querymatchoutoptions->ggemi,
                               greedyextension,
                               seedlen,
                               db_seqstart + dbstart_relative,
                               ulen,
                               query_seqstart + querystart_rel,
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
  gt_assert(db_seedpos_rel >= dbstart_relative + left_best_polished_point.row);
  coords->uoffset = db_seedpos_rel - left_best_polished_point.row -
                                     dbstart_relative;
  coords->ulen = seedlen + left_best_polished_point.row +
                 right_best_polished_point.row;
  leftcolumn = left_best_polished_point.alignedlen -
               left_best_polished_point.row;
  rightcolumn = right_best_polished_point.alignedlen -
                right_best_polished_point.row;
  gt_assert(query_seedpos_rel >= leftcolumn + querystart_rel);
  coords->voffset = query_seedpos_rel - leftcolumn - querystart_rel;
  coords->vlen = seedlen + leftcolumn + rightcolumn;
  coords->sumdist = left_best_polished_point.distance +
                    right_best_polished_point.distance;
  coords->sum_max_mismatches = left_best_polished_point.max_mismatches +
                               right_best_polished_point.max_mismatches;
  gt_eoplist_reverse_end(querymatchoutoptions->eoplist,0);
  if (verify_alignment)
  {
    gt_querymtch_alignment_verification(querymatchoutoptions,
                                        db_seqstart + dbstart_relative,
                                        coords->ulen,
                                        query_seqstart + querystart_rel,
                                        coords->vlen,
                                        coords->sumdist);
  }
  gt_eoplist_set_seedoffset(querymatchoutoptions->eoplist,
                            db_seedpos_rel - dbstart_relative,
                            seedlen);
}

void gt_frontprune2eoplist(GtQuerymatchoutoptions *querymatchoutoptions,
                           const GtSeqorEncseq *dbes,
                           GtUword dbstart,
                           GtUword dblen,
                           const GtSeqorEncseq *queryes,
                           GtReadmode query_readmode,
                           GtUword query_seqstart,
                           GtUword query_seqlen,
                           GtUword querystart,
                           GtUword querylen,
                           bool verify_alignment)
{
  GtFtPolished_point right_best_polished_point = {0,0,0,0,0};
  GtUword pol_size;
  GtSeqpaircoordinates *coords;
  const bool greedyextension = true, rightextension = true;

  gt_assert(querymatchoutoptions != NULL &&
            querymatchoutoptions->pol_info != NULL);
  pol_size = GT_MULT2(querymatchoutoptions->pol_info->cut_depth);
  gt_eoplist_reset(querymatchoutoptions->eoplist);
  gt_assert(dblen > 0 && querylen > 0);
  gt_align_front_prune_edist(rightextension,
                             &right_best_polished_point,
                             querymatchoutoptions->front_trace,
                             dbes,
                             queryes,
                             query_readmode,
                             query_seqstart,
                             query_seqlen,
                             querymatchoutoptions->ggemi,
                             greedyextension,
                             0,
                             dbstart,
                             dblen,
                             query_seqstart + querystart,
                             querylen);
  gt_assert(querymatchoutoptions->front_trace != NULL);
  front_trace2eoplist(querymatchoutoptions->always_polished_ends,
                      querymatchoutoptions->eoplist,
                      querymatchoutoptions->front_trace,
                      &right_best_polished_point,
                      pol_size,
                      querymatchoutoptions->pol_info->match_score,
                      querymatchoutoptions->pol_info->difference_score,
                      NULL,
                      dblen,
                      NULL,
                      querylen);
  gt_eoplist_reverse_end(querymatchoutoptions->eoplist,0);
  front_trace_reset(querymatchoutoptions->front_trace,dblen+querylen);
  coords = &querymatchoutoptions->correction_info;
  coords->uoffset = 0;
  coords->ulen = right_best_polished_point.row;
  coords->voffset = 0;
  coords->vlen = right_best_polished_point.alignedlen -
                 right_best_polished_point.row;
  coords->sumdist = right_best_polished_point.distance;
  coords->sum_max_mismatches = right_best_polished_point.max_mismatches;
  if (verify_alignment)
  {
    gt_querymtch_alignment_verification(querymatchoutoptions,
                                        dbstart,
                                        coords->ulen,
                                        query_seqstart + querystart,
                                        coords->vlen,
                                        coords->sumdist);
  }
}

static void gt_querymatchoutoptions_set_sequences(GtQuerymatchoutoptions
                                             *querymatchoutoptions,
                                           GtUword dbstart_relative,
                                           GtUword dblen,
                                           GtUword querystart,
                                           GtUword querylen,
                                           bool withcorrection)
{
  gt_assert(querymatchoutoptions != NULL);

  if (withcorrection)
  {
    gt_eoplist_set_sequences(querymatchoutoptions->eoplist,
                             querymatchoutoptions->useqbuffer +
                               querymatchoutoptions->correction_info.uoffset,
                             dbstart_relative +
                               querymatchoutoptions->correction_info.uoffset,
                             querymatchoutoptions->correction_info.ulen,
                             querymatchoutoptions->vseqbuffer +
                               querymatchoutoptions->correction_info.voffset,
                             querystart +
                               querymatchoutoptions->correction_info.voffset,
                             querymatchoutoptions->correction_info.vlen);
  } else
  {
    gt_eoplist_set_sequences(querymatchoutoptions->eoplist,
                             querymatchoutoptions->useqbuffer,
                             dbstart_relative,
                             dblen,
                             querymatchoutoptions->vseqbuffer,
                             querystart,
                             querylen);
  }
}

void gt_querymatchoutoptions_extract_seq(GtQuerymatchoutoptions
                                           *querymatchoutoptions,
                                         const GtSeqorEncseq *dbes,
                                         GtUword dbstart_relative,
                                         GtUword dbstart,
                                         GtUword dblen,
                                         GtReadmode query_readmode,
                                         const GtSeqorEncseq *queryes,
                                         GtUword querystart,
                                         GtUword abs_querystart_fwdstrand,
                                         GtUword querylen,
                                         bool withcorrection)
{
  gt_assert(querymatchoutoptions != NULL);
  if (querymatchoutoptions->characters == NULL)
  {
    if (dbes->encseq != NULL)
    {
      querymatchoutoptions->characters
        = gt_encseq_alphabetcharacters(dbes->encseq);
      querymatchoutoptions->wildcardshow
        = gt_alphabet_wildcard_show(gt_encseq_alphabet(dbes->encseq));
    } else
    {
      querymatchoutoptions->characters = dbes->characters;
      querymatchoutoptions->wildcardshow = dbes->wildcardshow;
    }
  }
  if (dbes->encseq != NULL)
  {
    if (querymatchoutoptions->db_esr_for_align_show == NULL)
    {
      querymatchoutoptions->db_esr_for_align_show
        = gt_encseq_create_reader_with_readmode(dbes->encseq,
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
                            querymatchoutoptions->db_esr_for_align_show,
                            dbes->encseq,
                            querymatchoutoptions->useqbuffer,
                            dbstart,
                            dbstart + dblen - 1);
  } else
  {
    querymatchoutoptions->useqbuffer = (GtUchar *) (dbes->seq + dbstart);
  }
  if ((queryes->encseq != NULL || query_readmode != GT_READMODE_FORWARD) &&
      querylen > querymatchoutoptions->vseqbuffer_size)
  {
    querymatchoutoptions->vseqbuffer
      = gt_realloc(querymatchoutoptions->vseqbuffer,
                   sizeof *querymatchoutoptions->vseqbuffer * querylen);
    querymatchoutoptions->vseqbuffer_size = querylen;
  }
  if (queryes->encseq != NULL)
  {
    if (querymatchoutoptions->query_esr_for_align_show == NULL)
    {
      querymatchoutoptions->query_esr_for_align_show
        = gt_encseq_create_reader_with_readmode(queryes->encseq,
                                                GT_READMODE_FORWARD,
                                                0);
    }
    gt_encseq_extract_encoded_with_reader(
                            querymatchoutoptions->query_esr_for_align_show,
                            queryes->encseq,
                            querymatchoutoptions->vseqbuffer,
                            abs_querystart_fwdstrand,
                            abs_querystart_fwdstrand + querylen - 1);
  } else
  {
    if (query_readmode != GT_READMODE_FORWARD)
    {
      memcpy(querymatchoutoptions->vseqbuffer,
             queryes->seq + abs_querystart_fwdstrand,
             querylen * sizeof *querymatchoutoptions->vseqbuffer);
    } else
    {
      querymatchoutoptions->vseqbuffer
        = (GtUchar *) (queryes->seq + abs_querystart_fwdstrand);
    }
  }
  if (query_readmode == GT_READMODE_REVERSE)
  {
    gt_inplace_reverse(querymatchoutoptions->vseqbuffer,querylen);
  } else
  {
    if (query_readmode == GT_READMODE_REVCOMPL)
    {
      gt_inplace_reverse_complement(querymatchoutoptions->vseqbuffer,querylen);
    } else
    {
      gt_assert(query_readmode == GT_READMODE_FORWARD);
    }
  }
  gt_querymatchoutoptions_set_sequences(querymatchoutoptions,
                                        dbstart_relative,
                                        dblen,
                                        querystart,
                                        querylen,
                                        withcorrection);
}

void gt_querymatchoutoptions_cigar_show(const GtQuerymatchoutoptions
                                              *querymatchoutoptions,
                                        bool distinguish_mismatch_match,
                                        FILE *fp)
{
  gt_assert(querymatchoutoptions != NULL &&
            querymatchoutoptions->eoplist != NULL);
  gt_eoplist_reader_reset(querymatchoutoptions->eoplist_reader,
                          querymatchoutoptions->eoplist,true);
  gt_eoplist_show_cigar(querymatchoutoptions->eoplist_reader,
                        distinguish_mismatch_match,fp);
}

void gt_querymatchoutoptions_trace_show(const GtQuerymatchoutoptions
                                              *querymatchoutoptions,
                                        bool dtrace,
                                        FILE *fp)
{
  GtEoplistSegment segment;
  bool first = true;

  gt_assert(querymatchoutoptions != NULL);
  gt_eoplist_reader_reset(querymatchoutoptions->eoplist_reader,
                          querymatchoutoptions->eoplist,true);
  while (gt_eoplist_reader_next_segment(&segment,
                                        querymatchoutoptions->eoplist_reader,
                                        querymatchoutoptions->trace_delta))
  {
    if (!first)
    {
      fputc(',',fp);
    } else
    {
      first = false;
    }
    fprintf(fp,"%d",dtrace ? ((int) querymatchoutoptions->trace_delta -
                              (int) segment.aligned_v)
                           : (int) segment.aligned_v);
  }
}

void gt_querymatchoutoptions_alignment_show(const GtQuerymatchoutoptions
                                              *querymatchoutoptions,
                                            GtUword subject_seqlength,
                                            GtUword query_reference,
                                            GtUword one_off,
                                            GtUword distance,
                                            bool distinguish_mismatch_match,
                                            bool verify_alignment,
                                            bool subject_first,
                                            bool alignment_show_forward,
                                            bool show_complement_characters,
                                            FILE *fp)
{
  gt_assert(querymatchoutoptions != NULL);
  if (distance > 0)
  {
    if (verify_alignment)
    {
      gt_assert(querymatchoutoptions->eoplist_reader_verify != NULL);
      gt_eoplist_verify(querymatchoutoptions->eoplist,
                        querymatchoutoptions->eoplist_reader_verify,
                        distance);
    }
    gt_eoplist_format_generic(fp,
                              querymatchoutoptions->eoplist,
                              querymatchoutoptions->eoplist_reader,
                              querymatchoutoptions->characters,
                              subject_seqlength,
                              query_reference,
                              one_off,
                              distinguish_mismatch_match,
                              subject_first,
                              alignment_show_forward,
                              show_complement_characters,
                              querymatchoutoptions->wildcardshow);
  } else
  {
    gt_eoplist_format_exact(fp,
                            querymatchoutoptions->eoplist,
                            querymatchoutoptions->eoplist_reader,
                            subject_seqlength,
                            query_reference,
                            one_off,
                            subject_first,
                            alignment_show_forward,
                            show_complement_characters,
                            querymatchoutoptions->characters);
  }
}

const GtSeqpaircoordinates *gt_querymatchoutoptions_correction_get(
              const GtQuerymatchoutoptions *querymatchoutoptions)
{
  gt_assert(querymatchoutoptions != NULL);
  return &querymatchoutoptions->correction_info;
}

GtEoplist *gt_querymatchoutoptions_eoplist(const GtQuerymatchoutoptions
                                             *querymatchoutoptions)
{
  if (querymatchoutoptions != NULL)
  {
    return querymatchoutoptions->eoplist;
  }
  return NULL;
}
