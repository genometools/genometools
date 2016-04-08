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
#include "core/unused_api.h"
#include "core/types_api.h"
#include "core/minmax.h"
#include "extended/alignment.h"
#include "extended/linearalign.h"
#include "extended/linspace_management.h"
#include "revcompl.h"
#include "seed-extend.h"
#include "ft-polish.h"
#include "querymatch-align.h"

static void eoplist_reverse_order(uint8_t *start,uint8_t *end)
{
  if (start < end)
  {
    uint8_t *fwd, *bwd;

    for (fwd = start, bwd = end; fwd < bwd; fwd++, bwd--)
    {
      uint8_t tmp = *fwd;
      *fwd = *bwd;
      *bwd = tmp;
    }
  }
}

static void converteoplist2alignment(GtAlignment *alignment,
                                     const GtArrayuint8_t *eoplist)
{
  uint8_t *ptr;

  for (ptr = eoplist->spaceuint8_t;
       ptr < eoplist->spaceuint8_t + eoplist->nextfreeuint8_t; ptr++)
  {
    if (*ptr == FT_EOPCODE_DELETION)
    {
      gt_alignment_add_deletion(alignment);
    } else
    {
      if (*ptr == FT_EOPCODE_INSERTION)
      {
        gt_alignment_add_insertion(alignment);
      } else
      {
        gt_alignment_add_replacement_multi(alignment,(GtUword) ((*ptr) + 1));
      }
    }
  }
}

void gt_querymatch_showeoplist(const GtArrayuint8_t *eoplist)
{
  uint8_t *ptr;

  gt_assert(eoplist != NULL);
  if (eoplist->nextfreeuint8_t == 0)
  {
    printf("[]\n");
    return;
  }
  for (ptr = eoplist->spaceuint8_t;
       ptr < eoplist->spaceuint8_t + eoplist->nextfreeuint8_t;
       ptr++)
  {
    if (*ptr == FT_EOP_DELETION)
    {
      printf("D");
    } else
    {
      if (*ptr == FT_EOP_INSERTION)
      {
        printf("I");
      } else
      {
        printf("R %d",(int) ((*ptr) + 1));
      }
    }
    if (ptr < eoplist->spaceuint8_t + eoplist->nextfreeuint8_t - 1)
    {
      printf(",");
    } else
    {
      printf("\n");
    }
  }
}

struct GtQuerymatchoutoptions
{
  GtUword totallength,
          alignmentwidth,
          useqbuffer_size,
          vseqbuffer_size;
  GtArrayuint8_t eoplist;
  GtFronttrace *front_trace;
  GtGreedyextendmatchinfo *ggemi;
  GtUchar *useqbuffer, *vseqbuffer, *alignment_show_buffer;
  const GtUchar *characters;
  GtAlignment *alignment;
  GtEncseqReader *esr_for_align_show;
  GtUchar wildcardshow;
  GtSeqpaircoordinates correction_info;
  GtLinspaceManagement *linspace_spacemanager;
  GtScoreHandler *linspace_scorehandler;
  bool always_polished_ends,
       generatealignment,
       showeoplist;
  Polishing_info *pol_info;
};

GtQuerymatchoutoptions *gt_querymatchoutoptions_new(bool generatealignment,
                                                    bool showeoplist,
                                                    GtUword alignmentwidth)
{
  GtQuerymatchoutoptions *querymatchoutoptions
    = gt_malloc(sizeof *querymatchoutoptions);

  querymatchoutoptions->generatealignment
    = generatealignment || alignmentwidth > 0 || showeoplist;
  querymatchoutoptions->showeoplist = showeoplist;
  querymatchoutoptions->alignmentwidth = alignmentwidth;
  querymatchoutoptions->front_trace = NULL;
  querymatchoutoptions->ggemi = NULL;
  querymatchoutoptions->useqbuffer = NULL;
  querymatchoutoptions->useqbuffer_size = 0;
  querymatchoutoptions->vseqbuffer = NULL;
  querymatchoutoptions->vseqbuffer_size = 0;
  GT_INITARRAY(&querymatchoutoptions->eoplist,uint8_t);
  querymatchoutoptions->totallength = GT_UWORD_MAX; /* not yet known */
  querymatchoutoptions->esr_for_align_show = NULL;
  querymatchoutoptions->characters = NULL;
  querymatchoutoptions->alignment = NULL;
  querymatchoutoptions->linspace_spacemanager = NULL;
  querymatchoutoptions->linspace_scorehandler = NULL;
  querymatchoutoptions->always_polished_ends = true;
  querymatchoutoptions->alignment_show_buffer = NULL;
  if (generatealignment)
  {
    querymatchoutoptions->alignment = gt_alignment_new();
  }
  if (alignmentwidth > 0)
  {
    querymatchoutoptions->alignment_show_buffer
      = gt_alignment_buffer_new(alignmentwidth);
    querymatchoutoptions->linspace_spacemanager = gt_linspace_management_new();
    querymatchoutoptions->linspace_scorehandler = gt_scorehandler_new(0,1,0,1);
  }
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
    if (querymatchoutoptions->generatealignment)
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
    if (querymatchoutoptions->alignment != NULL)
    {
      if (always_polished_ends)
      {
        gt_alignment_polished_ends(querymatchoutoptions->alignment,
                                   querymatchoutoptions->pol_info,true);
      }
      if (gt_querymatch_seed_display(display_flag))
      {
        gt_alignment_seed_display_set(querymatchoutoptions->alignment);
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
    gt_alignment_delete(querymatchoutoptions->alignment);
    GT_FREEARRAY(&querymatchoutoptions->eoplist,uint8_t);
    gt_alignment_buffer_delete(querymatchoutoptions->alignment_show_buffer);
    gt_encseq_reader_delete(querymatchoutoptions->esr_for_align_show);
    gt_linspace_management_delete(querymatchoutoptions->linspace_spacemanager);
    gt_scorehandler_delete(querymatchoutoptions->linspace_scorehandler);
    polishing_info_delete(querymatchoutoptions->pol_info);
    gt_free(querymatchoutoptions);
  }
}

static bool seededmatch2eoplist(GtQuerymatchoutoptions *querymatchoutoptions,
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
                                bool greedyextension)
{
  GtUword ulen, vlen, ustart, vstart;
  bool alignment_succeeded = true;
  Polished_point right_best_polished_point = {0,0,0};
  Polished_point left_best_polished_point = {0,0,0};
  GtUword pol_size;
  GtSeqpaircoordinates *coords;

  gt_assert(querymatchoutoptions != NULL &&
            querymatchoutoptions->pol_info != NULL);
  pol_size = GT_MULT2(querymatchoutoptions->pol_info->cut_depth);
  querymatchoutoptions->eoplist.nextfreeuint8_t = 0;
  if (querymatchoutoptions->totallength == GT_UWORD_MAX)
  {
    querymatchoutoptions->totallength = gt_encseq_total_length(encseq);
  }
  ustart = seedpos1 + seedlen;
  vstart = seedpos2 + seedlen;
  gt_assert(dbstart + dblen >= ustart);
  ulen = dbstart + dblen - ustart;
  gt_assert(abs_querystart + querylen >= vstart);
  vlen = abs_querystart + querylen - vstart;
  if (ulen > 0 && vlen > 0)
  {
    if (gt_align_front_prune_edist(true,
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
                                   vlen) == ulen + vlen + 1)
    {
      alignment_succeeded = false;
    } else
    {
      if (querymatchoutoptions->front_trace != NULL)
      {
        front_trace2eoplist(querymatchoutoptions->always_polished_ends,
                            &querymatchoutoptions->eoplist,
                            querymatchoutoptions->front_trace,
                            &right_best_polished_point,
                            pol_size,
                            querymatchoutoptions->pol_info->match_score,
                            querymatchoutoptions->pol_info->difference_score,
                            NULL,
                            ulen,
                            NULL,
                            vlen);
      }
    }
    if (querymatchoutoptions->front_trace != NULL)
    {
      gt_assert(querymatchoutoptions->generatealignment);
      front_trace_reset(querymatchoutoptions->front_trace,ulen+vlen);
    }
  }
  if (alignment_succeeded)
  {
    if (querymatchoutoptions->generatealignment)
    {
      front_trace_multireplacement(&querymatchoutoptions->eoplist,seedlen);
    }
    if (seedpos1 > dbstart && seedpos2 > abs_querystart)
    {
      ulen = seedpos1 - dbstart;
      vlen = seedpos2 - abs_querystart;
      if (gt_align_front_prune_edist(false,
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
                                     vlen) == ulen + vlen + 1)
      {
        alignment_succeeded = false;
      } else
      {
        if (querymatchoutoptions->front_trace != NULL)
        {
          GtUword eoplistlen = querymatchoutoptions->eoplist.nextfreeuint8_t;

          front_trace2eoplist(querymatchoutoptions->always_polished_ends,
                              &querymatchoutoptions->eoplist,
                              querymatchoutoptions->front_trace,
                              &left_best_polished_point,
                              pol_size,
                              querymatchoutoptions->pol_info->match_score,
                              querymatchoutoptions->pol_info->difference_score,
                              NULL,
                              ulen,
                              NULL,
                              vlen);
          eoplist_reverse_order(querymatchoutoptions->eoplist.spaceuint8_t +
                                 eoplistlen,
                                querymatchoutoptions->eoplist.spaceuint8_t +
                                 querymatchoutoptions->eoplist.nextfreeuint8_t -
                                 1);
        }
      }
      if (querymatchoutoptions->front_trace != NULL)
      {
        front_trace_reset(querymatchoutoptions->front_trace,ulen+vlen);
      }
    }
  }
  coords = &querymatchoutoptions->correction_info;
  if (alignment_succeeded)
  {
    GtUword leftcolumn, rightcolumn;

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
  } else
  {
    coords->uoffset = coords->voffset = 0;
    coords->ulen = dblen;
    coords->vlen = querylen;
    coords->sumdist = 0;
  }
  return alignment_succeeded;
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
                                               bool greedyextension)
{
  bool seededalignment = false;

  gt_assert(querymatchoutoptions != NULL);
  if (querymatchoutoptions->alignmentwidth > 0 &&
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
    if (seededmatch2eoplist(querymatchoutoptions,
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
                            greedyextension))
    {
      if (querymatchoutoptions->generatealignment)
      {
        gt_alignment_reset(querymatchoutoptions->alignment);
        converteoplist2alignment(querymatchoutoptions->alignment,
                                 &querymatchoutoptions->eoplist);
      }
      if (querymatchoutoptions->alignmentwidth > 0)
      {
        gt_assert(dbstart <= seedpos1);
        gt_alignment_set_seedoffset(querymatchoutoptions->alignment,
                                    seedpos1 - dbstart,
                                    seedlen);
        gt_alignment_set_seqs(querymatchoutoptions->alignment,
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
      if (querymatchoutoptions->generatealignment)
      {
#ifndef NDEBUG
        GtUword linedist;
#endif
        gt_assert(!greedyextension);
        gt_alignment_set_seqs(querymatchoutoptions->alignment,
                              querymatchoutoptions->useqbuffer,
                              dblen,
                              querymatchoutoptions->vseqbuffer,
                              querylen);
#ifndef NDEBUG
        linedist =
#else
        (void)
#endif
        gt_linearalign_compute_generic(
                              querymatchoutoptions->linspace_spacemanager,
                              querymatchoutoptions->linspace_scorehandler,
                              querymatchoutoptions->alignment,
                              querymatchoutoptions->useqbuffer,
                              0,
                              dblen,
                              querymatchoutoptions->vseqbuffer,
                              0,
                              querylen);
        gt_assert(linedist <= edist);
      }
    }
  } else
  {
    if (querymatchoutoptions->alignmentwidth > 0)
    {
      gt_alignment_set_seqs(querymatchoutoptions->alignment,
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
                                            GT_UNUSED bool verify_alignment,
                                            FILE *fp)
{
  if (querymatchoutoptions != NULL)
  {
    if (querymatchoutoptions->alignmentwidth > 0)
    {
      if (distance > 0)
      {
        gt_alignment_show_generic(querymatchoutoptions->alignment_show_buffer,
                                  false,
                                  querymatchoutoptions->alignment,
                                  fp,
                                  (unsigned int)
                                  querymatchoutoptions->alignmentwidth,
                                  querymatchoutoptions->characters,
                                  querymatchoutoptions->wildcardshow);
      } else
      {
        gt_alignment_exact_show(querymatchoutoptions->alignment_show_buffer,
                                querymatchoutoptions->alignment,
                                fp,
                                querymatchoutoptions->alignmentwidth,
                                querymatchoutoptions->characters);
      }
      if (verify_alignment)
      {
        (void) gt_alignment_check_edist(querymatchoutoptions->alignment,
                                        distance,NULL);
      }
    }
    if (querymatchoutoptions->showeoplist)
    {
      if (distance > 0)
      {
        gt_alignment_show_multieop_list(querymatchoutoptions->alignment, fp);
      } else
      {
        fprintf(fp, "[]\n");
      }
    }
    if (querymatchoutoptions->alignmentwidth > 0 ||
        querymatchoutoptions->showeoplist)
    {
      gt_alignment_reset(querymatchoutoptions->alignment);
    }
  }
}

const GtSeqpaircoordinates *gt_querymatchoutoptions_correction_get(
              const GtQuerymatchoutoptions *querymatchoutoptions)
{
  gt_assert(querymatchoutoptions != NULL);
  return &querymatchoutoptions->correction_info;
}

const GtAlignment *gt_querymatchoutoptions_alignment_get(
              const GtQuerymatchoutoptions *querymatchoutoptions)
{
  gt_assert(querymatchoutoptions != NULL);
  return querymatchoutoptions->alignment;
}
