/*
  Copyright (c) 2017 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2017 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/ma_api.h"
#include "core/assert_api.h"
#include "core/unused_api.h"
#include "core/str_api.h"
#include "core/minmax.h"
#include "core/intbits.h"
#include "core/chardef.h"
#include "core/score_matrix.h"
#include "match/diagband-struct.h"

/* called with real bpos */
GtUword gt_diagband_struct_num_diagbands(GtUword amaxlen,GtUword bmaxlen,
                                         GtUword logdiagbandwidth)
{
  return 1 + ((amaxlen + bmaxlen) >> logdiagbandwidth);
}

typedef uint32_t GtDiagbandseedScore;

struct GtDiagbandStruct
{
  GtUword amaxlen, bmaxlen, logdiagbandwidth, num_diagbands, used_diagbands,
          reset_from_matches, reset_with_memset;
  GtDiagbandseedScore *bcov;
  GtDiagbandseedPosition *lastpos;
  bool bpos_sorted;
};

static GtUword gt_diagbandseed_diagonalband(const GtDiagbandStruct
                                              *diagband_struct,
                                            GtDiagbandseedPosition apos,
                                            GtDiagbandseedPosition bpos)
{
  if (diagband_struct->bpos_sorted)
  {
    return GT_DIAGBANDSEED_DIAGONAL(diagband_struct->amaxlen,apos,bpos)
           >> diagband_struct->logdiagbandwidth;
  }
  return GT_DIAGBANDSEED_DIAGONAL(diagband_struct->bmaxlen,bpos,apos)
         >> diagband_struct->logdiagbandwidth;
}

bool gt_diagband_struct_empty(const GtDiagbandStruct *diagband_struct)
{
  return diagband_struct->used_diagbands == 0 ? true : false;
}

GtDiagbandStruct *gt_diagband_struct_new(GtUword amaxlen,GtUword bmaxlen,
                                         GtUword logdiagbandwidth)
{
  GtDiagbandStruct *diagband_struct = gt_malloc(sizeof *diagband_struct);

  diagband_struct->used_diagbands = 0;
  diagband_struct->num_diagbands
    = gt_diagband_struct_num_diagbands(amaxlen,bmaxlen,logdiagbandwidth);
  diagband_struct->bpos_sorted = true;
  diagband_struct->amaxlen = amaxlen;
  diagband_struct->bmaxlen = bmaxlen;
  diagband_struct->logdiagbandwidth = logdiagbandwidth;
  /* bcov[0] and bcov[num_diagbands+1] remain zero as boundaries */
  diagband_struct->bcov = gt_calloc(diagband_struct->num_diagbands + 2,
                                    sizeof *diagband_struct->bcov);
  diagband_struct->bcov++; /* so we need not increment the index when
                               accessing bcov */
  diagband_struct->lastpos
    = gt_calloc(diagband_struct->num_diagbands,
                sizeof *diagband_struct->lastpos);
  diagband_struct->reset_from_matches = 0;
  diagband_struct->reset_with_memset = 0;
  return diagband_struct;
}

void gt_diagband_struct_bpos_sorted_set(GtDiagbandStruct *diagband_struct,
                                        bool value)
{
  gt_assert(diagband_struct != NULL);
  diagband_struct->bpos_sorted = value;
}

void gt_diagbandseed_maxlen_update(GtDiagbandStruct *diagband_struct,
                                   GtUword amaxlen,GtUword bmaxlen)
{
  gt_assert(diagband_struct != NULL);
  diagband_struct->amaxlen = amaxlen;
  diagband_struct->bmaxlen = bmaxlen;
}

void gt_diagband_struct_reset_counts(const GtDiagbandStruct *diagband_struct,
                                     FILE *stream)
{
  fprintf(stream,"# number of resets of all used diagonal bands: " GT_WU,
          diagband_struct->reset_with_memset +
          diagband_struct->reset_from_matches);
  if (diagband_struct->reset_from_matches > 0)
  {
    fprintf(stream,"; resets from matches: " GT_WU,
            diagband_struct->reset_from_matches);
  }
  fprintf(stream,"\n");
}

void gt_diagband_struct_delete(GtDiagbandStruct *diagband_struct)
{
  if (diagband_struct != NULL)
  {
    diagband_struct->bcov--; /* need to recover original base adress */
    gt_free(diagband_struct->bcov);
    gt_free(diagband_struct->lastpos);
    gt_free(diagband_struct);
  }
}

/* for a given match of length <matchlength> ending a positions <apos> and
   <bpos> the sequence from A and from B, respectively, update the
   diagonal band b-coverage, which is the number of positions on B covered by
   the match. If previous matches have been added to the band before, then
   the positions overlapping with these on the B-sequence are not counted.*/

static void gt_diagband_struct_single_update(GtDiagbandStruct *diagband_struct,
                                             GtDiagbandseedPosition apos,
                                             GtDiagbandseedPosition bpos,
                                             GtDiagbandseedPosition matchlength)
{
  GtUword diagband_idx, keypos;

  gt_assert(diagband_struct != NULL);
  diagband_idx = gt_diagbandseed_diagonalband(diagband_struct, apos, bpos);
  keypos = diagband_struct->bpos_sorted ? bpos : apos;
  gt_assert(diagband_idx < diagband_struct->num_diagbands);
  if (diagband_struct->lastpos[diagband_idx] == 0 /* first matches */||
      /* match with end position keypos begins strictly after previous match */
      diagband_struct->lastpos[diagband_idx] + matchlength <= keypos)
  {
    /* no overlap */
    diagband_struct->lastpos[diagband_idx] = keypos;
    if (diagband_struct->bcov[diagband_idx] == 0)
    {
      diagband_struct->used_diagbands++;
    }
    diagband_struct->bcov[diagband_idx] += matchlength;
  } else
  {
    /* overlap: add positions after last counted position */
    if (diagband_struct->lastpos[diagband_idx] < keypos)
    {
      const GtUword addlength = keypos - diagband_struct->lastpos[diagband_idx];

      diagband_struct->lastpos[diagband_idx] = keypos;
      if (diagband_struct->bcov[diagband_idx] == 0)
      {
        diagband_struct->used_diagbands++;
      }
      diagband_struct->bcov[diagband_idx] += addlength;
    }
  }
}

static GtUword gt_diagband_struct_dband_coverage(
                                    const GtDiagbandStruct *diagband_struct,
                                    GtUword diagband_idx)
{
  gt_assert(diagband_struct != NULL);
  return (GtUword) MAX(diagband_struct->bcov[diagband_idx + 1],
                       diagband_struct->bcov[diagband_idx - 1])
         + (GtUword) diagband_struct->bcov[diagband_idx];
}

GtUword gt_diagband_struct_coverage(const GtDiagbandStruct *diagband_struct,
                                    GtDiagbandseedPosition apos,
                                    GtDiagbandseedPosition bpos)
{
  GtUword diagband_idx;

  gt_assert(diagband_struct != NULL);
  diagband_idx = gt_diagbandseed_diagonalband(diagband_struct,apos,bpos);
  return gt_diagband_struct_dband_coverage(diagband_struct,diagband_idx);
}

void gt_diagband_struct_mem_multi_update(GtDiagbandStruct *diagband_struct,
                                         const GtDiagbandseedMaximalmatch
                                           *memstore,
                                         GtUword numofmatches)
{
  GtUword idx;

  gt_assert(memstore != NULL);
  for (idx = 0; idx < numofmatches; idx++)
  {
    gt_diagband_struct_single_update(diagband_struct,
                                     memstore[idx].apos,
                                     memstore[idx].bpos,
                                     memstore[idx].len);
  }
}

void gt_diagband_struct_seed_multi_update(GtDiagbandStruct *diagband_struct,
                                          const GtSeedpairPositions *seedstore,
                                          GtUword segment_length,
                                          GtUword seedlength)
{
  GtUword idx;

  gt_assert(seedstore != NULL);
  for (idx = 0; idx < segment_length; idx++)
  {
    gt_assert (idx == 0 ||
               (diagband_struct->bpos_sorted &&
                seedstore[idx-1].bpos <= seedstore[idx].bpos) ||
               (!diagband_struct->bpos_sorted &&
                seedstore[idx-1].apos <= seedstore[idx].apos));
    gt_diagband_struct_single_update(diagband_struct,
                                     seedstore[idx].apos,
                                     seedstore[idx].bpos,
                                     seedlength);
  }
}

void gt_diagband_struct_reset(GtDiagbandStruct *diagband_struct,
                              const GtSeedpairPositions *seedstore,
                              const GtDiagbandseedMaximalmatch *memstore,
                              GtUword segment_length)
{
  gt_assert(diagband_struct != NULL);
  if (diagband_struct->used_diagbands * 3 >= diagband_struct->num_diagbands)
  { /* >= 33% of diagbands are used */
    memset(diagband_struct->bcov,0,
           sizeof *diagband_struct->bcov * diagband_struct->num_diagbands);
    memset(diagband_struct->lastpos,0,
           sizeof *diagband_struct->lastpos * diagband_struct->num_diagbands);
    diagband_struct->reset_with_memset++;
  } else
  {
    GtUword idx;

    if (seedstore != NULL)
    {
      for (idx = 0; idx < segment_length; idx++)
      {
        const GtUword diagband_idx
          = gt_diagbandseed_diagonalband(diagband_struct,
                                         seedstore[idx].apos,
                                         seedstore[idx].bpos);
        diagband_struct->bcov[diagband_idx] = 0;
        diagband_struct->lastpos[diagband_idx] = 0;
      }
    } else
    {
      gt_assert(memstore != NULL);
      for (idx = 0; idx < segment_length; idx++)
      {
        const GtUword diagband_idx
          = gt_diagbandseed_diagonalband(diagband_struct,
                                         memstore[idx].apos,
                                         memstore[idx].bpos);
        diagband_struct->bcov[diagband_idx] = 0;
        diagband_struct->lastpos[diagband_idx] = 0;
      }
    }
    diagband_struct->reset_from_matches++;
  }
  diagband_struct->used_diagbands = 0;
}

typedef int GtDiagbandScore;

struct GtDiagbandStatistics
{
  bool forward,
       compute_total_bcov,
       compute_total_score_seqpair,
       compute_has_two_seeds_on_same_diagonal;
  GtUword total_bcov,
          total_score_show_min;
  GtBitsequence *diagband_track;
  GtDiagbandseedPosition *previouspos;
  const GtScoreMatrix *score_matrix;
  GtUchar *a_buffer, *b_buffer;
  GtUword a_len, b_len;
  GtDiagbandScore score_threshold;
};

GtDiagbandStatistics *gt_diagband_statistics_new(const GtStr
                                                   *diagband_distance_arg,
                                                 bool forward)
{
  const char *arg = gt_str_get(diagband_distance_arg);
  GtDiagbandStatistics *diagband_statistics
    = gt_malloc(sizeof *diagband_statistics);

  diagband_statistics->forward = forward;
  diagband_statistics->compute_total_bcov = false;
  diagband_statistics->compute_total_score_seqpair = false;
  diagband_statistics->compute_has_two_seeds_on_same_diagonal = false;
  diagband_statistics->total_score_show_min = 0;
  if (strcmp(arg,"total_bcov") == 0)
  {
    diagband_statistics->compute_total_bcov = true;
  } else
  {
    if (strcmp(arg,"total_score_seqpair") == 0)
    {
      diagband_statistics->compute_total_score_seqpair = true;
    } else
    {
      if (strcmp(arg,"two_seeds") == 0)
      {
        diagband_statistics->compute_has_two_seeds_on_same_diagonal = true;
      } else
      {
        gt_assert(false);
      }
    }
  }
  diagband_statistics->total_bcov = 0;
  diagband_statistics->diagband_track = NULL;
  diagband_statistics->previouspos = NULL;
  diagband_statistics->score_matrix = NULL;
  diagband_statistics->a_buffer = NULL;
  diagband_statistics->b_buffer = NULL;
  diagband_statistics->score_threshold = INT_MIN;
  diagband_statistics->a_len = 0;
  diagband_statistics->b_len = 0;
  return diagband_statistics;
}

void gt_diagband_statistics_total_score_show_min_set(
              GtDiagbandStatistics *diagband_statistics,
              GtUword total_score_show_min)
{
  gt_assert(diagband_statistics != NULL);
  diagband_statistics->total_score_show_min = total_score_show_min;
}

void gt_diagband_statistics_score_matrix_set(
                GtDiagbandStatistics *diagband_statistics,
                const GtScoreMatrix *score_matrix,
                GtDiagbandScore score_threshold)
{
  gt_assert(diagband_statistics != NULL);
  diagband_statistics->score_matrix = score_matrix;
  diagband_statistics->score_threshold = score_threshold;
}

void gt_diagband_statistics_delete(GtDiagbandStatistics *diagband_statistics)
{
  if (diagband_statistics != NULL)
  {
    gt_free(diagband_statistics->a_buffer);
    gt_free(diagband_statistics->b_buffer);
    gt_free(diagband_statistics->diagband_track);
    gt_free(diagband_statistics->previouspos);
    gt_free(diagband_statistics);
  }
}

void gt_diagband_statistics_display(const GtDiagbandStatistics
                                           *diagband_statistics)
{
  gt_assert(diagband_statistics != NULL);
  if (diagband_statistics->compute_total_bcov)
  {
    printf("# forward=%s, total_diagband_bcov=" GT_WU "\n",
            diagband_statistics->forward ? "true" : "false",
            diagband_statistics->total_bcov);
  } else
  {
    gt_assert(diagband_statistics->compute_total_score_seqpair ||
              diagband_statistics->compute_has_two_seeds_on_same_diagonal);
  }
}

static void gt_diagband_statistics_bcov_add(
                                      GtDiagbandStatistics *diagband_statistics,
                                      const GtDiagbandStruct *diagband_struct,
                                      GtUword diagband_idx)
{
  if (!GT_ISIBITSET(diagband_statistics->diagband_track,diagband_idx))
  {
    diagband_statistics->total_bcov += diagband_struct->bcov[diagband_idx];
    GT_SETIBIT(diagband_statistics->diagband_track,diagband_idx);
  }
}

#ifndef NDEBUG
static GtDiagbandScore gt_diagband_statistics_eval_score(
                                             const GtScoreMatrix *score_matrix,
                                             const GtUchar *a_encoded,
                                             const GtUchar *b_encoded,
                                             GtUword length)
{
  GtDiagbandScore score_sum = 0;
  GtUword idx;

  for (idx = 0; idx < length; idx++)
  {
    if (ISSPECIAL(a_encoded[idx]) || ISSPECIAL(b_encoded[idx]))
    {
      return 0;
    }
    score_sum += gt_score_matrix_get_score(score_matrix,a_encoded[idx],
                                           b_encoded[idx]);
  }
  return score_sum;
}

static GtDiagbandScore gt_diagband_statistics_total_score_in_seqpair(
                                   const GtScoreMatrix *score_matrix,
                                   bool same_seq,
                                   unsigned int q_value,
                                   GtDiagbandScore score_threshold,
                                   const GtUchar *a_encoded,
                                   GtUword alen,
                                   const GtUchar *b_encoded,
                                   GtUword blen)
{
  const GtUchar *aptr, *bptr;
  GtDiagbandScore total_score = 0;

  if (alen >= q_value && blen >= q_value)
  {
    for (aptr = a_encoded; aptr <= a_encoded + alen - q_value; aptr++)
    {
      if (same_seq)
      {
        bptr = b_encoded + (aptr - a_encoded) + q_value;
      } else
      {
        bptr = b_encoded;
      }
      for (/* Nothing */; bptr <= b_encoded + blen - q_value;
           bptr++)
      {
        gt_assert(aptr + q_value - 1 <= a_encoded + alen - 1);
        gt_assert(bptr + q_value - 1 <= b_encoded + blen - 1);
        const GtDiagbandScore this_score
          = gt_diagband_statistics_eval_score(score_matrix,
                                              aptr,
                                              bptr,
                                              q_value);
        if (this_score >= score_threshold)
        {
          /*printf(GT_WU " " GT_WU " with score %d\n",
                 (GtUword) (aptr - a_encoded),
                 (GtUword) (bptr - b_encoded),this_score);*/
          total_score += this_score;
        }
      }
    }
  }
  return total_score;
}
#endif

void gt_diagband_statistics_header(const GtDiagbandStatistics
                                     *diagband_statistics)
{
  if (diagband_statistics->compute_has_two_seeds_on_same_diagonal)
  {
    printf("# Fields: dbnum, querynum\n");
  }
}

void gt_diagband_statistics_add(void *v_diagband_statistics,
                                GT_UNUSED bool bpos_sorted,
                                /* remove GT_UNUSED once arguments are used */
                                const GtEncseq *aencseq,
                                const GtEncseq *bencseq,
                                GtUword aseqnum,
                                GtUword bseqnum,
                                const GtDiagbandStruct *diagband_struct,
                                const GtDiagbandseedMaximalmatch *memstore,
                                unsigned int seedlength,
                                const GtSeedpairPositions *seedstore,
                                const uint8_t *segment_scores,
                                GtUword segment_length)
{
  GtUword idx;
  GtDiagbandStatistics *diagband_statistics
    = (GtDiagbandStatistics *) v_diagband_statistics;

  if (diagband_statistics->compute_has_two_seeds_on_same_diagonal ||
      diagband_statistics->compute_total_bcov)
  {
    if (diagband_statistics->diagband_track == NULL)
    {
      GT_INITBITTAB(diagband_statistics->diagband_track,
                    diagband_struct->num_diagbands);
    } else
    {
      GT_CLEARBITTAB(diagband_statistics->diagband_track,
                     diagband_struct->num_diagbands);
    }
  }
  if (diagband_statistics->compute_has_two_seeds_on_same_diagonal)
  {
    if (diagband_statistics->previouspos == NULL)
    {
      diagband_statistics->previouspos
        = gt_calloc(diagband_struct->num_diagbands,
                    sizeof *diagband_statistics->previouspos);
    } else
    {
      memset(diagband_statistics->previouspos,0,
             sizeof *diagband_statistics->previouspos *
             diagband_struct->num_diagbands);
    }
  }
  if (diagband_statistics->score_matrix != NULL)
  {
    GtUword startpos;

    if (diagband_statistics->a_buffer == NULL)
    {
      diagband_statistics->a_buffer
        = gt_malloc(sizeof *diagband_statistics->a_buffer *
                    gt_encseq_max_seq_length(aencseq));
    }
    startpos = gt_encseq_seqstartpos(aencseq, aseqnum);
    diagband_statistics->a_len = gt_encseq_seqlength(aencseq, aseqnum);
    gt_encseq_extract_encoded(aencseq, diagband_statistics->a_buffer,startpos,
                              startpos + diagband_statistics->a_len - 1);
    if (diagband_statistics->b_buffer == NULL)
    {
      diagband_statistics->b_buffer
        = gt_malloc(sizeof *diagband_statistics->b_buffer *
                    gt_encseq_max_seq_length(bencseq));
    }
    startpos = gt_encseq_seqstartpos(bencseq, bseqnum);
    diagband_statistics->b_len = gt_encseq_seqlength(bencseq, bseqnum);
    gt_encseq_extract_encoded(bencseq, diagband_statistics->b_buffer,startpos,
                              startpos + diagband_statistics->b_len - 1);
  }
  if (seedstore != NULL)
  {
    GtUword total_score_seqpair = 0, count_two_seeds_on_same_diagonal = 0;

    for (idx = 0; idx < segment_length; idx++)
    {
      const GtUword diagband_idx
        = gt_diagbandseed_diagonalband(diagband_struct,
                                       seedstore[idx].apos,
                                       seedstore[idx].bpos);
      if (diagband_statistics->compute_total_bcov)
      {
        gt_diagband_statistics_bcov_add(diagband_statistics,diagband_struct,
                                        diagband_idx);
      } else
      {
        if (diagband_statistics->compute_total_score_seqpair)
        {
          if (segment_scores != NULL)
          {
            total_score_seqpair += (GtUword) segment_scores[idx];
          }
        } else
        {
          const GtDiagbandseedPosition bpos = seedstore[idx].bpos;

          gt_assert(diagband_statistics->
                      compute_has_two_seeds_on_same_diagonal &&
                    diagband_idx < diagband_struct->num_diagbands &&
                    diagband_statistics->previouspos[diagband_idx] <= bpos);
          /* there was a previous match and current match with end position
             bpos begins strictly after previous match */

#ifdef OUTDIAGONALHITS
          if (diagband_statistics->previouspos[diagband_idx] == 0)
          {
            printf("aseq=" GT_WU ",bseq=" GT_WU ",diag= " GT_WU
                   "apos=%u,bpos=%u is first\n",
                   aseqnum,bseqnum,diagband_idx,seedstore[idx].apos,bpos);
          }
#endif
          if (diagband_statistics->previouspos[diagband_idx] > 0 &&
              diagband_statistics->previouspos[diagband_idx] + seedlength
              <= bpos)
          {
            if (!GT_ISIBITSET(diagband_statistics->diagband_track,diagband_idx))
            {
              /* Event has not been found on this diagonal */
              count_two_seeds_on_same_diagonal++;
              GT_SETIBIT(diagband_statistics->diagband_track,diagband_idx);
            }
#ifdef OUTDIAGONALHITS
            printf("aseq=" GT_WU ",bseq=" GT_WU ",diag= " GT_WU
                   "apos=%u,bpos=%u\n",
                   aseqnum,bseqnum,diagband_idx,seedstore[idx].apos,bpos);
#endif
          }
          diagband_statistics->previouspos[diagband_idx] = bpos;
        }
      }
    }
#ifndef NDEBUG
    if (segment_scores != NULL &&
        diagband_statistics->compute_total_score_seqpair &&
        diagband_statistics->score_matrix != NULL)
    {
      GtDiagbandScore bf_total_score_seqpair;

      bf_total_score_seqpair = gt_diagband_statistics_total_score_in_seqpair(
                                   diagband_statistics->score_matrix,
                                   (aencseq == bencseq && aseqnum == bseqnum)
                                      ? true
                                      : false,
                                   seedlength,
                                   diagband_statistics->score_threshold,
                                   diagband_statistics->a_buffer,
                                   diagband_statistics->a_len,
                                   diagband_statistics->b_buffer,
                                   diagband_statistics->b_len);
      gt_assert(bf_total_score_seqpair == total_score_seqpair);
    }
#endif
    if (total_score_seqpair >= diagband_statistics->total_score_show_min)
    {
      printf(GT_WU " " GT_WU " " GT_WU "\n",aseqnum,bseqnum,
                                            total_score_seqpair);
    }
    if (count_two_seeds_on_same_diagonal > 0)
    {
      printf(GT_WU "\t" GT_WU "\n",aseqnum,bseqnum);
    }
  } else
  {
    gt_assert(memstore != NULL &&
              !diagband_statistics->compute_total_score_seqpair &&
              !diagband_statistics->compute_has_two_seeds_on_same_diagonal);
    for (idx = 0; idx < segment_length; idx++)
    {
      const GtUword diagband_idx
        = gt_diagbandseed_diagonalband(diagband_struct,
                                       memstore[idx].apos,
                                       memstore[idx].bpos);
      if (diagband_statistics->compute_total_bcov)
      {
        gt_diagband_statistics_bcov_add(diagband_statistics,diagband_struct,
                                        diagband_idx);
      }
    }
  }
}
