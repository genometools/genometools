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
#include "core/minmax.h"
#include "match/diagband-struct.h"

/* called with real bpos */
#define GT_DIAGBANDSEED_DIAGONALBAND(AMAXLEN,LOGDIAGBANDWIDTH,APOS,BPOS)\
        (GT_DIAGBANDSEED_DIAGONAL(AMAXLEN,APOS,BPOS) >> (LOGDIAGBANDWIDTH))

typedef uint32_t GtDiagbandseedScore;

struct GtDiagbandStruct
{
  GtUword amaxlen, logdiagbandwidth, num_diagbands, used_diagbands;
  GtDiagbandseedScore *score;
  GtDiagbandseedPosition *lastpos;
};

GtUword gt_diagband_struct_num_diagbands(GtUword amaxlen,GtUword bmaxlen,
                                         GtUword logdiagbandwidth)
{
  return 1 + ((amaxlen + bmaxlen) >> logdiagbandwidth);
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
  diagband_struct->amaxlen = amaxlen;
  diagband_struct->logdiagbandwidth = logdiagbandwidth;
  /* diagband_score[0] and diagband_score[num_diagbands+1] remain zero as
     boundaries */
  diagband_struct->score = gt_calloc(diagband_struct->num_diagbands + 2,
                                     sizeof *diagband_struct->score);
  diagband_struct->score++; /* so we need not increment the index when
                               accessing score */
  diagband_struct->lastpos
    = gt_calloc(diagband_struct->num_diagbands,
                sizeof *diagband_struct->lastpos);
  return diagband_struct;
}

void gt_diagband_struct_delete(GtDiagbandStruct *diagband_struct)
{
  if (diagband_struct != NULL)
  {
    diagband_struct->score--; /* need to recover original base adress */
    gt_free(diagband_struct->score);
    gt_free(diagband_struct->lastpos);
    gt_free(diagband_struct);
  }
}

void gt_diagband_struct_single_update(GtDiagbandStruct *diagband_struct,
                                      GtUword apos,
                                      GtUword bpos,
                                      GtUword matchlength)
{
  GtUword diagband_idx;

  gt_assert(diagband_struct != NULL);
  diagband_idx = GT_DIAGBANDSEED_DIAGONALBAND(diagband_struct->amaxlen,
                                              diagband_struct->logdiagbandwidth,
                                              apos,
                                              bpos);
  gt_assert(diagband_idx < diagband_struct->num_diagbands);
  if (bpos >= diagband_struct->lastpos[diagband_idx] + matchlength)
  {
    /* no overlap */
    diagband_struct->lastpos[diagband_idx] = bpos;
    if (diagband_struct->score[diagband_idx] == 0)
    {
      diagband_struct->used_diagbands++;
    }
    diagband_struct->score[diagband_idx] += matchlength;
  } else
  {
    /* overlap: add positions after last counted position */
    if (diagband_struct->lastpos[diagband_idx] < bpos)
    {
      const GtUword addlength = bpos - diagband_struct->lastpos[diagband_idx];

      diagband_struct->lastpos[diagband_idx] = bpos;
      if (diagband_struct->score[diagband_idx] == 0)
      {
        diagband_struct->used_diagbands++;
      }
      diagband_struct->score[diagband_idx] += addlength;
    }
  }
}

GtUword gt_diagband_struct_coverage(const GtDiagbandStruct *diagband_struct,
                                    GtDiagbandseedPosition apos,
                                    GtDiagbandseedPosition bpos)
{
  GtUword diagband_idx;

  gt_assert(diagband_struct != NULL);
  diagband_idx = GT_DIAGBANDSEED_DIAGONALBAND(diagband_struct->amaxlen,
                                              diagband_struct->logdiagbandwidth,
                                              apos, bpos);
  return (GtUword) MAX(diagband_struct->score[diagband_idx + 1],
                       diagband_struct->score[diagband_idx - 1])
         + (GtUword) diagband_struct->score[diagband_idx];
}

void gt_diagband_struct_multi_update(
                         GtDiagbandStruct *diagband_struct,
                         const GtArrayGtDiagbandseedMaximalmatch *memstore)
{
  const GtDiagbandseedMaximalmatch *memstore_ptr;

  gt_assert(memstore != NULL);
  for (memstore_ptr = memstore->spaceGtDiagbandseedMaximalmatch;
       memstore_ptr < memstore->spaceGtDiagbandseedMaximalmatch +
                      memstore->nextfreeGtDiagbandseedMaximalmatch;
       memstore_ptr++)
  {
    gt_diagband_struct_single_update(diagband_struct,
                                     memstore_ptr->apos,
                                     memstore_ptr->bpos,
                                     memstore_ptr->len);
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
    memset(diagband_struct->score,0,
           sizeof *diagband_struct->score * diagband_struct->num_diagbands);
    memset(diagband_struct->lastpos,0,
           sizeof *diagband_struct->lastpos * diagband_struct->num_diagbands);
  } else
  {
    GtUword idx;

    if (seedstore != NULL)
    {
      for (idx = 0; idx < segment_length; idx++)
      {
        const GtUword diagband_idx
          = GT_DIAGBANDSEED_DIAGONALBAND(diagband_struct->amaxlen,
                                         diagband_struct->logdiagbandwidth,
                                         seedstore[idx].apos,
                                         seedstore[idx].bpos);
        diagband_struct->score[diagband_idx] = 0;
        diagband_struct->lastpos[diagband_idx] = 0;
      }
    } else
    {
      gt_assert(memstore != NULL);
      for (idx = 0; idx < segment_length; idx++)
      {
        const GtUword diagband_idx
          = GT_DIAGBANDSEED_DIAGONALBAND(diagband_struct->amaxlen,
                                         diagband_struct->logdiagbandwidth,
                                         memstore[idx].apos,
                                         memstore[idx].bpos);
        diagband_struct->score[diagband_idx] = 0;
        diagband_struct->lastpos[diagband_idx] = 0;
      }
    }
  }
  diagband_struct->used_diagbands = 0;
}
