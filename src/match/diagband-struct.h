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

#ifndef DIAGBAND_STRUCT_H
#define DIAGBAND_STRUCT_H
#include <inttypes.h>
#include "core/types_api.h"

#define GT_DIAGBANDSEED_DIAGONAL(AMAXLEN,APOS,BPOS)\
        ((AMAXLEN) + (GtUword) (BPOS) - (GtUword) (APOS))

typedef struct GtDiagbandStruct GtDiagbandStruct;

GtUword gt_diagband_struct_num_diagbands(GtUword amaxlen,GtUword bmaxlen,
                                         GtUword logdiagbandwidth);

bool gt_diagband_struct_empty(const GtDiagbandStruct *diagband_struct);

GtDiagbandStruct *gt_diagband_struct_new(GtUword amaxlen,GtUword bmaxlen,
                                         GtUword logdiagbandwidth);

void gt_diagband_struct_delete(GtDiagbandStruct *diagband_struct);

void gt_diagband_struct_single_update(GtDiagbandStruct *diagband_struct,
                                      GtUword apos,
                                      GtUword bpos,
                                      GtUword matchlength);

typedef uint32_t GtDiagbandseedPosition;

GtUword gt_diagband_struct_coverage(const GtDiagbandStruct *diagband_struct,
                                    GtDiagbandseedPosition apos,
                                    GtDiagbandseedPosition bpos);

typedef struct
{
  GtDiagbandseedPosition apos, bpos, len;
} GtDiagbandseedMaximalmatch;

typedef struct
{
  GtDiagbandseedMaximalmatch *spaceGtDiagbandseedMaximalmatch;
  GtUword allocatedGtDiagbandseedMaximalmatch,
          nextfreeGtDiagbandseedMaximalmatch;
} GtArrayGtDiagbandseedMaximalmatch;

void gt_diagband_struct_multi_update(GtDiagbandStruct *diagband_struct,
                                     const GtArrayGtDiagbandseedMaximalmatch
                                       *memstore);

typedef struct
{
  GtDiagbandseedPosition apos, /* secondary key */
                         bpos; /* primary key */
} GtSeedpairPositions;

void gt_diagband_struct_statistics(const GtDiagbandStruct *diagband_struct,
                                   FILE *stream);

void gt_diagband_struct_reset(GtDiagbandStruct *diagband_struct,
                              const GtSeedpairPositions *seedstore,
                              const GtDiagbandseedMaximalmatch *memstore,
                              GtUword segment_length);

typedef struct GtDiagbandStatistics GtDiagbandStatistics;

GtDiagbandStatistics *gt_diagband_statistics_new(GT_UNUSED const GtStr
                                                   *diagband_statistics_arg);

void gt_diagband_statistics_delete(GtDiagbandStatistics *diagband_statistics);

void gt_diagband_statistics_display(const GtDiagbandStatistics
                                      *diagband_statistics);

void gt_diagband_statistics_add(void *v_diagband_statistics,
                                GtUword aseqnum,
                                GtUword bseqnum,
                                GtDiagbandStruct *diagband_struct,
                                const GtDiagbandseedMaximalmatch *memstore,
                                unsigned int seedlength,
                                const GtSeedpairPositions *seedstore,
                                GtUword segment_length);

#endif
