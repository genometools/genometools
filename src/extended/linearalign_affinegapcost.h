/*
  Copyright (c) 2015 Annika <annika.seidel@studium.uni-hamburg.de>
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

#ifndef LINEARALIGN_AFFINEGAPCOST_H
#define LINEARALIGN_AFFINEGAPCOST_H

#include "core/unused_api.h"
#include "core/types_api.h"
#include "extended/affinealign.h"
#include "extended/alignment.h"

typedef struct {
  GtUword idx;
  AffineAlignEdge edge;
} Rnode;

typedef struct {
  Rnode val_R, val_D, val_I;
} Rtabentry;

void gt_checkaffinelinearspace(GT_UNUSED bool forward,
                               const GtUchar *useq,
                               GtUword ulen,
                               const GtUchar *vseq,
                               GtUword vlen);

void gt_checkaffinelinearspace_local(GT_UNUSED bool forward,
                                     const GtUchar *useq,
                                     GtUword ulen,
                                     const GtUchar *vseq,
                                     GtUword vlen);

void gt_computeaffinelinearspace(GtAlignment *align,
                                 const GtUchar *useq,
                                 GtUword ustart,
                                 GtUword ulen,
                                 const GtUchar *vseq,
                                 GtUword vstart,
                                 GtUword vlen,
                                 GtUword matchcost,
                                 GtUword mismatchcost,
                                 GtUword gap_opening,
                                 GtUword gap_extension);

void gt_computeaffinelinearspace_local(GtAlignment *align,
                                       const GtUchar *useq,
                                       GtUword ustart,
                                       GtUword ulen,
                                       const GtUchar *vseq,
                                       GtUword vstart,
                                       GtUword vlen,
                                       GtWord matchscore,
                                       GtWord mismatchscore,
                                       GtWord gap_opening,
                                       GtWord gap_extension);

AffineAlignEdge minAdditionalCosts(const AffinealignDPentry *entry,
                                   const AffineAlignEdge edge,
                                   GtUword gap_opening);

inline void firstAtabRtabentry(AffinealignDPentry *Atabcolumn,
                               GtUword gap_opening,
                               AffineAlignEdge edge);
#endif
