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
#include "extended/linspaceManagement.h"
#include "extended/scorehandler.h"

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

/* global alignment with affine gapcosts in linear space*/
GtUword gt_computeaffinelinearspace_generic(LinspaceManagement *spacemanager,
                                            GtScoreHandler *scorehandler,
                                            GtAlignment *align,
                                            const GtUchar *useq,
                                            GtUword ustart,
                                            GtUword ulen,
                                            const GtUchar *vseq,
                                            GtUword vstart,
                                            GtUword vlen);

/* global alignment with constant affine gapcosts in linear space,
 * only useful for DNA sequences*/
GtUword gt_computeaffinelinearspace(LinspaceManagement *spacemanager,
                                    GtAlignment *align,
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

/* local alignment with linear gapcosts in linear space */
GtWord gt_computeaffinelinearspace_local_generic(
                                              LinspaceManagement *spacemanager,
                                                 GtScoreHandler *scorehandler,
                                                 GtAlignment *align,
                                                 const GtUchar *useq,
                                                 GtUword ustart,
                                                 GtUword ulen,
                                                 const GtUchar *vseq,
                                                 GtUword vstart,
                                                 GtUword vlen);

/* local alignment with constant linear gapcosts in linear space,
 * only useful for DNA sequences */
GtWord gt_computeaffinelinearspace_local(LinspaceManagement *spacemanager,
                                         GtAlignment *align,
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

inline AffineAlignEdge set_edge(GtWord Rdist, GtWord Ddist, GtWord Idist);
#endif
