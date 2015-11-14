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

#ifndef RECONSTRUCTALIGNMENT_H
#define RECONSTRUCTALIGNMENT_H
#include "core/types_api.h"
#include "extended/affinealign.h"
#include "extended/alignment.h"
#include "extended/diagonalbandalign.h"
#include "extended/diagonalbandalign_affinegapcost.h"
#include "extended/maxcoordvalue.h"
#include "extended/scorehandler.h"

GtUword construct_trivial_deletion_alignment(GtAlignment *align,
                                             GtUword len,
                                             GtUword gapcost);

GtUword construct_trivial_insertion_alignment(GtAlignment *align,
                                              GtUword len,
                                              GtUword gapcost);

/* reconstruct alignment from square space table ED
 * use this function for global alignment with linear gapcosts*/
void reconstructalignment_from_EDtab(GtAlignment *align,
                                     GtUword * const *E,
                                     const GtUchar *useq,
                                     GtUword ustart,
                                     GtUword ulen,
                                     const GtUchar *vseq,
                                     GtUword vstart,
                                     GtUword vlen,
                                     const GtScoreHandler *scorehandler);

/* reconstruct alignment from square space table Ltab
 * use this function for lcoal alignment with linear gapscores*/
void reconstructalignment_from_Ltab(GtAlignment *align,
                                    GtWord **Ltabcolumn,
                                    Gtmaxcoordvalue *max,
                                    const GtUchar *useq,
                                    GtUword ustart,
                                    GtUword ulen,
                                    const GtUchar *vseq,
                                    GtUword vstart,
                                    GtUword vlen,
                                    const GtScoreHandler *scorehandler);

/* reconstruct alignment from crosspoint table, realting to midcolumn,
 * use this function for global or local alignment with linear or affine
 * gapcosts in diagonalband */
void reconstructalignment_from_Ctab(GtAlignment *align,
                                    const GtUword *Ctab,
                                    const GtUchar *useq,
                                    GtUword ustart,
                                    const GtUchar *vseq,
                                    GtUword vstart,
                                    GtUword vlen,
                                    const GtScoreHandler *scorehandler);

/* reconstruct alignment from crosspoints, crosspoints relating to diagonalband
 * use this function for alignment with linear gapcosts in diagonalband */
void reconstructalignment_from_Dtab(GtAlignment *align,
                                    const Diagentry *Dtab,GtUword ulen,
                                    GtUword vlen);

/* reconstruct alignment from crosspoints (affine gapcosts),
 * crosspointsrelating to diagonalband, use this function for alignment with
 * affine gapcosts in diagonalband */
void reconstructalignment_from_affineDtab(GtAlignment *align,
                                          const AffineDiagentry *Dtab,
                                          AffineAlignEdge edge,
                                          const GtUchar *useq, GtUword ulen,
                                          const GtUchar *vseq, GtUword vlen);
#endif
