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

#ifndef DIAGONALBANDALIGN_AFFINEGAPCOST_H
#define DIAGONALBANDALIGN_AFFINEGAPCOST_H
#include "core/error.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "extended/alignment.h"
#include "extended/linspaceManagement.h"

typedef struct {
  Diagentry val_R, val_D, val_I;
} AffineDiagentry;

void gt_checkdiagonalbandaffinealign(GT_UNUSED bool forward,
                                     const GtUchar *useq,
                                     GtUword ulen,
                                     const GtUchar *vseq,
                                     GtUword vlen);

/* creating alignment with diagonalband in linear space
 * with affine gapcosts (DNA or protein) */
void gt_computediagonalbandaffinealign_generic(LinspaceManagement *spacemanager,
                                              GtScoreHandler *scorehandler,
                                              GtAlignment *align,
                                              const GtUchar *useq,
                                              GtUword ustart, GtUword ulen,
                                              const GtUchar *vseq,
                                              GtUword vstart, GtUword vlen,
                                              GtWord left_dist,
                                              GtWord right_dist);

/* creating alignment with diagonalband in linear space
 * with constant affine gapcosts, only useful for DNA sequences */
void gt_computediagonalbandaffinealign(LinspaceManagement *spacemanager,
                                          GtAlignment *align,
                                          const GtUchar *useq,
                                          GtUword ustart, GtUword ulen,
                                          const GtUchar *vseq,
                                          GtUword vstart, GtUword vlen,
                                          GtWord left_dist,
                                          GtWord right_dist,
                                          GtUword matchcost,
                                          GtUword mismatchcost,
                                          GtUword gap_opening,
                                          GtUword gap_extension);

/* calculate only distance with diagonalband in square space with
 * affine gapcosts */
GtWord diagonalband_square_space_affine(const GtUchar *useq,
                                         GtUword ustart,
                                         GtUword ulen,
                                         const GtUchar *vseq,
                                         GtUword vstart,
                                         GtUword vlen,
                                         GtWord left_dist,
                                         GtWord right_dist,
                                         GtScoreHandler *scorehandler);

/* calculate alignment with diagonalband in square space
 * with affine gapcosts, to use it in linear context you
 * have to generate an spacemanager before, in any other case it can be NULL */
GtWord diagonalbandalignment_in_square_space_affine_generic(
                                                   LinspaceManagement *space,
                                                   GtScoreHandler *scorehandler,
                                                   GtAlignment *align,
                                                   const GtUchar *useq,
                                                   GtUword ustart,
                                                   GtUword ulen,
                                                   const GtUchar *vseq,
                                                   GtUword vstart,
                                                   GtUword vlen,
                                                   GtWord left_dist,
                                                   GtWord right_dist);

/* same with constant cost values, only useful for DNA sequences */
GtWord diagonalbandalignment_in_square_space_affine(LinspaceManagement *space,
                                                    GtAlignment *align,
                                                    const GtUchar *useq,
                                                    GtUword ustart,
                                                    GtUword ulen,
                                                    const GtUchar *vseq,
                                                    GtUword vstart,
                                                    GtUword vlen,
                                                    GtWord left_dist,
                                                    GtWord right_dist,
                                                    GtUword matchcost,
                                                    GtUword mismatchcost,
                                                    GtUword gap_opening,
                                                    GtUword gap_extension);
#endif
