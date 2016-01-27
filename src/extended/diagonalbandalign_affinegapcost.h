/*
  Copyright (c) 2015 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
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
#include "extended/linspace_management.h"

/* <AffineDiagentry> objects are entries of DP-matrix with a diagonal band in
   context of affine gap costs. */
typedef struct {
  GtDiagAlignentry val_R, val_D, val_I;
} GtAffineDiagAlignentry;

/* Computes a global alignment within a diagonal band with affine gapcosts in
   linear space. Use of this function requires an initialised <spacemanager>,
   the target alignment <align> and input sequences <useq> and <vseq>, with the
   regions to align given by their start positions <ustart> and <vstart> and
   lengths <ulen> and <vlen>. The cost values are specified by an initialised
   <scorehandler>. <left_dist> and <right_dist> give lower and upper bound of
   a diagonal band in which DP-matrix is valid. Returns affine cost
   value of calculated global alignment. */
void   gt_diagonalbandalign_affinegapcost_compute_generic(
                                            GtLinspaceManagement *spacemanager,
                                            const GtScoreHandler *scorehandler,
                                            GtAlignment *align,
                                            const GtUchar *useq,
                                            GtUword ustart, GtUword ulen,
                                            const GtUchar *vseq,
                                            GtUword vstart, GtUword vlen,
                                            GtWord left_dist,
                                            GtWord right_dist);

/* Computes a global alignment within a diaognal band with affine gapcosts in
   linear space and constant cost values. Use of this function requires an
   initialised <spacemanager>, the target alignment <align> and input sequences
   <useq> and <vseq>, with the regions to align given by their start positions
   <ustart> and <vstart> and lengths <ulen> and <vlen>. The cost values are
   specified by <matchcost>, <mismatchcost> and <gapcost>. <left_dist> and
   <right_dist> give lower and upper bound of a diagonal band in which DP-matrix
   is valid. Returns affine cost value of calculated global alignment. */
void   gt_diagonalbandalign_affinegapcost_compute(
                                            GtLinspaceManagement *spacemanager,
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

/* Computes only the distance of a global alignment within a diagonal band with
   affine gapcosts in square space. Use of this function requires the input
   sequences <useq> and <vseq>, with the regions to align given by their start
   positions <ustart> and <vstart> and lengths <ulen> and <vlen>. <left_dist>
   and <right_dist> give lower and upper bound of a diagonal band in which
   DP-matrix is valid. The cost values are specified by <scorehandler>. Returns
   cost value of global alignment. */
GtWord gt_diagonalbandalign_affinegapcost_square_space_distance_only(
                                                           const GtUchar *useq,
                                                           GtUword ustart,
                                                           GtUword ulen,
                                                           const GtUchar *vseq,
                                                           GtUword vstart,
                                                           GtUword vlen,
                                                           GtWord left_dist,
                                                           GtWord right_dist,
                                                           const GtScoreHandler
                                                           *scorehandler);

/* Computes a global alignment within a diagonal band with linear gapcosts in
   square space. Use of this function requires an initialised <scorehandler>
   with cost values, the target alignment <align> and input sequences <useq> and
   <vseq>, with the regions to align given by their start positions <ustart> and
   <vstart> and lengths <ulen> and <vlen>. <left_dist> and <right_dist> give
   lower and upper bound of a diagonal band in which DP-matrix is valid.
   <spacemanager> is required to use this function in linear space context, in
   any other case it can be NULL. <scorehandler> manages linear gap costs.
   Returns cost value of global alignment. */
GtWord gt_diagonalbandalign_affinegapcost_in_square_space_generic(
                                             GtLinspaceManagement *space,
                                             const GtScoreHandler *scorehandler,
                                             GtAlignment *align,
                                             const GtUchar *useq,
                                             GtUword ustart,
                                             GtUword ulen,
                                             const GtUchar *vseq,
                                             GtUword vstart,
                                             GtUword vlen,
                                             GtWord left_dist,
                                             GtWord right_dist);

/* Computes a global alignment within a diagonal band with linear gapcosts in
   square space and constant cost values. Use of this function requires the
   target alignment <align> and input sequences <useq> and <vseq>, with the
   regions to align given by their start positions <ustart> and <vstart> and
   lengths <ulen> and <vlen>. <left_dist> and <right_dist> give lower and upper
   bound of a diagonal band in which DP-matrix is valid. The cost values
   are specified by <matchcost>, <mismatchcost> and <gapcost>.
   <spacemanager> is required to use this function in linear space context, in
   any other case it can be NULL. <scorehandler> manages linear gap costs.
   Returns cost value of global alignment. */
GtWord gt_diagonalbandalign_affinegapcost_in_square_space(GtLinspaceManagement
                                                        *space,
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

void   gt_diagonalbandalign_affinegapcost_check(GT_UNUSED bool forward,
                                                const GtUchar *useq,
                                                GtUword ulen,
                                                const GtUchar *vseq,
                                                GtUword vlen);
#endif
