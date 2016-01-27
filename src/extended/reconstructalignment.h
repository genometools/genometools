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

#ifndef RECONSTRUCTALIGNMENT_H
#define RECONSTRUCTALIGNMENT_H
#include "core/types_api.h"
#include "extended/affinealign.h"
#include "extended/alignment.h"
#include "extended/diagonalbandalign.h"
#include "extended/diagonalbandalign_affinegapcost.h"
#include "extended/maxcoordvalue.h"
#include "extended/scorehandler.h"

/* Add an alignment of length <len> to <align> containing only deletion
   operations. Returns linear cost value of constructed alignment using
   specified <gapcost>. */
GtUword gt_reconstructalignment_trivial_deletion(GtAlignment *align,
                                                 GtUword len,
                                                 GtUword gapcost);

/* Add an alignment of length <len> to <align> containing only insertion
   operations. Returns linear cost value of constructed alignment using
   specified <gapcost>. */
GtUword gt_reconstructalignment_trivial_insertion(GtAlignment *align,
                                                  GtUword len,
                                                  GtUword gapcost);

/* Reconstruct an object <align> from square space table <E>, which describes
   a global alignment between two sequences  <useq> and <vseq>, with the
   regions to align given by their start positions <ustart> and <vstart> and
   lengths <ulen> and <vlen>. <scorehandler> specifies the cost values with
   which the <E>-matrix was created. */
void    gt_reconstructalignment_from_EDtab(GtAlignment *align,
                                           GtUword * const *E,
                                           const GtUchar *useq,
                                           GtUword ustart,
                                           GtUword ulen,
                                           const GtUchar *vseq,
                                           GtUword vstart,
                                           GtUword vlen,
                                           const GtScoreHandler *scorehandler);

/* Reconstruct an object <align> from square space table <Ltabcolumn>, which
   describes a local alignment between two sequences  <useq> and <vseq>, with
   the regions to align given by their start positions <ustart> and <vstart> and
   lengths <ulen> and <vlen>. <scorehandler> specifies the score values with
   which the <Ltabcolumn>-matrix was created. */
void    gt_reconstructalignment_from_Ltab(GtAlignment *align,
                                          GtWord **Ltabcolumn,
                                          GtMaxcoordvalue *max,
                                          const GtUchar *useq,
                                          GtUword ustart,
                                          GtUword ulen,
                                          const GtUchar *vseq,
                                          GtUword vstart,
                                          GtUword vlen,
                                          const GtScoreHandler *scorehandler);

/* Reconstruct an object <align> from crosspoints in <Ctab> between two
   sequences  <useq> and <vseq>, with the regions to align given by their start
   positions <ustart> and <vstart> and lengths <ulen> and <vlen>. The
   crosspoints are located on midcolumns. <scorehandler> specifies the cost
   values with which the <Ctab> was created. */
void    gt_reconstructalignment_from_Ctab(GtAlignment *align,
                                          const GtUword *Ctab,
                                          const GtUchar *useq,
                                          GtUword ustart,
                                          const GtUchar *vseq,
                                          GtUword vstart,
                                          GtUword vlen,
                                          const GtScoreHandler *scorehandler);

/* Reconstruct an object <align> from crosspoints in <Dtab> between two
   sequences  <useq> and <vseq>, with the regions to align given by their start
   positions <ustart> and <vstart> and lengths <ulen> and <vlen>. The
   crosspoints are located on diagonals and the evaluation is based on
   linear gap costs. */
void    gt_reconstructalignment_from_Dtab(GtAlignment *align,
                                          const GtDiagAlignentry *Dtab,
                                          GtUword ulen,
                                          GtUword vlen);

/* Reconstruct an object <align> from crosspoints in <Dtab> between two
   sequences  <useq> and <vseq>, with the regions to align given by their start
   positions <ustart> and <vstart> and lengths <ulen> and <vlen>. The
   crosspoints are located on diagonals and the and the evaluation is based on
   affine gap costs. */
void     gt_reconstructalignment_from_affineDtab(GtAlignment *align,
                                             const GtAffineDiagAlignentry *Dtab,
                                             GtAffineAlignEdge edge,
                                             const GtUchar *useq,
                                             GtUword ulen,
                                             const GtUchar *vseq,
                                             GtUword vlen);
#endif
