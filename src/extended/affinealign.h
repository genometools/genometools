/*
  Copyright (C) 2015 Annika Seidel, annika.seidel@studium.uni-hamburg.de
  Copyright (c) 2007 Gordon Gremme <gordon@gremme.org>
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

#ifndef AFFINEALIGN_H
#define AFFINEALIGN_H

#include "extended/linspaceManagement.h"
#include "extended/scorehandler.h"

#include "extended/alignment.h"

typedef enum {
  Affine_R,
  Affine_D,
  Affine_I,
  Affine_X /* unknown */
} AffineAlignEdge;

typedef struct {
  GtWord Rvalue, Dvalue, Ivalue, totalvalue;
  AffineAlignEdge Redge,
                  Dedge,
                  Iedge;
} AffinealignDPentry;

/* (globally) align u and v (affine gap costs) and return one optimal
   GtAlignment */
GtAlignment* gt_affinealign(const GtUchar *u, GtUword ulen,
                            const GtUchar *v, GtUword vlen,
                            GtUword matchcost, GtUword mismatchcost,
                            GtUword gap_opening_cost,
                            GtUword gap_extension_cost);

/* globall alignment (DNA or protein)*/
GtWord gt_affinealign_with_Management(LinspaceManagement *spacemanager,
                                    GtScoreHandler *scorehandler,
                                    GtAlignment *align,
                                    const GtUchar *u, GtUword ulen,
                                    const GtUchar *v, GtUword vlen);

GtWord affinealign_traceback(GtAlignment *a,
                           AffinealignDPentry * const *dptable,
                           GtUword i, GtUword j);

/* filling ctab to combine square calculating with linear calculating */
void affine_ctab_in_square_space(LinspaceManagement *spacemanager,
                                 GtScoreHandler *scorehandler,
                                 GtUword *Ctab,
                                 const GtUchar *useq,
                                 GtUword ustart,
                                 GtUword ulen,
                                 const GtUchar *vseq,
                                 GtUword vstart,
                                 GtUword vlen,
                                 GtUword rowoffset,
                                 AffineAlignEdge from_edge,
                                 AffineAlignEdge to_edge);

/* create an local alignment in square space, to use it in linear context you
 * have to generate an spacemanager before, in any other case it can be NULL,
 * (DNA or protein) */
GtWord affinealign_in_square_space_local_generic(LinspaceManagement *space,
                                                 GtScoreHandler *scorehandler,
                                                 GtAlignment *align,
                                                 const GtUchar *useq,
                                                 GtUword ustart,
                                                 GtUword ulen,
                                                 const GtUchar *vseq,
                                                 GtUword vstart,
                                                 GtUword vlen);

/* same with constant score values, use it only for DNA sequences! */
GtWord affinealign_in_square_space_local(LinspaceManagement *spacemanager,
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

#endif
