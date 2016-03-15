/*
  Copyright (C) 2015 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
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

#include "extended/alignment.h"
#include "extended/linspace_management.h"
#include "extended/scorehandler.h"

typedef enum {
  Affine_R,
  Affine_D,
  Affine_I,
  Affine_X /* unknown */
} GtAffineAlignEdge;

/* <AffinealignDPentry> objects describe the information of distance values and
   backtracing edges relating on last edit operation R,D,I. */
typedef struct {
  GtWord Rvalue, Dvalue, Ivalue, totalvalue;
  GtAffineAlignEdge Redge,
                    Dedge,
                    Iedge;
} GtAffinealignDPentry;

/* Computes a global alignment with affine gapcosts in square space
   and constant cost values. Use of this function requires input sequences
   <useq> and <vseq> and lengths <ulen> and <vlen>. The cost values are
   specified by <matchcost>, <mismatchcost>, <gap_opening_cost> and
   <gap_extension_cost>. Returns an object of the <GtAlignment> class. */
GtAlignment* gt_affinealign(const GtUchar *u, GtUword ulen,
                            const GtUchar *v, GtUword vlen,
                            GtUword matchcost, GtUword mismatchcost,
                            GtUword gap_opening_cost,
                            GtUword gap_extension_cost);

/* Computes a global alignment with affine gapcosts in square space. Use of this
   function requires an initialised <scorehandler> with cost values an
   initialised <spacemanager>, the target alignment <align> and input sequences
   <useq> and <vseq> and lengths <ulen> and <vlen>. Returns affine cost
   value of global alignment. */
GtWord       gt_affinealign_with_Management(GtLinspaceManagement *spacemanager,
                                            const GtScoreHandler *scorehandler,
                                            GtAlignment *align,
                                            const GtUchar *u, GtUword ulen,
                                            const GtUchar *v, GtUword vlen);

GtWord       gt_affinealign_traceback(GtAlignment *align,
                                      GtAffinealignDPentry * const *dptable,
                                      GtUword i, GtUword j);

/* Computes crosspoints for a global alignment with affine gapcosts in square
   space. Use of this function requires an initialised <spacemanager> an
   initialised <scorehandler> with cost values, the target crosspoint table
   <Ctab> and input sequences <useq> and <vseq>, with the regions to align given
   by their start positions <ustart> and <vstart> and lengths <ulen> and <vlen>.
   If this function is used in linear context, <rowoffset> is the offset value
   of the subproblem and <from_edge> and <to_edge> are the in- and outcoming
   edge of this subproblem. Otherwise set default values 0 and Affine_X. Returns
   affine distance value of global alignment. */
void         gt_affinealign_ctab(GtLinspaceManagement *spacemanager,
                                 const GtScoreHandler *scorehandler,
                                 GtUword *Ctab,
                                 const GtUchar *useq,
                                 GtUword ustart,
                                 GtUword ulen,
                                 const GtUchar *vseq,
                                 GtUword vstart,
                                 GtUword vlen,
                                 GtUword rowoffset,
                                 GtAffineAlignEdge from_edge,
                                 GtAffineAlignEdge to_edge);

/* Computes a local alignment with linear gapcosts in square space. Use of this
   function requires an initialised <scorehandler> with score values, the target
   alignment <align> and input sequences <useq> and <vseq>, with the regions to
   align given by their start positions <ustart> and <vstart> and lengths <ulen>
   and <vlen>. An initialised <spacemanager> is required to use this function in
   linear space context, in any other case it can be NULL. Returns score
   value of local alignment. */
GtWord       gt_affinealign_calculate_local_generic(GtLinspaceManagement
                                                    *spacemanager,
                                                    const GtScoreHandler
                                                    *scorehandler,
                                                    GtAlignment *align,
                                                    const GtUchar *useq,
                                                    GtUword ustart,
                                                    GtUword ulen,
                                                    const GtUchar *vseq,
                                                    GtUword vstart,
                                                    GtUword vlen);

/* Computes a local alignment with affine gapcosts in square space
   and constant score values. Use of this function requires the target alignment
   <align> and input sequences <useq> and <vseq>, with the regions to align
   given by their start positions <ustart> and <vstart> and lengths <ulen> and
   <vlen>. The score values are specified by <matchscore>, <mismatchscore>,
   <gap_opening> and <gap_extension>. An initialised <spacemanager> is required
   to use this function in linear space context, in any other case it can be
   NULL. Returns affine score value of local alignment. */
GtWord       gt_affinealign_calculate_local(GtLinspaceManagement
                                            *spacemanager,
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
