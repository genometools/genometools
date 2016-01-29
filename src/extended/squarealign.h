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

#ifndef SQUAREALIGN_H
#define SQUAREALIGN_H

#include "core/types_api.h"
#include "extended/alignment.h"
#include "extended/linspace_management.h"
#include "extended/scorehandler.h"

/* Computes a global alignment with linear gapcosts in square space. Use of this
   function requires an initialised <scorehandler> with cost values, the target
   alignment <align> and input sequences <useq> and <vseq>, with the regions to
   align given by their start positions <ustart> and <vstart> and lengths <ulen>
   and <vlen>. An initialised <spacemanager> is required to use this function in
   linear space context, in any other case it can be NULL. Returns distance
   value of global alignment. */
GtUword gt_squarealign_calculate_generic(GtLinspaceManagement *spacemanager,
                                         GtAlignment *align,
                                         const GtUchar *useq,
                                         GtUword ustart,
                                         GtUword ulen,
                                         const GtUchar *vseq,
                                         GtUword vstart,
                                         GtUword vlen,
                                         const GtScoreHandler *scorehandler);

/* Computes a global alignment with linear gapcosts in square space
   and constant cost values. Use of this function requires the target alignment
   <align> and input sequences <useq> and <vseq>, with the regions to align
   given by their start positions <ustart> and <vstart> and lengths <ulen> and
   <vlen>. The cost values are specified by <matchcost>, <mismatchcost> and
   <gapcost>. An initialised <spacemanager> is required to use this function in
   linear space context, in any other case it can be NULL. Returns distance
   value of global alignment. */
GtUword gt_squarealign_calculate(GtLinspaceManagement *spacemanager,
                                 GtAlignment *align,
                                 const GtUchar *useq,
                                 GtUword ustart,
                                 GtUword ulen,
                                 const GtUchar *vseq,
                                 GtUword vstart,
                                 GtUword vlen,
                                 GtUword matchcost,
                                 GtUword mismatchcost,
                                 GtUword gapcost);

/* Computes only the distance of a global alignment with linear gapcosts in
   square space. Use of this function requires an initialised <scorehandler>
   with cost values and input sequences <useq> and <vseq>, with the regions to
   align given by their start positions <ustart> and <vstart> and lengths <ulen>
   and <vlen>. */
GtUword gt_squarealign_global_distance_only(const GtUchar *useq,
                                            GtUword ustart,
                                            GtUword ulen,
                                            const GtUchar *vseq,
                                            GtUword vstart,
                                            GtUword vlen,
                                            const GtScoreHandler *scorehandler);

/* Shows a global alignment, which is computed with linear gapcosts in
   square space. Use of this function requires input sequences <useq> and
   <vseq>, with the regions to align given by their start positions <ustart> and
   <vstart> and lengths <ulen> and <vlen>. */
void    gt_squarealign_print_edit_alignment(const GtUchar *useq, GtUword ustart,
                                            GtUword ulen, const GtUchar *vseq,
                                            GtUword vstart, GtUword vlen);

/* Computes crosspoints for a global alignment with linear gapcosts in square
   space. Use of this function requires an initialised <spacemanager> an
   initialised <scorehandler> with cost values, the target crosspoint table
   <Ctab> and input sequences <useq> and <vseq>, with the regions to align given
   by their start positions <ustart> and <vstart> and lengths <ulen> and <vlen>.
   If this function is use in linear context, <rowoffset> is the offset value of
   the subproblem, else zero. Returns distance value of global alignment. */
GtUword gt_squarealign_ctab(GtLinspaceManagement *spacemanager,
                            const GtScoreHandler *scorehandler,
                            GtUword *Ctab,
                            const GtUchar *useq,
                            GtUword ustart,
                            GtUword ulen,
                            const GtUchar *vseq,
                            GtUword vstart,
                            GtUword vlen,
                            GtUword rowoffset);

/* Computes a local alignment with linear gapcosts in square space. Use of this
   function requires an initialised <scorehandler> with score values, the target
   alignment <align> and input sequences <useq> and <vseq>, with the regions to
   align given by their start positions <ustart> and <vstart> and lengths <ulen>
   and <vlen>. An initialised <spacemanager> is required to use this function in
   linear space context, in any other case it can be NULL. Returns score
   value of local alignment. */
GtWord  gt_squarealign_calculate_local_generic(GtLinspaceManagement
                                               *spacemanager,
                                               GtAlignment *align,
                                               const GtUchar *useq,
                                               GtUword ustart,
                                               GtUword ulen,
                                               const GtUchar *vseq,
                                               GtUword vstart,
                                               GtUword vlen,
                                               const GtScoreHandler
                                               *scorehandler);

/* Computes a local alignment with linear gapcosts in square space
   and constant score values. Use of this function requires the target alignment
   <align> and input sequences <useq> and <vseq>, with the regions to align
   given by their start positions <ustart> and <vstart> and lengths <ulen> and
   <vlen>. The score values are specified by <matchscore>, <mismatchscore> and
   <gapscore>. An initialised <spacemanager> is required to use this function in
   linear space context, in any other case it can be NULL. Returns score
   value of local alignment. */
GtWord  gt_squarealign_calculate_local(GtLinspaceManagement *spacemanager,
                                       GtAlignment *align,
                                       const GtUchar *useq,
                                       GtUword ustart,
                                       GtUword ulen,
                                       const GtUchar *vseq,
                                       GtUword vstart,
                                       GtUword vlen,
                                       GtWord matchscore,
                                       GtWord mismatchscore,
                                       GtWord gapscore);
#endif
