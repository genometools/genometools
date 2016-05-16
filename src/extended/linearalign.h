/*
  Copyright (c) 2015 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
  Copyright (C) 2015 Joerg Winkler <j.winkler@posteo.de>
  Copyright (c) 2006-2007 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2015 Center for Bioinformatics, University of Hamburg

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

#ifndef LINEARALIGN_H
#define LINEARALIGN_H

#include "core/unused_api.h"
#include "core/error.h"
#include "extended/linspace_management.h"
#include "extended/scorehandler.h"

/* Computes a global alignment with linear gapcosts in linear space. Use of this
   function requires an initialised <spacemanager>, the target alignment <align>
   and input sequences <useq> and <vseq>, with the regions to align given by
   their start positions <ustart> and <vstart> and lengths <ulen> and <vlen>.
   The cost values are specified by an initialised <scorehandler>. Returns
   distance value of calculated global alignment. */
GtUword gt_linearalign_compute_generic(GtLinspaceManagement *spacemanager,
                                       const GtScoreHandler *scorehandler,
                                       GtAlignment *align,
                                       const GtUchar *useq,
                                       GtUword ustart,
                                       GtUword ulen,
                                       const GtUchar *vseq,
                                       GtUword vstart,
                                       GtUword vlen);

/* Computes a global alignment with linear gapcosts in linear space
   and constant cost values. Use of this function requires an initialised
   <spacemanager>, the target alignment <align> and input sequences <useq> and
   <vseq>, with the regions to align given by their start positions <ustart> and
   <vstart> and lengths <ulen> and <vlen>. The cost values are specified by
   <matchcost>, <mismatchcost> and <gapcost>. Returns distance value of
   calculated global alignment. */
GtUword gt_linearalign_compute(GtLinspaceManagement *spacemanager,
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

/* Computes a local alignment with linear gapcosts in linear space. Use of this
   function requires an initialised <spacemanager>, the target alignment <align>
   and input sequences <useq> and <vseq>, with the regions to align given by
   their start positions <ustart> and <vstart> and lengths <ulen> and <vlen>.
   The score values are specified by an initialised <scorehandler>. Returns
   score value of calculated local alignment. */
GtWord  gt_linearalign_compute_local_generic(GtLinspaceManagement *spacemanager,
                                             const GtScoreHandler *scorehandler,
                                             GtAlignment *align,
                                             const GtUchar *useq,
                                             GtUword ustart,
                                             GtUword ulen,
                                             const GtUchar *vseq,
                                             GtUword vstart,
                                             GtUword vlen);

/* Computes a local alignment with linear gapcosts in linear space
   and constant score values. Use of this function requires an initialised
   <spacemanager>, the target alignment <align> and input sequences <useq> and
   <vseq>, with the regions to align given by their start positions <ustart> and
   <vstart> and lengths <ulen> and <vlen>. The score values are specified by
   <matchscore>, <mismatchscore> and <gapscore>. Returns score value of
   calculated local alignment. */
GtWord  gt_linearalign_compute_local(GtLinspaceManagement *spacemanager,
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

void    gt_linearalign_check(GT_UNUSED bool forward,
                             const GtUchar *useq,
                             GtUword ulen,
                             const GtUchar *vseq,
                             GtUword vlen);

void    gt_linearalign_check_local(GT_UNUSED bool forward,
                                   const GtUchar *useq,
                                   GtUword ulen,
                                   const GtUchar *vseq,
                                   GtUword vlen);
#endif
