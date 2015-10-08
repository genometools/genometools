/*
  Copyright (c) 2015 Annika <annika.seidel@studium.uni-hamburg.de>
  Copyright (C) 2015 Joerg Winkler, joerg.winkler@studium.uni-hamburg.de
  Copyright (c) 2006-2007 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include "extended/linspaceManagement.h"
#include "extended/scorehandler.h"

void gt_checklinearspace(GT_UNUSED bool forward,
                         const GtUchar *useq,
                         GtUword ulen,
                         const GtUchar *vseq,
                         GtUword vlen);

void gt_checklinearspace_local(GT_UNUSED bool forward,
                               const GtUchar *useq,
                               GtUword ulen,
                               const GtUchar *vseq,
                               GtUword vlen);

/* global alignment with linear gapcosts in linear space (DNA or protein)*/
GtUword gt_computelinearspace_generic(LinspaceManagement *spacemanager,
                                      GtScoreHandler *scorehandler,
                                      GtAlignment *align,
                                      const GtUchar *useq,
                                      GtUword ustart,
                                      GtUword ulen,
                                      const GtUchar *vseq,
                                      GtUword vstart,
                                      GtUword vlen);

/* global alignment with linear gapcosts in linear space
 * with constant cost values, only useful for DNA sequences */
GtUword gt_computelinearspace(LinspaceManagement *spacemanager,
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

/* local alignment with linear gapcosts in linear space (DNA or protein)*/
GtWord gt_computelinearspace_local_generic(LinspaceManagement *spacemanager,
                                           GtScoreHandler *scorehandler,
                                           GtAlignment *align,
                                           const GtUchar *useq,
                                           GtUword ustart,
                                           GtUword ulen,
                                           const GtUchar *vseq,
                                           GtUword vstart,
                                           GtUword vlen);

/* local alignment with linear gapcosts in linear space
 * with constant score values, only useful for DNA sequences */
GtWord gt_computelinearspace_local(LinspaceManagement *spacemanager,
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

/* edit distance of sequences useq and vseq in linear space*/
GtUword gt_calc_linearedist(const GtUchar *useq, GtUword ulen,
                            const GtUchar *vseq, GtUword vlen);

#endif
