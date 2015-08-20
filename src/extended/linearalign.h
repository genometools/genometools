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
#include "extended/alignment.h"

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

/*global alignment with linear gapcosts*/
GtAlignment *gt_computelinearspace( const GtUchar *useq,
                                     const GtUword ustart,
                                     const GtUword ulen,
                                     const GtUchar *vseq,
                                     const GtUword vstart,
                                     const GtUword vlen,
                                     const GtWord matchcost,
                                     const GtWord mismatchcost,
                                     const GtWord gapcost);

/*local alignment with linear gapcosts*/
GtAlignment * gt_computelinearspace_local(const GtUchar *useq,
                                          const GtUword ustart,
                                          const GtUword ulen,
                                          const GtUchar *vseq,
                                          const GtUword vstart,
                                          const GtUword vlen,
                                          const GtWord matchscore,
                                          const GtWord mismatchscore,
                                          const GtWord gapscore);

/*edit distance of sequences u and v*/
GtUword gt_calc_linearedist(const GtUchar *u, GtUword ulen,
                            const GtUchar *v, GtUword vlen);

void gt_print_edist_alignment(const GtUchar *useq, GtUword ustart,
                           GtUword ulen,
                           const GtUchar *vseq,GtUword vstart,
                           GtUword vlen);
#endif
