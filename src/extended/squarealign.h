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

#ifndef SQUAREALIGN_H
#define SQUAREALIGN_H
#include "core/types_api.h"
#include "extended/alignment.h"
#include "extended/linspaceManagement.h"
#include "extended/scorehandler.h"

/* create an global alignment in square space, to use it in linear context you
 * have to generate an spacemanager before, in any other case it can be NULL,
 * (DNA or protein) */
GtUword alignment_in_square_space_generic (LinspaceManagement *spacemanager,
                                           GtAlignment *align,
                                           const GtUchar *useq,
                                           GtUword ustart,
                                           GtUword ulen,
                                           const GtUchar *vseq,
                                           GtUword vstart,
                                           GtUword vlen,
                                           const GtScoreHandler *scorehandler);

/* same with constant cost values, only useful for DNA sequences */
GtUword alignment_in_square_space(LinspaceManagement *spacemanager,
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

/*distance only for global alignment */
GtUword distance_only_global_alignment(const GtUchar *useq,
                                       GtUword ustart,
                                       GtUword ulen,
                                       const GtUchar *vseq,
                                       GtUword vstart,
                                       GtUword vlen,
                                       const GtScoreHandler *scorehandler);

void gt_print_edist_alignment(const GtUchar *useq, GtUword ustart, GtUword ulen,
                             const GtUchar *vseq, GtUword vstart, GtUword vlen);

/* fill crosspointtable ctab for part of sequences useq and vseq in square
 * space, use it to combine square calculating with linear calculating */
GtUword ctab_in_square_space(LinspaceManagement *spacemanager,
                             const GtScoreHandler *scorehandler,
                             GtUword *Ctab,
                             const GtUchar *useq,
                             GtUword ustart,
                             GtUword ulen,
                             const GtUchar *vseq,
                             GtUword vstart,
                             GtUword vlen,
                             GtUword rowoffset);

/* create an local alignment in square space, to use it in linear context you
 * have to generate an spacemanager before, in any other case it can be NULL,
 * (DNA or protein) */
GtWord alignment_in_square_space_local_generic(LinspaceManagement *spacemanager,
                                               GtAlignment *align,
                                               const GtUchar *useq,
                                               GtUword ustart,
                                               GtUword ulen,
                                               const GtUchar *vseq,
                                               GtUword vstart,
                                               GtUword vlen,
                                            const GtScoreHandler *scorehandler);

/* same with constant cost values, only useful for DNA sequences */
GtWord alignment_in_square_space_local(LinspaceManagement *spacemanager,
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
