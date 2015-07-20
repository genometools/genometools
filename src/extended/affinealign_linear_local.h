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

#ifndef AFFINEALIGN_LINEAR_LOCAL_H
#define AFFINEALIGN_LINEAR_LOCAL_H

#include "core/types_api.h"

void gt_computeaffinelinearspace_local(const GtUchar *useq,
                                       const GtUword ustart,
                                       const GtUword ulen,
                                       const GtUchar *vseq,
                                       const GtUword vstart,
                                       const GtUword vlen,
                                       const GtWord matchscore,
                                       const GtWord mismatchscore,
                                       const GtWord gap_opening,
                                       const GtWord gap_extension,
                                       FILE *fp);

void gt_checkaffinelinearspace_local(GT_UNUSED bool forward,
                               const GtUchar *useq,
                               GtUword ulen,
                               const GtUchar *vseq,
                               GtUword vlen);
#endif
