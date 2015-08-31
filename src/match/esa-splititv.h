/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef ESA_SPLITITV_H
#define ESA_SPLITITV_H

#include "core/arraydef.h"
#include "core/types_api.h"

#include "core/encseq.h"
#include "match/sarr-def.h"
#include "splititv.h"

typedef struct
{
  GtUword left,
         right;
} Simplelcpinterval;

bool gt_lcpintervalfindcharchildintv(const GtEncseq *encseq,
                                  GtReadmode readmode,
                                  GtUword totallength,
                                  const ESASuffixptr *suftab,
                                  Simplelcpinterval *itv,
                                  GtUchar cc,
                                  GtUword offset,
                                  GtUword left,
                                  GtUword right);

void gt_lcpintervalsplitwithoutspecial(GtArrayBoundswithchar *bwci,
                                    const GtEncseq *encseq,
                                    GtReadmode readmode,
                                    GtUword totallength,
                                    const ESASuffixptr *suftab,
                                    GtUword parentoffset,
                                    GtUword parentleft,
                                    GtUword parentright);

GtUchar gt_lcpintervalextendlcp(const GtEncseq *encseq,
                           GtReadmode readmode,
                           const ESASuffixptr *suftab,
                           GtUword totallength,
                           GtUchar alphasize,
                           GtUword parentoffset,
                           GtUword parentleft,
                           GtUword parentright);

GtUword gt_findmaximalprefixinESA(const GtEncseq *encseq,
                                  GtReadmode readmode,
                                  GtUword totallength,
                                  const ESASuffixptr *suftab,
                                  const GtUchar *query,
                                  GtUword querylen);

#endif
