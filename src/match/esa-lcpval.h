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

#ifndef ESA_LCPVAL_H
#define ESA_LCPVAL_H

#include "core/readmode.h"
#include "core/encseq.h"
#include "match/sarr-def.h"

typedef struct Lcpvalueiterator Lcpvalueiterator;

Lcpvalueiterator *gt_newLcpvalueiterator(const GtEncseq *encseq,
                                      GtReadmode readmode);

unsigned long gt_nextLcpvalueiterator(Lcpvalueiterator *lvi,
                                      bool firstpage,
                                      const ESASuffixptr *suftabptr,
                                      unsigned long numberofsuffixes);

void gt_freeLcpvalueiterator(Lcpvalueiterator *lvi);

#endif
