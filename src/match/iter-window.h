/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef ITER_WINDOW_H
#define ITER_WINDOW_H

#include "core/encseq.h"

typedef struct Windowiterator Windowiterator;

Windowiterator *gt_windowiterator_new(const GtEncseq *encseq,
                                   unsigned long windowsize,
                                   unsigned long startpos,
                                   unsigned long endpos);

void gt_windowiterator_delete(Windowiterator *wit);

const GtUchar *gt_windowiterator_next(unsigned long *currentpos,
                                   unsigned long *firstpos,
                                   Windowiterator *wit);

#endif
