/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#ifndef LAYOUT_H
#define LAYOUT_H

#include "annotationsketch/block_api.h"
#include "annotationsketch/layout_api.h"

typedef int (*GtBlockOrderingFunc)(const GtBlock *b1, const GtBlock *b2,
                                   void *data);

/* Sets the ordering function in the layout which determines in what order
   the blocks (per track) are inserted into lines. The default is to
   sort using <gt_block_compare()>. */
void                   gt_layout_set_block_ordering_func(GtLayout*,
                                                         GtBlockOrderingFunc,
                                                         void *data);
/* Unsets any block ordering function, processing the blocks in the order they
   were inserted. */
void                   gt_layout_unset_block_ordering_func(GtLayout*);
/* Returns the interval shown in the layout. */
GtRange                gt_layout_get_range(const GtLayout*);
/* Returns the TextWidthCalculator object used in the layout. */
GtTextWidthCalculator* gt_layout_get_twc(const GtLayout*);

#endif
