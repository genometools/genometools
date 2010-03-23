/*
  Copyright (c) 2005-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef ARRAY_H
#define ARRAY_H

#include "core/array_api.h"
#include "core/error.h"

typedef int (*GtArrayProcessor)(void *elem, void *info, GtError*);

/* Compare the raw content of <array_a> with the content of <array_b>.
   <array_a> and <array_b> must have the same gt_array_size() and
   gt_array_elem_size(). */
int           gt_array_cmp(const GtArray  *array_a, const GtArray *array_b);
/* Compare the content of <array_a> with the content of <array_b> using the
   comparator function <cmpfunc>. If the elements of both arrays are equal
   w.r.t. <cmpfunc>, true is returned. If the array sizes or content w.r.t.
   <cmpfunc> are different, false is returned. */
bool          gt_array_equal(const GtArray *a, const GtArray *b,
                             GtCompare cmpfunc);
/* Iterate over all elements in <array> and call <array_processor> with them.
   <info> and <err> are passed to <array_processor>.
   If <array_processor> returns a value != 0, the iteration is stopped and the
   return value of <array_processor> is returned. */
int           gt_array_iterate(GtArray *array,
                               GtArrayProcessor array_processor,
                               void *info, GtError *err);
/* Similar to <array_iterate>, except that the <array> is traversed in reverse
   order. */
int           gt_array_iterate_reverse(GtArray *array,
                                       GtArrayProcessor array_processor,
                                       void *info, GtError *err);
void          gt_array_prepend_array(GtArray *dest, const GtArray *src);
int           gt_array_example(GtError*);
int           gt_array_unit_test(GtError*);

#endif
