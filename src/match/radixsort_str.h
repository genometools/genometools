/*
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#ifndef RADIXSORT_STR_H
#define RADIXSORT_STR_H
#include "core/intbits.h"
#include "sfx-lcpvalues.h"

/*
 * suffixes = array of suffixes to be sorted
 * width = number of suffixes to sort
 * depth = start sorting at given depth
 * sortmaxdepth = stop sorting at given depth (0 means infinite)
 *
 * realtotallength = "real" length (not the mirror logical length)
 * equallengthplus1 = sequence length, including the separator
 *
 * */

typedef struct GtRadixsortstringinfo GtRadixsortstringinfo;

unsigned long gt_radixsort_str_minwidth(void);

unsigned long gt_radixsort_str_maxwidth(const GtRadixsortstringinfo *rsi);

GtRadixsortstringinfo *gt_radixsort_str_new(const GtTwobitencoding
                                             *twobitencoding,
                                            unsigned long realtotallength,
                                            unsigned long equallengthplus1,
                                            unsigned long maxwidth);

void gt_radixsort_str_delete(GtRadixsortstringinfo *rsi);

void gt_radixsort_str_eqlen(GtRadixsortstringinfo *rsi,
                            unsigned long *suffixes,
                            GtLcpvalues *lcpvalues,
                            unsigned long subbucketleft,
                            unsigned long depth,
                            unsigned long sortmaxdepth,
                            unsigned long width);

#endif
