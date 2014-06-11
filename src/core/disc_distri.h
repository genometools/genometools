/*
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

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

#ifndef DISC_DISTRI_H
#define DISC_DISTRI_H

#include "core/error_api.h"
#include "core/file_api.h"
#include "core/types_api.h"

#include "core/disc_distri_api.h"

/* Callback function called during iteration over each item of the
   distribution, where <key> is the counted value and <value> is the count.
   Has to return 0 on success, non zero on failure. On failure <err> must be
   set!
 */
typedef int (*GtDiscDistriIterFuncErr)(GtUword key, GtUint64 value, void *data,
                                       GtError *err);

/* Iterate over all non-empty entries in <d>, calling <func> for each one,
   from the smallest to the largest key. The <data> pointer can be used to pass
   arbitrary data to <func>.
   Returns 0 on success, or the error value returned by <func> which will set
   <err> on failure.
 */
int gt_disc_distri_foreach_err(const GtDiscDistri *d,
                               GtDiscDistriIterFuncErr func,
                               void *data, GtError *err);

/* Same as foreach, but from the longest to the smallest key. */
int gt_disc_distri_foreach_in_reverse_order_err(const GtDiscDistri *d,
                                                GtDiscDistriIterFuncErr func,
                                                void *data, GtError *err);
#endif
