/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/countingsort.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/unused_api.h"

void gt_countingsort(void *out, const void *in, size_t elem_size,
                     unsigned long size, unsigned long max_elemvalue,
                     void *data, GtGetElemvalue get_elemvalue)
{
  unsigned long i, k, *c;
  gt_assert(out && in && elem_size && size && max_elemvalue && get_elemvalue);

  /* allocate count array */
  c = gt_calloc(sizeof (unsigned long), max_elemvalue + 1);

  /* count number of elements of a given value */
  for (i = 0; i < size; i++) {
    k = get_elemvalue((const char*) in + elem_size * i, data);
    gt_assert(k <= max_elemvalue);
    c[k]++;
  }

  /* compute running sum of the count array */
  for (i = 1; i <= max_elemvalue; i++)
    c[i] += c[i-1];

  /* sorting (stable) */
  for (i = size; i > 0; i--) {
    k = get_elemvalue((const char*) in + elem_size * (i-1), data);
    memcpy((char*) out + elem_size * (c[k] - 1),
           (const char*) in + elem_size * (i-1), elem_size);
    c[k]--;
  }

  gt_free(c);
}

unsigned long gt_countingsort_get_max(const void *in, size_t elem_size,
                                      unsigned long size, void *data,
                                      GtGetElemvalue get_elemvalue)
{
  unsigned long i, value, max_value = 0;
  for (i = 0; i < size; i++) {
    value = get_elemvalue((const char*) in + elem_size * i, data);
    if (value > max_value)
      max_value = value;
  }
  return max_value;
}

static unsigned long get_int(const void *elem, GT_UNUSED void *data)
{
  gt_assert(elem);
  return *(unsigned int*) elem;
}

int gt_countingsort_unit_test(GtError *err)
{
  unsigned int numbers[]        = { 1, 2, 1, 2, 0 }, numbers_out[5],
               sorted_numbers[] = { 0, 1, 1, 2, 2 };
  int had_err = 0;
  gt_error_check(err);
  gt_countingsort(numbers_out, numbers, sizeof (unsigned int), 5,
                  gt_countingsort_get_max(numbers, sizeof (unsigned int), 5,
                                          NULL, get_int),
                  NULL,  get_int);
  gt_ensure(had_err,
         !memcmp(sorted_numbers, numbers_out, sizeof (unsigned int) * 5));
  return had_err;
}
