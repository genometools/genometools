/*
  Copyright (c) 2005-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <assert.h>
#include "libgtcore/ensure.h"
#include "libgtcore/fptr.h"
#include "libgtcore/minmax.h"
#include "libgtcore/msort.h"
#include "libgtcore/range.h"
#include "libgtcore/safearith.h"
#include "libgtcore/undef.h"

int range_compare(Range range_a, Range range_b)
{
  assert(range_a.start <= range_a.end && range_b.start <= range_b.end);

  if ((range_a.start == range_b.start) && (range_a.end == range_b.end))
    return 0; /* range_a == range_b */

  if ((range_a.start < range_b.start) ||
      ((range_a.start == range_b.start) && (range_a.end < range_b.end)))
    return -1; /* range_a < range_b */

  return 1; /* range_a > range_b */
}

int range_compare_ptr(const Range *range_a, const Range *range_b)
{
  return range_compare(*range_a, *range_b);
}

int range_compare_with_delta(Range range_a, Range range_b, unsigned long delta)
{
  unsigned long start_min, start_max, end_min, end_max;

  assert(range_a.start <= range_a.end && range_b.start <= range_b.end);

  start_min = MIN(range_a.start, range_b.start);
  start_max = MAX(range_a.start, range_b.start);
  end_min   = MIN(range_a.end, range_b.end);
  end_max   = MAX(range_a.end, range_b.end);

  if (start_max - start_min <= delta && end_max - end_min <= delta)
    return 0; /* range_a == range_b */

  if ((range_a.start < range_b.start) ||
      ((range_a.start == range_b.start) && (range_a.end < range_b.end)))
    return -1; /* range_a < range_b */

  return 1; /* range_a > range_b */
}

int range_compare_by_length_ptr(const Range *range_a, const Range *range_b)
{
  unsigned long range_a_length, range_b_length;
  assert(range_a && range_b);
  range_a_length = range_length(*range_a);
  range_b_length = range_length(*range_b);
  if (range_a_length == range_b_length)
    return 0;
  if (range_a_length > range_b_length)
    return -1;
  return 1;
}

bool range_overlap(Range range_a, Range range_b)
{
  assert(range_a.start <= range_a.end && range_b.start <= range_b.end);

  if (range_a.start <= range_b.end && range_a.end >= range_b.start)
    return true;
  return false;
}

bool range_contains(Range range_a, Range range_b)
{
  assert(range_a.start <= range_a.end && range_b.start <= range_b.end);

  if (range_a.start <= range_b.start && range_a.end >= range_b.end)
    return true;
  return false;
}

bool range_within(Range range, unsigned long point)
{
  assert(range.start <= range.end);

  if (range.start <= point && range.end >= point)
    return true;
  return false;
}

Range range_join(Range range_a, Range range_b)
{
  Range r;

  assert(range_a.start <= range_a.end && range_b.start <= range_b.end);

  r.start = range_a.start < range_b.start ? range_a.start : range_b.start;
  r.end   = range_a.end   > range_b.end   ? range_a.end   : range_b.end;

  return r;
}

Range range_offset(Range range, long offset)
{
  Range transformed_range = { 0, 0 };
  assert(range.start <= range.end);
  safe_add(transformed_range.start, range.start, offset);
  safe_add(transformed_range.end, range.end, offset);
  assert(transformed_range.start <= transformed_range.end);
  return transformed_range;
}

Range range_reorder(Range range)
{
  Range ordered_range;
  if (range.start <= range.end)
    return range;
  ordered_range.start = range.end;
  ordered_range.end   = range.start;
  return ordered_range;
}

unsigned long range_length(Range range)
{
  assert(range.start <= range.end);
  return range.end - range.start + 1;
}

int range_unit_test(Error *err)
{
  static Range ranges_in[] = {  { 620432, 620536 }, { 620432, 620536 },
                                { 620957, 621056 }, { 620957, 621056 },
                                { 625234, 625253 }, { 625500, 625655 },
                                { 625533, 625655 }, { 625533, 625655 },
                                { 627618, 627729 }, { 627618, 627729 },
                                { 627618, 627729 }, { 662083, 662194 },
                                { 662083, 662194 }, { 662083, 662194 },
                                { 663032, 663166 }, { 663032, 663166 },
                                { 663032, 663166 }, { 664782, 664906 },
                                { 664782, 664906 }, { 664782, 664906 },
                                { 665748, 665823 }, { 665748, 665823 },
                                { 665748, 665823 }, { 666825, 666881 },
                                { 666825, 666881 }, { 667797, 667954 },
                                { 667845, 667954 }, { 667845, 667954 },
                                { 679175, 679280 }, { 679175, 679280 },
                                { 679175, 679280 }, { 680427, 680540 },
                                { 680427, 680540 }, { 680427, 680540 },
                                { 684144, 684293 }, { 684144, 684293 },
                                { 684144, 684293 }, { 724903, 724985 },
                                { 724903, 724985 }, { 727099, 727325 },
                                { 727099, 727325 }, { 732544, 732821 },
                                { 732544, 732821 }, { 750016, 750280 },
                                { 750016, 750280 }, { 769508, 769734 },
                                { 769508, 769734 } },
               ranges_out[] = { { 620432, 620536 }, { 620957, 621056 },
                                { 625234, 625253 }, { 625500, 625655 },
                                { 625533, 625655 }, { 627618, 627729 },
                                { 662083, 662194 }, { 663032, 663166 },
                                { 664782, 664906 }, { 665748, 665823 },
                                { 666825, 666881 }, { 667797, 667954 },
                                { 667845, 667954 }, { 679175, 679280 },
                                { 680427, 680540 }, { 684144, 684293 },
                                { 724903, 724985 }, { 727099, 727325 },
                                { 732544, 732821 }, { 750016, 750280 },
                                { 769508, 769734 }};
  unsigned long counts[] = { 2, 2, 1, 1, 2, 3, 3, 3, 3, 3, 2, 1, 2, 3, 3, 3, 2,
                             2, 2, 2, 2 };
  Array *ranges, *tmp_ranges, *ctr;
  unsigned long i;
  int had_err = 0;
  error_check(err);

  ensure(had_err, sizeof (ranges_out) / sizeof (ranges_out[0]) ==
                  sizeof (counts)     / sizeof (counts[0]));

  /* test ranges_uniq() */
  ranges = array_new(sizeof (Range));
  tmp_ranges = array_new(sizeof (Range));
  for (i = 0; i < sizeof (ranges_in) / sizeof (ranges_in[0]) && !had_err; i++)
    array_add(ranges, ranges_in[i]);
  ranges_uniq(tmp_ranges, ranges);
  ensure(had_err, array_size(ranges) ==
                  sizeof (ranges_in) / sizeof (ranges_in[0]));
  ensure(had_err, array_size(tmp_ranges) ==
                  sizeof (ranges_out) / sizeof (ranges_out[0]));
  for (i = 0; i < array_size(tmp_ranges) && !had_err; i++) {
    ensure(had_err,
           ranges_out[i].start == (*(Range*) array_get(tmp_ranges, i)).start);
    ensure(had_err,
           ranges_out[i].end == (*(Range*) array_get(tmp_ranges, i)).end);
  }

  /* test ranges_uniq_in_place() */
  array_reset(tmp_ranges);
  array_add_array(tmp_ranges, ranges);
  ranges_uniq_in_place(tmp_ranges);
  for (i = 0; i < array_size(tmp_ranges) && !had_err; i++) {
    ensure(had_err,
           ranges_out[i].start == (*(Range*) array_get(tmp_ranges, i)).start);
    ensure(had_err,
           ranges_out[i].end == (*(Range*) array_get(tmp_ranges, i)).end);
  }

  /* test ranges_uniq_count() */
  array_reset(tmp_ranges);
  ctr = ranges_uniq_count(tmp_ranges, ranges);
  ensure(had_err, array_size(tmp_ranges) == array_size(ctr));
  ensure(had_err, array_size(ctr) == sizeof (counts) / sizeof (counts[0]));
  for (i = 0; i < array_size(ctr) && !had_err; i++) {
    ensure(had_err, counts[i] == *(unsigned long*) array_get(ctr, i));
    ensure(had_err,
           ranges_out[i].start == (*(Range*) array_get(tmp_ranges, i)).start);
    ensure(had_err,
           ranges_out[i].end == (*(Range*) array_get(tmp_ranges, i)).end);
  }
  array_delete(ctr);

  /* test ranges_uniq_in_place_count() */
  ctr = ranges_uniq_in_place_count(ranges);
  ensure(had_err, array_size(ranges) == array_size(ctr));
  ensure(had_err, array_size(ctr) == sizeof (counts) / sizeof (counts[0]));
  for (i = 0; i < array_size(ctr) && !had_err; i++) {
    ensure(had_err, counts[i] == *(unsigned long*) array_get(ctr, i));
    ensure(had_err,
           ranges_out[i].start == (*(Range*) array_get(ranges, i)).start);
    ensure(had_err,
           ranges_out[i].end == (*(Range*) array_get(ranges, i)).end);
  }
  array_delete(ctr);

  /* test range_reorder() */
  if (!had_err) {
    Range range = { 1, 100 };
    range = range_reorder(range);
    ensure(had_err, range.start == 1 && range.end == 100);
    range.start = 100;
    range.end = 1;
    range = range_reorder(range);
    ensure(had_err, range.start == 1 && range.end == 100);
  }

  /* free */
  array_delete(ranges);
  array_delete(tmp_ranges);
  return had_err;
}

void ranges_sort(Array *ranges)
{
  assert(ranges);
  qsort(array_get_space(ranges), array_size(ranges), sizeof (Range),
        (Compare) range_compare_ptr);
}

void ranges_sort_by_length_stable(Array *ranges)
{
  assert(ranges);
  msort(array_get_space(ranges), array_size(ranges), sizeof (Range),
        (Compare) range_compare_by_length_ptr);
}

bool ranges_are_sorted(const Array *ranges)
{
  unsigned long i;

  assert(ranges);

  for (i = 1; i < array_size(ranges); i++) {
    if (range_compare(*(Range*) array_get(ranges, i-1),
                      *(Range*) array_get(ranges, i)) == 1) {
      return false;
    }
  }
  return true;
}

bool ranges_do_not_overlap(const Array *ranges)
{
  unsigned long i;

  assert(ranges && array_size(ranges));

  for (i = 1; i < array_size(ranges); i++) {
    if (range_overlap(*(Range*) array_get(ranges, i-1),
                      *(Range*) array_get(ranges, i))) {
      return false;
    }
  }
  return true;
}

bool ranges_are_sorted_and_do_not_overlap(const Array *ranges)
{
  return ranges_are_sorted(ranges) && ranges_do_not_overlap(ranges);
}

bool ranges_are_equal(const Array *ranges_1, const Array *ranges_2)
{
  unsigned long i;
  Range range_1, range_2;

  assert(ranges_are_sorted(ranges_1) && ranges_are_sorted(ranges_2));

  if (array_size(ranges_1) != array_size(ranges_2))
    return false;

  for (i = 0; i < array_size(ranges_1); i++) {
    range_1 = *(Range*) array_get(ranges_1, i);
    range_2 = *(Range*) array_get(ranges_2, i);
    if (range_compare(range_1, range_2))
      return false;
  }

  return true;
}

static Array* generic_ranges_uniq(Array *out_ranges, const Array *in_ranges,
                                  bool count)
{
  unsigned long i, *ctr_ptr, ctr = 1;
  Array *count_array = NULL;
  Range cur  = { UNDEF_ULONG, UNDEF_ULONG },
        prev = { UNDEF_ULONG, UNDEF_ULONG };
  assert(out_ranges && in_ranges);
  assert(ranges_are_sorted(in_ranges));
  if (count)
    count_array = array_new(sizeof (unsigned long));
  for (i = 0; i < array_size(in_ranges); i++) {
    cur = *(Range*) array_get(in_ranges, i);
    if (!i) {
      array_add(out_ranges, cur);
      if (count)
        array_add(count_array, ctr);
    }
    else {
      if (prev.start == cur.start && prev.end == cur.end) {
        if (count) {
          ctr_ptr = array_get_last(count_array);
          (*ctr_ptr)++;
        }
      }
      else {
        array_add(out_ranges, cur);
        if (count)
          array_add(count_array, ctr);
      }
    }
    prev = cur;
  }
  return count_array;
}

static Array* generic_ranges_uniq_in_place(Array *ranges, bool count)
{
  Array *out_ranges, *count_array;
  assert(ranges);
  out_ranges = array_new(sizeof (Range));
  count_array = generic_ranges_uniq(out_ranges, ranges, count);
  array_reset(ranges);
  array_add_array(ranges, out_ranges); /* XXX: could be more efficient
                                               with something like
                                               array_replace(ranges,
                                                             out_ranges) */
  array_delete(out_ranges);
  return count_array;
}

void ranges_uniq(Array *out_ranges, const Array *in_ranges)
{
  assert(out_ranges && in_ranges);
  (void) generic_ranges_uniq(out_ranges, in_ranges, false);
}

void ranges_uniq_in_place(Array *ranges)
{
  assert(ranges);
  (void) generic_ranges_uniq_in_place(ranges, false);
}

Array* ranges_uniq_count(Array *out_ranges, const Array *in_ranges)
{
  assert(out_ranges && in_ranges);
  return generic_ranges_uniq(out_ranges, in_ranges, true);
}

Array* ranges_uniq_in_place_count(Array *ranges)
{
  assert(ranges);
  return generic_ranges_uniq_in_place(ranges, true);
}
