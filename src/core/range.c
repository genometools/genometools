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

#include "core/ensure.h"
#include "core/minmax.h"
#include "core/msort.h"
#include "core/range.h"
#include "core/safearith.h"
#include "core/undef_api.h"

int gt_range_compare(const GtRange *range_a, const GtRange *range_b)
{
  gt_assert(range_a->start <= range_a->end && range_b->start <= range_b->end);

  if ((range_a->start == range_b->start) && (range_a->end == range_b->end))
    return 0; /* range_a == range_b */

  if ((range_a->start < range_b->start) ||
      ((range_a->start == range_b->start) && (range_a->end < range_b->end)))
    return -1; /* range_a < range_b */

  return 1; /* range_a > range_b */
}

int gt_range_compare_with_delta(const GtRange *range_a, const GtRange *range_b,
                                unsigned long delta)
{
  unsigned long start_min, start_max, end_min, end_max;

  gt_assert(range_a->start <= range_a->end && range_b->start <= range_b->end);

  start_min = MIN(range_a->start, range_b->start);
  start_max = MAX(range_a->start, range_b->start);
  end_min   = MIN(range_a->end, range_b->end);
  end_max   = MAX(range_a->end, range_b->end);

  if (start_max - start_min <= delta && end_max - end_min <= delta)
    return 0; /* range_a == range_b */

  if ((range_a->start < range_b->start) ||
      ((range_a->start == range_b->start) && (range_a->end < range_b->end)))
    return -1; /* range_a < range_b */

  return 1; /* range_a > range_b */
}

int gt_range_compare_by_length_ptr(const GtRange *range_a,
                                   const GtRange *range_b)
{
  unsigned long range_a_length, range_b_length;
  gt_assert(range_a && range_b);
  range_a_length = gt_range_length(range_a);
  range_b_length = gt_range_length(range_b);
  if (range_a_length == range_b_length)
    return 0;
  if (range_a_length > range_b_length)
    return -1;
  return 1;
}

bool gt_range_overlap(const GtRange *range_a, const GtRange *range_b)
{
  gt_assert(range_a->start <= range_a->end && range_b->start <= range_b->end);

  if (range_a->start <= range_b->end && range_a->end >= range_b->start)
    return true;
  return false;
}

bool gt_range_overlap_delta(const GtRange *range_a, const GtRange *range_b,
                            unsigned long delta)
{
  unsigned long range_a_length, range_b_length;
  gt_assert(range_a->start <= range_a->end && range_b->start <= range_b->end);

  range_a_length = range_a->end - range_a->start + 1;
  range_b_length = range_b->end - range_b->start + 1;

  if (range_a_length < delta || range_b_length < delta) {
    /* no overlap of delta possible */
    return false;
  }

  if (gt_range_overlap(range_a, range_b)) {
    if (range_a->start <= range_b->start) {
      if (range_a->end >= range_b->end) {
        /* ----A----
            ---B---  */
        if (range_b_length >= delta)
          return true;
      }
      else { /* range_a->end < range_b->end */
        /* ----A----
            ----B---- */
        if (range_a->end - range_b->start + 1 >= delta)
          return true;
      }
    }
    else { /* range_a->start > range_b->start */
      if (range_a->end <= range_b->end) {
        /*  ---A---
           ----B---- */
        if (range_a_length >= delta)
          return true;
      }
      else { /* range_a->end > range_b->end */
        /*  ----A----
           ----B----  */
        if (range_b->end - range_a->start + 1 >= delta)
          return true;
      }
    }
  }

  return false;
}

bool gt_range_contains(const GtRange *range_a, const GtRange *range_b)
{
  gt_assert(range_a->start <= range_a->end && range_b->start <= range_b->end);

  if (range_a->start <= range_b->start && range_a->end >= range_b->end)
    return true;
  return false;
}

bool gt_range_within(const GtRange *range, unsigned long point)
{
  gt_assert(range->start <= range->end);

  if (range->start <= point && range->end >= point)
    return true;
  return false;
}

GtRange gt_range_join(const GtRange *range_a, const GtRange *range_b)
{
  GtRange r;

  gt_assert(range_a->start <= range_a->end && range_b->start <= range_b->end);

  r.start = range_a->start < range_b->start ? range_a->start : range_b->start;
  r.end   = range_a->end   > range_b->end   ? range_a->end   : range_b->end;

  return r;
}

GtRange gt_range_offset(const GtRange *range, long offset)
{
  GtRange transformed_range = { 0, 0 };
  gt_assert(range->start <= range->end);
  gt_safe_add(transformed_range.start, range->start, offset);
  gt_safe_add(transformed_range.end, range->end, offset);
  gt_assert(transformed_range.start <= transformed_range.end);
  return transformed_range;
}

GtRange gt_range_reorder(GtRange range)
{
  GtRange ordered_range;
  if (range.start <= range.end)
    return range;
  ordered_range.start = range.end;
  ordered_range.end   = range.start;
  return ordered_range;
}

unsigned long gt_range_length(const GtRange *range)
{
  gt_assert(range->start <= range->end);
  return range->end - range->start + 1;
}

int gt_range_unit_test(GtError *err)
{
  static GtRange ranges_in[] = {  { 620432, 620536 }, { 620432, 620536 },
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
  GtArray *ranges, *tmp_ranges, *ctr;
  unsigned long i;
  int had_err = 0;
  gt_error_check(err);

  gt_ensure(had_err, sizeof (ranges_out) / sizeof (ranges_out[0]) ==
                  sizeof (counts)     / sizeof (counts[0]));

  /* test gt_ranges_uniq() */
  ranges = gt_array_new(sizeof (GtRange));
  tmp_ranges = gt_array_new(sizeof (GtRange));
  for (i = 0;
       i < sizeof (ranges_in) / sizeof (ranges_in[0]) && !had_err;
       i++)
    gt_array_add(ranges, ranges_in[i]);
  gt_ranges_uniq(tmp_ranges, ranges);
  gt_ensure(had_err, gt_array_size(ranges) ==
                  sizeof (ranges_in) / sizeof (ranges_in[0]));
  gt_ensure(had_err, gt_array_size(tmp_ranges) ==
                  sizeof (ranges_out) / sizeof (ranges_out[0]));
  for (i = 0; i < gt_array_size(tmp_ranges) && !had_err; i++) {
    gt_ensure(had_err, ranges_out[i].start ==
                    (*(GtRange*) gt_array_get(tmp_ranges, i)).start);
    gt_ensure(had_err, ranges_out[i].end ==
                    (*(GtRange*) gt_array_get(tmp_ranges, i)).end);
  }

  /* test gt_ranges_uniq_in_place() */
  gt_array_reset(tmp_ranges);
  gt_array_add_array(tmp_ranges, ranges);
  gt_ranges_uniq_in_place(tmp_ranges);
  for (i = 0; i < gt_array_size(tmp_ranges) && !had_err; i++) {
    gt_ensure(had_err, ranges_out[i].start ==
                    (*(GtRange*) gt_array_get(tmp_ranges, i)).start);
    gt_ensure(had_err, ranges_out[i].end ==
                    (*(GtRange*) gt_array_get(tmp_ranges, i)).end);
  }

  /* test gt_ranges_uniq_count() */
  gt_array_reset(tmp_ranges);
  ctr = gt_ranges_uniq_count(tmp_ranges, ranges);
  gt_ensure(had_err, gt_array_size(tmp_ranges) == gt_array_size(ctr));
  gt_ensure(had_err,
            gt_array_size(ctr) == sizeof (counts) / sizeof (counts[0]));
  for (i = 0; i < gt_array_size(ctr) && !had_err; i++) {
    gt_ensure(had_err, counts[i] == *(unsigned long*) gt_array_get(ctr, i));
    gt_ensure(had_err, ranges_out[i].start ==
                    (*(GtRange*) gt_array_get(tmp_ranges, i)).start);
    gt_ensure(had_err, ranges_out[i].end ==
                    (*(GtRange*) gt_array_get(tmp_ranges, i)).end);
  }
  gt_array_delete(ctr);

  /* test gt_ranges_uniq_in_place_count() */
  ctr = gt_ranges_uniq_in_place_count(ranges);
  gt_ensure(had_err, gt_array_size(ranges) == gt_array_size(ctr));
  gt_ensure(had_err,
            gt_array_size(ctr) == sizeof (counts) / sizeof (counts[0]));
  for (i = 0; i < gt_array_size(ctr) && !had_err; i++) {
    gt_ensure(had_err, counts[i] == *(unsigned long*) gt_array_get(ctr, i));
    gt_ensure(had_err,
           ranges_out[i].start == (*(GtRange*)
                                             gt_array_get(ranges, i)).start);
    gt_ensure(had_err,
           ranges_out[i].end == (*(GtRange*) gt_array_get(ranges, i)).end);
  }
  gt_array_delete(ctr);

  /* test gt_range_reorder() */
  if (!had_err) {
    GtRange range = { 1, 100 };
    range = gt_range_reorder(range);
    gt_ensure(had_err, range.start == 1 && range.end == 100);
    range.start = 100;
    range.end = 1;
    range = gt_range_reorder(range);
    gt_ensure(had_err, range.start == 1 && range.end == 100);
  }

  /* free */
  gt_array_delete(ranges);
  gt_array_delete(tmp_ranges);
  return had_err;
}

void gt_ranges_sort(GtArray *ranges)
{
  gt_assert(ranges);
  qsort(gt_array_get_space(ranges), gt_array_size(ranges), sizeof (GtRange),
        (GtCompare) gt_range_compare);
}

void gt_ranges_sort_by_length_stable(GtArray *ranges)
{
  gt_assert(ranges);
  gt_msort(gt_array_get_space(ranges), gt_array_size(ranges), sizeof (GtRange),
           (GtCompare) gt_range_compare_by_length_ptr);
}

bool gt_ranges_are_sorted(const GtArray *ranges)
{
  unsigned long i;

  gt_assert(ranges);

  for (i = 1; i < gt_array_size(ranges); i++) {
    if (gt_range_compare(gt_array_get(ranges, i-1),
                         gt_array_get(ranges, i)) == 1) {
      return false;
    }
  }
  return true;
}

bool gt_ranges_do_not_overlap(const GtArray *ranges)
{
  unsigned long i;

  gt_assert(ranges && gt_array_size(ranges));

  for (i = 1; i < gt_array_size(ranges); i++) {
    if (gt_range_overlap(gt_array_get(ranges, i-1), gt_array_get(ranges, i)))
      return false;
  }
  return true;
}

bool gt_ranges_are_sorted_and_do_not_overlap(const GtArray *ranges)
{
  return gt_ranges_are_sorted(ranges) && gt_ranges_do_not_overlap(ranges);
}

bool gt_ranges_are_equal(const GtArray *ranges_1, const GtArray *ranges_2)
{
  unsigned long i;

  gt_assert(gt_ranges_are_sorted(ranges_1) && gt_ranges_are_sorted(ranges_2));

  if (gt_array_size(ranges_1) != gt_array_size(ranges_2))
    return false;

  for (i = 0; i < gt_array_size(ranges_1); i++) {
    if (gt_range_compare(gt_array_get(ranges_1, i), gt_array_get(ranges_2, i)))
      return false;
  }

  return true;
}

static GtArray* generic_ranges_uniq(GtArray *out_ranges,
                                    const GtArray *in_ranges, bool count)
{
  unsigned long i, *ctr_ptr, ctr = 1;
  GtArray *count_array = NULL;
  GtRange cur  = { GT_UNDEF_ULONG, GT_UNDEF_ULONG },
        prev = { GT_UNDEF_ULONG, GT_UNDEF_ULONG };
  gt_assert(out_ranges && in_ranges);
  gt_assert(gt_ranges_are_sorted(in_ranges));
  if (count)
    count_array = gt_array_new(sizeof (unsigned long));
  for (i = 0; i < gt_array_size(in_ranges); i++) {
    cur = *(GtRange*) gt_array_get(in_ranges, i);
    if (!i) {
      gt_array_add(out_ranges, cur);
      if (count)
        gt_array_add(count_array, ctr);
    }
    else {
      if (prev.start == cur.start && prev.end == cur.end) {
        if (count) {
          ctr_ptr = gt_array_get_last(count_array);
          (*ctr_ptr)++;
        }
      }
      else {
        gt_array_add(out_ranges, cur);
        if (count)
          gt_array_add(count_array, ctr);
      }
    }
    prev = cur;
  }
  return count_array;
}

static GtArray* generic_ranges_uniq_in_place(GtArray *ranges, bool count)
{
  GtArray *out_ranges, *count_array;
  gt_assert(ranges);
  out_ranges = gt_array_new(sizeof (GtRange));
  count_array = generic_ranges_uniq(out_ranges, ranges, count);
  gt_array_reset(ranges);
  gt_array_add_array(ranges, out_ranges); /* XXX: could be more efficient
                                               with something like
                                               gt_array_replace(ranges,
                                                             out_ranges) */
  gt_array_delete(out_ranges);
  return count_array;
}

void gt_ranges_uniq(GtArray *out_ranges, const GtArray *in_ranges)
{
  gt_assert(out_ranges && in_ranges);
  (void) generic_ranges_uniq(out_ranges, in_ranges, false);
}

void gt_ranges_uniq_in_place(GtArray *ranges)
{
  gt_assert(ranges);
  (void) generic_ranges_uniq_in_place(ranges, false);
}

GtArray* gt_ranges_uniq_count(GtArray *out_ranges, const GtArray *in_ranges)
{
  gt_assert(out_ranges && in_ranges);
  return generic_ranges_uniq(out_ranges, in_ranges, true);
}

GtArray* gt_ranges_uniq_in_place_count(GtArray *ranges)
{
  gt_assert(ranges);
  return generic_ranges_uniq_in_place(ranges, true);
}

bool gt_ranges_are_consecutive(const GtArray *ranges)
{
  unsigned long i;
  for (i = 0; i < gt_array_size(ranges); i++) {
    gt_assert(((GtRange*) gt_array_get(ranges, i))->start <=
              ((GtRange*) gt_array_get(ranges, i))->end);
    if (i) {
      /* check if ranges are consecutive */
      if (((GtRange*) gt_array_get(ranges, i-1))->end >=
          ((GtRange*) gt_array_get(ranges, i))->start) {
        return false;
      }
    }
  }
  return true;
}

unsigned long gt_ranges_total_length(const GtArray *ranges)
{
  unsigned long i, totallen = 0;
  GtRange *range;
  gt_assert(ranges);
  for (i = 0; i < gt_array_size(ranges); i++) {
    range = gt_array_get(ranges, i);
    totallen += range->end - range->start + 1;
  }
  return totallen;
}

unsigned long gt_ranges_spanned_length(const GtArray *ranges)
{
  GtRange spanned_range;
  gt_assert(ranges);
  spanned_range.start = ((GtRange*) gt_array_get_first(ranges))->start;
  spanned_range.end   = ((GtRange*) gt_array_get_last(ranges))->end;
  return gt_range_length(&spanned_range);
}

void gt_ranges_copy_to_opposite_strand(GtArray *outranges,
                                       const GtArray *inranges,
                                       unsigned long gen_total_length,
                                       unsigned long gen_offset)
{
  GtRange range;
  unsigned long i;

  /* outranges are empty */
  gt_assert(!gt_array_size(outranges));
  /* inranges are not empty */
  gt_assert(gt_array_size(inranges));

  for (i = gt_array_size(inranges); i > 0; i--) {
    /* genomic offset is defined */
    gt_assert(gen_offset != GT_UNDEF_ULONG);
    range.start  = gen_total_length - 1
                  - (((GtRange*) gt_array_get(inranges, i-1))->end -
                     gen_offset)
                  + gen_offset;
    range.end = gen_total_length - 1
                  - (((GtRange*) gt_array_get(inranges, i-1))->start -
                     gen_offset)
                  + gen_offset;
    gt_array_add(outranges, range);
  }

  /* outranges has the same number of elements as inranges */
  gt_assert(gt_array_size(inranges) == gt_array_size(outranges));
}

bool gt_ranges_borders_are_in_region(GtArray *ranges, const GtRange *region)
{
  gt_assert(ranges && region);

  /* check region start */
  if (((GtRange*) gt_array_get_first(ranges))->start < region->start)
    return false;

  /* check region end */
  if (((GtRange*) gt_array_get_last(ranges))->end > region->end)
    return false;

  return true;
}

void gt_ranges_show(GtArray *ranges, GtFile *outfp)
{
  GtRange *range;
  unsigned long i;
  gt_assert(ranges);
  for (i = 0; i < gt_array_size(ranges); i++) {
    range = gt_array_get(ranges, i);
    gt_file_xprintf(outfp, "(%lu,%lu)", range->start, range->end);
  }
  gt_file_xfputc('\n', outfp);
}
