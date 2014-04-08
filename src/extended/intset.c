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

#include <inttypes.h>
#include <limits.h>

#include "core/assert_api.h"
#include "core/ensure.h"
#include "core/intbits.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/types_api.h"
#include "extended/intset.h"

#define BITS_FOR_SIZE(SIZE)     ((SIZE) * CHAR_BIT)
#define ELEM2SECTION(LOGVAL,X)  ((X) >> (LOGVAL))

GtIntsetType gt_intset_best_type(GtUword maxelement,
                                 GtUword num_of_elems)
{
  GtIntsetType best = GT_INTSET_8;
  size_t min, tmp;
  min = gt_intset_8_size(maxelement, num_of_elems);
  tmp = gt_intset_16_size(maxelement, num_of_elems);
  best = min < tmp ? best : GT_INTSET_16;
  min = min < tmp ? min : tmp;
  tmp = gt_intset_32_size(maxelement, num_of_elems);
  return min < tmp ? best : GT_INTSET_32;
}

typedef bool (*in_set_func)(const void *set, GtUword elem);

int gt_intset_unit_test_notinset(const void *intset, GtUword start, GtUword end,
                                 GtIntsetType type, GtError *err)
{
  int had_err = 0;
  GtUword test;
  in_set_func callback = NULL;
  switch (type) {
    case GT_INTSET_8:
      callback = gt_intset_8_is_member_fp;
      break;
    case GT_INTSET_16:
      callback = gt_intset_16_is_member_fp;
      break;
    case GT_INTSET_32:
      callback = gt_intset_32_is_member_fp;
      break;
  }
  for (test = start; test <= end; ++test) {
    gt_ensure(callback(intset, test) == false);
  }
  return had_err;
}

typedef GtUword (*seqnum_func)(const void *set, GtUword elem);

int gt_intset_unit_test_check_seqnum(const void *intset, GtUword start,
                                     GtUword end, GtUword num,
                                     GtIntsetType type, GtError *err)
{
  int had_err = 0;
  GtUword test;
  seqnum_func callback = NULL;
  switch (type) {
    case GT_INTSET_8:
      callback = gt_intset_8_pos2seqnum_fp;
      break;
    case GT_INTSET_16:
      callback = gt_intset_16_pos2seqnum_fp;
      break;
    case GT_INTSET_32:
      callback = gt_intset_32_pos2seqnum_fp;
      break;
  }
  for (test = start; test <= end; ++test) {
    gt_ensure(callback(intset, test) == num);
  }
  return had_err;
}

int gt_intset_unit_test(GT_UNUSED GtError *err) {
  int had_err = 0;
  GtIntset8 *is8 = NULL;
  GtIntset16 *is16 = NULL;
  GtIntset32 *is32 = NULL;
  GtUword num_of_elems = gt_rand_max(((GtUword) 1) << 10) + 1,
          *arr = gt_malloc(sizeof (*arr) * num_of_elems),
          stepsize = (num_of_elems <<4 / num_of_elems) >> 1,
          idx;
  size_t is8size, is16size, is32size;
  GtIntsetType type;

  gt_error_check(err);

  arr[0] = gt_rand_max(stepsize) + 1;
  for (idx = (GtUword) 1; idx < num_of_elems; ++idx) {
    arr[idx] = arr[idx - 1] + gt_rand_max(stepsize) + 1;
  }

  is8size = gt_intset_8_size(arr[num_of_elems - 1], num_of_elems);
  is16size = gt_intset_16_size(arr[num_of_elems - 1], num_of_elems);
  is32size = gt_intset_32_size(arr[num_of_elems - 1], num_of_elems);
  type = gt_intset_best_type(arr[num_of_elems - 1], num_of_elems);
  switch (type) {
    case GT_INTSET_8:
      gt_ensure(is8size < is16size && is8size < is32size);
      break;
    case GT_INTSET_16:
      gt_ensure(is16size <= is8size && is16size < is32size);
      break;
    case GT_INTSET_32:
      gt_ensure(is32size <= is8size && is32size <= is16size);
      break;
  }

  if (!had_err) {
    if (is8size < (size_t) UINT_MAX) {
      is8 = gt_intset_8_new(arr[num_of_elems - 1], num_of_elems);
      for (idx = 0; idx < num_of_elems; idx++) {
        gt_intset_8_add(is8, arr[idx]);
      }
      had_err = gt_intset_unit_test_notinset(is8, 0, arr[0] - 1,
                                             GT_INTSET_8, err);
      if (!had_err)
        had_err = gt_intset_unit_test_check_seqnum(is8, 0, arr[0] - 1, 0,
                                                   GT_INTSET_8, err);
      for (idx = (GtUword) 1; !had_err && idx < num_of_elems; idx++) {
        had_err = gt_intset_unit_test_notinset(is8, arr[idx - 1] + 1,
                                               arr[idx] - 1, GT_INTSET_8, err);
        if (!had_err)
          had_err = gt_intset_unit_test_check_seqnum(is8, arr[idx - 1] + 1,
                                                     arr[idx] - 1, idx,
                                                     GT_INTSET_8, err);
      }
      gt_intset_8_delete(is8);
    }

    if (!had_err && is16size < (size_t) UINT_MAX) {
      is16 = gt_intset_16_new(arr[num_of_elems - 1], num_of_elems);
      for (idx = 0; idx < num_of_elems; idx++) {
        gt_intset_16_add(is16, arr[idx]);
      }
      had_err = gt_intset_unit_test_notinset(is16, 0, arr[0] - 1,
                                             GT_INTSET_16, err);
      if (!had_err)
        had_err = gt_intset_unit_test_check_seqnum(is16, 0, arr[0] - 1, 0,
                                                   GT_INTSET_16, err);
      for (idx = (GtUword) 1; !had_err && idx < num_of_elems; idx++) {
        had_err = gt_intset_unit_test_notinset(is16, arr[idx - 1] + 1,
                                               arr[idx] - 1, GT_INTSET_16, err);
        if (!had_err)
          had_err = gt_intset_unit_test_check_seqnum(is16, arr[idx - 1] + 1,
                                                     arr[idx] - 1, idx,
                                                     GT_INTSET_16, err);
      }
      gt_intset_16_delete(is16);
    }
    if (!had_err && is32size < (size_t) UINT_MAX) {
      is32 = gt_intset_32_new(arr[num_of_elems - 1], num_of_elems);
      for (idx = 0; idx < num_of_elems; idx++) {
        gt_intset_32_add(is32, arr[idx]);
      }
      had_err = gt_intset_unit_test_notinset(is32, 0, arr[0] - 1,
                                             GT_INTSET_32, err);
      if (!had_err)
        had_err = gt_intset_unit_test_check_seqnum(is32, 0, arr[0] - 1, 0,
                                                   GT_INTSET_32, err);
      for (idx = (GtUword) 1; !had_err && idx < num_of_elems; idx++) {
        had_err = gt_intset_unit_test_notinset(is32, arr[idx - 1] + 1,
                                               arr[idx] - 1, GT_INTSET_32, err);
        if (!had_err)
          had_err = gt_intset_unit_test_check_seqnum(is32, arr[idx - 1] + 1,
                                                     arr[idx] - 1, idx,
                                                     GT_INTSET_32, err);
      }
      gt_intset_32_delete(is32);
    }
  }
  gt_free(arr);
  return had_err;
}
