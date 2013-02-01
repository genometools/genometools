/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include "core/array.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/str.h"
#include "extended/reverse_api.h"
#include "extended/splicedseq.h"

struct Splicedseq {
  GtStr *splicedseq;
  GtArray *positionmapping;
  bool forward;
};

Splicedseq* gt_splicedseq_new(void)
{
  Splicedseq *ss = gt_malloc(sizeof (Splicedseq));
  ss->splicedseq = gt_str_new();
  ss->positionmapping = gt_array_new(sizeof (unsigned long));
  ss->forward = true;
  return ss;
}

void gt_splicedseq_add(Splicedseq *ss, unsigned long start, unsigned long end,
                       const char *original_sequence)
{
  unsigned long i;
  gt_assert(ss && start <= end && original_sequence);
  gt_str_append_cstr_nt(ss->splicedseq, original_sequence,
                        end - start + 1);
  /* make sure elements are added in ascending order */
  gt_assert(!gt_array_size(ss->positionmapping) ||
            start > *(unsigned long*) gt_array_get_last(ss->positionmapping));
  for (i = start; i <= end; i++)
    gt_array_add(ss->positionmapping, i);
}

char* gt_splicedseq_get(const Splicedseq *ss)
{
  return gt_str_get(ss->splicedseq);
}

bool gt_splicedseq_pos_is_border(const Splicedseq *ss, unsigned long pos)
{
  gt_assert(ss &&
         gt_str_length(ss->splicedseq) == gt_array_size(ss->positionmapping));
  gt_assert(pos < gt_str_length(ss->splicedseq)); /* legal position */
  if (ss->forward && pos + 1 < gt_array_size(ss->positionmapping) &&
      *(unsigned long*) gt_array_get(ss->positionmapping, pos) + 1 !=
      *(unsigned long*) gt_array_get(ss->positionmapping, pos+1)) {
    return true;
  }
  if (!ss->forward && pos &&
      *(unsigned long*) gt_array_get(ss->positionmapping, pos-1) - 1 !=
      *(unsigned long*) gt_array_get(ss->positionmapping, pos)) {
    return true;
  }
  return false;
}

unsigned long gt_splicedseq_map(const Splicedseq *ss, unsigned long pos)
{
  gt_assert(ss &&
         gt_str_length(ss->splicedseq) == gt_array_size(ss->positionmapping));
  gt_assert(pos < gt_str_length(ss->splicedseq)); /* legal position */
  return *(unsigned long*) gt_array_get(ss->positionmapping, pos);
}

unsigned long gt_splicedseq_length(const Splicedseq *ss)
{
  gt_assert(ss);
  return gt_str_length(ss->splicedseq);
}

int gt_splicedseq_reverse(Splicedseq *ss, GtError *err)
{
  int had_err;
  gt_error_check(err);
  gt_assert(ss);
  had_err = gt_reverse_complement(gt_str_get(ss->splicedseq),
                               gt_str_length(ss->splicedseq), err);
  if (!had_err) {
    gt_array_reverse(ss->positionmapping);
    ss->forward = !ss->forward;
  }
  return had_err;
}

void gt_splicedseq_reset(Splicedseq *ss)
{
  gt_assert(ss);
  gt_str_reset(ss->splicedseq);
  gt_array_reset(ss->positionmapping);
  ss->forward = true;
}

static int check_splicedseq(Splicedseq *ss, GtError *err)
{                       /*0123456789*/
  static char *origseq = "aaccaagtga", *splicedseq = "ccgtg";
  int had_err = 0;
  gt_error_check(err);
  gt_splicedseq_add(ss, 2, 3, origseq + 2);
  gt_splicedseq_add(ss, 6, 8, origseq + 6);
  gt_ensure(had_err, strcmp(gt_splicedseq_get(ss), splicedseq) == 0);
  gt_ensure(had_err, !gt_splicedseq_pos_is_border(ss, 0));
  gt_ensure(had_err,  gt_splicedseq_pos_is_border(ss, 1));
  gt_ensure(had_err, !gt_splicedseq_pos_is_border(ss, 2));
  gt_ensure(had_err, !gt_splicedseq_pos_is_border(ss, 3));
  gt_ensure(had_err, !gt_splicedseq_pos_is_border(ss, 4));
  return had_err;
}

int gt_splicedseq_unit_test(GtError *err)
{
  Splicedseq *ss;
  int had_err = 0;
  gt_error_check(err);
  ss = gt_splicedseq_new();
  had_err = check_splicedseq(ss, err);
  if (!had_err) {
    gt_splicedseq_reset(ss);
    had_err = check_splicedseq(ss, err);
  }
  gt_splicedseq_delete(ss);
  return had_err;
}

void gt_splicedseq_delete(Splicedseq *ss)
{
  if (!ss) return;
  gt_str_delete(ss->splicedseq);
  gt_array_delete(ss->positionmapping);
  gt_free(ss);
}
