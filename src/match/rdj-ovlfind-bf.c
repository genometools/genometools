/*
  Copyright (c) 2010-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg

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

#include <string.h> /* memcmp */
#include "core/ma.h"
#include "core/assert_api.h"
#include "core/minmax.h"
#include "match/rdj-ovlfind-bf.h"
/* for unit test: */
#include "core/array_api.h"
#include "core/ensure.h"

static inline bool contfind_bf(const char *a, unsigned long alen,
    const char *b, unsigned long blen)
{
  unsigned long offset;
  gt_assert(alen < blen);
  for (offset = 0; offset <= blen-alen; offset++)
    if (memcmp(a, b+offset, (size_t)alen) == 0) return true;
  return false;
}

static inline void spmfind_bf(const char *a, unsigned long alen,
    const char *b, unsigned long blen, unsigned long minlen,
    bool find_nonmaximal, bool self_comparison,
    void(*proc)(unsigned long, bool, void*), bool direction,
    void* procdata)
{
  unsigned long len, from;
  from = MIN(alen, blen);
  if (self_comparison)
    from -= 1;
  for (len = from; len >= minlen; len--)
    if (memcmp(a+alen-len, b, (size_t)len) == 0)
    {
      proc(len, direction, procdata);
      if (!find_nonmaximal) break;
    }
}

GtContfind gt_ovlfind_bf(const char *u, unsigned long u_length,
    const char *v, unsigned long v_length, GtOvlfindMode m,
    unsigned long min_length, bool find_nonmaximal,
    void(*spmproc) (unsigned long /* overlap length */,
        bool /* true if suffix of u == prefix of v;
            false if prefix of u == suffix of v */,
        void* /* spmprocdata */), void* spmprocdata)
{
  GtContfind retval = GT_CONTFIND_OFF;
  bool self_comparison;

  gt_assert(u != NULL);
  gt_assert(u_length > 0);
  self_comparison = (v == NULL);
  gt_assert(self_comparison || v_length > 0);
  gt_assert(!self_comparison || m == GT_OVLFIND_SPM || m == GT_OVLFIND_ALL);

  if (self_comparison && m == GT_OVLFIND_ALL)
    retval = GT_CONTFIND_EQ;
  if (m != GT_OVLFIND_SPM && !self_comparison)
  {
    if (u_length == v_length)
      retval = memcmp(u, v, (size_t)u_length) == 0
               ? GT_CONTFIND_EQ : GT_CONTFIND_NO;
    else if (u_length < v_length)
      retval = contfind_bf(u, u_length, v, v_length)
               ? GT_CONTFIND_U : GT_CONTFIND_NO;
    else /* u_length > v_length */
      retval = contfind_bf(v, v_length, u, u_length)
               ? GT_CONTFIND_V : GT_CONTFIND_NO;
    if (m == GT_OVLFIND_PROPER_SPM && retval != GT_CONTFIND_NO)
      return retval;
  }
  if (m != GT_OVLFIND_CNT)
  {
    if (self_comparison)
    {
      spmfind_bf(u, u_length, u, u_length, min_length, find_nonmaximal,
                 self_comparison, spmproc, true, spmprocdata);
    }
    else
    {
      spmfind_bf(u, u_length, v, v_length, min_length, find_nonmaximal,
                 self_comparison, spmproc, true, spmprocdata);
      spmfind_bf(v, v_length, u, u_length, min_length, find_nonmaximal,
                 self_comparison, spmproc, false, spmprocdata);
    }
  }
  return retval;
}

/*--------------------------   UNIT TEST   --------------------------*/

struct GtOvlfindBfResult { bool u_suffix; unsigned long length; };

static
void ovlfind_bf_test_save(unsigned long length, bool u_suffix, void *a)
{
  struct GtOvlfindBfResult r = {u_suffix, length};
  gt_array_add((GtArray*)a,r);
}

#define GT_OVLFIND_BF_EXPECT_RESULT(N,U_SUF,LEN)        \
        if (!had_err) r = gt_array_get(a, (N));         \
        gt_ensure(had_err, r->u_suffix == (U_SUF));        \
        gt_ensure(had_err, r->length == (LEN))

int gt_ovlfind_bf_unit_test(GtError *err)
{
  int had_err = 0;
  GtArray *a;
  struct GtOvlfindBfResult *r;
  GtContfind retval;

  /*@i1@*/ gt_error_check(err);
  a = gt_array_new(sizeof (struct GtOvlfindBfResult));

  /* u suffix == v prefix */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_bf("aacgcacctg", 10UL, "acctgatttc", 10UL,
                           GT_OVLFIND_PROPER_SPM, 1UL, false,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_NO);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    GT_OVLFIND_BF_EXPECT_RESULT(0UL, true, 5UL);
  }
  /* v suffix == u prefix */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_bf("atccgtgacgtg", 12UL, "aagaagaatccg", 12UL,
                           GT_OVLFIND_ALL, 1UL, false,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_NO);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    GT_OVLFIND_BF_EXPECT_RESULT(0UL, false, 5UL);
  }
  /* no overlap */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_bf("aac", 3UL, "tgc", 3UL,
                           GT_OVLFIND_PROPER_SPM, 1UL, false,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_NO);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
  }
  /* u suffix of v */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_bf("acagc", 5UL, "gtacagc", 7UL,
                           GT_OVLFIND_SPM, 1UL, false,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
    GT_OVLFIND_BF_EXPECT_RESULT(0UL, false, 5UL);

    gt_array_reset(a);
    retval = gt_ovlfind_bf("acagc", 5UL, "gtacagc", 7UL,
                           GT_OVLFIND_PROPER_SPM, 1UL, false,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_U);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
    gt_array_reset(a);

    retval = gt_ovlfind_bf("acagc", 5UL, "gtacagc", 7UL,
                           GT_OVLFIND_CNT, 1UL, false,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
    gt_ensure(had_err, retval == GT_CONTFIND_U);

    gt_array_reset(a);
    retval = gt_ovlfind_bf("acagc", 5UL, "gtacagc", 7UL,
                           GT_OVLFIND_ALL, 1UL, false,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    gt_ensure(had_err, retval == GT_CONTFIND_U);
    GT_OVLFIND_BF_EXPECT_RESULT(0UL, false, 5UL);
  }
  /* v suffix of u */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_bf("gtacagc", 7UL, "acagc", 5UL,
                           GT_OVLFIND_SPM, 1UL, false,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
    GT_OVLFIND_BF_EXPECT_RESULT(0UL, true, 5UL);

    gt_array_reset(a);
    retval = gt_ovlfind_bf("gtacagc", 7UL, "acagc", 5UL,
                           GT_OVLFIND_PROPER_SPM, 1UL, false,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_V);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
  }
  /* u prefix of v */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_bf("ctat", 4UL, "ctatacagg", 9UL,
                           GT_OVLFIND_SPM, 1UL, false,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    GT_OVLFIND_BF_EXPECT_RESULT(0UL, true, 4UL);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);

    gt_array_reset(a);
    retval = gt_ovlfind_bf("ctat", 4UL, "ctatacagg", 9UL,
                           GT_OVLFIND_PROPER_SPM, 1UL, false,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
    gt_ensure(had_err, retval == GT_CONTFIND_U);
  }
  /* v prefix of u */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_bf("ctatacagg", 9UL, "ctat", 4UL,
                           GT_OVLFIND_SPM, 1UL, false,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    GT_OVLFIND_BF_EXPECT_RESULT(0UL, false, 4UL);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);

    gt_array_reset(a);
    retval = gt_ovlfind_bf("ctatacagg", 9UL, "ctat", 4UL,
                           GT_OVLFIND_PROPER_SPM, 1UL, false,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
    gt_ensure(had_err, retval == GT_CONTFIND_V);
  }
  /* identical sequences */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_bf("acagc", 5UL, "acagc", 5UL,
                           GT_OVLFIND_SPM, 1UL, false,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 2UL);
    GT_OVLFIND_BF_EXPECT_RESULT(0UL, true, 5UL);
    GT_OVLFIND_BF_EXPECT_RESULT(1UL, false, 5UL);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);

    gt_array_reset(a);
    retval = gt_ovlfind_bf("acagc", 5UL, "acagc", 5UL,
                           GT_OVLFIND_PROPER_SPM, 1UL, false,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
    gt_ensure(had_err, retval == GT_CONTFIND_EQ);
  }
  /* find_nonmaximal */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_bf("aacagtagtagt", 12UL, "agtagtagttaa", 12UL,
                           GT_OVLFIND_SPM, 1UL, false,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 2UL);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
    GT_OVLFIND_BF_EXPECT_RESULT(0UL, true, 9UL);
    GT_OVLFIND_BF_EXPECT_RESULT(1UL, false, 2UL);

    gt_array_reset(a);
    retval = gt_ovlfind_bf("aacagtagtagt", 12UL, "agtagtagttaa", 12UL,
                           GT_OVLFIND_SPM, 1UL, true,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 5UL);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
  }
  /* min_length */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_bf("aggaccagtagt", 12UL, "agtagttactac", 12UL,
                           GT_OVLFIND_SPM, 1UL, true,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 2UL);

    gt_array_reset(a);
    retval = gt_ovlfind_bf("aggaccagtagt", 12UL, "agtagttactac", 12UL,
                           GT_OVLFIND_SPM, 4UL, true,
                           ovlfind_bf_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
  }
  gt_array_delete(a);
  return had_err;
}
