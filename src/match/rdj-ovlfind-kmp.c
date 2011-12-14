/*
  Copyright (c) 2009-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2009-2011 Center for Bioinformatics, University of Hamburg

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
#include "core/ma.h"
#include "core/minmax.h"
#include "core/assert_api.h"
#include "match/rdj-ovlfind-kmp.h"
/* for unit test: */
#include "core/array_api.h"
#include "core/ensure.h"

gt_kmp_t* gt_kmp_preproc(const char *seq, unsigned long seqlen)
{
  gt_kmp_t* pi;
  gt_kmp_t k;
  unsigned long q;

  pi = gt_malloc(sizeof (gt_kmp_t) * seqlen);
  pi[0] = 0;
  k = 0;
  for (q = 1UL; q < seqlen; q++)
  {
    while ((k > 0) && (seq[k] != seq[q])) k = pi[k-1];
    if (seq[k] == seq[q]) k++;
    gt_assert(k < GT_KMP_MAX);
    pi[q] = k;
  }
  return pi;
}

static inline gt_kmp_t findnextmatch(const char *a, const char *b,
    unsigned long from, unsigned long to, const gt_kmp_t *pi)
{
  gt_kmp_t q;
  unsigned long i;
  q = 0;
  for (i = from; i <= to; i++)
  {
    while ((q > 0) && (a[i] != b[q]))
      q = pi[q - 1];
    if (a[i] == b[q])
      q++;
  }
  return q;
}

static inline void spmfind_kmp(const char *a, unsigned long alen,
    const char *b, unsigned long blen, const gt_kmp_t *pi,
    unsigned long min_length, bool find_submaximal, bool self_comparison,
    void(*proc)(unsigned long, bool, void*), bool direction,
    void* procdata)
{
  gt_kmp_t q;
  unsigned long max_matchlen;
  max_matchlen = MIN(alen, blen);
  if (self_comparison)
    max_matchlen -= 1;
  q = findnextmatch(a, b, alen - max_matchlen, alen - 1, pi);
  if ((unsigned long)q >= min_length)
    proc((unsigned long)q, direction, procdata);
  if (find_submaximal)
    while (q > 0)
    {
      q = findnextmatch(a, b, alen - q + 1, alen - 1, pi);
      if ((unsigned long)q >= min_length)
        proc((unsigned long)q, direction, procdata);
    }
}

static inline bool contfind_kmp(const char *a, unsigned long alen,
    const gt_kmp_t *pi, const char *b, unsigned long blen)
{
  gt_kmp_t q;
  unsigned long i;
  q = 0;
  gt_assert(alen < blen);
  for (i = 0; i < blen; i++)
  {
    while ((q > 0) && (a[q] != b[i]))
      q = pi[q - 1];
    if (a[q] == b[i])
      q++;
    if ((unsigned long)q == alen)
      return true;
  }
  return false;
}

GtContfind gt_ovlfind_kmp(const char *u, unsigned long u_length,
    const gt_kmp_t *u_pi, const char *v, unsigned long v_length,
    const gt_kmp_t *v_pi, GtOvlfindMode m, unsigned long min_length,
    bool find_nonmaximal, void(*spmproc) (unsigned long /* overlap length */,
    bool /* true if suffix of u == prefix of v, false if prefix of
    u == suffix of v */, void* /* spmprocdata */), void* spmprocdata)
{
  GtContfind retval = GT_CONTFIND_OFF;
  bool self_comparison;

  gt_assert(u != NULL);
  gt_assert(u_length > 0);
  gt_assert(u_pi != NULL);
  self_comparison = (v == NULL);
  gt_assert(self_comparison || v_length > 0);
  gt_assert(self_comparison || v_pi != NULL);
  gt_assert(!self_comparison || m == GT_OVLFIND_SPM || m == GT_OVLFIND_ALL);

  if (self_comparison && m == GT_OVLFIND_ALL)
    retval = GT_CONTFIND_EQ;
  if (m != GT_OVLFIND_SPM && !self_comparison)
  {
    if (u_length == v_length)
      retval = memcmp(u, v, (size_t)u_length) == 0
               ? GT_CONTFIND_EQ : GT_CONTFIND_NO;
    else if (u_length < v_length)
      retval = contfind_kmp(u, u_length, u_pi, v, v_length)
               ? GT_CONTFIND_U : GT_CONTFIND_NO;
    else /* u_length > v_length */
      retval = contfind_kmp(v, v_length, v_pi, u, u_length)
               ? GT_CONTFIND_V : GT_CONTFIND_NO;
    if (m == GT_OVLFIND_PROPER_SPM && retval != GT_CONTFIND_NO)
      return retval;
  }
  if (m != GT_OVLFIND_CNT)
  {
    if (self_comparison)
    {
      spmfind_kmp(u, u_length, u, u_length, u_pi, min_length, find_nonmaximal,
                 self_comparison, spmproc, true, spmprocdata);
    }
    else
    {
      spmfind_kmp(u, u_length, v, v_length, v_pi, min_length, find_nonmaximal,
                 self_comparison, spmproc, true, spmprocdata);
      spmfind_kmp(v, v_length, u, u_length, u_pi, min_length, find_nonmaximal,
                 self_comparison, spmproc, false, spmprocdata);
    }
  }
  return retval;
}

/*--------------------------   UNIT TEST   --------------------------*/

int gt_kmp_preproc_unit_test(GtError *err)
{
  int had_err = 0;
  int i;
  gt_kmp_t *pi, expected_pi[10] = {(gt_kmp_t)0, (gt_kmp_t)0, (gt_kmp_t)1,
    (gt_kmp_t)2, (gt_kmp_t)3, (gt_kmp_t)4, (gt_kmp_t)5, (gt_kmp_t)6,
    (gt_kmp_t)0, (gt_kmp_t)1};
  pi = gt_kmp_preproc("ababababca", 10UL);
  for (i = 0; i < 10; i++)
    gt_ensure(had_err, pi[i] == expected_pi[i]);
  gt_free(pi);
  return had_err;
}

struct GtOvlfindKmpResult { bool u_suffix; unsigned long length; };

static
void ovlfind_kmp_test_save(unsigned long length, bool u_suffix, void *a)
{
  struct GtOvlfindKmpResult r = {u_suffix, length};
  gt_array_add((GtArray*)a,r);
}

#define GT_OVLFIND_KMP_EXPECT_RESULT(N,U_SUF,LEN) \
        if (!had_err) r = gt_array_get(a, (N)); \
        gt_ensure(had_err, r->u_suffix == (U_SUF)); \
        gt_ensure(had_err, r->length == (LEN))

int gt_ovlfind_kmp_unit_test(GtError *err)
{
  int had_err = 0;
  GtArray *a;
  struct GtOvlfindKmpResult *r;
  GtContfind retval;
  gt_kmp_t *u_pi, *v_pi;

  /*@i1@*/ gt_error_check(err);

  had_err = gt_kmp_preproc_unit_test(err);
  if (had_err != 0)
    return had_err;

  a = gt_array_new(sizeof (struct GtOvlfindKmpResult));

  /* u suffix == v prefix */
  if (!had_err)
  {
    gt_array_reset(a);
    u_pi = gt_kmp_preproc("aacgcacctg", 10UL);
    v_pi = gt_kmp_preproc("acctgatttc", 10UL);
    retval = gt_ovlfind_kmp("aacgcacctg", 10UL, u_pi, "acctgatttc", 10UL, v_pi,
        GT_OVLFIND_PROPER_SPM, 1UL, false, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_NO);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    GT_OVLFIND_KMP_EXPECT_RESULT(0UL, true, 5UL);
    gt_free(u_pi);
    gt_free(v_pi);
  }
  /* v suffix == u prefix */
  if (!had_err)
  {
    gt_array_reset(a);
    u_pi = gt_kmp_preproc("atccgtgacgtg", 12UL);
    v_pi = gt_kmp_preproc("aagaagaatccg", 12UL);
    retval = gt_ovlfind_kmp("atccgtgacgtg", 12UL, u_pi, "aagaagaatccg", 12UL,
        v_pi, GT_OVLFIND_ALL, 1UL, false, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_NO);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    GT_OVLFIND_KMP_EXPECT_RESULT(0UL, false, 5UL);
    gt_free(u_pi);
    gt_free(v_pi);
  }
  /* no overlap */
  if (!had_err)
  {
    gt_array_reset(a);
    u_pi = gt_kmp_preproc("aac", 3UL);
    v_pi = gt_kmp_preproc("tgc", 3UL);
    retval = gt_ovlfind_kmp("aac", 3UL, u_pi, "tgc", 3UL, v_pi,
        GT_OVLFIND_PROPER_SPM, 1UL, false, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_NO);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
    gt_free(u_pi);
    gt_free(v_pi);
  }
  /* u suffix of v */
  if (!had_err)
  {
    gt_array_reset(a);
    u_pi = gt_kmp_preproc("acagc", 5UL);
    v_pi = gt_kmp_preproc("gtacagc", 7UL);
    retval = gt_ovlfind_kmp("acagc", 5UL, u_pi, "gtacagc", 7UL, v_pi,
        GT_OVLFIND_SPM, 1UL, false, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
    GT_OVLFIND_KMP_EXPECT_RESULT(0UL, false, 5UL);

    gt_array_reset(a);
    retval = gt_ovlfind_kmp("acagc", 5UL, u_pi, "gtacagc", 7UL, v_pi,
        GT_OVLFIND_PROPER_SPM, 1UL, false, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_U);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
    gt_array_reset(a);

    retval = gt_ovlfind_kmp("acagc", 5UL, u_pi, "gtacagc", 7UL, v_pi,
        GT_OVLFIND_CNT, 1UL, false, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
    gt_ensure(had_err, retval == GT_CONTFIND_U);

    gt_array_reset(a);
    retval = gt_ovlfind_kmp("acagc", 5UL, u_pi, "gtacagc", 7UL, v_pi,
        GT_OVLFIND_ALL, 1UL, false, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    gt_ensure(had_err, retval == GT_CONTFIND_U);
    GT_OVLFIND_KMP_EXPECT_RESULT(0UL, false, 5UL);
    gt_free(u_pi);
    gt_free(v_pi);
  }
  /* v suffix of u */
  if (!had_err)
  {
    gt_array_reset(a);
    u_pi = gt_kmp_preproc("gtacagc", 7UL);
    v_pi = gt_kmp_preproc("acagc", 5UL);

    retval = gt_ovlfind_kmp("gtacagc", 7UL, u_pi, "acagc", 5UL, v_pi,
        GT_OVLFIND_SPM, 1UL, false, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    GT_OVLFIND_KMP_EXPECT_RESULT(0UL, true, 5UL);

    gt_array_reset(a);
    retval = gt_ovlfind_kmp("gtacagc", 7UL, u_pi, "acagc", 5UL, v_pi,
        GT_OVLFIND_PROPER_SPM, 1UL, false, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_V);
    gt_ensure(had_err, gt_array_size(a) == 0UL);

    gt_free(u_pi);
    gt_free(v_pi);
  }
  /* u prefix of v */
  if (!had_err)
  {
    gt_array_reset(a);
    u_pi = gt_kmp_preproc("ctat", 4UL);
    v_pi = gt_kmp_preproc("ctatacagg", 9UL);
    retval = gt_ovlfind_kmp("ctat", 4UL, u_pi, "ctatacagg", 9UL, v_pi,
        GT_OVLFIND_SPM, 1UL, false, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    GT_OVLFIND_KMP_EXPECT_RESULT(0UL, true, 4UL);

    gt_array_reset(a);
    retval = gt_ovlfind_kmp("ctat", 4UL, u_pi, "ctatacagg", 9UL, v_pi,
        GT_OVLFIND_PROPER_SPM, 1UL, false, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_U);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
    gt_free(u_pi);
    gt_free(v_pi);
  }
  /* v prefix of u */
  if (!had_err)
  {
    gt_array_reset(a);
    u_pi = gt_kmp_preproc("ctatacagg", 9UL);
    v_pi = gt_kmp_preproc("ctat", 4UL);
    retval = gt_ovlfind_kmp("ctatacagg", 9UL, u_pi, "ctat", 4UL, v_pi,
        GT_OVLFIND_SPM, 1UL, false, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    GT_OVLFIND_KMP_EXPECT_RESULT(0UL, false, 4UL);

    gt_array_reset(a);
    retval = gt_ovlfind_kmp("ctatacagg", 9UL, u_pi, "ctat", 4UL, v_pi,
        GT_OVLFIND_PROPER_SPM, 1UL, false, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_V);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
    gt_free(u_pi);
    gt_free(v_pi);
  }
  /* identical sequences */
  if (!had_err)
  {
    gt_array_reset(a);
    u_pi = gt_kmp_preproc("acagc", 5UL);
    retval = gt_ovlfind_kmp("acagc", 5UL, u_pi, "acagc", 5UL, u_pi,
        GT_OVLFIND_SPM, 1UL, false, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
    gt_ensure(had_err, gt_array_size(a) == 2UL);
    GT_OVLFIND_KMP_EXPECT_RESULT(0UL, true, 5UL);
    GT_OVLFIND_KMP_EXPECT_RESULT(1UL, false, 5UL);

    gt_array_reset(a);
    retval = gt_ovlfind_kmp("acagc", 5UL, u_pi, "acagc", 5UL, u_pi,
        GT_OVLFIND_PROPER_SPM, 1UL, false, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_EQ);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
    gt_free(u_pi);
  }
  /* find_nonmaximal */
  if (!had_err)
  {
    gt_array_reset(a);
    u_pi = gt_kmp_preproc("aacagtagtagt", 12UL);
    v_pi = gt_kmp_preproc("agtagtagttaa", 12UL);
    retval = gt_ovlfind_kmp("aacagtagtagt", 12UL, u_pi, "agtagtagttaa", 12UL,
        v_pi, GT_OVLFIND_SPM, 1UL, false, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
    gt_ensure(had_err, gt_array_size(a) == 2UL);
    GT_OVLFIND_KMP_EXPECT_RESULT(0UL, true, 9UL);
    GT_OVLFIND_KMP_EXPECT_RESULT(1UL, false, 2UL);
    gt_array_reset(a);
    retval = gt_ovlfind_kmp("aacagtagtagt", 12UL, u_pi, "agtagtagttaa", 12UL,
        v_pi, GT_OVLFIND_SPM, 1UL, true, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 5UL);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
    gt_free(u_pi);
    gt_free(v_pi);
  }
  /* min_length */
  if (!had_err)
  {
    gt_array_reset(a);
    u_pi = gt_kmp_preproc("aggaccagtagt", 12UL);
    v_pi = gt_kmp_preproc("agtagttactac", 12UL);
    retval = gt_ovlfind_kmp("aggaccagtagt", 12UL, u_pi, "agtagttactac", 12UL,
        v_pi, GT_OVLFIND_SPM, 1UL, true, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
    gt_ensure(had_err, gt_array_size(a) == 2UL);
    gt_array_reset(a);
    retval = gt_ovlfind_kmp("aggaccagtagt", 12UL, u_pi, "agtagttactac", 12UL,
        v_pi, GT_OVLFIND_SPM, 4UL, true, ovlfind_kmp_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
    gt_free(u_pi);
    gt_free(v_pi);
  }
  gt_array_delete(a);
  return had_err;
}
