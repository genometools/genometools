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

#include "core/ma.h"
#include "core/minmax.h"
#include "core/unused_api.h"
#include "core/log.h"
#include "match/rdj-ovlfind-dp.h"
/* unit test: */
#include <string.h>
#include "core/ensure.h"
#include "core/array_api.h"

typedef struct
{
  unsigned long edist, start_i, start_j;
} GtOvlfindDpCell;

/* true if u == v (approximatively) */
static inline
bool eqlfind_dp(unsigned long max_edist, GtOvlfindDpCell *mn_cell)
{
  return mn_cell->edist <= max_edist &&
         mn_cell->start_i == 0 &&
         mn_cell->start_j == 0;
}

/* true if u is approximatively contained in v */
static inline
bool u_contfind_dp(unsigned long max_edist, unsigned long n,
                  GtOvlfindDpCell **lastrow)
{
  unsigned long j;
  for (j = 0; j <= n; j++)
    if ((*lastrow)[j].edist <= max_edist &&
        (*lastrow)[j].start_i == 0)
      return true;
  return false;
}

/* true if v is approximatively contained in u */
static inline
bool v_contfind_dp(unsigned long max_edist, unsigned long m,
                  GtOvlfindDpCell **lastcol)
{
  unsigned long i;
  for (i = 0; i <= m; i++)
    if ((*lastcol)[i].edist <= max_edist &&
        (*lastcol)[i].start_j == 0)
      return true;
  return false;
}

/* find approximate suffix-prefix matches, with the suffix coming from u */
static inline void u_smpfind_dp(unsigned long min_length,
    unsigned long max_edist, bool find_submaximal,
    unsigned long m, unsigned long n, GtOvlfindDpCell **lastrow,
    void (*proc)(unsigned long, unsigned long, unsigned long, bool,
    void*), void* procdata)
{
  unsigned long j, suffix_length;

  gt_assert(m > 0);
  gt_assert(n > 0);
  gt_assert(lastrow != NULL);
  for (j = n; j >= min_length; j--)
  {
    if ((*lastrow)[j].edist <= max_edist &&
        (*lastrow)[j].start_j == 0)
    {
      suffix_length = m - (*lastrow)[j].start_i;
      if (suffix_length >= min_length)
        proc(suffix_length, j /* prefix length */,
             (*lastrow)[j].edist, true, procdata);
      if (!find_submaximal) break;
    }
  }
}

/* find approximate suffix-prefix matches, with the suffix coming from v */
static inline
void v_smpfind_dp(unsigned long min_length, unsigned long max_edist,
                  bool find_submaximal, unsigned long m,
                  unsigned long n, GtOvlfindDpCell **lastcol,
                  void (*proc)(unsigned long, unsigned long,
                               unsigned long, bool, void*),
                  void* procdata)
{
  unsigned long i, prefix_length;

  gt_assert(m > 0);
  gt_assert(n > 0);
  gt_assert(lastcol != NULL);
  for (i = m; i >= min_length; i--)
  {
    if ((*lastcol)[i].edist <= max_edist &&
        (*lastcol)[i].start_i == 0)
    {
      prefix_length = n - (*lastcol)[i].start_j;
      if (prefix_length >= min_length)
        proc(i /* suffix length */, prefix_length,
             (*lastcol)[i].edist, false, procdata);
      if (!find_submaximal) break;
    }
  }
}

static inline
void run_dp(const char *u, unsigned long m, const char *v, unsigned long n,
            GtOvlfindDpCell **row, GtOvlfindDpCell **col)
{
  unsigned long i, j, cost_r, cost_d, cost_i;
  GtOvlfindDpCell cell_r, next_cell_r;

  /* first row */
  for (i = 0; i <= m; i++)
  {
    (*col)[i].start_i = i;
    (*col)[i].start_j = 0;
    (*col)[i].edist = 0;
  }
  (*row)[0] = (*col)[m];
  /* following (*row)s */
  for (j = 1UL; j <= n; j++)
  {
    (*col)[0].start_j = j;
    cell_r.start_i = 0;
    cell_r.start_j = j - 1;
    cell_r.edist = 0;
    for (i = 1UL; i <= m; i++)
    {
      next_cell_r = (*col)[i];
      cost_i = (*col)[i].edist + 1; /* cost(e --> v[j]) */
      cost_d = (*col)[i - 1].edist + 1; /* cost(u[i] --> e)*/
      cost_r = cell_r.edist + ((u[i - 1] == v[j - 1]) ? 0 : 1);
      if (cost_r < cost_i && cost_r < cost_d)
      {
        /* replacement */
        (*col)[i].start_i = cell_r.start_i;
        (*col)[i].start_j = cell_r.start_j;
        (*col)[i].edist = cost_r;
      }
      else if (cost_i < cost_d)
      {
        /* insertion */
        /* start_i and start_j already correct */
        (*col)[i].edist = cost_i;
      }
      else
      {
        /* deletion */
        (*col)[i].start_i = (*col)[i-1].start_i;
        (*col)[i].start_j = (*col)[i-1].start_j;
        (*col)[i].edist = cost_d;
      }
      cell_r = next_cell_r;
    }
    (*row)[j] = (*col)[m];
  }
}

GT_UNUSED static void print_row_and_col(GtOvlfindDpCell *row, unsigned long m,
    GtOvlfindDpCell *col, unsigned long n)
{
  unsigned long i, j;

  for (i = 0; i <= m; i++)
    gt_log_log("col[%lu]=((%lu,%lu),%lu)\n", i, col[i].start_i,
        col[i].start_j, col[i].edist);

  for (j = 0; j <= n; j++)
    gt_log_log("row[%lu]=((%lu,%lu),%lu)\n", j, row[j].start_i,
        row[j].start_j, row[j].edist);
}

GtContfind gt_ovlfind_dp(const char *u, unsigned long m,
                         const char *v, unsigned long n,
                         double max_error, GtOvlfindMode mode,
                         unsigned long min_length, bool find_submaximal,
                         void (*smpproc)
                           (unsigned long /* length on u */,
                            unsigned long /* length on v */,
                            unsigned long /* unit edit distance */,
                            bool /* true if suffix comes from u,
                                    false if suffix comes from v */,
                            void* /* procdata */),
                         void* smpprocdata)
{
  GtOvlfindDpCell *row, *col;
  unsigned long max_edist;
  GtContfind retval = GT_CONTFIND_OFF;
  bool self_comparison;

  gt_assert(u != NULL);
  gt_assert(m > 0);
  self_comparison = (v == NULL);
  gt_assert(self_comparison || n > 0);
  gt_assert(!self_comparison || mode == GT_OVLFIND_SPM ||
      mode == GT_OVLFIND_ALL);

  if (self_comparison)
  {
    max_edist = (unsigned long)(max_error * m);
    row = gt_malloc(sizeof (GtOvlfindDpCell) * (m + 1));
    col = gt_malloc(sizeof (GtOvlfindDpCell) * (m + 1));

    run_dp(u, m, u, m, &row, &col);
  }
  else
  {
    max_edist = (unsigned long)(max_error * MAX(m,n));
    row = gt_malloc(sizeof (GtOvlfindDpCell) * (n + 1));
    col = gt_malloc(sizeof (GtOvlfindDpCell) * (m + 1));

    run_dp(u, m, v, n, &row, &col);
    /* decomment to show row and col to debug log: */
    /* print_row_and_col(row, m, col, n); */
  }

  if (self_comparison && mode == GT_OVLFIND_ALL)
    retval = GT_CONTFIND_EQ;
  if (mode != GT_OVLFIND_SPM && !self_comparison)
  {
    if (eqlfind_dp(max_edist, row+n))
      retval = GT_CONTFIND_EQ;
    else if (v_contfind_dp(max_edist, m, &col))
      retval = GT_CONTFIND_V;
    else if (u_contfind_dp(max_edist, n, &row))
      retval = GT_CONTFIND_U;
    else
      retval = GT_CONTFIND_NO;
  }
  if (mode == GT_OVLFIND_SPM || mode == GT_OVLFIND_ALL ||
      (mode == GT_OVLFIND_PROPER_SPM && retval == GT_CONTFIND_NO))
  {
    if (self_comparison)
    {
      u_smpfind_dp(min_length, max_edist, find_submaximal, m, m-1, &row,
                   smpproc, smpprocdata);
    }
    else
    {
      u_smpfind_dp(min_length, max_edist, find_submaximal, m, n, &row,
                   smpproc, smpprocdata);
      v_smpfind_dp(min_length, max_edist, find_submaximal, m, n, &col,
                   smpproc, smpprocdata);
    }
  }

  gt_free(row);
  gt_free(col);
  return retval;
}

/*--------------------------   UNIT TEST   --------------------------*/

struct GtOvlfindDpResult { unsigned long ulen, vlen, edist; bool dir; };

static
void ovlfind_dp_test_save(unsigned long ulen, unsigned long vlen,
                          unsigned long edist, bool dir, void *a)
{
  struct GtOvlfindDpResult r = {ulen, vlen, edist, dir};
  gt_array_add((GtArray*)a,r);
}

#define GT_OVLFIND_DP_EXPECT_RESULT(N,ULEN,VLEN,EDST,DIR)        \
        if (!had_err) r = gt_array_get(a, (N));                  \
        gt_ensure(had_err, r->ulen == (ULEN));                      \
        gt_ensure(had_err, r->vlen == (VLEN));                      \
        gt_ensure(had_err, r->edist == (EDST));                     \
        gt_ensure(had_err, r->dir == (DIR))

int gt_ovlfind_dp_unit_test(GtError *err)
{
  int had_err = 0;
  GtArray *a;
  struct GtOvlfindDpResult *r;
  GtContfind retval;

  /*@i1@*/ gt_error_check(err);
  a = gt_array_new(sizeof (struct GtOvlfindDpResult));

  /* u suffix == v prefix */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_dp("aacgcacctg", 10UL, "acctgatttc", 10UL, 0.0,
                           GT_OVLFIND_PROPER_SPM, 1UL, false,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_NO);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    GT_OVLFIND_DP_EXPECT_RESULT(0UL, 5UL, 5UL, 0UL, true);
  }
  /* v suffix == u prefix */
  if (!had_err)
  {
    gt_array_reset(a);

    retval = gt_ovlfind_dp("atccgtgacgtg", 12UL, "aagaagaatccg", 12UL, 0.0,
                           GT_OVLFIND_ALL, 1UL, false,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_NO);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    GT_OVLFIND_DP_EXPECT_RESULT(0UL, 5UL, 5UL, 0UL, false);
  }
  /* no overlap */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_dp("aac", 3UL, "tgc", 3UL, 0.0,
                           GT_OVLFIND_PROPER_SPM, 1UL, false,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_NO);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
  }
  /* u suffix of v */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_dp("acagc", 5UL, "gtacagc", 7UL, 0.0,
                           GT_OVLFIND_SPM, 1UL, false,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    GT_OVLFIND_DP_EXPECT_RESULT(0UL, 5UL, 5UL, 0UL, false);

    gt_array_reset(a);
    retval = gt_ovlfind_dp("acagc", 5UL, "gtacagc", 7UL, 0.0,
                           GT_OVLFIND_PROPER_SPM, 1UL, false,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_U);
    gt_ensure(had_err, gt_array_size(a) == 0UL);

    gt_array_reset(a);
    retval = gt_ovlfind_dp("acagc", 5UL, "gtacagc", 7UL, 0.0,
                           GT_OVLFIND_CNT, 1UL, false,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
    gt_ensure(had_err, retval == GT_CONTFIND_U);

    gt_array_reset(a);
    retval = gt_ovlfind_dp("acagc", 5UL, "gtacagc", 7UL, 0.0,
                           GT_OVLFIND_ALL, 1UL, false,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    gt_ensure(had_err, retval == GT_CONTFIND_U);
    GT_OVLFIND_DP_EXPECT_RESULT(0UL, 5UL, 5UL, 0UL, false);
  }
  /* v suffix of u */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_dp("gtacagc", 7UL, "acagc", 5UL, 0.0,
                           GT_OVLFIND_SPM, 1UL, false,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    GT_OVLFIND_DP_EXPECT_RESULT(0UL, 5UL, 5UL, 0UL, true);

    gt_array_reset(a);
    retval = gt_ovlfind_dp("gtacagc", 7UL, "acagc", 5UL, 0.0,
                           GT_OVLFIND_PROPER_SPM, 1UL, false,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_V);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
  }
  /* u prefix of v */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_dp("ctat", 4UL, "ctatacagg", 9UL, 0.0,
                           GT_OVLFIND_SPM, 1UL, false,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    GT_OVLFIND_DP_EXPECT_RESULT(0UL, 4UL, 4UL, 0UL, true);

    gt_array_reset(a);
    retval = gt_ovlfind_dp("ctat", 4UL, "ctatacagg", 9UL, 0.0,
                           GT_OVLFIND_PROPER_SPM, 1UL, false,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_U);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
  }
  /* v prefix of u */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_dp("ctatacagg", 9UL, "ctat", 4UL, 0.0,
                           GT_OVLFIND_SPM, 1UL, false,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    GT_OVLFIND_DP_EXPECT_RESULT(0UL, 4UL, 4UL, 0UL, false);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);

    gt_array_reset(a);
    retval = gt_ovlfind_dp("ctatacagg", 9UL, "ctat", 4UL, 0.0,
                           GT_OVLFIND_PROPER_SPM, 1UL, false,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_V);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
  }
  /* identical sequences */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_dp("acagc", 5UL, "acagc", 5UL, 0.0,
                           GT_OVLFIND_SPM, 1UL, false,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 2UL);
    GT_OVLFIND_DP_EXPECT_RESULT(0UL, 5UL, 5UL, 0UL, true);
    GT_OVLFIND_DP_EXPECT_RESULT(1UL, 5UL, 5UL, 0UL, false);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);

    gt_array_reset(a);
    retval = gt_ovlfind_dp("acagc", 5UL, "acagc", 5UL, 0.0,
                           GT_OVLFIND_PROPER_SPM, 1UL, false,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
    gt_ensure(had_err, retval == GT_CONTFIND_EQ);
  }
  /* find_nonmaximal */
  if (!had_err)
  {
    gt_array_reset(a);

    retval = gt_ovlfind_dp("aacagtagtagt", 12UL, "agtagtagttaa", 12UL, 0.0,
                           GT_OVLFIND_SPM, 1UL, false,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
    gt_ensure(had_err, gt_array_size(a) == 2UL);
    GT_OVLFIND_DP_EXPECT_RESULT(0UL, 9UL, 9UL, 0UL, true);
    GT_OVLFIND_DP_EXPECT_RESULT(1UL, 2UL, 2UL, 0UL, false);

    gt_array_reset(a);
    retval = gt_ovlfind_dp("aacagtagtagt", 12UL, "agtagtagttaa", 12UL, 0.0,
                           GT_OVLFIND_SPM, 1UL, true,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 5UL);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
  }
  /* min_length */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_dp("aggaccagtagt", 12UL, "agtagttactac", 12UL, 0.0,
                           GT_OVLFIND_SPM, 1UL, true,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 2UL);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);

    gt_array_reset(a);
    retval = gt_ovlfind_dp("aggaccagtagt", 12UL, "agtagttactac", 12UL, 0.0,
                           GT_OVLFIND_SPM, 4UL, true,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
  }
  /* max_error */
  if (!had_err)
  {
    gt_array_reset(a);
    retval = gt_ovlfind_dp("gaccactagt", 10UL, "agtagttact", 10UL, 0.0,
                           GT_OVLFIND_SPM, 6UL, true,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 0UL);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);

    gt_array_reset(a);
    retval = gt_ovlfind_dp("gaccactagt", 10UL, "agtagttact", 10UL, 0.1,
                           GT_OVLFIND_SPM, 6UL, true,
                           ovlfind_dp_test_save, a);
    gt_ensure(had_err, gt_array_size(a) == 1UL);
    gt_ensure(had_err, retval == GT_CONTFIND_OFF);
  }
  gt_array_delete(a);
  return had_err;
}
