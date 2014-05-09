/*
  Copyright (c) 2007      David Ellinghaus <d.ellinghaus@ikmb.uni-kiel.de>
  Copyright (c) 2012-2013 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2007-2014 Center for Bioinformatics, University of Hamburg

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

#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/ensure.h"
#include "core/log_api.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "extended/alignment.h"
#include "match/seqabstract.h"
#include "match/xdrop.h"

typedef struct
{
  int mis,
      ins,
      del,
      gcd;  /* greatest common divisor */
} GtXdropArbitrarydistances;

typedef struct
{
  GtWord row;
  unsigned char direction; /* one of the bits GT_XDROP_REPLACEMENTBIT,
                                              GT_XDROP_DELETIONBIT,
                                              GT_XDROP_INSERTIONBIT */
} GtXdropfrontvalue;

GT_DECLAREARRAYSTRUCT(GtXdropfrontvalue);

#define GT_XDROP_FRONTIDX(D,K)    ((GtUword) (D) * (D) + (D) + (K))

/*
  For each entry in the DP-matrix we store a single byte, and
  use the three rightmost bits to mark which edge in the edit distance
  graph to trace back.
*/

#define GT_XDROP_REPLACEMENTBIT   ((unsigned char) 1)
#define GT_XDROP_DELETIONBIT      (((unsigned char) 1) << 1)
#define GT_XDROP_INSERTIONBIT     (((unsigned char) 1) << 2)

/*
 The following function shows the matrix of the calculated fronts.
 */
/* CAUTION: fronts, that run over the matrix boundaries are not shown in
   the printed matrix.
 */
void gt_showfrontvalues(const GtArrayGtXdropfrontvalue *fronts,
                       GtWord distance,
                       const unsigned char *useqptr,
                       const unsigned char *vseqptr,
                       GtWord ulen,
                       GtWord vlen)
{
  GtUword l, max_l;
  GtWord i, j, k, d, filled = 0, integermax = MAX(ulen,vlen),
       integermin = -integermax;

  printf("frontvalues:\n");
  printf("        ");
  printf("%-3c ", vseqptr[0]);
  /* print vseq */
  max_l = MIN(fronts->nextfreeGtXdropfrontvalue,
              (GtUword) GT_XDROP_FRONTIDX(distance, distance));
  for (i = 1L; i < vlen; i++)
    printf("%-3c ", vseqptr[i]);

  for (i = 0; i <= ulen; i++) {
    printf("\n");
    /* print useq */
    if (i != 0)
      printf("%-3c ", useqptr[i - 1]);
    else
      printf("    ");

    for (j = 0; j <= vlen; j++) {
      d = distance + 1;
      k = i - j;
      for (l = 0; l < max_l; l++) {
        if (fronts->spaceGtXdropfrontvalue[l].row == integermin)
          continue;

        if (i == fronts->spaceGtXdropfrontvalue[l].row) {
          for (d = 0; d <= distance; d++) {
            if (k >= -d && k <= d && l == GT_XDROP_FRONTIDX(d, i - j)) {
#ifndef S_SPLINT_S
              printf("%-3"GT_WDS" ", d);
#else
              printf("%-3ld ", d);
#endif
              l = fronts->nextfreeGtXdropfrontvalue;
              filled++;
              break;
            }
          }
        }
      }
      if (d > distance)
        printf(".   ");

    }
  }
  printf("\n%.2f percent of matrix filled\n",
           (double) filled * 100.00 / ((ulen + 1) * (vlen + 1)));
}

static void
gt_calculatedistancesfromscores(const GtXdropArbitraryscores *arbitscores,
                                GtXdropArbitrarydistances *dist)
{
  int mat, mis, ins, del;

  /* if mat is odd double all scores */
  gt_assert(arbitscores->mat > 0);
  if (GT_MOD2((unsigned int) arbitscores->mat) > 0)
  {
    mat = arbitscores->mat * 2;
    mis = arbitscores->mis * 2;
    ins = arbitscores->ins * 2;
    del = arbitscores->del * 2;
  } else
  {
    mat = arbitscores->mat;
    mis = arbitscores->mis;
    ins = arbitscores->ins;
    del = arbitscores->del;
  }
  gt_assert(mat >= mis && mat/2 >= ins && mat/2 >= del);
  dist->gcd =
    (int) gt_gcd_uint(gt_gcd_uint((unsigned int) (mat-mis),
                                  (unsigned int) (mat/2-ins)),
                      (unsigned int) (mat/2-del));
  dist->mis = (mat - mis) / dist->gcd;
  dist->ins = (mat/2 - ins) / dist->gcd;
  dist->del = (mat/2 - del) / dist->gcd;
}

struct GtXdropresources
{
  const GtXdropArbitraryscores *arbitscores;
  GtXdropArbitrarydistances arbitdistances;
  GtArrayGtXdropfrontvalue fronts;
  GtArrayGtXdropscore big_t;
};

void gt_xdrop_resources_reset(GtXdropresources *res)
{
  res->fronts.nextfreeGtXdropfrontvalue = 0;
  res->big_t.nextfreeGtXdropscore = 0;
}

GtXdropresources *gt_xdrop_resources_new(const GtXdropArbitraryscores *scores)
{
  GtXdropresources *res = gt_malloc(sizeof *res);

  res->arbitscores = scores;
  GT_INITARRAY (&res->fronts, GtXdropfrontvalue);
  GT_INITARRAY (&res->big_t, GtXdropscore);
  gt_calculatedistancesfromscores(scores,&res->arbitdistances);
  return res;
}

void gt_xdrop_resources_delete(GtXdropresources *res)
{
  if (res != NULL)
  {
    GT_FREEARRAY (&res->fronts, GtXdropfrontvalue);
    GT_FREEARRAY (&res->big_t, GtXdropscore);
    gt_free(res);
  }
}

#define GT_XDROP_EVAL(K,D)\
        ((K) * res->arbitscores->mat/2 - (D) * res->arbitdistances.gcd)

#define GT_XDROP_SETDBACK(XDROPBELOWSCORE)\
        (XDROPBELOWSCORE + res->arbitscores->mat/2)/res->arbitdistances.gcd + 1

static GtWord gt_xdrop_frontvalue_get(const GtXdropresources *res, GtWord d,
                                      GtWord k)
{
  const GtUword frontidx = GT_XDROP_FRONTIDX(d, k);
  return res->fronts.spaceGtXdropfrontvalue[frontidx].row;
}

static void gt_xdrop_frontvalue_set(GtXdropresources *res, GtWord d, GtWord k,
                                    GtXdropfrontvalue value)
{
  const GtUword frontidx = GT_XDROP_FRONTIDX(d, k);

  if (frontidx >= res->fronts.allocatedGtXdropfrontvalue)
  {
    res->fronts.allocatedGtXdropfrontvalue = frontidx + 32UL;
    res->fronts.spaceGtXdropfrontvalue
      = gt_realloc_mem(res->fronts.spaceGtXdropfrontvalue,
                       sizeof (*res->fronts.spaceGtXdropfrontvalue) *
                       res->fronts.allocatedGtXdropfrontvalue,
                       __FILE__, __LINE__);
  }
  res->fronts.spaceGtXdropfrontvalue[frontidx] = value;
  res->fronts.nextfreeGtXdropfrontvalue = frontidx + 1;
}

void gt_evalxdroparbitscoresextend(bool forward,
                                   GtXdropbest *xdropbest,
                                   GtXdropresources *res,
                                   const GtSeqabstract *useq,
                                   const GtSeqabstract *vseq,
                                   GtXdropscore xdropbelowscore)
{
  const GtWord ulen = (GtWord) gt_seqabstract_length(useq),
               vlen = (GtWord) gt_seqabstract_length(vseq),
               end_k =
                 (GtWord) ulen - vlen, /* diagonal of endpoint (ulen, vlen) */
               integermax = (GtWord) MAX(ulen, vlen),
               integermin = -integermax,
               dback = GT_XDROP_SETDBACK(xdropbelowscore);
  GtWord idx,
         lbound,    /* diagonal lower bound */
         ubound,    /* diagonal upper bound */
         currd = 0, /* distance */
         k;         /* lbound - 1 <= k <= ubound + 1*/
  /*The following function calculates the maximal allowed number of
    generations with all front values equal minus infinity.*/
  const int allowedMININFINITYINTgenerations = MAX(MAX(res->arbitdistances.mis,
                                                       res->arbitdistances.ins),
                                                   res->arbitdistances.del) - 1;
  int currentMININFINITYINTgeneration = 0;
  GtXdropfrontvalue tmpfront;
  GtXdropscore bigt_tmp;        /* best score T' seen already */
  bool alwaysMININFINITYINT = true;

  gt_assert(ulen != 0 && vlen != 0);

  res->big_t.nextfreeGtXdropscore = 0;
  res->fronts.nextfreeGtXdropfrontvalue = 0;
  /* phase 0 */
  idx =  (GtWord) gt_seqabstract_lcp(forward, useq, vseq,
                                     forward ? 0 : (GtUword) ulen - 1,
                                     forward ? 0 : (GtUword) vlen - 1);
  /* alignment already finished */
  if (idx >= ulen || idx >= vlen) {
    lbound =  1L;
    ubound = -1L;
  }
  else {
    lbound = 0;
    ubound = 0;
  }
  tmpfront.row = (GtWord) idx;
  tmpfront.direction = (GtUchar) 0;   /* no predecessor */
  gt_xdrop_frontvalue_set(res,0,0,tmpfront);
  xdropbest->score = bigt_tmp = GT_XDROP_EVAL(idx + idx, 0);
  gt_assert(idx >= 0);
  xdropbest->ivalue = xdropbest->jvalue = (GtUword) idx;
  xdropbest->best_d = currd;
  xdropbest->best_k = 0;
  GT_STOREINARRAY (&res->big_t, GtXdropscore, 10, bigt_tmp);

  /* phase d > 0 */
  while (lbound <= ubound) {
    currd++;
    /* calculate fronts */
    for (k = lbound - 1; k <= ubound + 1; k++) {
      GtWord i = integermin, row;
      /* case 1 : DELETION-EDGE  */
      if (lbound < k &&
          currd - res->arbitdistances.del >= 0 &&
          -(currd - res->arbitdistances.del) <= k - 1 &&
          k - 1 <= currd - res->arbitdistances.del) {
        i = gt_xdrop_frontvalue_get(res, currd - res->arbitdistances.del, k-1)
            + 1;
        tmpfront.direction = GT_XDROP_DELETIONBIT;
      }
      /* case 2: REPLACEMENT-EDGE */
      if (lbound <= k &&
          k <= ubound &&
          currd - res->arbitdistances.mis >= 0 &&
          -(currd - res->arbitdistances.mis) <= k &&
          k <= currd - res->arbitdistances.mis) {
        row = gt_xdrop_frontvalue_get(res, currd - res->arbitdistances.mis, k)
              + 1;
        /* test, if case 1 has happened. */
        if (!(tmpfront.direction & GT_XDROP_DELETIONBIT) || row > i) {
          i = row;
          tmpfront.direction = GT_XDROP_REPLACEMENTBIT;
        }
      }
      /* case 3: INSERTION-EDGE */
      if (k < ubound &&
          currd - res->arbitdistances.ins >= 0 &&
          -(currd - res->arbitdistances.ins) <= k + 1 &&
           k + 1 <= currd - res->arbitdistances.ins) {
        row =
          gt_xdrop_frontvalue_get(res, currd - res->arbitdistances.ins, k+1);
        if (!(tmpfront.direction & (GT_XDROP_DELETIONBIT |
                                    GT_XDROP_REPLACEMENTBIT)) || row > i) {
          i = row;
          tmpfront.direction = GT_XDROP_INSERTIONBIT;
        }
      }
      /* if i = MINUSINFINITYINY or MINUSINFINITYINY + 1 */
      if (i < 0) {
        if (tmpfront.direction == (GtUchar) 0)
          alwaysMININFINITYINT = false;

        tmpfront.row = integermin;
      }
      else {
        GtWord j = i - k;
        const GtWord previousd = currd - dback;

        /* alignment score smaller than T - X */
        if (previousd > 0 &&
            res->big_t.spaceGtXdropscore != NULL &&
            GT_XDROP_EVAL (i + j, currd) <
            res->big_t.spaceGtXdropscore[previousd] - xdropbelowscore) {
          tmpfront.row = integermin;
        }
        else {
          if (k <= -currd || k >= currd ||
              (gt_xdrop_frontvalue_get(res, currd-1, k) < i &&
               i <= MIN(ulen, vlen + k))) {
            if (ulen > i && vlen > j) {
              GtUword lcp;
              gt_assert(forward || (ulen - 1 >= (GtWord) i &&
                                    vlen - 1 >= (GtWord) j));
              lcp = gt_seqabstract_lcp(forward, useq, vseq,
                                       (GtUword) (forward ? i : ulen - i - 1),
                                       (GtUword) (forward ? j : vlen - j - 1));
              i += lcp;
              j += lcp;
            }
            alwaysMININFINITYINT = false;
            tmpfront.row = i;
            if (GT_XDROP_EVAL(i + j, currd) > bigt_tmp) {
              xdropbest->score = bigt_tmp = GT_XDROP_EVAL(i + j, currd);
              gt_assert(i >= 0 && j >= 0);
              xdropbest->ivalue = (GtUword) i;
              xdropbest->jvalue = (GtUword) j;
              xdropbest->best_d = currd;
              xdropbest->best_k = k;
            }
          }
          else {
            alwaysMININFINITYINT = false;
            tmpfront.row = gt_xdrop_frontvalue_get(res,currd-1,k);
          }
        }
      }
      gt_xdrop_frontvalue_set(res, currd, k, tmpfront);
    }
    /* if all front values are integermin, alignment prematurely finished if
       allowedMININFINITYINTgenerations exceeded (full front has already ended
       at currd - currentMININFINITYINTgeneration). */
    if (alwaysMININFINITYINT) {
      currentMININFINITYINTgeneration++;
      if (currentMININFINITYINTgeneration > allowedMININFINITYINTgenerations)
        break;

    }
    else {
      currentMININFINITYINTgeneration = 0;
      alwaysMININFINITYINT = true;
    }
    GT_STOREINARRAY (&res->big_t, GtXdropscore, 10, bigt_tmp);
    /* fill out of bounds values of integermin
       needed for gt_showfrontvalues function */
    for (k = -currd; k < lbound - 1; k++) {
      tmpfront.row = integermin;
      gt_xdrop_frontvalue_set(res,currd,k,tmpfront);
    }
    for (k = ubound + 2; k <= currd; k++) {
      tmpfront.row = integermin;
      gt_xdrop_frontvalue_set(res,currd,k,tmpfront);
    }
    /* alignment finished */
    if (-currd <= end_k && end_k <= currd &&
        gt_xdrop_frontvalue_get(res,currd,end_k) == ulen)
      break;

    /* pruning lower bound
       lbound may decrease by one or increase/stays the same
       l <- min{k:R(d,k) > -inf} */
    for (k = lbound - 1; k <= ubound + 1; k++) {
      if (gt_xdrop_frontvalue_get(res,currd,k) > integermin) {
        lbound = k;
        break;
      }
    }
    /* pruning upper bound
       ubound may increase by one or decrease/stays the same
       u <- max{k:R(d,k) > -inf} */
    for (k = ubound + 1; k >= lbound - 1; k--) {
      if (gt_xdrop_frontvalue_get(res,currd,k) > integermin) {
        ubound = k;
        break;
      }
    }
    /* handling boundaries lower bound */
    for (k = 0; k >= lbound; k--) {
      if (gt_xdrop_frontvalue_get(res,currd,k) == vlen + k) {
        lbound = k;
        break;
      }
    }
    /* handling boundaries upper bound */
    for (k = 0; k <= ubound; k++) {
      if (gt_xdrop_frontvalue_get(res,currd,k) == ulen) {
        ubound = k;
        break;
      }
    }
  }
}

GtMultieoplist * gt_xdrop_backtrack(GtXdropresources *res,
                                    GtXdropbest *best)
{
  GtMultieoplist *meops = gt_multieoplist_new();
  GtUword idx, i;
  GtWord k = best->best_k,
       d = best->best_d,
       old_row = (GtWord) best->ivalue;
  GtXdropfrontvalue *fronts = res->fronts.spaceGtXdropfrontvalue,
                    currfront;
  gt_assert(best->ivalue && best->jvalue);

  idx = GT_XDROP_FRONTIDX(d, k);
  currfront = fronts[idx];

  while (d > 0) {
    if (currfront.direction == GT_XDROP_INSERTIONBIT) {
      d -= res->arbitdistances.ins; k++;
      idx = GT_XDROP_FRONTIDX(d, k);
      currfront = fronts[idx];
      for (i = 0; i < (GtUword) old_row - (currfront.row); ++i) {
        gt_multieoplist_add_match(meops);
      }
      gt_multieoplist_add_insertion(meops);
    }
    else if (currfront.direction == GT_XDROP_DELETIONBIT) {
      d -= res->arbitdistances.del; k--;
      idx = GT_XDROP_FRONTIDX(d, k);
      currfront = fronts[idx];
      for (i = 0; i < (GtUword) old_row - (currfront.row + 1); ++i) {
        gt_multieoplist_add_match(meops);
      }
      gt_multieoplist_add_deletion(meops);
    }
    else if (currfront.direction == GT_XDROP_REPLACEMENTBIT) {
      d -= res->arbitdistances.mis;
      idx = GT_XDROP_FRONTIDX(d, k);
      currfront = fronts[idx];
      for (i = 0; i < (GtUword) old_row - (currfront.row + 1); ++i) {
        gt_multieoplist_add_match(meops);
      }
      gt_multieoplist_add_mismatch(meops);
    }
    else {
      gt_assert(false && "this should not be reached");
    }

    gt_assert(currfront.row >= 0 && old_row >= currfront.row);
    old_row = currfront.row;
  }
  while (old_row > 0) {
    gt_multieoplist_add_match(meops);
    old_row--;
  }
  gt_assert(d == 0);
  return meops;
}

#define GT_XDROP_NUM_OF_TESTS 8
int gt_xdrop_unit_test(GT_UNUSED GtError *err)
{
  int had_err = 0, i, j, s;
  const GtUchar *strings[GT_XDROP_NUM_OF_TESTS] =
    {(const GtUchar*) "TTTTTTTTTTTTTTTAAAGGGTTTCCCAAAGGGTTTCCCTTTTTTTTTTTTTTT",
     (const GtUchar*) "TTTTTTTTTTTTTTTTTTTGGGGCCCCAAAATTTTTTTTTTTTTTT",
     (const GtUchar*) "TTTTTTTTTTTTTTTNNNNTTTTGGGGCCCCAAAATTTTTTTTTTTTTTT",
     (const GtUchar*) "TTTTTTTTTTTTTTTAAAGGGTTTCGCAAAGGGTTTCCCTTTTTTTTTTTTTTT",
     (const GtUchar*) "TTTTTTTTTTTTTTTAAAGGGTTTCCAAAGGGTTTCCCCTTTTTTTTTTTTTTT",
     (const GtUchar*) "TTTTTTTTTTTTTTTAAAGGGTTTCCTCAAAGGGTTTCCTTTTTTTTTTTTTTT",
     (const GtUchar*) "TTTTTTTTTTTTTTTAAACAGATCACCCGCTTTTTTTTTTTTTTTT",
     (const GtUchar*) "TTTTTTTTTTTTTTTAAACGGGTTTCTCAAAGGGTTCCCTTTTTTTTTTTTTTT"};
  GtUword lengths[GT_XDROP_NUM_OF_TESTS] =
  {54UL, 46UL, 50UL, 54UL, 54UL, 54UL, 46UL, 54UL},
    eval_scores[GT_XDROP_NUM_OF_TESTS *
      GT_XDROP_NUM_OF_TESTS *
      GT_XDROP_NUM_OF_TESTS] =
      {0, 13UL, 0, 1UL, 4UL, 1UL, 0, 7UL,
        13UL, 0, 0, 14UL, 15UL, 14UL, 0, 15UL,
        0, 0, 0, 0, 0, 0, 0, 0,
        1UL, 14UL, 0, 0, 1UL, 2UL, 0, 1UL,
        4UL, 15UL, 0, 1UL, 0, 8UL, 0, 1UL,
        1UL, 14UL, 0, 2UL, 8UL, 0, 0, 4UL,
        0, 0, 0, 0, 0, 0, 0, 0,
        7UL, 15UL, 0, 1UL, 1UL, 4UL, 0, 0,

        0, 13UL, 0, 1UL, 4UL, 5UL, 14UL, 7UL,
        13UL, 0, 0, 14UL, 15UL, 14UL, 12UL, 15UL,
        0, 0, 0, 20UL, 0, 19UL, 17UL, 0,
        1UL, 14UL, 20UL, 0, 5UL, 6UL, 15UL, 8UL,
        4UL, 15UL, 0, 5UL, 0, 8UL, 15UL, 10UL,
        5UL, 14UL, 19UL, 6UL, 8UL, 0, 14UL, 4UL,
        14UL, 12UL, 17UL, 15UL, 15UL, 14UL, 0, 14UL,
        7UL, 15UL, 0, 8UL, 10UL, 4UL, 14UL, 0,

        0, 13UL, 19UL, 1UL, 2UL, 2UL, 13UL, 3UL,
        13UL, 0, 9UL, 14UL, 14UL, 13UL, 12UL, 14UL,
        17UL, 4UL, 0, 18UL, 19UL, 16UL, 16UL, 18UL,
        1UL, 14UL, 18UL, 0, 2UL, 3UL, 13UL, 3UL,
        2UL, 14UL, 18UL, 2UL, 0, 4UL, 13UL, 4UL,
        2UL, 13UL, 19UL, 3UL, 4UL, 0, 13UL, 3UL,
        14UL, 12UL, 17UL, 13UL, 13UL, 14UL, 0, 14UL,
        3UL, 14UL, 18UL, 3UL, 4UL, 3UL, 13UL, 0,

        0, 13UL, 17UL, 1UL, 2UL, 2UL, 14UL, 3UL,
        13UL, 0, 4UL, 14UL, 15UL, 13UL, 12UL, 14UL,
        19UL, 9UL, 0, 18UL, 18UL, 19UL, 17UL, 18UL,
        1UL, 14UL, 18UL, 0, 2UL, 3UL, 13UL, 3UL,
        2UL, 14UL, 19UL, 2UL, 0, 4UL, 13UL, 4UL,
        2UL, 13UL, 16UL, 3UL, 4UL, 0, 14UL, 3UL,
        13UL, 12UL, 16UL, 13UL, 13UL, 13UL, 0, 13UL,
        3UL, 14UL, 18UL, 3UL, 4UL, 3UL, 14UL, 0,

        0, 0, 0, 1UL, 1UL, 1UL, 0, 1UL,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        1UL, 0, 0, 0, 1UL, 0, 0, 1UL,
        1UL, 0, 0, 1UL, 0, 0, 0, 1UL,
        1UL, 0, 0, 0, 0, 0, 0, 1UL,
        0, 0, 0, 0, 0, 0, 0, 0,
        1UL, 0, 0, 1UL, 1UL, 1UL, 0, 0,

        0, 0, 0, 1UL, 1UL, 1UL, 0, 1UL,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        1UL, 0, 0, 0, 1UL, 2UL, 0, 1UL,
        1UL, 0, 0, 1UL, 0, 0, 0, 1UL,
        1UL, 0, 0, 2UL, 0, 0, 0, 1UL,
        0, 0, 0, 0, 0, 0, 0, 0,
        1UL, 0, 0, 1UL, 1UL, 1UL, 0, 0,

        0, 13UL, 17UL, 1UL, 2UL, 2UL, 13UL, 3UL,
        13UL, 0, 4UL, 14UL, 14UL, 13UL, 12UL, 14UL,
        17UL, 4UL, 0, 18UL, 19UL, 16UL, 16UL, 19UL,
        1UL, 14UL, 18UL, 0, 2UL, 3UL, 13UL, 3UL,
        2UL, 14UL, 19UL, 2UL, 0, 4UL, 13UL, 4UL,
        2UL, 13UL, 16UL, 3UL, 4UL, 0, 13UL, 3UL,
        13UL, 12UL, 16UL, 13UL, 13UL, 13UL, 0, 13UL,
        3UL, 14UL, 19UL, 3UL, 4UL, 3UL, 13UL, 0,

        0, 13UL, 0, 1UL, 2UL, 2UL, 5UL, 3UL,
        13UL, 0, 0, 14UL, 15UL, 13UL, 0, 14UL,
        0, 0, 0, 0, 0, 0, 0, 0,
        1UL, 14UL, 0, 0, 2UL, 3UL, 5UL, 3UL,
        2UL, 15UL, 0, 2UL, 0, 4UL, 5UL, 4UL,
        2UL, 13UL, 0, 3UL, 4UL, 0, 6UL, 3UL,
        5UL, 0, 0, 5UL, 5UL, 6UL, 0, 5UL,
        3UL, 14UL, 0, 3UL, 4UL, 3UL, 5UL, 0};
  GtSeqabstract *useq, *vseq;

  GtXdropArbitraryscores score[GT_XDROP_NUM_OF_TESTS] = {{2, -2, -2, -2},
                                                         {2, -1, -1, -1},
                                                         {2, -1, -5, -2},
                                                         {2, -1, -2, -5},
                                                         {3, -2, -3, -3},
                                                         {3, -1, -1, -1},
                                                         {4, -1, -3, -3},
                                                         {10, -3, -8, -8}};
  GtXdropresources *resources;
  GtXdropbest best;
  GtXdropscore dropscore = (GtXdropscore) 12;
  GtMultieoplist *edit_ops = NULL;
  GtAlignment *alignment;

  gt_error_check(err);

  for (s = 0; s < GT_XDROP_NUM_OF_TESTS; ++s) {
    resources = gt_xdrop_resources_new(&score[s]);
    for (i = 0; i < GT_XDROP_NUM_OF_TESTS && !had_err; ++i) {
      for (j = 0; j < GT_XDROP_NUM_OF_TESTS; ++j) {
        useq = gt_seqabstract_new_gtuchar(strings[i], lengths[i], 0);
        vseq = gt_seqabstract_new_gtuchar(strings[j], lengths[j], 0);
        gt_evalxdroparbitscoresextend(true, &best, resources, useq, vseq,
                                      dropscore);

        edit_ops = gt_xdrop_backtrack(resources, &best);
        gt_ensure(edit_ops != NULL);
        alignment = gt_alignment_new_with_seqs(strings[i], best.ivalue,
                                               strings[j], best.jvalue);
        gt_alignment_set_multieop_list(alignment, edit_ops);
        gt_ensure(eval_scores[s*64+i*8+j] == gt_alignment_eval(alignment));

        gt_multieoplist_delete(edit_ops);
        gt_alignment_delete(alignment);
        if (i == j) {
          gt_evalxdroparbitscoresextend(false, &best, resources, useq, vseq,
                                        dropscore);

          edit_ops = gt_xdrop_backtrack(resources, &best);
          alignment = gt_alignment_new_with_seqs(strings[i], best.ivalue,
                                                 strings[j], best.jvalue);
          gt_alignment_set_multieop_list(alignment, edit_ops);
          gt_ensure(eval_scores[s*64+i*8+j] == gt_alignment_eval(alignment));
          gt_multieoplist_delete(edit_ops);
          gt_alignment_delete(alignment);
        }
        gt_seqabstract_delete(useq);
        gt_seqabstract_delete(vseq);
      }
    }
    gt_xdrop_resources_delete(resources);
  }

  return had_err;
}
