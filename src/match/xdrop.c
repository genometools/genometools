/*
  Copyright (c) 2007 David Ellinghaus <d.ellinghaus@ikmb.uni-kiel.de>
  Copyright (c) 2012-2013 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2013 Center for Bioinformatics, University of Hamburg

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
#include "core/minmax.h"
#include "xdrop.h"

typedef struct
{
  int mis,
      ins,
      del,
      gcd;  /* greatest common divisor */
} GtXdropArbitrarydistances;

typedef struct
{
  long row;
  unsigned char direction; /* one of the bits GT_XDROP_REPLACEMENTBIT,
                                              GT_XDROP_DELETIONBIT,
                                              GT_XDROP_INSERTIONBIT */
} GtXdropfrontvalue;

GT_DECLAREARRAYSTRUCT(GtXdropfrontvalue);

#define GT_XDROP_FRONTIDX(D,K)    ((unsigned long) (D) * (D) + (D) + (K))

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
                       long distance,
                       const unsigned char *useqptr,
                       const unsigned char *vseqptr,
                       long ulen,
                       long vlen)
{
  unsigned long l;
  long i, j, k, d = distance + 1, filled = 0, integermax = MAX(ulen,vlen),
       integermin = -integermax;

  printf("frontvalues:\n");
  printf("        ");
  printf("%-3c ", vseqptr[0]);
  /* print vseq */
  for (i = 1L; i < vlen; i++)
  {
    printf("%-3c ", vseqptr[i]);
  }
  for (i = 0; i <= ulen; i++)
  {
    printf("\n");
    /* print useq */
    if (i != 0)
    {
      printf("%-3c ", useqptr[i - 1]);
    } else
    {
      printf("    ");
    }
    for (j = 0; j <= vlen; j++)
    {
      d = distance + 1;
      k = i - j;
      for (l = 0; l < fronts->nextfreeGtXdropfrontvalue; l++)
      {
        if (fronts->spaceGtXdropfrontvalue[l].row == integermin)
        {
          continue;
        }
        if (i == fronts->spaceGtXdropfrontvalue[l].row)
        {
          for (d = 0; d <= distance; d++)
          {
            if (k >= -d && k <= d && l == GT_XDROP_FRONTIDX(d, i - j))
            {
              printf("%-3ld ", d);
              l = fronts->nextfreeGtXdropfrontvalue;
              filled++;
              break;
            }
          }
        }
      }
      if (d > distance)
      {
        printf("U   ");
      }
    }
  }
  printf("\n%.2f percent of matrix filled\n",
           (double) filled * 100.00 / ((ulen + 1) * (vlen + 1)));
}

static unsigned int gt_xdrop_gcd(unsigned int m, unsigned int n)
{
  unsigned int r;

  if (m < n)
  {
    r = m;
    m = n;
    n = r;
  }
  do
  {
    gt_assert(m >= n);
    r = m % n;
    m = n;
    n = r;
  } while (r != 0);
  return m;
}

/*
 The following function calculates the distance from the given scores.
 */
static void gt_calculatedistancesfromscores(
                                    const GtXdropArbitraryscores *arbitscores,
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
  dist->gcd
    = (int) gt_xdrop_gcd(gt_xdrop_gcd((unsigned int) (mat-mis),
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

static long gt_xdrop_frontvalue_get(const GtXdropresources *res,long d,long k)
{
  const unsigned long frontidx = GT_XDROP_FRONTIDX(d,k);
  return res->fronts.spaceGtXdropfrontvalue[frontidx].row;
}

static void gt_xdrop_frontvalue_set(GtXdropresources *res,long d,long k,
                                    GtXdropfrontvalue value)
{
  const unsigned long frontidx = GT_XDROP_FRONTIDX(d,k);

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
                                   unsigned long uoffset,
                                   unsigned long voffset,
                                   GtXdropscore xdropbelowscore)
{
  const long ulen = (long) gt_seqabstract_length_get(useq),
             vlen = (long) gt_seqabstract_length_get(vseq),
             end_k = ulen - vlen, /* diagonal of endpoint (ulen, vlen) */
             integermax = MAX(ulen, vlen),
             integermin = -integermax,
             dback = GT_XDROP_SETDBACK(xdropbelowscore);
  long idx,
       lbound,    /* diagonal lower bound */
       ubound,    /* diagonal upper bound */
       currd = 0; /* distance */
  /*The following function calculates the maximal allowed number of
    generations with all front values equal minus infinity.*/
  const int allowedMININFINITYINTgenerations
              = MAX(MAX(res->arbitdistances.mis,res->arbitdistances.ins),
                    res->arbitdistances.del) - 1;
  int currentMININFINITYINTgeneration = 0;
  GtXdropfrontvalue tmpfront;
  GtXdropscore bigt_tmp;        /* best score T' seen already */
  bool alwaysMININFINITYINT = true;

  res->big_t.nextfreeGtXdropscore = 0;
  res->fronts.nextfreeGtXdropfrontvalue = 0;
  /* phase 0 */
  idx = (long) gt_seqabstract_lcp(forward,
                                  useq,
                                  vseq,
                                  forward ? uoffset : uoffset - 1,
                                  forward ? voffset : voffset - 1,
                                  (unsigned long) MIN(ulen,vlen));
  /* alignment already finished */
  if (idx >= ulen || idx >= vlen)
  {
    lbound =  1L;
    ubound = -1L;
  } else
  {
    lbound = 0;
    ubound = 0;
  }
  tmpfront.row = idx;
  tmpfront.direction = (GtUchar) 0;   /* no predecessor */
  gt_xdrop_frontvalue_set(res,0,0,tmpfront);
  xdropbest->score = bigt_tmp = GT_XDROP_EVAL(idx + idx, 0);
  gt_assert(idx >= 0);
  xdropbest->ivalue = xdropbest->jvalue = (unsigned long) idx;
  GT_STOREINARRAY (&res->big_t, GtXdropscore, 10, bigt_tmp);
  /* phase d > 0 */
  while (lbound <= ubound)
  {
    long k;

    currd++;
    /* calculate fronts */
    for (k = lbound - 1; k <= ubound + 1; k++)
    {
      long i = integermin;
      /* case 1 : DELETION-EDGE  */
      if (lbound < k &&
          currd - res->arbitdistances.del >= 0 &&
          -(currd - res->arbitdistances.del) <= k - 1 &&
          k - 1 <= currd - res->arbitdistances.del)
      {
        i = gt_xdrop_frontvalue_get(res,currd - res->arbitdistances.del,k-1)
            + 1;
        tmpfront.direction = GT_XDROP_DELETIONBIT;
      }
      /* case 2: REPLACEMENT-EDGE */
      if (lbound <= k &&
          k <= ubound &&
          currd - res->arbitdistances.mis >= 0 &&
          -(currd - res->arbitdistances.mis) <= k &&
          k <= currd - res->arbitdistances.mis)
      {
        /* test, if case 1 has happened. */
        long row = gt_xdrop_frontvalue_get(res,
                                           currd - res->arbitdistances.mis,k)
                   + 1;
        if (!(tmpfront.direction & GT_XDROP_DELETIONBIT) || row > i)
        {
          i = row;
          tmpfront.direction = GT_XDROP_REPLACEMENTBIT;
        }
      }
      /* case 3: INSERTION-EDGE */
      if (k < ubound &&
          currd - res->arbitdistances.ins >= 0 &&
          -(currd - res->arbitdistances.ins) <= k + 1 &&
           k + 1 <= currd - res->arbitdistances.ins)
      {
        long row = gt_xdrop_frontvalue_get(res,
                                           currd - res->arbitdistances.ins,k+1);
        if (!(tmpfront.direction & (GT_XDROP_DELETIONBIT |
                                    GT_XDROP_REPLACEMENTBIT)) || row > i)
        {
          i = row;
          tmpfront.direction = GT_XDROP_INSERTIONBIT;
        }
      }
      /* if i = MINUSINFINITYINY or MINUSINFINITYINY + 1 */
      if (i < 0)
      {
        if (tmpfront.direction == (GtUchar) 0)
        {
          alwaysMININFINITYINT = false;
        }
        tmpfront.row = integermin;
      } else
      {
        long j = i - k;
        const long previousd = currd - dback;

        /* alignment score smaller than T - X */
        if (previousd > 0 &&
            res->big_t.spaceGtXdropscore != NULL &&
            GT_XDROP_EVAL (i + j, currd) <
            res->big_t.spaceGtXdropscore[previousd] - xdropbelowscore)
        {
          tmpfront.row = integermin;
        } else
        {
          if (k <= -currd || k >= currd ||
              (gt_xdrop_frontvalue_get(res,currd-1,k) < i &&
               i <= MIN(ulen,vlen + k)))
          {
            if (ulen > i && vlen > j)
            {
              unsigned long lcp;
              gt_assert(forward || (uoffset > (unsigned long) i &&
                                    voffset > (unsigned long) j));
              lcp = gt_seqabstract_lcp(forward,
                                       useq,
                                       vseq,
                                       forward ? uoffset + i : uoffset - i - 1,
                                       forward ? voffset + j : voffset - j - 1,
                                       (unsigned long) MIN(ulen - i,vlen - j));
              i += lcp;
              j += lcp;
            }
            alwaysMININFINITYINT = false;
            tmpfront.row = i;
            if (GT_XDROP_EVAL(i + j, currd) > bigt_tmp)
            {
              xdropbest->score = bigt_tmp = GT_XDROP_EVAL(i + j,currd);
              gt_assert(i >= 0 && j >= 0);
              xdropbest->ivalue = (unsigned long) i;
              xdropbest->jvalue = (unsigned long) j;
            }
          } else
          {
            alwaysMININFINITYINT = false;
            tmpfront.row = gt_xdrop_frontvalue_get(res,currd-1,k);
          }
        }
      }
      gt_xdrop_frontvalue_set(res,currd,k,tmpfront);
    }
    /* if all front values are integermin, aligment prematurely finished if
       allowedMININFINITYINTgenerations exceeded (full front has already ended
       at currd - currentMININFINITYINTgeneration). */
    if (alwaysMININFINITYINT)
    {
      currentMININFINITYINTgeneration++;
      if (currentMININFINITYINTgeneration > allowedMININFINITYINTgenerations)
      {
        currd -= currentMININFINITYINTgeneration;
        break;
      }
    } else
    {
      currentMININFINITYINTgeneration = 0;
      alwaysMININFINITYINT = true;
    }
    GT_STOREINARRAY (&res->big_t, GtXdropscore, 10, bigt_tmp);
    /* fill out of bounds values of integermin
       needed for gt_showfrontvalues function */
    for (k = -currd; k < lbound - 1; k++)
    {
      tmpfront.row = integermin;
      gt_xdrop_frontvalue_set(res,currd,k,tmpfront);
    }
    for (k = ubound + 2; k <= currd; k++)
    {
      tmpfront.row = integermin;
      gt_xdrop_frontvalue_set(res,currd,k,tmpfront);
    }
    /* alignment finished */
    if (-currd <= end_k && end_k <= currd &&
        gt_xdrop_frontvalue_get(res,currd,end_k) == ulen)
    {
      break;
    }
    /* pruning lower bound
       lbound may decrease by one or increase/stays the same */
    for (k = lbound - 1; k <= ubound + 1; k++)
    {
      if (gt_xdrop_frontvalue_get(res,currd,k) > integermin)
      {
        lbound = k;
        break;
      }
    }
    /* pruning upper bound
       ubound may increase by one or decrease/stays the same */
    for (k = ubound + 1; k >= lbound - 1; k--)
    {
      if (gt_xdrop_frontvalue_get(res,currd,k) > integermin)
      {
        ubound = k;
        break;
      }
    }
    /* handling boundaries lower bound */
    for (k = 0; k >= lbound; k--)
    {
      if (gt_xdrop_frontvalue_get(res,currd,k) == vlen + k)
      {
        lbound = k;
        break;
      }
    }
    /* handling boundaries upper bound */
    for (k = 0; k <= ubound; k++)
    {
      if (gt_xdrop_frontvalue_get(res,currd,k) == ulen)
      {
        ubound = k;
        break;
      }
    }
  }
}
