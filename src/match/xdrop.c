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
  long dptabrow;
  unsigned char dptabdirection; /* one of the bits REPLACEMENTBIT,
                                                   DELETIONBIT,
                                                   INSERTIONBIT */
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
        if (fronts->spaceGtXdropfrontvalue[l].dptabrow == integermin)
        {
          continue;
        }
        if (i == fronts->spaceGtXdropfrontvalue[l].dptabrow)
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

/*
 POS+1 necessary, since at POS = 0 and nextfree##TYPE = 0
 nothing will be allocated in addition!
 */

#define STOREINARRAYFRONTS(POS,VAL)\
        GT_CHECKARRAYSPACEMULTI(&res->fronts,GtXdropfrontvalue,POS+1);\
        res->fronts.spaceGtXdropfrontvalue[POS] = VAL

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

/*
 The following function calculates the maximal allowed number of
 generations with all front values equal minus infinity.
 */

static int gt_calculateallowedMININFINITYINTgenerations(
                   const GtXdropArbitrarydistances *arbitdistances)
{
  return MAX(MAX(arbitdistances->mis,arbitdistances->ins),
             arbitdistances->del) - 1;
}

#define LOOKAHEAD 10

#define GT_XDROP_EVAL(K,D)\
        ((K) * res->arbitscores->mat/2 - (D) * res->arbitdistances.gcd)
#define GT_XDROP_SETDBACK(XDROPBELOWSCORE)\
        (XDROPBELOWSCORE + res->arbitscores->mat/2)/res->arbitdistances.gcd + 1

/*
 The following macro checks for wildcard and separator symbols.
 */

#undef GT_XDROP_CHECKINCREMENT
#ifdef GT_XDROP_CHECKINCREMENT
#define GT_XDROP_UINCREMENT(I)\
          characcess++;\
          if (forward)\
          {\
            if (lastaccess_useq != ULONG_MAX)\
            {\
              if (lastaccess_useq > uoffset + (I))\
              {\
                unsigned long goback = lastaccess_useq - (uoffset + (I));\
                if (maxgoback < goback)\
                {\
                  maxgoback = goback;\
                }\
              } else\
              {\
                if (lastaccess_useq + LOOKAHEAD >= uoffset + (I))\
                {\
                  consecutive++;\
                }\
              }\
            }\
            lastaccess_useq = uoffset + (I);\
          }

#define GT_XDROP_VINCREMENT(J)\
          characcess++;\
          if (forward)\
          {\
            if (lastaccess_vseq != ULONG_MAX)\
            {\
              if (lastaccess_vseq > voffset + (J))\
              {\
                unsigned long goback = lastaccess_vseq - (voffset + (J));\
                if (maxgoback < goback)\
                {\
                  maxgoback = goback;\
                }\
              } else\
              {\
                if (lastaccess_vseq + LOOKAHEAD >= voffset + (J))\
                {\
                  consecutive++;\
                }\
              }\
            }\
            lastaccess_vseq = voffset + (J);\
          }
#else
#define GT_XDROP_UINCREMENT(I) /* Nothing */
#define GT_XDROP_VINCREMENT(J) /* Nothing */
#endif

#define GT_XDROP_SEQACC(SEQ,OFF,IDX)\
        gt_seqabstract_encoded_char(SEQ,forward ? (OFF) + (IDX)\
                                                : (OFF) - 1UL - (IDX))

#define GT_XDROP_COMPARESYMBOLSSEP(I,J)\
        {\
          GtUchar a, b;\
          a = GT_XDROP_SEQACC(useq,uoffset,I);\
          GT_XDROP_UINCREMENT(I)\
          if (a == (GtUchar) SEPARATOR)\
          {\
            ulen = I;\
            break;\
          }\
          b = GT_XDROP_SEQACC(vseq,voffset,J);\
          GT_XDROP_VINCREMENT(J)\
          if (b == (GtUchar) SEPARATOR)\
          {\
            vlen = J;\
            break;\
          }\
          if (a != b || a == (GtUchar) WILDCARD)\
          {\
            break;\
          }\
        }

static long gt_xdrop_frontvalue(const GtXdropresources *res,long d,long k)
{
  return res->fronts.spaceGtXdropfrontvalue[GT_XDROP_FRONTIDX(d,k)].dptabrow;
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
  long integermax,
      integermin,
      ulen,
      vlen,
      i,
      j = 0,
      k = 0,         /* diagonal */
      end_k,         /* diagonal of endpoint (ulen, vlen) */
      lbound = 0,    /* diagonal lower bound */
      ubound = 0,    /* diagonal upper bound */
      lboundtmp = 0,
      uboundtmp = 0,
      tmprow,
      d = 0,         /* distance */
      d_pre = 0,     /* previous distance */
      dback;
  int allowedMININFINITYINTgenerations,
      currentMININFINITYINTgeneration = 0;
  GtXdropfrontvalue tmpfront;
  GtXdropscore bigt_tmp;        /* best score T' seen already */
  bool alwaysMININFINITYINT = true;
#ifdef GT_XDROP_CHECKINCREMENT
  unsigned long lastaccess_useq = ULONG_MAX,
                lastaccess_vseq = ULONG_MAX,
                maxgoback = 0;
  unsigned long characcess = 0, consecutive = 0;
#endif

#undef DEBUG_POSTEL
#ifdef DEBUG_POSTEL
  printf("%s(forward=%s,useq.len=%lu,uoffset=%lu,vseq.len=%lu,voffset=%lu)\n",
          __func__,forward ? "true" : "false",
          gt_seqabstract_length_get(useq),
          uoffset,
          gt_seqabstract_length_get(vseq),
          voffset);
#endif
  ulen = (long) gt_seqabstract_length_get(useq);
  vlen = (long) gt_seqabstract_length_get(vseq);
  end_k = ulen - vlen;              /* diagonal of endpoint (ulen, vlen) */
  allowedMININFINITYINTgenerations
    = gt_calculateallowedMININFINITYINTgenerations(&res->arbitdistances);
  dback = GT_XDROP_SETDBACK(xdropbelowscore);
  res->big_t.nextfreeGtXdropscore = 0;
  res->fronts.nextfreeGtXdropfrontvalue = 0;
  integermax = MAX(ulen, vlen);
  integermin = -integermax;
  /* phase 0 */
  for (i = 0; i < MIN(ulen,vlen); i++)
  {
    GT_XDROP_COMPARESYMBOLSSEP(i,i);
  }
  /* alignment already finished */
  if (i >= ulen || i >= vlen)
  {
    lbound =  1L;
    ubound = -1L;
  }
  tmpfront.dptabrow = i;
  tmpfront.dptabdirection = (GtUchar) 0;   /* no predecessor */
  STOREINARRAYFRONTS(0, tmpfront);
  xdropbest->score = bigt_tmp = GT_XDROP_EVAL(i + i, 0);
  gt_assert(i >= 0);
  xdropbest->ivalue = xdropbest->jvalue = (unsigned long) i;
  GT_STOREINARRAY (&res->big_t, GtXdropscore, 10, bigt_tmp);
  /* phase d > 0 */
  while (lbound <= ubound)
  {
    d++;
    d_pre = d - dback;
    /* calculate fronts */
    for (k = lbound - 1; k <= ubound + 1; k++)
    {
      i = integermin;
      /* case 1 : DELETION-EDGE  */
      if (lbound < k &&
          d - res->arbitdistances.del >= 0 &&
          -(d - res->arbitdistances.del) <= k - 1 &&
          k - 1 <= d - res->arbitdistances.del)
      {
        i = gt_xdrop_frontvalue(res,d - res->arbitdistances.del,k-1) + 1;
        tmpfront.dptabdirection = GT_XDROP_DELETIONBIT;
      }
      /* case 2: REPLACEMENT-EDGE */
      if (lbound <= k &&
          k <= ubound &&
          d - res->arbitdistances.mis >= 0 &&
          -(d - res->arbitdistances.mis) <= k &&
          k <= d - res->arbitdistances.mis)
      {
        /* test, if case 1 has happened. */
        tmprow = gt_xdrop_frontvalue(res,d - res->arbitdistances.mis,k) + 1;
        if (!(tmpfront.dptabdirection & GT_XDROP_DELETIONBIT) ||
            tmprow > i)
        {
          i = tmprow;
          tmpfront.dptabdirection = GT_XDROP_REPLACEMENTBIT;
        }
      }
      /* case 3: INSERTION-EDGE */
      if (k < ubound &&
          d - res->arbitdistances.ins >= 0 &&
          -(d - res->arbitdistances.ins) <= k + 1 &&
           k + 1 <= d - res->arbitdistances.ins)
      {
        tmprow = gt_xdrop_frontvalue(res,d - res->arbitdistances.ins,k+1);
        if (!(tmpfront.dptabdirection & (GT_XDROP_DELETIONBIT |
                                         GT_XDROP_REPLACEMENTBIT)) ||
            tmprow > i)
        {
          i = tmprow;
          tmpfront.dptabdirection = GT_XDROP_INSERTIONBIT;
        }
      }
      j = i - k;
      /* if i = MINUSINFINITYINY or MINUSINFINITYINY + 1 */
      if (i < 0)
      {
        if (tmpfront.dptabdirection == (GtUchar) 0)
        {
          alwaysMININFINITYINT = false;
        }
        tmpfront.dptabrow = integermin;
      } else
      {
        /* alignment score smaller than T - X */
        if (d_pre > 0 &&
            res->big_t.spaceGtXdropscore != NULL &&
            GT_XDROP_EVAL (i + j, d) <
            res->big_t.spaceGtXdropscore[d_pre] - xdropbelowscore)
        {
          tmpfront.dptabrow = integermin;
        } else
        {
          if (k <= -d || k >= d || /* not correct boundaries for
                                      ACCESTOFRONT(d-1,k) */
              (gt_xdrop_frontvalue(res,d-1,k) < i && i <= MIN(ulen,vlen + k)))
          {
            while (i < ulen && j < vlen)
            {
              GT_XDROP_COMPARESYMBOLSSEP(i,j);
              i++;
              j++;
            }
            alwaysMININFINITYINT = false;
            tmpfront.dptabrow = i;
            if (GT_XDROP_EVAL(i + j, d) > bigt_tmp)
            {
              xdropbest->score = bigt_tmp = GT_XDROP_EVAL(i + j, d);
              gt_assert(i >= 0 && j >= 0);
              xdropbest->ivalue = (unsigned long) i;
              xdropbest->jvalue = (unsigned long) j;
            }
          } else
          {
            alwaysMININFINITYINT = false;
            tmpfront.dptabrow = gt_xdrop_frontvalue(res,d-1,k);
          }
        }
      }
      STOREINARRAYFRONTS(GT_XDROP_FRONTIDX(d,k),tmpfront);
      /* delete value for test for */
      /* INSERTIONBIT/REPLACEMENTBIT/DELETIONBIT above
      tmpfront.dptabdirection = (GtUchar) 0; */
    }
    /* if all front values are integermin,
       aligment prematurely finished if allowedMININFINITYINTgenerations
       exceeded
       (full front has already ended at d - currentMININFINITYINTgeneration). */
    if (alwaysMININFINITYINT)
    {
      currentMININFINITYINTgeneration++;
      if (currentMININFINITYINTgeneration > allowedMININFINITYINTgenerations)
      {
        d -= currentMININFINITYINTgeneration;
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
    for (k = -d; k < lbound - 1; k++)
    {
      tmpfront.dptabrow = integermin;
      STOREINARRAYFRONTS(GT_XDROP_FRONTIDX(d,k), tmpfront);
    }
    for (k = ubound + 2; k <= d; k++)
    {
      tmpfront.dptabrow = integermin;
      STOREINARRAYFRONTS(GT_XDROP_FRONTIDX(d,k), tmpfront);
    }
    /* alignment finished */
    if (-d <= end_k && end_k <= d && gt_xdrop_frontvalue(res,d,end_k) == ulen)
    {
      break;
    }
    /* pruning lower bound
       lbound may decrease by one or increase/stays the same */
    for (k = lbound - 1; k <= ubound + 1; k++)
    {
      if (gt_xdrop_frontvalue(res,d,k) > integermin)
      {
        lbound = k;
        break;
      }
    }
    /* pruning upper bound
       ubound may increase by one or decrease/stays the same */
    for (k = ubound + 1; k >= lbound - 1; k--)
    {
      if (gt_xdrop_frontvalue(res,d,k) > integermin)
      {
        ubound = k;
        break;
      }
    }
    /* handling boundaries lower bound */
    lboundtmp = lbound;
    for (k = 0; k >= lbound; k--)
    {
      if (gt_xdrop_frontvalue(res,d,k) == vlen + k)
      {
        lboundtmp = k;
        break;
      }
    }
    /* handling boundaries upper bound */
    uboundtmp = ubound;
    for (k = 0; k <= ubound; k++)
    {
      if (gt_xdrop_frontvalue(res,d,k) == ulen)
      {
        uboundtmp = k;
        break;
      }
    }
    lbound = MAX(lbound, lboundtmp);
    ubound = MIN(ubound, uboundtmp);
  }
#ifdef GT_XDROP_CHECKINCREMENT
  if (maxgoback > 0)
  {
    printf("maxgoback = %lu\n",maxgoback);
  }
  printf("characcess=%lu,consecutive=%lu(%.2f)\n",characcess,consecutive,
                                       100.0 * (double) consecutive/characcess);
#endif
}
