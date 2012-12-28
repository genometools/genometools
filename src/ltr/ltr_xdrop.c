/*
  Copyright (c) 2007 David Ellinghaus <d.ellinghaus@ikmb.uni-kiel.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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
#include "ltr_xdrop.h"

typedef struct
{
  int mis,
      ins,
      del;
} GtXdropArbitrarydistances;

#define GT_XDROP_MINUSINFINITYINT ((int)integermin)
#define ACCESSTOFRONT(D,K) ((unsigned long) (D) * (D) + (D) + (K))

/*
  For each entry in the DP-matrix we store a single byte, and
  use the three rightmost bits to mark which edge in the edit distance
  graph to trace back.
*/

#define MYREPLACEMENTBIT   ((unsigned char) 1)          /* replacement */
#define MYDELETIONBIT      (((unsigned char) 1) << 1)   /* deletion */
#define MYINSERTIONBIT     (((unsigned char) 1) << 2)   /* insertion */

/*
 The following function shows the matrix of the calculated fronts.
 */
/* CAUTION: fronts, that run over the matrix boundaries are not shown in
   the printed matrix.
 */
int gt_showfrontvalues(GtArrayGtXdropfrontvalue * fronts,
                       int distance,
                       unsigned char *useq,
                       unsigned char *vseq,
                       int ulen,
                       int vlen)
{
  unsigned long l;
  int i, j, k, d = distance + 1, filled = 0, integermax = MAX (ulen,vlen),
      integermin = -integermax;

  printf("frontvalues:\n");
  printf("        ");
  printf("%-3c ", vseq[0]);
  /* print vseq */
  for (i = 1; i < vlen; i++)
  {
    printf("%-3c ", vseq[i]);
  }
  for (i = 0; i <= ulen; i++)
  {
    printf("\n");
    /* print useq */
    if (i != 0)
    {
      printf("%-3c ", useq[i - 1]);
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
        if (fronts->spaceGtXdropfrontvalue[l].dptabrow ==
            GT_XDROP_MINUSINFINITYINT)
        {
          continue;
        }
        if (i == fronts->spaceGtXdropfrontvalue[l].dptabrow)
        {
          for (d = 0; d <= distance; d++)
          {
            if (k >= -d && k <= d && l == ACCESSTOFRONT(d, i - j))
            {
              printf("%-3d ", d);
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
           filled * 100.00 / ((ulen + 1) * (vlen + 1)));
  return 0;
}

/* If gt_showfrontvalues = "True", then
    (A)->nextfree##TYPE = POS+1; must be set:

 ** The other case, saves more space

   #define STOREINARRAYFRONTS(A,POS,TYPE,VAL)\
           GT_CHECKARRAYSPACEMULTI(A,TYPE,POS+1);\
           (A)->space##TYPE[POS] = VAL;\
           (A)->nextfree##TYPE = POS+1
*/

/*
 POS+1 necessary, since at POS = 0 and nextfree##TYPE = 0
 nothing will be allocated in addition!
 */
#define STOREINARRAYFRONTS(A,POS,TYPE,VAL)\
        GT_CHECKARRAYSPACEMULTI(A,TYPE,POS+1);\
        (A)->space##TYPE[POS] = VAL

/*
 The following macro swaps two values.
 */
#define SWAPMAYBE\
        if (m < n)\
        {\
          r = m;\
          m = n;\
          n = r;\
        }
/*
 The following macro determines the greatest common divisor.
 */
#define GGT\
        do\
        {\
          r = m % n;\
          m = n;\
          n = r;\
        } while (r != 0)

/*
 The following function calculates the distance from the given scores.
 */
static void gt_calculatedistancesfromscores(
                                    GtXdropArbitraryscores *arbitscores,
                                    GtXdropArbitrarydistances *arbitdistances)
{
  unsigned int m, n, r;
  int mat, mis, ins, del;

  /* if mat is odd double all scores */
  if (GT_MOD2((unsigned int)arbitscores->mat))
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

  m = (unsigned int)(mat - mis);
  n = (unsigned int)(mat/2 - ins);
  SWAPMAYBE
  GGT;

  n = (unsigned int)(mat/2 - del);
  SWAPMAYBE
  GGT;

  arbitscores->gcd = (int) m;

  arbitdistances->mis = (mat - mis) / arbitscores->gcd;
  arbitdistances->ins = (mat/2 - ins) / arbitscores->gcd;
  arbitdistances->del = (mat/2 - del) / arbitscores->gcd;

}

/*
 The following function calculates the maximal allowed number of
 generations with all front values equal minus infinity.
 */
static void gt_calculateallowedMININFINITYINTgenerations(
                   int *allowedMININFINITYINTgenerations,
                   GtXdropArbitrarydistances *arbitdistances)
{
  *allowedMININFINITYINTgenerations = MAX(arbitdistances->mis,
                                          arbitdistances->ins);
  *allowedMININFINITYINTgenerations = MAX(*allowedMININFINITYINTgenerations,
                                          arbitdistances->del);
  (*allowedMININFINITYINTgenerations)--;
}

#define GT_XDROP_EVAL(K,D)\
        ((K) * arbitscores->mat/2 - (D) * arbitscores->gcd)
#define GT_XDROP_SETDBACK(XDROPBELOWSCORE)\
        (XDROPBELOWSCORE + arbitscores->mat/2) / arbitscores->gcd + 1

/*
 The following macro checks for wildcard and separator symbols.
 */

#define GT_XDROP_COMPARESYMBOLSSEP(VARA,VARB,I,J)\
        GT_XDROP_USEQ(VARA,I);\
        if (VARA == (GtUchar) SEPARATOR)\
        {\
          ulen = I;\
          break;\
        }\
        GT_XDROP_VSEQ(VARB,J);\
        if (VARB == (GtUchar) SEPARATOR)\
        {\
          vlen = J;\
          break;\
        }\
        if (VARA != VARB || VARA == (GtUchar) WILDCARD)\
        {\
          break;\
        }

#define GT_XDROP_EVALXDROPARBITSCORES GT_XDROP_EVALXDROPARBITSCORESRIGHT
#define GT_XDROP_USEQ(VAR,I) VAR = gt_encseq_get_encoded_char(str_useq,\
                                                              useq + (I),\
                                                           GT_READMODE_FORWARD)
#define GT_XDROP_VSEQ(VAR,J) VAR = gt_encseq_get_encoded_char(str_vseq,\
                                                              vseq + (J),\
                                                           GT_READMODE_FORWARD)

#include "myxdrop.gen"

#undef GT_XDROP_EVALXDROPARBITSCORES
#undef GT_XDROP_USEQ
#undef GT_XDROP_VSEQ

/*
  Now we redefine the macros to compute the left to right
  we use macros to abstract from the differences. The
  following 4 macros are defined for the right to left extension beginning
  with the last character of the strings under consideration.
*/

#define GT_XDROP_EVALXDROPARBITSCORES GT_XDROP_EVALXDROPARBITSCORESLEFT
#define GT_XDROP_USEQ(VAR,I) VAR = gt_encseq_get_encoded_char(str_useq,\
                                                     useq - 1UL - (I),\
                                                     GT_READMODE_FORWARD)
#define GT_XDROP_VSEQ(VAR,J) VAR = gt_encseq_get_encoded_char(str_vseq,\
                                                     vseq - 1UL - (J),\
                                                     GT_READMODE_FORWARD)

#include "myxdrop.gen"
