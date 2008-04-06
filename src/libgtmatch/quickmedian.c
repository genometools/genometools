/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#include <assert.h>
#include "divmodmul.h"

/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */

typedef int Elem;

#define ElemDEREF(V)     *(V)
#define ElemGREATER(A,B) (ElemDEREF(A) > ElemDEREF(B))
#define ElemSWAP(A,B)    {\
                           register Elem tmp = ElemDEREF(A);\
                                         ElemDEREF(A) = ElemDEREF(B);\
                                         ElemDEREF(B) = tmp;\
                         }

Elem *quickmedian (Elem *arr,unsigned long n)
{
  Elem *low, *high, *median, *middle, *ll, *hh;

  assert(n > 0);
  low = arr;
  high = arr + n - 1;
  median = low + DIV2(high - low + 1);
  for (;;)
  {
   if (high <= low)                   /* One element only */
    {
      return median;
    }
    if (high == low + 1)
    {                                  /* Two elements only */
      if (ElemGREATER(low,high))
      {
        ElemSWAP (low, high);
      }
      return median;
    }

    /* Find median of low, middle and high items; swap into position low */
    middle = low + DIV2(high - low + 1);
    if (ElemGREATER(middle,high))
    {
      ElemSWAP (middle, high);
    }
    if (ElemGREATER(low,high))
    {
      ElemSWAP (low, high);
    }
    if (ElemGREATER(middle,low))
    {
      ElemSWAP (middle, low);
    }
    /* Swap low item (now in position middle) into position (low+1) */
    ElemSWAP (middle, low + 1);

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;)
    {
      do
      {
        ll++;
      } while (ElemGREATER(low,ll));
      do
      {
        hh--;
      } while  (ElemGREATER(hh,low));
      if (hh < ll)
      {
        break;
      }
      ElemSWAP (ll, hh);
    }

    /* Swap middle item (in position low) back into correct position */
    ElemSWAP (low, hh);

    /* Re-set active partition */
    if (hh <= median)
    {
      low = ll;
    }
    if (hh >= median)
    {
      high = hh - 1;
    }
  }
}
