/*
  Modifications for integration with genometools
  2008 Thomas Jahns <Thomas.Jahns@gmx.net>

  The advertising clause 3. was removed due to the corresponding
  revoke by William Hoskins on July 22, 1999.
  <ftp://ftp.cs.berkeley.edu/pub/4bsd/README.Impt.License.Change>
*/
/*-
 * Copyright (c) 1992, 1993
 *        The Regents of the University of California.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 4. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#include <stdlib.h>

#include "core/fptr_api.h"
#include "core/minmax.h"
#include "core/types_api.h"

static inline void *
med3(void *a, void *b, void *c, GtCompareWithData cmp, void *data);
static inline void       swapfunc(char *, char *, int, int);

/*
 * Qsort routine from Bentley & McIlroy's "Engineering a Sort Function".
 */
#define swapcode(TYPE, parmi, parmj, n) do {            \
    GtWord i = (n) / sizeof (TYPE);                     \
    TYPE *pi = (TYPE *) (parmi);                        \
    TYPE *pj = (TYPE *) (parmj);                        \
    do {                                                \
      TYPE t = *pi;                                     \
      *pi++ = *pj;                                      \
      *pj++ = t;                                        \
    } while (--i > 0);                                  \
  } while (0)

#define SWAPINIT(a, es) \
  swaptype = ((char *)a - (char *)0) % sizeof (GtWord) || \
  es % sizeof (GtWord) ? 2 : es == sizeof (GtWord)? 0 : 1;

static inline void
swapfunc(char *a, char *b, int n, int swaptype)
{
  if (swaptype <= 1)
    swapcode(GtWord, a, b, n);
  else
    swapcode(char, a, b, n);
}

#define swap(a, b)                                  \
  do {                                              \
    if (swaptype == 0) {                            \
      GtWord t = *(GtWord *)(a);                    \
      *(GtWord *)(a) = *(GtWord *)(b);              \
      *(GtWord *)(b) = t;                           \
    } else                                          \
      swapfunc(a, b, es, swaptype);                 \
  } while (0)

#define vecswap(a, b, n)       if ((n) > 0) swapfunc(a, b, n, swaptype)

static inline void *
med3(void *a, void *b, void *c, GtCompareWithData cmp, void *data)
{
  return cmp(a, b, data) < 0 ?
    (cmp(b, c, data) < 0 ? b : (cmp(a, c, data) < 0 ? c : a ))
    :(cmp(b, c, data) > 0 ? b : (cmp(a, c, data) < 0 ? a : c ));
}

void
gt_qsort_r(void *a, size_t n, size_t es, void *data, GtCompareWithData cmp)
{
  char *pa, *pb, *pc, *pd, *pl, *pm, *pn;
  int d, r, swaptype, swap_cnt;

  while (1)
  {
    SWAPINIT(a, es);
    swap_cnt = 0;
    if (n < 7) {
      for (pm = (char *)a + es; pm < (char *)a + n * es; pm += es)
        for (pl = pm;
             pl > (char *)a && cmp(pl - es, pl, data) > 0;
             pl -= es)
          swap(pl, pl - es);
      return;
    }
    pm = (char *)a + (n / 2) * es;
    if (n > 7) {
      pl = a;
      pn = (char *)a + (n - 1) * es;
      if (n > 40) {
        d = (n / 8) * es;
        pl = med3(pl, pl + d, pl + 2 * d, cmp, data);
        pm = med3(pm - d, pm, pm + d, cmp, data);
        pn = med3(pn - 2 * d, pn - d, pn, cmp, data);
      }
      pm = med3(pl, pm, pn, cmp, data);
    }
    swap(a, pm);
    pa = pb = (char *)a + es;

    pc = pd = (char *)a + (n - 1) * es;
    for (;;) {
      while (pb <= pc && (r = cmp(pb, a, data)) <= 0) {
        if (r == 0) {
          swap_cnt = 1;
          swap(pa, pb);
          pa += es;
        }
        pb += es;
      }
      while (pb <= pc && (r = cmp(pc, a, data)) >= 0) {
        if (r == 0) {
          swap_cnt = 1;
          swap(pc, pd);
          pd -= es;
        }
        pc -= es;
      }
      if (pb > pc)
        break;
      swap(pb, pc);
      swap_cnt = 1;
      pb += es;
      pc -= es;
    }
    if (swap_cnt == 0) {  /* Switch to insertion sort */
      for (pm = (char *)a + es; pm < (char *)a + n * es; pm += es)
        for (pl = pm;
             pl > (char *)a && cmp(pl - es, pl, data) > 0;
             pl -= es)
          swap(pl, pl - es);
      return;
    }

    pn = (char *)a + n * es;
    r = MIN(pa - (char *)a, pb - pa);
    vecswap(a, pb - r, r);
    r = MIN(pd - pc, pn - pd - es);
    vecswap(pb, pn - r, r);
    if ((r = pb - pa) > es)
      gt_qsort_r(a, r / es, es, data, cmp);
    if ((r = pd - pc) > es) {
      /* Iterate rather than recurse to save stack space */
      a = pn - r;
      n = r / es;
    }
    else
      break;
  }
/*            qsort(pn - r, r / es, es, cmp);*/
}
