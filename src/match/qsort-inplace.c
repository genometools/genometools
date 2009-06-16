#include <stdbool.h>
#include "core/unused_api.h"
#include "core/minmax.h"
#include "match/divmodmul.h"

/*
 * Qsort routine from Bentley & McIlroy's "Engineering a Sort Function".
 */

#define SWAP(A,B)\
        if ((A) != (B))\
        {\
          tmp = *(A);\
          *(A) = *(B);\
          *(B) = tmp;\
        }

#define VECSWAP(A,B,N)\
        aptr = A;\
        bptr = B;\
        while ((N)-- > 0)\
        {\
          tmp = *aptr;\
          *aptr++ = *bptr;\
          *bptr++ = tmp;\
        }

typedef unsigned int Sorttype;

static int qsortcmp (Sorttype * a,Sorttype * b,GT_UNUSED void *data)
{
  if ((*a) < (*b))
  {
    return -1;
  }
  if ((*a) > (*b))
  {
    return 1;
  }
  return 0;
}

static inline Sorttype *med3 (Sorttype * a, Sorttype * b, Sorttype * c,
                              void *data)
{
  return qsortcmp (a, b, data) < 0 ?
    (qsortcmp (b, c, data) <
     0 ? b : (qsortcmp (a, c, data) <
              0 ? c : a)) : (qsortcmp (b, c,
                                              data) >
                             0 ? b
                             : (qsortcmp (a, c, data) <
                                0 ? a : c));
}

void gt_inlined_qsort_r (Sorttype * a,unsigned long n,void *data)
{
  Sorttype tmp, *pa, *pb, *pc, *pd, *pl, *pm, *pn, *aptr, *bptr;
  unsigned long d, minval;
  int r;
  bool swapped;

  while (1)
  {
    swapped = false;
    if (n < 7UL)
    {
      for (pm = a + 1; pm < a + n; pm++)
      {
        for (pl = pm;
             pl > a && qsortcmp (pl - 1, pl, data) > 0;
             pl--)
        {
          SWAP (pl, pl - 1);
        }
      }
      return;
    }
    pm = a + DIV2 (n);
    if (n > 7UL)
    {
      pl = a;
      pn = a + n - 1;
      if (n > 40UL)
      {
        d = DIV8 (n);
        pl = med3 (pl, pl + d, pl + MULT2 (d), data);
        pm = med3 (pm - d, pm, pm + d, data);
        pn = med3 (pn - MULT2 (d), pn - d, pn, data);
      }
      pm = med3 (pl, pm, pn, data);
    }
    SWAP (a, pm);
    pa = pb = a + 1;
    pc = pd = a + n - 1;
    while (1)
    {
      while (pb <= pc
             && (r = qsortcmp (pb, a, data)) <= 0)
      {
        if (r == 0)
        {
          swapped = true;
          SWAP (pa, pb);
          pa++;
        }
        pb++;
      }
      while (pb <= pc
             && (r = qsortcmp (pc, a, data)) >= 0)
      {
        if (r == 0)
        {
          swapped = true;
          SWAP (pc, pd);
          pd--;
        }
        pc--;
      }
      if (pb > pc)
      {
        break;
      }
      SWAP (pb, pc);
      swapped = true;
      pb++;
      pc--;
    }
    if (!swapped)
    {                                  /* Switch to insertion sort */
      for (pm = a + 1; pm < a + n; pm++)
      {
        for (pl = pm;
             pl > a && qsortcmp (pl - 1, pl, data) > 0;
             pl--)
        {
          SWAP (pl, pl - 1);
        }
      }
      return;
    }
    pn = a + n;
    minval = MIN ((unsigned long) (pa - a), (unsigned long) (pb - pa));
    VECSWAP (a, pb - minval, minval);
    minval = MIN ((unsigned long) (pd - pc), (unsigned long) (pn - pd - 1));
    VECSWAP (pb, pn - minval, minval);
    if ((minval = (unsigned long) (pb - pa)) > 1UL)
    {
      gt_inlined_qsort_r (a, minval, data);
    }
    if ((minval = (unsigned long) (pd - pc)) > 1UL)
    {
      /*
         Iterate rather than recurse to save stack space
       */
      a = pn - minval;
      n = minval;
    } else
    {
      break;
    }
  }
}
