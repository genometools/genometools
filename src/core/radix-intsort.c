/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include "core/assert_api.h"
#include "core/stack-inlined.h"
#include "core/types_api.h"

#define GT_RADIX_ACCESS_KEY(MASK,SHIFT,PTR) ((*(PTR) >> (SHIFT)) & (MASK))
#define GT_RADIX_ACCESS_KEY_PAIR(MASK,SHIFT,PTR)\
        (((PTR)->a >> (SHIFT)) & (MASK))
#define GT_RADIX_ACCESS_UINT8(SHIFT,SP)\
        GT_RADIX_ACCESS_KEY(UINT8_MAX,SHIFT,SP)

#ifdef SKDEBUG
static void showbytewise(unsigned long value)
{
  int shift;

  for (shift = 56; shift >= 0; shift-=8)
  {
    printf("%lu ",GT_RADIX_ACCESS_UINT8(shift,&value));
    if (shift > 0)
    {
      printf(" ");
    }
  }
}
#endif

/* if sorting for table large than UINT_MAX is required, set the following
   type to unsigned long */

typedef unsigned int Countbasetype;

static void gt_radixsort_GtUlong_linear_phase(size_t shift,
                                              size_t maxvalue,
                                              Countbasetype *count,
                                              GtUlong *source,
                                              GtUlong *dest,
                                              unsigned long len)
{
  Countbasetype *cptr, idx;
  GtUlong *sptr;

  /* count occurences of every byte value */
  gt_assert(len <= (unsigned long) UINT_MAX);
  for (cptr = count; cptr <= count + maxvalue; cptr++)
  {
    *cptr = 0;
  }
  for (sptr = source; sptr < source + len; sptr++)
  {
    count[GT_RADIX_ACCESS_KEY(maxvalue,shift,sptr)]++;
  }

  /* compute partial sums */
  for (cptr = count+1; cptr <= count + maxvalue; cptr++)
  {
    *cptr += *(cptr-1);
  }

  /* fill dest with the right values in the right place */
  for (sptr = source + len - 1; sptr >= source; sptr--)
  {
    idx = --count[GT_RADIX_ACCESS_KEY(maxvalue,shift,sptr)];
    dest[idx] = *sptr;
  }
}

void gt_radixsort_GtUlong_linear(bool smalltables,GtUlong *source,
                                 GtUlong *temp,unsigned long len)
{
  unsigned int iter;
  Countbasetype *count;
  size_t basesize, maxvalue;

  gt_assert(source != NULL && temp != NULL);
  if (smalltables)
  {
    basesize = sizeof (uint8_t);
    maxvalue = UINT8_MAX;
  } else
  {
    basesize = sizeof (uint16_t);
    maxvalue = UINT16_MAX;
  }
  count = gt_malloc(sizeof(*count) * (maxvalue+1));
  for (iter = 0; iter <(unsigned int) (sizeof(unsigned long)/basesize);
       iter++)
  {
    GtUlong *ptr;

    /*printf("phase %lu\n",(unsigned long) (iter * CHAR_BIT * basesize));*/
    gt_radixsort_GtUlong_linear_phase (iter * CHAR_BIT * basesize, maxvalue,
                                       count, source, temp, len);
    ptr = source;
    source = temp;
    temp = ptr;
  }
  gt_free(count);
}

void gt_merge_sorted_inplace(GtUlong *leftpart,GtUlong *rightpart,
                             unsigned long len)
{
  GtUlong tmp, *left = leftpart, *right = rightpart;

  while(left < right && right < rightpart + len)
  {
    if (*left > *right)
    {
      tmp = *left;
      *left = *right;
      *right = tmp;
      left++;
      right++;
    } else
    {
      left++;
    }
  }
}

static void gt_radixsort_GtUlongPair_linear_phase(size_t shift,
                                                  size_t maxvalue,
                                                  Countbasetype *count,
                                                  GtUlongPair *source,
                                                  GtUlongPair *dest,
                                                  unsigned long len)
{
  Countbasetype *cptr, idx;
  GtUlongPair *sptr;

  /* count occurences of every byte value */
  gt_assert(len <= (unsigned long) UINT_MAX);
  for (cptr = count; cptr <= count + maxvalue; cptr++)
  {
    *cptr = 0;
  }
  for (sptr = source; sptr < source + len; sptr++)
  {
    count[GT_RADIX_ACCESS_KEY_PAIR(maxvalue,shift,sptr)]++;
  }

  /* compute partial sums */
  for (cptr = count+1; cptr <= count + maxvalue; cptr++)
  {
    *cptr += *(cptr-1);
  }

  /* fill dest with the right values in the right place */
  for (sptr = source + len - 1; sptr >= source; sptr--)
  {
    idx = --count[GT_RADIX_ACCESS_KEY_PAIR(maxvalue,shift,sptr)];
    dest[idx] = *sptr;
  }
}

void gt_radixsort_GtUlongPair_linear(bool smalltables,GtUlongPair *source,
                                     GtUlongPair *temp,unsigned long len)
{
  unsigned int iter;
  Countbasetype *count;
  size_t basesize, maxvalue;

  gt_assert(source != NULL && temp != NULL);
  if (smalltables)
  {
    basesize = sizeof (uint8_t);
    maxvalue = UINT8_MAX;
  } else
  {
    basesize = sizeof (uint16_t);
    maxvalue = UINT16_MAX;
  }
  count = gt_malloc(sizeof(*count) * (maxvalue+1));
  for (iter = 0; iter <(unsigned int) (sizeof(unsigned long)/basesize);
       iter++)
  {
    GtUlongPair *ptr;

    /*printf("phase %lu\n",(unsigned long) (iter * CHAR_BIT * basesize));*/
    gt_radixsort_GtUlongPair_linear_phase (iter * CHAR_BIT * basesize, maxvalue,
                                           count, source, temp, len);
    ptr = source;
    source = temp;
    temp = ptr;
  }
  gt_free(count);
}

static void gt_radix_phase_GtUlong_recursive(size_t offset,
                                             GtUlong *source,
                                             GtUlong *dest,
                                             unsigned long len)
{
  unsigned long idx, s, c, *sp, *cp;
  const size_t maxoffset = sizeof (unsigned long) - 1;
  unsigned long count[UINT8_MAX+1] = {0};
  const size_t shift = (maxoffset - offset) * CHAR_BIT;

  /* count occurences of every byte value */
  for (sp = source; sp < source+len; sp++)
  {
    count[GT_RADIX_ACCESS_UINT8(shift,sp)]++;
  }
  /* compute partial sums */
  for (s = 0, cp = count; cp <= count + UINT8_MAX; cp++)
  {
    c = *cp;
    *cp = s;
    s += c;
  }
  /* fill dest with the right values in the right place */
  for (sp = source; sp < source+len; sp++)
  {
    dest[count[GT_RADIX_ACCESS_UINT8(shift,sp)]++] = *sp;
  }
  memcpy(source,dest,(size_t) sizeof (*source) * len);
  if (offset < maxoffset)
  {
    for (idx = 0; idx <= UINT8_MAX; idx++)
    {
      unsigned long newleft = (idx == 0) ? 0 : count[idx-1];
      /* |newleft .. count[idx]-1| = count[idx]-1-newleft+1
                                   = count[idx]-newleft > 1
      => count[idx] > newleft + 1 */
      if (newleft+1 < count[idx])
      {
        gt_radix_phase_GtUlong_recursive(offset+1,
                                         source+newleft,
                                         dest+newleft,
                                         count[idx]-newleft);
      }
    }
  }
}

void gt_radixsort_GtUlong_recursive(GtUlong *source, GtUlong *dest,
                                    unsigned long len)
{
  gt_radix_phase_GtUlong_recursive(0,source,dest,len);
}

typedef struct
{
  GtUlong *left;
  unsigned long len;
  uint8_t shift;
} GtRadixsort_stackelem;

GT_STACK_DECLARESTRUCT(GtRadixsort_stackelem,512);

static void gt_radixsort_GtUlong_initstack(GtStackGtRadixsort_stackelem *stack,
                                           GtUlong *source,
                                           GtUlong *dest,
                                           unsigned long len)
{
  GtRadixsort_stackelem tmpelem;
  unsigned long idx, s, c, *sp, *cp, newleft, count[UINT16_MAX+1];
  const size_t mask = UINT16_MAX, shift = (size_t) 48;

  GT_STACK_INIT(stack,64UL);
  for (idx=0; idx<=UINT16_MAX; idx++)
  {
    count[idx] = 0;
  }
  for (sp = source; sp < source + len; sp++)
  {
    count[GT_RADIX_ACCESS_KEY(mask,shift,sp)]++;
  }
  for (s = 0, cp = count; cp <= count + UINT16_MAX; cp++)
  {
    c = *cp;
    *cp = s;
    s += c;
  }
  /* fill dest with the right values in the right place */
  for (sp = source; sp < source + len; sp++)
  {
    dest[count[GT_RADIX_ACCESS_KEY(mask,shift,sp)]++] = *sp;
  }
  memcpy(source,dest,(size_t) sizeof (*source) * len);
  for (idx = 0; idx <= UINT16_MAX; idx++)
  {
    newleft = (idx == 0) ? 0 : count[idx-1];
    /* |newleft .. count[idx]-1| = count[idx]-1-newleft+1
                                 = count[idx]-newleft > 1
        => count[idx] > newleft + 1
    */
    if (newleft+1 < count[idx])
    {
      tmpelem.shift = (uint8_t) 40;
      tmpelem.left = source + newleft;
      tmpelem.len = count[idx] - newleft;
      GT_STACK_PUSH(stack,tmpelem);
    }
  }
}

void gt_radixsort_GtUlong_divide(GtUlong *source, GtUlong *dest,
                                 unsigned long len)
{
  GtStackGtRadixsort_stackelem stack;
  GtRadixsort_stackelem tmpelem, current;
  unsigned long idx, s, c, *sp, *cp, newleft, count[UINT8_MAX+1] = {0};
  const bool simple = false;

  if (simple)
  {
    GT_STACK_INIT(&stack,64UL);
    tmpelem.shift = (sizeof (unsigned long) - 1) * CHAR_BIT;
    tmpelem.left = source;
    tmpelem.len = len;
    GT_STACK_PUSH(&stack,tmpelem);
  } else
  {
    gt_radixsort_GtUlong_initstack(&stack, source, dest, len);
  }
  while (!GT_STACK_ISEMPTY(&stack))
  {
    current = GT_STACK_POP(&stack);
    /* count occurences of every byte value */
    for (sp = current.left; sp < current.left+current.len; sp++)
    {
      count[GT_RADIX_ACCESS_UINT8(current.shift,sp)]++;
    }
    /* compute partial sums */
    for (s = 0, cp = count; cp <= count + UINT8_MAX; cp++)
    {
      c = *cp;
      *cp = s;
      s += c;
    }
    /* fill dest with the right values in the right place */
    for (sp = current.left; sp < current.left+current.len; sp++)
    {
      dest[count[GT_RADIX_ACCESS_UINT8(current.shift,sp)]++] = *sp;
    }
    memcpy(current.left,dest,(size_t) sizeof (*source) * current.len);
    if (current.shift > 0)
    {
      for (idx = 0; idx <= UINT8_MAX; idx++)
      {
        newleft = (idx == 0) ? 0 : count[idx-1];
        /* |newleft .. count[idx]-1| = count[idx]-1-newleft+1
                                     = count[idx]-newleft > 1
        => count[idx] > newleft + 1 */
        if (newleft+1 < count[idx])
        {
          tmpelem.shift = current.shift - CHAR_BIT;
          tmpelem.left = current.left + newleft;
          tmpelem.len = count[idx] - newleft;
          GT_STACK_PUSH(&stack,tmpelem);
        }
      }
    }
    memset(count,0,(size_t) sizeof (*count) * (UINT8_MAX+1));
  }
  GT_STACK_DELETE(&stack);
}
