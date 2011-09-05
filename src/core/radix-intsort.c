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
#include "core/endianess_api.h"
#include "core/stack-inlined.h"
#include "core/types_api.h"

/* replaced byte with offset to avoid * 8 operation in loop */

/* be careful only works for litte endian byte order */

static void gt_radix_phase_GtUlong(unsigned int offset,
                                   GtUlong *source,
                                   GtUlong *dest,
                                   unsigned long len)
{
  unsigned long idx, s, c, *sp, *cp, count[256] = {0};
  const size_t increment = sizeof (*source);
  uint8_t *bp;

  /* count occurences of every byte value */
  bp = ((uint8_t *) source) + offset;
  for (idx = 0; idx < len; idx++, bp += increment)
  {
    count[*bp]++;
  }

  /* compute partial sums */
  for (s = 0, cp = count, idx = 0; idx < 256UL; idx++, cp++)
  {
    c = *cp;
    *cp = s;
    s += c;
  }

  /* fill dest with the right values in the right place */
  bp = ((uint8_t *) source) + offset;
  for (sp = source; sp < source + len; bp += increment, sp++)
  {
    dest[count[*bp]++] = *sp;
  }
}

void gt_radixsort_GtUlong_linear(GtUlong *source, GtUlong *temp,
                                 unsigned long len)
{
  gt_assert(gt_is_little_endian());
  /* allocate heap memory to avoid the need of additional parameter */
  gt_assert(temp != NULL && source != NULL);
  gt_radix_phase_GtUlong(0, source, temp, len);
  gt_radix_phase_GtUlong(1U, temp, source, len);
  gt_radix_phase_GtUlong(2U, source, temp, len);
  gt_radix_phase_GtUlong(3U, temp, source, len);
#ifdef _LP64
  gt_radix_phase_GtUlong(4U, source, temp, len);
  gt_radix_phase_GtUlong(5U, temp, source, len);
  gt_radix_phase_GtUlong(6U, source, temp, len);
  gt_radix_phase_GtUlong(7U, temp, source, len);
#endif
}

#define GT_RADIX_ACCESS_CHAR(SHIFT,SP) ((*(SP) >> (SHIFT)) & UINT8_MAX)

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
    count[GT_RADIX_ACCESS_CHAR(shift,sp)]++;
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
    dest[count[GT_RADIX_ACCESS_CHAR(shift,sp)]++] = *sp;
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

void gt_radixsort_GtUlong_divide(GtUlong *source, GtUlong *dest,
                                 unsigned long len)
{
  GtStackGtRadixsort_stackelem stack;
  GtRadixsort_stackelem tmpelem, current;
  unsigned long idx, s, c, *sp, *cp, newleft, count[UINT8_MAX+1] = {0};

  GT_STACK_INIT(&stack,64UL);
  tmpelem.shift = (sizeof (unsigned long) - 1) * CHAR_BIT;
  tmpelem.left = source;
  tmpelem.len = len;
  GT_STACK_PUSH(&stack,tmpelem);
  while (!GT_STACK_ISEMPTY(&stack))
  {
    current = GT_STACK_POP(&stack);
    /* count occurences of every byte value */
    for (sp = current.left; sp < current.left+current.len; sp++)
    {
      count[GT_RADIX_ACCESS_CHAR(current.shift,sp)]++;
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
      dest[count[GT_RADIX_ACCESS_CHAR(current.shift,sp)]++] = *sp;
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

/* assume that the first element in GtUlongPair is the sort key */

static void gt_radix_phase_GtUlongPair(unsigned int offset,
                                       GtUlongPair *source,
                                       GtUlongPair *dest,
                                       unsigned long len)
{
  unsigned long idx, *cp, s, c, count[256] = {0};
  GtUlongPair *sp;
  const size_t increment = sizeof (*source);
  uint8_t *bp;

  gt_assert(gt_is_little_endian());
  /* count occurences of every byte value */
  bp = ((uint8_t *) source) + offset;
  for (idx = 0; idx < len; idx++, bp += increment)
  {
    count[*bp]++;
  }

  /* compute partial sums */
  for (s = 0, cp = count, idx = 0; idx < 256UL; idx++, cp++)
  {
    c = *cp;
    *cp = s;
    s += c;
  }

  /* fill dest with the right values in the right place */
  bp = ((uint8_t *) source) + offset;
  for (sp = source; sp < source + len; bp += increment, sp++)
  {
    dest[count[*bp]++] = *sp;
  }
}

void gt_radixsort_GtUlongPair(GtUlongPair *source, GtUlongPair *temp,
                              unsigned long len)
{
  /* allocate heap memory to avoid the need of additional parameter */
  gt_assert(temp != NULL && source != NULL);
  gt_radix_phase_GtUlongPair(0, source, temp, len);
  gt_radix_phase_GtUlongPair(1U, temp, source, len);
  gt_radix_phase_GtUlongPair(2U, source, temp, len);
  gt_radix_phase_GtUlongPair(3U, temp, source, len);
#ifdef _LP64
  gt_radix_phase_GtUlongPair(4U, source, temp, len);
  gt_radix_phase_GtUlongPair(5U, temp, source, len);
  gt_radix_phase_GtUlongPair(6U, source, temp, len);
  gt_radix_phase_GtUlongPair(7U, temp, source, len);
#endif
}
