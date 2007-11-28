/* This module has been derived from the file msort.c of the glibc 2.3.2
   by Gordon Gremme <gremme@zbh.uni-hamburg.de>.

   Copyright (C) 1992,95-97,99,2000,01,02 Free Software Foundation, Inc.
   Written by Mike Haertel, September 1988.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.
*/

#include <string.h>
#include "libgtcore/ma.h"
#include "libgtcore/msort.h"
#include "libgtcore/xansi.h"

static void msort_r_withbuf(void *base, size_t numofelems, size_t size,
                            void *cmpinfo,
                            int (*cmpfunc)(void *, const void *, const void *),
                            char *buf)
{
  char *tmpbuf;
  char *b1, *b2;
  size_t n1, n2;

  if (numofelems <= (size_t) 1)
    return;

  n1 = numofelems >> 1; /* /2 */
  n2 = numofelems - n1;
  b1 = base;
  b2 = (char *) base + (n1 * size);

  msort_r_withbuf(b1, n1, size, cmpinfo, cmpfunc, buf);
  msort_r_withbuf(b2, n2, size, cmpinfo, cmpfunc, buf);

  tmpbuf = buf;

  while (n1 > 0 && n2 > 0) {
    if ((*cmpfunc)(cmpinfo, b1, b2) <= 0) {
      memcpy(tmpbuf, b1, size);
      tmpbuf += size;
      b1 += size;
      n1--;
    }
    else {
      memcpy(tmpbuf, b2, size);
      tmpbuf += size;
      b2 += size;
      n2--;
    }
  }

  if (n1 > 0)
    memcpy (tmpbuf, b1, n1 * size);

  memcpy (base, buf, (numofelems - n2) * size);
}

void msort_r(void *base, size_t nmemb, size_t size, void *comparinfo,
             int (*compar)(void *, const void *, const void *))
{
  void *buf;
  buf = ma_malloc(size * nmemb);
  msort_r_withbuf (base, nmemb, size, comparinfo, compar, buf);
  ma_free(buf);
}

int non_r_cmpfunc(void *compar, const void *a, const void *b)
{
  int (*cmpfunc)(const void *, const void *) = compar, rval;
  rval = (*cmpfunc)(a, b);
  return rval;
}

void msort(void *base, size_t nmemb, size_t size,
           int (*compar)(const void *, const void *))
{
  msort_r(base, nmemb, size, compar, non_r_cmpfunc);
}
