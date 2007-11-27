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

#include <stdio.h>
#include <stdlib.h>
#include "libgtcore/array.h"

int array_compare(const Array *a,const Array *b,
                  int(*compar)(const void *, const void *))
{
  unsigned long idx, size_a, size_b;
  int cmp;

  size_a = array_size(a);
  size_b = array_size(b);
  if (size_a < size_b)
  {
    fprintf(stderr,"array_size(a) = %lu < %lu = array_size(b)\n",
                  size_a,
                  size_b);
    return -1;
  }
  if (size_a > size_b)
  {
    fprintf(stderr,"array_size(a) = %lu > %lu = array_size(b)\n",
                  size_a,
                  size_b);
    return 1;
  }
  for (idx=0; idx < size_a; idx++)
  {
    cmp = compar(array_get(a,idx),array_get(b,idx));
    if (cmp != 0)
    {
      return cmp;
    }
  }
  return 0;
}
