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

#include "core/types_api.h"
#include "core/ma_api.h"
#include "squarededist.h"

static unsigned long gt_squarededistunit2 (const GtUchar *u, unsigned long m,
                                        const GtUchar *v, unsigned long n)
{
  unsigned long val, we, nw, *ecol, *ecolptr;
  const GtUchar *uptr, *vptr;

  ecol = gt_malloc(sizeof *ecol * (m+1));
  for (*ecol = 0, ecolptr = ecol+1, uptr = u; uptr < u + m; ecolptr++, uptr++)
  {
    *ecolptr = *(ecolptr-1) + 1;
  }
  for (vptr = v; vptr < v + n; vptr++)
  {
    nw = *ecol;
    *ecol = nw + 1;
    for (ecolptr = ecol+1, uptr = u; uptr < u + m; ecolptr++, uptr++)
    {
      we = *ecolptr;
      *ecolptr = *(ecolptr-1) + 1;
      if (*uptr == *vptr)
      {
        val = nw;
      } else
      {
        val = nw + 1;
      }
      if (val < *ecolptr)
      {
        *ecolptr = val;
      }
      if ((val = we + 1) < *ecolptr)
      {
        *ecolptr = val;
      }
      nw = we;
    }
  }
  val = *(ecolptr-1);
  gt_free(ecol);
  return val;
}

unsigned long gt_squarededistunit (const GtUchar *u, unsigned long m,
                                const GtUchar *v, unsigned long n)
{
  if (m < n)
  {
    return gt_squarededistunit2(u,m,v,n);
  }
  return gt_squarededistunit2(v,n,u,m);
}
