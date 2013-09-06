/*
  Copyright (c) 2006-2009 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "core/minmax.h"
#include "extended/linearedist.h"

static void fillDPtable(GtUword *dptable,
                        const char *u, GtUword n,
                        const char *v, GtUword m)
{
  GtUword i, j , nw, we;
  for (i = 0; i <= n; i++)
    dptable[i] = i;
  for (j = 1; j <= m; j++) {
    nw = dptable[0];
    dptable[0] = j;
    for (i = 1; i <= n; i++) {
      we = dptable[i];
      dptable[i] = nw + (u[i-1] == v[j-1] ? 0 : 1); /* replacement */
      if (dptable[i-1] + 1 < dptable[i]) /* deletion */
        dptable[i] = dptable[i-1] + 1;
      if (we + 1 < dptable[i]) /* insertion */
        dptable[i] = we + 1;
      nw = we;
    }
  }
}

GtUword gt_calc_linearedist(const char *u, GtUword n,
                                  const char *v, GtUword m)
{
  GtUword *dptable, edist;
  dptable = gt_malloc(sizeof (GtUword) * (MIN(n,m) + 1));
  fillDPtable(dptable, n <= m ? u : v, MIN(n,m), n <= m ? v : u, MAX(n,m));
  edist = dptable[MIN(n,m)];
  gt_free(dptable);
  return edist;
}
