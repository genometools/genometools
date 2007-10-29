/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtext/matchcount.h"

void matchcount(const char *u, int m, const char *v, int n, int k,
                ProcMC proc_match_count)
{
  int diagonal, i, j, mc;

  /* compute match count for every diagonal */
  for (diagonal = -(m - k); diagonal <= n - k; diagonal++) {
    /* set i and j subject to diagonal */
    if (diagonal < 0) {
      i = -diagonal;
      j = 0;
    }
    else if (diagonal == 0) {
      i = 0;
      j = 0;
    }
    else { /* diagonal > 0 */
      i = 0;
      j = diagonal;
    }

    for (mc = 0; i < m && j < n; i++, j++) {
      /* computation */
      if (i < k  || j < k) {
        if (u[i] == v[j])
          mc++;
      }
      else {
        if (u[i-k] == v[j-k])
          mc--;
        if (u[i] == v[j])
          mc++;
      }

      /* process match count */
      if (i >= k-1 && j >= k-1) {
        assert(mc >= 0); /* match count is always >= 0 */
        proc_match_count(i-(k-1), j-(k-1), mc);
      }
    }
  }
}
