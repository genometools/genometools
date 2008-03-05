/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <assert.h>
#include <limits.h>
#include "libgtcore/bittab.h"
#include "libgtext/string_matching.h"

void string_matching_bmh(const char *s, unsigned long n,
                         const char *p, unsigned long m,
                         ProcessMatch process_match, void *data)
{
  unsigned long i, j, pos, d[UCHAR_MAX];
  assert(s && p);
  if (m > n || !m || !n) /* no match possible */
    return;
  /* preprocessing */
  for (i = 0; i < UCHAR_MAX; i++)
    d[i] = m;
  for (j = 0; j < m-1; j++)
    d[(unsigned) p[j]] = m-(j+1);
  /* searching */
  pos = 0;
  while (pos <= n-m) {
    j = m;
    while (j > 0 && s[pos+j-1] == p[j-1])
      j--;
    if (j == 0 && process_match) {
      if (process_match(pos, data))
        return;
    }
    pos += d[(unsigned) s[pos+m-1]];
  }
}

void string_matching_shift_and(const char *s, unsigned long n,
                               const char *p, unsigned long m,
                               ProcessMatch process_match, void *data)
{
  Bittab *D, *B[UCHAR_MAX] = { NULL };
  unsigned long i, j;
  assert(s && p);
  if (m > n || !m || !n) /* no match possible */
    return;
  /* preprocessing */
  for (j = 0; j < m; j++) {
    if (!B[(unsigned) p[j]])
      B[(unsigned) p[j]] = bittab_new(m);
    bittab_set_bit(B[(unsigned) p[j]], j);
  }
  /* searching */
  D = bittab_new(m);
  for (i = 0; i < n; i++) {
    bittab_shift_left_equal(D);
    bittab_set_bit(D, 0);
    if (B[(unsigned) s[i]])
      bittab_and_equal(D, B[(unsigned) s[i]]);
    if (bittab_bit_is_set(D, m - 1) && process_match) {
      if (process_match(i - m + 1, data))
        break;
    }
  }
  /* free */
  for (i = 0; i < UCHAR_MAX; i++)
    bittab_delete(B[i]);
  bittab_delete(D);
}
