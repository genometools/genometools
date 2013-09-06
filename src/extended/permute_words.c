/*
  Copyright (c) 2006 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/assert_api.h"
#include "extended/permute_words.h"

void gt_permute_word_init(char *w, GtUword l)
{
  memset(w, 0, l);
}

unsigned int gt_permute_word_next(char *w, GtUword l,
                               GtUword alphabet_size)
{
  GtUword i;
  unsigned int next = 1;

  gt_assert(w && l && alphabet_size);

  for (i = 0; i < l; i++) {
    if (w[i] < alphabet_size - 1) {
      w[i]++;
      break;
    }
    else {
      w[i] = 0;
      if (i == l - 1)
        next = 0;
    }
  }

  return next;
}
