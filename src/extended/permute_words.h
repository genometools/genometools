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

#ifndef PERMUTE_WORDS_H
#define PERMUTE_WORDS_H

#include "core/types_api.h"

/* module for word permutation */

/* initializes the word w of length l with 0's */
void gt_permute_word_init(char *w, GtUword l);

/* computes the next permutation and returns True if one exists,
   False otherwise */
unsigned int gt_permute_word_next(char *w, GtUword l,
                               GtUword alphabet_size);

#if 0
  a typical use:

  gt_permute_word_init(w, l);
  do {
    /* ... */
  } while (gt_permute_word_next(w, l, alphabet_size));

#endif

#endif
