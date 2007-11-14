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
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libgtcore/undef.h"
#include "libgtext/multiset_matching.h"

void multiset_matching(unsigned char *multiset_string,
                       unsigned long multiset_size, unsigned char *text,
                       unsigned long text_length, void *data,
                       void (*procmatchfunc)(unsigned long pos, void *data))
{
  unsigned long i;
  long multiset[UCHAR_MAX],
       exhausted_characters = 0,
       alphabet_size = 0; /* with respect to the multiset */

  /* init multiset */
  for (i = 0; i < UCHAR_MAX; i++)
    multiset[i] = UNDEF_CHAR;

  /* construct the multiset and determine the alphabet size */
  for (i = 0; i < multiset_size; i++) {
    if (multiset[multiset_string[i]] == UNDEF_CHAR) {
      multiset[multiset_string[i]] = 0;
      alphabet_size++;
    }
    multiset[multiset_string[i]]++;
  }

  /* matching (sliding window) */
  for (i = 0; i < text_length; i++) {
    /* undo first character if necessary */
    if (i >= multiset_size) {
      if (multiset[text[i - multiset_size]] == 0)
        exhausted_characters--;
      if (multiset[text[i - multiset_size]] != UNDEF_CHAR)
        multiset[text[i - multiset_size]]++;
      if (multiset[text[i - multiset_size]] == 0)
        exhausted_characters++;
    }

    /* match next character */
    if (multiset[text[i]] == 0)
      exhausted_characters--;
    if (multiset[text[i]] != UNDEF_CHAR)
      multiset[text[i]]--;
    if (multiset[text[i]] == 0)
      exhausted_characters++;

    assert(exhausted_characters <= alphabet_size);
    if (exhausted_characters == alphabet_size) /* match found */
      procmatchfunc(i - multiset_size + 1, data);
  }
}
