/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <string.h>
#include "libgtext/permute_words.h"

void permute_word_init(char *w, unsigned long l)
{
  memset(w, 0, l);
}

unsigned int permute_word_next(char *w, unsigned long l,
                               unsigned long alphabet_size)
{
  unsigned long i;
  unsigned int next = 1;

  assert(w && l && alphabet_size);

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
