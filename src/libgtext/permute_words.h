/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef PERMUTE_WORDS_H
#define PERMUTE_WORDS_H

/* module for word permutation */

/* initializes the word w of length l with 0's */
void permute_word_init(char *w, unsigned long l);

/* computes the next permutation and returns True if one exists,
   False otherwise */
unsigned int permute_word_next(char *w, unsigned long l,
                               unsigned long alphabet_size);

#if 0
  a typical use:

  permute_word_init(w, l);
  do {
    /* ... */
  } while (permute_word_next(w, l, alphabet_size));

#endif

#endif
