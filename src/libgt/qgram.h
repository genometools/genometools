/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef QGRAM_H
#define QGRAM_H

#include <libgt/alpha.h>
#include <libgt/array.h>

/* encodes a word w of length q over an alphabet of size alpha_size as a unique
   number */
unsigned long qgram_encode(const char *w, unsigned long q,
                           unsigned long alphabet_size);

/* computes the next encoding */
unsigned long qgram_step(unsigned long current_code, char previous, char next,
                         unsigned long alphabet_size,
                         unsigned long
                         alpha_size_raised_to_the_power_of_q_minus_1);

/* computes all q-grams of the given encoded_seq (over an alphabet of size
   alpha_size) and stores them in the array qgrams. */
void          qgram_compute(Array *qgrams, const char *encoded_seq,
                            unsigned long seqlen, unsigned long alpha_size,
                            unsigned int q, Env*);

#endif
