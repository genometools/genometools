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

#ifndef QGRAM_H
#define QGRAM_H

#include "libgtcore/alpha.h"
#include "libgtcore/array.h"

/* Encodes a word <w> of length <q> over an alphabet of size <alpha_size> as a
   unique number. */
unsigned long qgram_encode(const char *w, unsigned long q,
                           unsigned long alphabet_size);

/* Computes the next encoding. */
unsigned long qgram_step(unsigned long current_code, char previous, char next,
                         unsigned long alphabet_size,
                         unsigned long
                         alpha_size_raised_to_the_power_of_q_minus_1);

/* Computes all q-grams of the given <encoded_seq> (over an alphabet of size
   <alpha_size>) and stores them in the array <qgrams>. */
void          qgram_compute(Array *qgrams, const char *encoded_seq,
                            unsigned long seqlen, unsigned long alpha_size,
                            unsigned int q);

/* Decode <code> (which encodes a qgram of length <q> over alphabet <alpha) and
   store the resulting string in <qgram>. */
void          qgram_decode(char *qgram, unsigned long code, unsigned long q,
                           Alpha *alpha);

#endif
