/*
  Copyright (c) 2017 Julian Elvers <julian.elvers@gmail.com>
  Copyright (c) 2017 Center for Bioinformatics, University of Hamburg

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

#ifndef K_ENV_ALPHABET_H
#define K_ENV_ALPHABET_H

#include <stdbool.h>
#include "core/alphabet_api.h"
#include "core/codetype.h"

/* Holds all necessary information regarding the used alphabet, including
   the <GtAlphabet> itself. */
typedef struct GtKenvAlphabet GtKenvAlphabet;

/* Returns a new <GtKenvAlphabet> object representing a dna alphabet,
   if <use_dna> is true, a protein alphabet otherwise. Also prepares
   the needed power values up to alph_size^(q_val - 1). */
GtKenvAlphabet *gt_kenv_alphabet_new(bool use_dna, unsigned int q_val);

/* Deletes the given <GtKenvAlphabet> object. */
void gt_kenv_alphabet_delete(GtKenvAlphabet *kenv_alph);

/* Returns the <GtAlphabet> stored in <GtKenvAlphabet>. */
GtAlphabet *gt_kenv_alphabet_get_alphabet(const GtKenvAlphabet *kenv_alph);

/* Returns the <GtCodetype> Array representing the power values up to
   alph_size^(q_val - 1). */
const GtCodetype *gt_kenv_alphabet_get_alph_powers(
                                            const GtKenvAlphabet *kenv_alph);

/* Returns the size, resp. num_of_chars, of the given <GtKenvAlphabet>. */
unsigned int gt_kenv_alphabet_get_alph_size(const GtKenvAlphabet *kenv_alph);

/* Returns the state of the <GtKenvAlphabet>. True, if <kenv_alph> contains
   a dna alphabet, false if it describes a protein alphabet. */
bool gt_kenv_alphabet_is_dna(const GtKenvAlphabet *kenv_alph);

/* Returns the array of symbols from the <GtAlphabet> inside <kenv_alph> such
   that the index of the character equals its encoding. The array contains
   all possible characters, regardless whether they are encodable or not. */
const GtUchar* gt_kenv_alphabet_get_symbolmap(const GtKenvAlphabet *kenv_alph);

/* Returns an array of the characters in the <GtAlphabet> inside <kenv_alph>
   such that the position of the character equals its encoding. */
const GtUchar* gt_kenv_alphabet_get_characters(const GtKenvAlphabet *kenv_alph);

#endif
