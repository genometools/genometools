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

#ifndef ALPHA_H
#define ALPHA_H

#include <stdbool.h>
#include <stdio.h>

/* the alphabet class */
typedef struct Alpha Alpha; /* XXX: Alpha -> Alphabet */

Alpha*       alpha_new(void);
Alpha*       alpha_new_dna(void);
Alpha*       alpha_new_protein(void);
Alpha*       alpha_guess(const char *seq, unsigned long seqlen);
Alpha*       alpha_ref(Alpha*);
/* add the mapping of all given characters to the alphabet, the first character
   is the result of subsequent alpha_decode() calls  */
void         alpha_add_mapping(Alpha*, const char *characters);
char         alpha_decode(const Alpha*, unsigned int);
unsigned int alpha_encode(const Alpha*, char);
void         alpha_decode_seq(const Alpha*, char *out, char *in,
                              unsigned long length); /* in can be == out */
void         alpha_encode_seq(const Alpha*, char *out, char *in,
                              unsigned long length); /* in can be == out */
bool         alpha_char_is_valid(const Alpha*, char);
bool         alpha_is_compatible_with_alpha(const Alpha*, const Alpha*);
unsigned int alpha_size(const Alpha*);
void         alpha_delete(Alpha*);

#endif
