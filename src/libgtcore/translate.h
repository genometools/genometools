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

#ifndef TRANSLATE_H
#define TRANSLATE_H

#include "libgtcore/str.h"

void translate_dna(Str*, const char*, unsigned long, unsigned int frame);

/* Translate <dna_sequence> of length <seqlen> in all three reading frames.
   The translations are stored in <frame1>, <frame2>, and <frame3>.
   The necessary space for the frames is allocated by translate_all_frames(),
   it is the responsibility of the caller to free it.
   All characters different from 'acgtu' (in lower and upper case) are mapped to
   't' before the translation is performed. */
void translate_all_frames(char **frame1, char **frame2, char **frame3,
                          const char *dna_sequence, unsigned long seqlen);

#endif
