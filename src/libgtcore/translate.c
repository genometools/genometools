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
#include "libgtcore/codon.h"
#include "libgtcore/translate.h"

void translate_dna(Str *protein, const char *dnaseq, unsigned long dnalen,
                   unsigned int frame, Env *env)
{
  const char *dnaptr;
  assert(protein && !str_length(protein) && dnaseq && frame < 3);
  /* translate the DNA in forward direction */
  for (dnaptr = dnaseq + frame; dnaptr < dnaseq + dnalen - 2; dnaptr += 3) {
    str_append_char(protein, codon2amino(*dnaptr, *(dnaptr+1), *(dnaptr+2)),
                    env);
  }
}
