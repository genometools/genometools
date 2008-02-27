/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef SFX_ENUMCODES_H
#define SFX_ENUMCODES_H

#include "intcode-def.h"
#include "encseq-def.h"
#include "seqpos-def.h"

typedef struct Enumcodeatposition Enumcodeatposition;

typedef struct
{
  unsigned int maxprefixindex;
  Seqpos position;
} Specialcontext;

Enumcodeatposition *newEnumcodeatposition(const Encodedsequence *encseq,
                                          Readmode readmode,
                                          unsigned int prefixlength,
                                          unsigned int numofchars);

bool nextEnumcodeatposition(Specialcontext *specialcontext,
                            Enumcodeatposition *ecp);

void freeEnumcodeatposition(Enumcodeatposition **ecp);

Codetype computefilledqgramcode(const Enumcodeatposition *ecp,
                                unsigned int prefixindex,
                                Seqpos pos);

bool computefilledqgramcodestopatmax(Codetype *code,
                                     const Enumcodeatposition *ecp,
                                     unsigned int prefixindex,
                                     Seqpos pos,
                                     Codetype stopcode);
#endif
