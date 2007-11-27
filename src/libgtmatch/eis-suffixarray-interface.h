/*
  Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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

#ifndef EIS_SUFFIXARRAY_INTERFACE_H
#define EIS_SUFFIXARRAY_INTERFACE_H

#include "libgtmatch/eis-mrangealphabet.h"

struct fileReadState
{
  FILE *fp;
  MRAEnc *alphabet;
};

extern int
saReadBWT(void *state, Symbol *dest, size_t readLen, Env *env);

extern int
saGetOrigSeqSym(void *state, Symbol *dest, Seqpos pos, size_t len);

extern int
saReadSeqpos(void *src, Seqpos *dest, size_t len, Env *env);

extern DefinedSeqpos
reportSALongest(void *state);

extern MRAEnc *
newMRAEncFromSA(void *state, Env *env);

#endif
