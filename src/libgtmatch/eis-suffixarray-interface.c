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

#include "libgtcore/chardef.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/eis-mrangealphabet.h"
#include "libgtmatch/eis-suffixarray-interface.h"

extern int
saReadBWT(void *state, Symbol *dest, size_t len, Error *err)
{
  struct suffixarrayReadState *saRState;
  assert(state);
  saRState = state;
  return MRAEncReadAndTransform(saRState->alphabet,
                                saRState->sa->bwttabstream.fp, len, dest);
}

extern int
saGetOrigSeqSym(void *state, Symbol *dest, Seqpos pos, size_t len)
{
  const Suffixarray *sa;
  size_t i;
  assert(state);
  sa = state;
  assert(sa->encseq);
  for (i = 0; i < len; ++i)
    dest[i] = getencodedchar(sa->encseq, pos + i, sa->readmode);
  return len;
}

DECLAREREADFUNCTION(Seqpos)

extern int
saReadSeqpos(void *src, Seqpos *dest, size_t len, Error *err)
{
  int rv = 1;
  size_t i;
  Suffixarray *sa = src;
  assert(sa);
  for (i = 0; i  < len; ++i)
    if ((rv = readnextSeqposfromstream(dest + i, &sa->suftabstream, err)) != 1)
      break;
  return ((rv == 1)? i : 0);
}

extern DefinedSeqpos
reportSALongest(void *state)
{
  Suffixarray *sa = state;
  assert(sa);
  return sa->longest;
}

extern MRAEnc *
newMRAEncFromSA(const void *state)
{
  MRAEnc *alphabet;
  const Suffixarray *sa = state;
  assert(state);
  alphabet = MRAEncGTAlphaNew(sa->alpha);
  MRAEncAddSymbolToRange(alphabet, SEPARATOR, 1);
  return alphabet;
}
