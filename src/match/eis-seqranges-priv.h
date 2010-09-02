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

#ifndef EIS_SEQRANGES_PRIV_H
#define EIS_SEQRANGES_PRIV_H

#include "match/eis-bitpackseqpos.h"
#include "match/eis-seqranges.h"

typedef uint16_t regionLength;

struct seqRange
{
  unsigned long startPos;
  BitElem symLenStr[sizeof (unsigned long) / sizeof (BitElem)];
  /**< holds the bits required for the symbol stored and the length
   * of the range */
};

enum {
  symLenStrBits = sizeof (unsigned long) / sizeof (BitElem) * bitElemBits,
};

struct seqRangeList
{
  size_t numRangesStorable, numRanges;
  struct seqRange *ranges;
  unsigned long *partialSymSums;
  const MRAEnc *alphabet;
  unsigned symBits;
  unsigned long maxRangeLen;
};

static inline unsigned long
seqRangeLen(const struct seqRange *p, unsigned symBits)
{
  return gt_bsGetUlong(p->symLenStr, symBits, symLenStrBits - symBits);
}

static inline Symbol
seqRangeSym(const struct seqRange *p, unsigned symBits)
{
  return gt_bsGetSymbol(p->symLenStr, 0, symBits);
}

static inline void
seqRangeSetLen(struct seqRange *p, unsigned long len, unsigned symBits)
{
  return gt_bsStoreUlong(p->symLenStr, symBits, symLenStrBits - symBits, len);
}

static inline void
seqRangeSetSym(struct seqRange *p, Symbol sym, unsigned symBits)
{
  return gt_bsStoreSymbol(p->symLenStr, 0, symBits, sym);
}

static inline int
posIsInSeqRange(struct seqRange *p, unsigned long pos, unsigned symBits)
{
  return (pos >= p->startPos && pos < p->startPos + seqRangeLen(p, symBits));
}

#endif
