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

#ifndef EIS_SEQRANGES_H
#define EIS_SEQRANGES_H

#include <stdlib.h>

#include "libgtcore/env.h"
#include "libgtmatch/seqpos-def.h"
#include "libgtmatch/eis-mrangealphabet.h"

typedef uint16_t regionLength;

struct seqRange
{
  Seqpos startPos;
  regionLength len;
  Symbol sym;
};

static inline int
posIsInSeqRange(struct seqRange *p, Seqpos pos)
{
  return (pos >= p->startPos && pos < p->startPos + p->len);
}

enum {
  MAX_SEQRANGE_LEN = UINT16_MAX,
};

enum SRLFeatures {
  SRL_NO_FEATURES = 0,
  SRL_PARTIAL_SYMBOL_SUMS = 1 << 0,
};

struct seqRangeList
{
  size_t numRangesStorable, numRanges;
  struct seqRange *ranges;
  Seqpos *partialSymSums;
  const MRAEnc *alphabet;
};

typedef size_t seqRangeListSearchHint;

extern struct seqRangeList *
newSeqRangeList(size_t rangesStartNum, const MRAEnc *alphabet,
                enum SRLFeatures features, Env *env);

extern void
SRLCompact(struct seqRangeList *rangeList, Env *env);

extern void
deleteSeqRangeList(struct seqRangeList *rangeList, Env *env);

extern void
SRLAppendNewRange(struct seqRangeList *rangeList, Seqpos pos, Seqpos len,
                  Symbol sym, Env *env);

extern void
SRLAddPosition(struct seqRangeList *rangeList, Seqpos pos,
               Symbol sym, Env *env);

extern void
SRLInitListSearchHint(struct seqRangeList *rangeList,
                      seqRangeListSearchHint *hint);

struct seqRange *
SRLFindPositionNext(struct seqRangeList *rangeList, Seqpos pos,
                    seqRangeListSearchHint *hint);

extern int
SRLOverlapsPosition(struct seqRangeList *rangeList, Seqpos pos,
                    seqRangeListSearchHint *hint, Symbol *symAtPos);

extern void
SRLSymbolsInSeqRegion(struct seqRangeList *rangeList, Seqpos start,
                      Seqpos end, Seqpos *occStore,
                      seqRangeListSearchHint *hint);

extern Seqpos
SRLSymbolCountInSeqRegion(struct seqRangeList *rangeList, Seqpos start,
                          Seqpos end, Symbol sym, seqRangeListSearchHint *hint);

extern Seqpos
SRLAllSymbolsCountInSeqRegion(struct seqRangeList *rangeList, Seqpos start,
                              Seqpos end, seqRangeListSearchHint *hint);

extern void
SRLapplyRangesToSubString(struct seqRangeList *rangeList, MRAEnc *alphabet,
                          Symbol *subString, Seqpos start, Seqpos len,
                          Seqpos subStringOffset, seqRangeListSearchHint *hint);

extern int
SRLSaveToStream(struct seqRangeList *rangeList, FILE *fp);

extern struct seqRangeList *
SRLReadFromStream(FILE *fp, const MRAEnc *alphabet,
                  enum SRLFeatures features, Env *env);

#endif
