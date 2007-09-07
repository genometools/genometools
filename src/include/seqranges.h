
/*
** Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
**  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**  
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**  
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
**  
*/
#ifndef SEQRANGES_H_INCLUDED
#define SEQRANGES_H_INCLUDED

#include <stdlib.h>

#include <libgtcore/env.h>
#include <libgtmatch/seqpos-def.h>

#include "mrangealphabet.h"

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

extern int
SRLSaveToStream(struct seqRangeList *rangeList, FILE *fp);

extern struct seqRangeList *
SRLReadFromStream(FILE *fp, const MRAEnc *alphabet,
                  enum SRLFeatures features, Env *env);

#endif /* SEQRANGES_H_INCLUDED */
