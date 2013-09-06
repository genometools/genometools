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

#ifndef EIS_ENCIDXSEQ_PRIV_H
#define EIS_ENCIDXSEQ_PRIV_H

#include <stdio.h>

#include "match/eis-mrangealphabet.h"
#include "match/eis-seqranges.h"
#include "match/eis-encidxseq.h"

struct encIdxSeqClass
{
  void (*delete)(EISeq *seq);
  GtUword (*rank)(EISeq *seq, Symbol sym, GtUword pos,
                 union EISHint *hint);
  GtUlongPair (*posPairRank)(EISeq *seq, Symbol tSym, GtUword posA,
                                   GtUword posB, union EISHint *hint);
  void (*rangeRank)(struct encIdxSeq *eSeqIdx, unsigned range,
                    GtUword pos, GtUword *rankCounts,
                    union EISHint *hint);
  void (*posPairRangeRank)(struct encIdxSeq *eSeqIdx, unsigned range,
                           GtUword posA, GtUword posB,
                           GtUword *rankCounts, union EISHint *hint);
  GtUword (*select)(EISeq *seq, Symbol sym, GtUword count,
                   union EISHint *hint);
  Symbol (*get)(EISeq *seq, GtUword pos, EISHint hint);
  union EISHint *(*newHint)(const EISeq *seq);
  void (*deleteHint)(EISeq *seq, EISHint hint);
  const MRAEnc *(*getAlphabet)(const EISeq *seq);
  void (*expose)(EISeq *seq, GtUword pos, int persistent,
                 struct extBitsRetrieval *retval, union EISHint *hint);
  FILE *(*seekToHeader)(const EISeq *seq, uint16_t headerID,
                        uint32_t *lenRet);
  int (*printPosDiags)(const EISeq *seq, GtUword pos, FILE *fp,
                       EISHint hint);
  int (*printExtPosDiags)(const EISeq *seq, GtUword pos, FILE *fp,
                          EISHint hint);
};

struct encIdxSeq
{
  const struct encIdxSeqClass *classInfo;
  MRAEnc *alphabet;
  GtUword seqLen;
};

struct seqCache
{
  size_t numEntries;
  GtUword *cachedPos;
  void **entriesPtr;
  void *entries;
};

struct blockEncIdxSeqHint
{
  struct seqCache sBlockCache;
  seqRangeListSearchHint rangeHint;
};

union EISHint
{
  struct blockEncIdxSeqHint bcHint;
};

unsigned
gt_blockEncIdxSeqSegmentLen(const struct blockEncParams *params);

#endif
