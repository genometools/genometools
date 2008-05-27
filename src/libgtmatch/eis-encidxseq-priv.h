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

#include "libgtmatch/eis-mrangealphabet.h"
#include "libgtmatch/eis-seqranges.h"
#include "libgtmatch/eis-encidxseq.h"

struct encIdxSeqClass
{
  void (*delete)(EISeq *seq);
  Seqpos (*rank)(EISeq *seq, Symbol sym, Seqpos pos,
                 union EISHint *hint);
  struct SeqposPair (*posPairRank)(EISeq *seq, Symbol tSym, Seqpos posA,
                                   Seqpos posB, union EISHint *hint);
  void (*rangeRank)(struct encIdxSeq *eSeqIdx, unsigned range, Seqpos pos,
                    Seqpos *rankCounts, union EISHint *hint);
  void (*posPairRangeRank)(struct encIdxSeq *eSeqIdx, unsigned range,
                           Seqpos posA, Seqpos posB, Seqpos *rankCounts,
                           union EISHint *hint);
  Seqpos (*select)(EISeq *seq, Symbol sym, Seqpos count,
                   union EISHint *hint);
  Symbol (*get)(EISeq *seq, Seqpos pos, EISHint hint);
  union EISHint *(*newHint)(const EISeq *seq);
  void (*deleteHint)(EISeq *seq, EISHint hint);
  const MRAEnc *(*getAlphabet)(const EISeq *seq);
  void (*expose)(EISeq *seq, Seqpos pos, int persistent,
                 struct extBitsRetrieval *retval, union EISHint *hint);
  FILE *(*seekToHeader)(const EISeq *seq, uint16_t headerID,
                        uint32_t *lenRet);
  int (*printPosDiags)(const EISeq *seq, Seqpos pos, FILE *fp, EISHint hint);
  int (*printExtPosDiags)(const EISeq *seq, Seqpos pos, FILE *fp, EISHint hint);
};

struct encIdxSeq
{
  const struct encIdxSeqClass *classInfo;
  MRAEnc *alphabet;
  Seqpos seqLen;
};

struct seqCache
{
  size_t numEntries;
  Seqpos *cachedPos;
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

extern unsigned
blockEncIdxSeqSegmentLen(const struct blockEncParams *params);

#endif
