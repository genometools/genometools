/*
  Copyright (c) 2008 Thomas Jahns <Thomas.Jahns@gmx.net>

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

#ifndef EIS_BWTSEQ_CONTEXT_SIOP_H
#define EIS_BWTSEQ_CONTEXT_SIOP_H

#include <stdio.h>
#ifndef S_SPLINT_S
#include <unistd.h>
#endif

#include "match/eis-bitpackseqpos.h"
#include "match/eis-bwtseq.h"
#include "match/eis-bwtseq-context.h"
#include "match/eis-bwtseq-context-priv.h"
#include "match/eis-encidxseq.h"

static inline unsigned long
numMapEntries(unsigned long seqLen, unsigned short mapIntervalLog2)
{
  return (seqLen + (1<<mapIntervalLog2) - 1)>> mapIntervalLog2;
}

static inline struct SeqMark
BWTSeqCRNextMark(const BWTSeqContextRetriever *bwtSeqCR, unsigned long pos)
{
  unsigned long seqLen, mapMask;
  struct SeqMark nextMark;
  gt_assert(bwtSeqCR);
  seqLen = BWTSeqLength(bwtSeqCR->bwtSeq);
  gt_assert(pos < seqLen);
  mapMask = bwtSeqCR->mapMask;
  /* rounds up to next multiple of mapInterval */
  nextMark.textPos = (pos + mapMask) & ~mapMask;
  if (nextMark.textPos < seqLen)
  {
    nextMark.bwtPos
      = gt_bsGetUlong(bwtSeqCR->revMap,
                    (nextMark.textPos >> bwtSeqCR->mapIntervalLog2)
                    * bwtSeqCR->bitsPerUlong, bwtSeqCR->bitsPerUlong);
  }
  else
  {
    nextMark.textPos = seqLen - 1;
    nextMark.bwtPos = BWTSeqTerminatorPos(bwtSeqCR->bwtSeq);
  }
  return nextMark;
}

#endif
