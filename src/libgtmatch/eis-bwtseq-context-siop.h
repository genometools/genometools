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
#include <unistd.h>

#include "libgtmatch/eis-bitpackseqpos.h"
#include "libgtmatch/eis-bwtseq.h"
#include "libgtmatch/eis-bwtseq-context.h"
#include "libgtmatch/eis-bwtseq-context-priv.h"
#include "libgtmatch/eis-encidxseq.h"

static inline Seqpos
numMapEntries(Seqpos seqLen, unsigned short mapIntervalLog2)
{
  return (seqLen + (1<<mapIntervalLog2) - 1)>> mapIntervalLog2;
}

static inline struct SeqMark
BWTSeqCRNextMark(const BWTSeqContextRetriever *bwtSeqCR, Seqpos pos)
{
  Seqpos seqLen, mapMask;
  struct SeqMark nextMark;
  assert(bwtSeqCR);
  seqLen = BWTSeqLength(bwtSeqCR->bwtSeq);
  assert(pos < seqLen);
  mapMask = bwtSeqCR->mapMask;
  /* rounds up to next multiple of mapInterval */
  nextMark.textPos = (pos + mapMask) & ~mapMask;
  if (nextMark.textPos < seqLen)
  {
    nextMark.bwtPos
      = bsGetSeqpos(bwtSeqCR->revMap,
                    (nextMark.textPos >> bwtSeqCR->mapIntervalLog2)
                    * bwtSeqCR->bitsPerSeqpos, bwtSeqCR->bitsPerSeqpos);
  }
  else
  {
    nextMark.textPos = seqLen - 1;
    nextMark.bwtPos = BWTSeqTerminatorPos(bwtSeqCR->bwtSeq);
  }
  return nextMark;
}

#endif
