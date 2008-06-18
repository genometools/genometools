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
#ifndef EIS_BWTSEQ_CONTEXT_PRIV_H
#define EIS_BWTSEQ_CONTEXT_PRIV_H

#include "libgtmatch/seqpos-def.h"
#include "libgtmatch/eis-bwtseq-context.h"

struct BWTSeqContextRetriever
{
  Seqpos mapInterval, mapMask;
  const BWTSeq *bwtSeq;
  BitString revMapMMap,         /**< holds the page aligned result
                                 * from mmap, later used in munmap  */
    revMap;                     /**< maps positions in the original
                                 * sequence to positions in the BWT */
  uint16_t mapIntervalLog2;
  uint16_t bitsPerSeqpos;
};

static inline Seqpos
numMapEntries(Seqpos seqLen, unsigned short mapIntervalLog2);

#endif
