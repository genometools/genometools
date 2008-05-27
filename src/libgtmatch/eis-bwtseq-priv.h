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
#ifndef EIS_BWTSEQ_PRIV_H
#define EIS_BWTSEQ_PRIV_H

#include "libgtcore/chardef.h"
#include "libgtmatch/eis-bwtseq.h"
#include "libgtmatch/eis-bwtseq-extinfo.h"
#include "libgtmatch/eis-encidxseq.h"

enum {
  bwtTerminatorSym = SEPARATOR - 3,
};

struct BWTSeq
{
  struct encIdxSeq *seqIdx;
  MRAEnc *alphabet;
  size_t alphabetSize;
  EISHint hint;
  unsigned locateSampleInterval; /**< no sampling if 0 */
  Symbol bwtTerminatorFallback;  /**< the terminator symbol has been
                                  * flattened into this symbol for
                                  * storage reasons */
  AlphabetRangeID bwtTerminatorFallbackRange;
  Seqpos rot0Pos;
  Seqpos *count;
  int featureToggles;
  unsigned bitsPerOrigRank;
  enum rangeSortMode *rangeSort;
};

struct BWTSeqExactMatchesIterator
{
  struct matchBound bounds;
  Seqpos nextMatchBWTPos;
  struct extBitsRetrieval extBits;
};

extern Seqpos
BWTSeqLocateMatch(const BWTSeq *bwtSeq, Seqpos pos,
                  struct extBitsRetrieval *extBits);

extern BWTSeq *
newBWTSeq(EISeq *seqIdx, MRAEnc *alphabet,
          const enum rangeSortMode *defaultRangeSort);

#endif
