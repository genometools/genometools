/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef EIS_VOIDITF_H
#define EIS_VOIDITF_H

#include "seqpos-def.h"
#include "sarr-def.h"

Seqpos bwtseqfirstmatch(const void *voidbwtSeq,Seqpos bound);

/** Iterator for positions defined by a lower and upper bound
 */

typedef struct Bwtseqpositioniterator Bwtseqpositioniterator;

Bwtseqpositioniterator *newBwtseqpositioniterator(const void *bwtSeq,
                                                  Seqpos lowerbound,
                                                  Seqpos upperbound);

bool nextBwtseqpositioniterator(Seqpos *pos,Bwtseqpositioniterator *bspi);

void freeBwtseqpositioniterator(Bwtseqpositioniterator **bspi);

typedef struct Bwtseqcontextiterator Bwtseqcontextiterator;

Bwtseqcontextiterator *newBwtseqcontextiterator(const void *voidbwtSeq,
                                                Seqpos bound);

Uchar nextBwtseqcontextiterator(Bwtseqcontextiterator *bsci);

void freeBwtseqcontextiterator(Bwtseqcontextiterator **bsci);

Uchar bwtseqintervalextendlcp(const void *voidBwtSeq,
                              const Lcpinterval *itv,
                              Uchar alphasize);

void bwtrangesplitwithoutspecial(Seqpos *rangeOccs,
                                 const void *voidBwtSeq,
                                 const Lcpinterval *parent);

void *loadvoidBWTSeqForSA(const Str *indexname,
                          const Suffixarray *suffixarray,
                          Seqpos totallength,
                          Error *err);

void deletevoidBWTSeq(void *packedindex);

unsigned long voidpackedindexuniqueforward(const void *genericindex,
                                           UNUSED unsigned long offset,
                                           UNUSED Seqpos left,
                                           UNUSED Seqpos right,
                                           UNUSED Seqpos *witnessposition,
                                           const Uchar *qstart,
                                           const Uchar *qend);

unsigned long voidpackedindexmstatsforward(const void *genericindex,
                                           UNUSED unsigned long offset,
                                           UNUSED Seqpos left,
                                           UNUSED Seqpos right,
                                           Seqpos *witnessposition,
                                           const Uchar *qstart,
                                           const Uchar *qend);

#endif
