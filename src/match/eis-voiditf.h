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

#include "core/unused_api.h"
#include "seqpos-def.h"
#include "sarr-def.h"
#include "splititv.h"
#include "procmatch.h"

Seqpos bwtseqfirstmatch(const void *voidbwtseq,Seqpos bound);

/** Iterator for positions defined by a lower and upper bound
 */

typedef struct Bwtseqpositioniterator Bwtseqpositioniterator;

Bwtseqpositioniterator *newBwtseqpositioniterator(const void *voidbwtseq,
                                                  Seqpos lowerbound,
                                                  Seqpos upperbound);

bool nextBwtseqpositioniterator(Seqpos *pos,Bwtseqpositioniterator *bspi);

bool nextBwtseqpositionwithoutSEPiterator(Seqpos *pos,
                                          Bwtseqpositioniterator *bspi);

GtUchar bwtseqgetsymbol(Seqpos bound,const void *voidbwtseq);

void freeBwtseqpositioniterator(Bwtseqpositioniterator **bspi);

typedef struct Bwtseqcontextiterator Bwtseqcontextiterator;

Bwtseqcontextiterator *newBwtseqcontextiterator(const void *voidbwtseq,
                                                Seqpos bound);

GtUchar nextBwtseqcontextiterator(Seqpos *bound,Bwtseqcontextiterator *bsci);

void freeBwtseqcontextiterator(Bwtseqcontextiterator **bsci);

void bwtrangesplitwithoutspecial(ArrayBoundswithchar *bwci,
                                 Seqpos *rangeOccs,
                                 const void *voidbwtseq,
                                 Seqpos lbound,
                                 Seqpos ubound);

void *loadvoidBWTSeqForSA(const GtStr *indexname,
                          const Suffixarray *suffixarray,
                          Seqpos totallength,
                          bool withpckbt,
                          GtError *err);

void deletevoidBWTSeq(void *packedindex);

unsigned long voidpackedindexuniqueforward(const void *voidbwtseq,
                                           GT_UNUSED unsigned long offset,
                                           GT_UNUSED Seqpos left,
                                           GT_UNUSED Seqpos right,
                                           GT_UNUSED Seqpos *witnessposition,
                                           const GtUchar *qstart,
                                           const GtUchar *qend);

unsigned long voidpackedindexmstatsforward(const void *voidbwtseq,
                                           GT_UNUSED unsigned long offset,
                                           GT_UNUSED Seqpos left,
                                           GT_UNUSED Seqpos right,
                                           Seqpos *witnessposition,
                                           const GtUchar *qstart,
                                           const GtUchar *qend);

bool pck_exactpatternmatching(const void *voidbwtseq,
                              const GtUchar *pattern,
                              unsigned long patternlength,
                              Seqpos totallength,
                              const GtUchar *dbsubstring,
                              Processmatch processmatch,
                              void *processmatchinfo);

Seqpos voidpackedfindfirstmatchconvert(const void *voidbwtseq,
                                       Seqpos witnessbound,
                                       unsigned long matchlength);

typedef struct
{
  Seqpos lowerbound, upperbound;
} Mbtab;

unsigned long bwtrangesplitallwithoutspecial(Mbtab *mbtab,
                                             Seqpos *rangeOccs,
                                             const void *voidbwtseq,
                                             Seqpos lbound,
                                             Seqpos ubound);

unsigned int bwtseq2maxdepth(const void *voidbwtseq);

const Mbtab **bwtseq2mbtab(const void *voidbwtseq);

#endif
