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

#include "sarr-def.h"
#include "splititv.h"
#include "procmatch.h"

unsigned long bwtseqfirstmatch(const void *voidbwtseq,unsigned long bound);

/** Iterator for positions defined by a lower and upper bound
 */

typedef struct Bwtseqpositioniterator Bwtseqpositioniterator;

Bwtseqpositioniterator *newBwtseqpositioniterator(const void *voidbwtseq,
                                                  unsigned long lowerbound,
                                                  unsigned long upperbound);

bool nextBwtseqpositioniterator(unsigned long *pos,
                                Bwtseqpositioniterator *bspi);

bool nextBwtseqpositionwithoutSEPiterator(unsigned long *pos,
                                          Bwtseqpositioniterator *bspi);

GtUchar bwtseqgetsymbol(unsigned long bound,const void *voidbwtseq);

void freeBwtseqpositioniterator(Bwtseqpositioniterator **bspi);

typedef struct Bwtseqcontextiterator Bwtseqcontextiterator;

Bwtseqcontextiterator *newBwtseqcontextiterator(const void *voidbwtseq,
                                                unsigned long bound);

GtUchar nextBwtseqcontextiterator(unsigned long *bound,
                                  Bwtseqcontextiterator *bsci);

void freeBwtseqcontextiterator(Bwtseqcontextiterator **bsci);

void bwtrangesplitwithoutspecial(GtArrayBoundswithchar *bwci,
                                 unsigned long *rangeOccs,
                                 const void *voidbwtseq,
                                 unsigned long lbound,
                                 unsigned long ubound);

void *loadvoidBWTSeqForSA(const GtStr *indexname,
                          const Suffixarray *suffixarray,
                          unsigned long totallength,
                          bool withpckbt,
                          GtError *err);

void deletevoidBWTSeq(void *packedindex);

unsigned long voidpackedindexuniqueforward(const void *voidbwtseq,
                                       GT_UNUSED unsigned long offset,
                                       GT_UNUSED unsigned long left,
                                       GT_UNUSED unsigned long right,
                                       GT_UNUSED unsigned long *witnessposition,
                                       const GtUchar *qstart,
                                       const GtUchar *qend);

unsigned long voidpackedindexmstatsforward(const void *voidbwtseq,
                                           GT_UNUSED unsigned long offset,
                                           GT_UNUSED unsigned long left,
                                           GT_UNUSED unsigned long right,
                                           unsigned long *witnessposition,
                                           const GtUchar *qstart,
                                           const GtUchar *qend);

bool pck_exactpatternmatching(const void *voidbwtseq,
                              const GtUchar *pattern,
                              unsigned long patternlength,
                              unsigned long totallength,
                              const GtUchar *dbsubstring,
                              Processmatch processmatch,
                              void *processmatchinfo);

unsigned long voidpackedfindfirstmatchconvert(const void *voidbwtseq,
                                       unsigned long witnessbound,
                                       unsigned long matchlength);

typedef struct
{
  unsigned long lowerbound, upperbound;
} Mbtab;

unsigned long bwtrangesplitallwithoutspecial(Mbtab *mbtab,
                                             unsigned long *rangeOccs,
                                             const void *voidbwtseq,
                                             unsigned long lbound,
                                             unsigned long ubound);

unsigned int bwtseq2maxdepth(const void *voidbwtseq);

const Mbtab **bwtseq2mbtab(const void *voidbwtseq);

#endif
