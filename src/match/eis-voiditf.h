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
#include "core/encseq_api.h"

#include "splititv.h"
#include "procmatch.h"

/* The following type is just used to stronlgy type the functions
   using the FMindex (rather than using void pointers) */

typedef struct FMindex FMindex;

unsigned long gt_bwtseqfirstmatch(const FMindex *voidbwtseq,
                                  unsigned long bound);

/** Iterator for positions defined by a lower and upper bound
 */

typedef struct Bwtseqpositioniterator Bwtseqpositioniterator;

Bwtseqpositioniterator *gt_newBwtseqpositioniterator(const FMindex *voidbwtseq,
                                                  unsigned long lowerbound,
                                                  unsigned long upperbound);

bool gt_nextBwtseqpositioniterator(unsigned long *pos,
                                Bwtseqpositioniterator *bspi);

bool gt_nextBwtseqpositionwithoutSEPiterator(unsigned long *pos,
                                          Bwtseqpositioniterator *bspi);

GtUchar gt_bwtseqgetsymbol(unsigned long bound,const FMindex *voidbwtseq);

void gt_freeBwtseqpositioniterator(Bwtseqpositioniterator **bspi);

typedef struct Bwtseqcontextiterator Bwtseqcontextiterator;

Bwtseqcontextiterator *gt_newBwtseqcontextiterator(const FMindex *voidbwtseq,
                                                unsigned long bound);

GtUchar gt_nextBwtseqcontextiterator(unsigned long *bound,
                                  Bwtseqcontextiterator *bsci);

void gt_freeBwtseqcontextiterator(Bwtseqcontextiterator **bsci);

void gt_bwtrangesplitwithoutspecial(GtArrayBoundswithchar *bwci,
                                 unsigned long *rangeOccs,
                                 const FMindex *voidbwtseq,
                                 unsigned long lbound,
                                 unsigned long ubound);

FMindex *gt_loadvoidBWTSeqForSA(const char *indexname,
                                const GtAlphabet *gtalphabet,
                                unsigned long totallength,
                                bool withpckbt,
                                GtError *err);

void gt_deletevoidBWTSeq(FMindex *packedindex);

unsigned long gt_voidpackedindexuniqueforward(const void *voidbwtseq,
                                              GT_UNUSED unsigned long offset,
                                              GT_UNUSED unsigned long left,
                                              GT_UNUSED unsigned long right,
                                              GT_UNUSED unsigned long
                                                        *witnessposition,
                                              const GtUchar *qstart,
                                              const GtUchar *qend);

unsigned long gt_voidpackedindexmstatsforward(const void *voidbwtseq,
                                           GT_UNUSED unsigned long offset,
                                           GT_UNUSED unsigned long left,
                                           GT_UNUSED unsigned long right,
                                           unsigned long *witnessposition,
                                           const GtUchar *qstart,
                                           const GtUchar *qend);

bool gt_pck_exactpatternmatching(const FMindex *voidbwtseq,
                              const GtUchar *pattern,
                              unsigned long patternlength,
                              unsigned long totallength,
                              const GtUchar *dbsubstring,
                              Processmatch processmatch,
                              void *processmatchinfo);

unsigned long gt_voidpackedfindfirstmatchconvert(const FMindex *voidbwtseq,
                                       unsigned long witnessbound,
                                       unsigned long matchlength);

typedef struct
{
  unsigned long lowerbound, upperbound;
} Mbtab;

unsigned long gt_bwtrangesplitallwithoutspecial(Mbtab *mbtab,
                                             unsigned long *rangeOccs,
                                             const FMindex *voidbwtseq,
                                             unsigned long lbound,
                                             unsigned long ubound);

unsigned int gt_bwtseq2maxdepth(const FMindex *voidbwtseq);

const Mbtab **gt_bwtseq2mbtab(const FMindex *voidbwtseq);

#endif
