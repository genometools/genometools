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

#include "eis-bwtseq.h"
#include "eis-bwtseq-construct.h"
#include "eis-voiditf.h"
#include "splititv.h"
#include "stamp.h"

Seqpos bwtseqfirstmatch(const void *voidbwtSeq,Seqpos bound)
{
  struct extBitsRetrieval extBits;
  Seqpos pos;

  initExtBitsRetrieval(&extBits);
  pos = BWTSeqLocateMatch((const BWTSeq *) voidbwtSeq,bound,&extBits);
  destructExtBitsRetrieval(&extBits);
  return pos;
}

struct Bwtseqpositioniterator
{
  struct extBitsRetrieval extBits;
  const BWTSeq *bwtSeq;
  Seqpos currentbound, upperbound;
};

Bwtseqpositioniterator *newBwtseqpositioniterator(const void *voidbwtSeq,
                                                  Seqpos lowerbound,
                                                  Seqpos upperbound)
{
  Bwtseqpositioniterator *bspi;

  bspi = ma_malloc(sizeof (*bspi));
  initExtBitsRetrieval(&bspi->extBits);
  bspi->bwtSeq = (const BWTSeq *) voidbwtSeq;
  bspi->currentbound = lowerbound;
  bspi->upperbound = upperbound;
  return bspi;
}

bool nextBwtseqpositioniterator(Seqpos *pos,Bwtseqpositioniterator *bspi)
{
  if (bspi->currentbound < bspi->upperbound)
  {
    *pos = BWTSeqLocateMatch(bspi->bwtSeq,bspi->currentbound,&bspi->extBits);
    bspi->currentbound++;
    return true;
  }
  return false;
}

void freeBwtseqpositioniterator(Bwtseqpositioniterator **bspi)
{
  destructExtBitsRetrieval(&(*bspi)->extBits);
  ma_free(*bspi);
  *bspi = NULL;
}

struct Bwtseqcontextiterator
{
  struct extBitsRetrieval extBits;
  const BWTSeq *bwtSeq;
  Seqpos bound;
};

Bwtseqcontextiterator *newBwtseqcontextiterator(const void *voidbwtSeq,
                                                Seqpos bound)
{
  Bwtseqcontextiterator *bsci;

  bsci = ma_malloc(sizeof (*bsci));
  initExtBitsRetrieval(&bsci->extBits);
  bsci->bwtSeq = (const BWTSeq *) voidbwtSeq;
  bsci->bound = bound;
  return bsci;
}

Uchar nextBwtseqcontextiterator(Bwtseqcontextiterator *bsci)
{
  Uchar cc;

  if (bsci->bound != BWTSeqTerminatorPos(bsci->bwtSeq))
  {
    cc = BWTSeqGetSym(bsci->bwtSeq, bsci->bound);
  } else
  {
    cc = SEPARATOR;
  }
  bsci->bound = BWTSeqLFMap(bsci->bwtSeq, bsci->bound, &bsci->extBits);
  return cc;
}

void freeBwtseqcontextiterator(Bwtseqcontextiterator **bsci)
{
  destructExtBitsRetrieval(&(*bsci)->extBits);
  ma_free(*bsci);
  *bsci = NULL;
}

void *loadvoidBWTSeqForSA(const Str *indexname,
                          const Suffixarray *suffixarray,
                          Seqpos totallength,
                          Error *err)
{
  return loadBWTSeqForSA(indexname,
                         BWT_ON_BLOCK_ENC,
                         BWTDEFOPT_MULTI_QUERY,
                         suffixarray,
                         totallength+1, err);
}

void bwtrangesplitwithoutspecial(ArrayBoundswithchar *bwci,
                                 Seqpos *rangeOccs,
                                 unsigned long alphasize,
                                 const void *voidBwtSeq,
                                 const Lcpinterval *parent)
{
  unsigned long idx;
  const BWTSeq *bwtseq = (const BWTSeq *) voidBwtSeq;

  bwci->nextfreeBoundswithchar = 0;
  BWTSeqPosPairRangeOcc(bwtseq, 0, parent->left, parent->right,rangeOccs);
  for (idx = 0; idx < alphasize; idx++)
  {
    if (rangeOccs[idx] < rangeOccs[alphasize+idx])
    {
      bwci->spaceBoundswithchar[bwci->nextfreeBoundswithchar].inchar = idx;
      bwci->spaceBoundswithchar[bwci->nextfreeBoundswithchar].lbound
        = bwtseq->count[idx] + rangeOccs[idx];
      bwci->spaceBoundswithchar[bwci->nextfreeBoundswithchar++].rbound
        = bwtseq->count[idx] + rangeOccs[alphasize+idx];
    }
  }
}

void deletevoidBWTSeq(void *packedindex)
{
  deleteBWTSeq((BWTSeq *) packedindex);
}

unsigned long voidpackedindexuniqueforward(const void *genericindex,
                                           UNUSED unsigned long offset,
                                           UNUSED Seqpos left,
                                           UNUSED Seqpos right,
                                           UNUSED Seqpos *witnessposition,
                                           const Uchar *qstart,
                                           const Uchar *qend)
{
  return packedindexuniqueforward((const BWTSeq *) genericindex,
                                  qstart,
                                  qend);
}

unsigned long voidpackedindexmstatsforward(const void *genericindex,
                                           UNUSED unsigned long offset,
                                           UNUSED Seqpos left,
                                           UNUSED Seqpos right,
                                           Seqpos *witnessposition,
                                           const Uchar *qstart,
                                           const Uchar *qend)
{
  return packedindexmstatsforward((const BWTSeq *) genericindex,
                                  witnessposition,
                                  qstart,
                                  qend);
}

void pck_exactpatternmatching(const void *genericindex,
                              const Uchar *pattern,
                              unsigned long patternlength,
                              Seqpos totallength,
                              void (*processmatch)(void *,bool,Seqpos,
                                                   Seqpos,Seqpos),
                              void *processmatchinfo)
{
  BWTSeqExactMatchesIterator *bsemi;
  Seqpos dbstartpos;

  bsemi = newEMIterator((const BWTSeq *) genericindex,
                        pattern,(size_t) patternlength, true);
  assert(bsemi != NULL);
  while (EMIGetNextMatch(bsemi,&dbstartpos,(const BWTSeq *) genericindex))
  {
    processmatch(processmatchinfo,false,totallength,
                 dbstartpos + patternlength,(Seqpos) patternlength);
  }
  if (bsemi != NULL)
  {
    deleteEMIterator(bsemi);
    bsemi = NULL;
  }
}
