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
#include "eis-iterpos.h"

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

  bspi = ma_malloc(sizeof (Bwtseqpositioniterator));
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
  /* how do I determine that we are at position totallength? */
  cc = BWTSeqGetSym(bsci->bwtSeq, bsci->bound);
  bsci->bound = BWTSeqLFMap(bsci->bwtSeq, bsci->bound, &bsci->extBits);
  return cc;
}

void freeBwtseqcontextiterator(Bwtseqcontextiterator **bsci)
{
  destructExtBitsRetrieval(&(*bsci)->extBits);
  ma_free(*bsci);
  *bsci = NULL;
}

Uchar bwtseqintervalextendlcp(const void *voidBwtSeq,
                              const Lcpinterval *itv,
                              Uchar alphasize)
{
  Uchar ccl, ccr;

  ccl = BWTSeqGetSym((const BWTSeq *) voidBwtSeq,itv->left);
  ccr = BWTSeqGetSym((const BWTSeq *) voidBwtSeq,itv->right);
  if (ccl != ccr || ISSPECIAL(ccl))
  {
    return alphasize;
  }
  assert(ccl < alphasize);
  return ccl;
}
