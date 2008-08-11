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
#include "divmodmul.h"
#include "splititv.h"

Seqpos bwtseqfirstmatch(const void *voidbwtseq,Seqpos bound)
{
  struct extBitsRetrieval extBits;
  Seqpos pos;

  initExtBitsRetrieval(&extBits);
  pos = BWTSeqLocateMatch((const BWTSeq *) voidbwtseq,bound,&extBits);
  destructExtBitsRetrieval(&extBits);
  return pos;
}

struct Bwtseqpositioniterator
{
  struct extBitsRetrieval extBits;
  const BWTSeq *bwtseq;
  Seqpos currentbound, upperbound;
};

Bwtseqpositioniterator *newBwtseqpositioniterator(const void *voidbwtseq,
                                                  Seqpos lowerbound,
                                                  Seqpos upperbound)
{
  Bwtseqpositioniterator *bspi;

  bspi = ma_malloc(sizeof (*bspi));
  initExtBitsRetrieval(&bspi->extBits);
  bspi->bwtseq = (const BWTSeq *) voidbwtseq;
  bspi->currentbound = lowerbound;
  bspi->upperbound = upperbound;
  return bspi;
}

bool nextBwtseqpositioniterator(Seqpos *pos,Bwtseqpositioniterator *bspi)
{
  if (bspi->currentbound < bspi->upperbound)
  {
    *pos = BWTSeqLocateMatch(bspi->bwtseq,bspi->currentbound,&bspi->extBits);
    bspi->currentbound++;
    return true;
  }
  return false;
}

bool nextBwtseqpositionwithoutSEPiterator(Seqpos *pos,
                                          Bwtseqpositioniterator *bspi)
{
  while (bspi->currentbound < bspi->upperbound)
  {
    Uchar cc;

    if (bspi->currentbound != BWTSeqTerminatorPos(bspi->bwtseq))
    {
      cc = BWTSeqGetSym(bspi->bwtseq, bspi->currentbound);
    } else
    {
      cc = SEPARATOR;
    }
    if (cc != SEPARATOR)
    {
      *pos = BWTSeqLocateMatch(bspi->bwtseq,bspi->currentbound,&bspi->extBits);
      bspi->currentbound++;
      return true;
    }
    bspi->currentbound++;
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
  const BWTSeq *bwtseq;
  Seqpos bound;
};

Bwtseqcontextiterator *newBwtseqcontextiterator(const void *voidbwtseq,
                                                Seqpos bound)
{
  Bwtseqcontextiterator *bsci;

  bsci = ma_malloc(sizeof (*bsci));
  initExtBitsRetrieval(&bsci->extBits);
  bsci->bwtseq = (const BWTSeq *) voidbwtseq;
  bsci->bound = bound;
  return bsci;
}

Uchar bwtseqgetsymbol(Seqpos bound,const void *voidbwtseq)
{
  if (bound != BWTSeqTerminatorPos(voidbwtseq))
  {
    return BWTSeqGetSym(voidbwtseq, bound);
  }
  return SEPARATOR;
}

Uchar nextBwtseqcontextiterator(Seqpos *bound,Bwtseqcontextiterator *bsci)
{
  Uchar cc;

  if (bsci->bound != BWTSeqTerminatorPos(bsci->bwtseq))
  {
    cc = BWTSeqGetSym(bsci->bwtseq, bsci->bound);
  } else
  {
    cc = SEPARATOR;
  }
  *bound = bsci->bound = BWTSeqLFMap(bsci->bwtseq, bsci->bound, &bsci->extBits);
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
                                 const void *voidBwtSeq,
                                 Seqpos lbound,
                                 Seqpos ubound)
{
  unsigned long idx;
  const BWTSeq *bwtseq = (const BWTSeq *) voidBwtSeq;
  AlphabetRangeSize rangesize
    = MRAEncGetRangeSize(EISGetAlphabet(bwtseq->seqIdx),0);

  bwci->nextfreeBoundswithchar = 0;
  BWTSeqPosPairRangeOcc(bwtseq, 0, lbound, ubound,rangeOccs);
  for (idx = 0; idx < rangesize; idx++)
  {
    if (rangeOccs[idx] < rangeOccs[rangesize+idx])
    {
      bwci->spaceBoundswithchar[bwci->nextfreeBoundswithchar].inchar = idx;
      bwci->spaceBoundswithchar[bwci->nextfreeBoundswithchar].lbound
        = bwtseq->count[idx] + rangeOccs[idx];
      bwci->spaceBoundswithchar[bwci->nextfreeBoundswithchar++].rbound
        = bwtseq->count[idx] + rangeOccs[rangesize+idx];
    }
  }
}

/*
void bwtrangewithspecial(UNUSED ArrayBoundswithchar *bwci,
                         Seqpos *rangeOccs,
                         UNUSED unsigned long alphasize,
                         const void *voidBwtSeq,
                         const Lcpinterval *parent)
{
  const BWTSeq *bwtseq = (const BWTSeq *) voidBwtSeq;

  AlphabetRangeSize rangesize
    = MRAEncGetRangeSize(EISGetAlphabet(bwtseq->seqIdx),1);
  assert(rangesize < (AlphabetRangeSize) 4);
  BWTSeqPosPairRangeOcc(bwtseq, 1, parent->left, parent->right,rangeOccs);
    inchar = WILDCARD
    bwtcode = MRAEncMapSymbol(EISGetAlphabet(bwtseq->seqIdx),WILDCARD);
    idx = bwtcode -  MRAEncGetRangeBase(EISGetAlphabet(bwtseq->seqIdx),1);
    if (rangeOccs[idx] < rangeOccs[rangesize+idx])
    {
      pos = BWTSeqLocateMatch((const BWTSeq *) voidbwtseq,bound,&extBits);
    }
  }
}
*/

void deletevoidBWTSeq(void *packedindex)
{
  deleteBWTSeq((BWTSeq *) packedindex);
}

unsigned long voidpackedindexuniqueforward(const void *voidbwtseq,
                                           UNUSED unsigned long offset,
                                           UNUSED Seqpos left,
                                           UNUSED Seqpos right,
                                           UNUSED Seqpos *witnessposition,
                                           const Uchar *qstart,
                                           const Uchar *qend)
{
  return packedindexuniqueforward((const BWTSeq *) voidbwtseq,
                                  qstart,
                                  qend);
}

Seqpos voidpackedfindfirstmatchconvert(const void *voidbwtseq,
                                       Seqpos witnessbound,
                                       unsigned long matchlength)
{
  const BWTSeq *bwtseq = (const BWTSeq *) voidbwtseq;
  Seqpos startpos;

  startpos = bwtseqfirstmatch(voidbwtseq,witnessbound);
  assert((bwtseq->seqIdx->seqLen-1) >= (startpos + matchlength));
  return (bwtseq->seqIdx->seqLen - 1) - (startpos + matchlength);
}

unsigned long voidpackedindexmstatsforward(const void *voidbwtseq,
                                           UNUSED unsigned long offset,
                                           UNUSED Seqpos left,
                                           UNUSED Seqpos right,
                                           Seqpos *witnessposition,
                                           const Uchar *qstart,
                                           const Uchar *qend)
{
  const BWTSeq *bwtseq = (const BWTSeq *) voidbwtseq;
  unsigned long matchlength;

  matchlength = packedindexmstatsforward(bwtseq,witnessposition,qstart,qend);
  if (matchlength > 0 && witnessposition != NULL)
  {
    *witnessposition = voidpackedfindfirstmatchconvert(voidbwtseq,
                                                       *witnessposition,
                                                       matchlength);
  }
  return matchlength;
}

void pck_exactpatternmatching(const void *voidbwtseq,
                              const Uchar *pattern,
                              unsigned long patternlength,
                              Seqpos totallength,
                              void (*processmatch)(void *,bool,Seqpos,
                                                   Seqpos,Seqpos,unsigned long),
                              void *processmatchinfo)
{
  BWTSeqExactMatchesIterator *bsemi;
  Seqpos dbstartpos;

  bsemi = newEMIterator((const BWTSeq *) voidbwtseq,
                        pattern,(size_t) patternlength, true);
  assert(bsemi != NULL);
  while (EMIGetNextMatch(bsemi,&dbstartpos,(const BWTSeq *) voidbwtseq))
  {
    processmatch(processmatchinfo,false,totallength,
                 dbstartpos + patternlength,(Seqpos) patternlength,
                 patternlength);
  }
  if (bsemi != NULL)
  {
    deleteEMIterator(bsemi);
    bsemi = NULL;
  }
}

typedef struct
{
  Seqpos lowerbound,
         upperbound;
  unsigned long depth;
} Boundsatdepth;

DECLAREARRAYSTRUCT(Boundsatdepth);

void pck_precomputebounds(Matchbound *boundsarray,
                          unsigned long numofbounds,
                          const void *voidbwtseq,
                          unsigned int alphasize,
                          Seqpos totallength,
                          unsigned long maxdepth)
{
  const BWTSeq *bwtseq = (const BWTSeq *) voidbwtseq;
  ArrayBoundsatdepth stack;
  Boundsatdepth *stackptr, parent, child;
  unsigned long idx;
  Seqpos *rangeOccs;
  Matchbound *bptr;
  AlphabetRangeSize rangesize;

  rangesize = MRAEncGetRangeSize(EISGetAlphabet(bwtseq->seqIdx),0);
  INITARRAY(&stack,Boundsatdepth);
  GETNEXTFREEINARRAY(stackptr,&stack,Boundsatdepth,128);
  stackptr->lowerbound = 0;
  stackptr->upperbound = totallength+1;
  stackptr->depth = 0;
  rangeOccs = ma_malloc(sizeof(*rangeOccs) * MULT2(alphasize));
  bptr = boundsarray;
  while (stack.nextfreeBoundsatdepth > 0)
  {
    parent = stack.spaceBoundsatdepth[--stack.nextfreeBoundsatdepth];
    BWTSeqPosPairRangeOcc(bwtseq, 0,
                          parent.lowerbound,parent.upperbound,rangeOccs);
    for (idx = 0; idx < rangesize; idx++)
    {
      if (rangeOccs[idx] < rangeOccs[rangesize+idx])
      {
        child.lowerbound = bwtseq->count[idx] + rangeOccs[idx];
        child.upperbound = bwtseq->count[idx] + rangeOccs[rangesize+idx];
      } else
      {
        child.lowerbound = child.upperbound = 0;
      }
      child.depth = parent.depth + 1;
      if (child.depth == maxdepth)
      {
        assert(bptr < boundsarray + numofbounds);
        bptr->lowerbound = child.lowerbound;
        bptr->upperbound = child.upperbound;
        bptr++;
      } else
      {
        GETNEXTFREEINARRAY(stackptr,&stack,Boundsatdepth,128);
        *stackptr = child;
      }
    }
  }
  FREEARRAY(&stack,Boundsatdepth);
  ma_free(rangeOccs);
}
