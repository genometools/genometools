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

#include "initbasepower.pr"

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
                         UNUSED unsigned long numofchars,
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
  Codetype code;
} Boundsatdepth;

DECLAREARRAYSTRUCT(Boundsatdepth);

typedef struct
{
  Seqpos lowerbound, upperbound;
} Matchbound;

struct Pckbuckettable
{
  Matchbound **mbtab;
  unsigned int maxdepth;
  unsigned long numofvalues;
  Codetype *basepower, maxnumofvalues;
};

static Pckbuckettable *allocandinitpckbuckettable(unsigned int numofchars,
                                                  unsigned int maxdepth)
{
  Matchbound *cptr;
  unsigned int idx;
  Pckbuckettable *pckbt;

  pckbt = ma_malloc(sizeof(Pckbuckettable));
  pckbt->basepower = initbasepower(numofchars,maxdepth);
  pckbt->maxdepth = maxdepth;
  pckbt->maxnumofvalues = pckbt->numofvalues = 0;
  for (idx=0; idx <= maxdepth; idx++)
  {
    pckbt->maxnumofvalues += pckbt->basepower[idx];
  }
  pckbt->mbtab = ma_malloc(sizeof(Matchbound *) * (maxdepth+1));
  pckbt->mbtab[0] = ma_malloc(sizeof(Matchbound) * pckbt->maxnumofvalues);
  for (cptr = pckbt->mbtab[0];
       cptr < pckbt->mbtab[0] + pckbt->maxnumofvalues; cptr++)
  {
    cptr->lowerbound = cptr->upperbound = 0;
  }
  for (idx=0; idx<maxdepth; idx++)
  {
    pckbt->mbtab[idx+1] = pckbt->mbtab[idx] + pckbt->basepower[idx];
  }
  return pckbt;
}

void pckbuckettable_free(Pckbuckettable *pckbt)
{
  ma_free(pckbt->mbtab[0]);
  ma_free(pckbt->mbtab);
  ma_free(pckbt->basepower);
  ma_free(pckbt);
}

static void storeBoundsatdepth(Pckbuckettable *pckbt,const Boundsatdepth *bd)
{
  assert(bd->depth <= pckbt->maxdepth);
  assert(bd->code <= pckbt->basepower[bd->depth]);
  assert(pckbt->mbtab[bd->depth][bd->code].lowerbound == 0 &&
         pckbt->mbtab[bd->depth][bd->code].upperbound == 0);
  assert(pckbt->numofvalues < pckbt->maxnumofvalues);
  pckbt->numofvalues++;
  pckbt->mbtab[bd->depth][bd->code].lowerbound = bd->lowerbound;
  pckbt->mbtab[bd->depth][bd->code].upperbound = bd->upperbound;
}

Pckbuckettable *pckbuckettable_new(const void *voidbwtseq,
                                   unsigned int numofchars,
                                   Seqpos totallength,
                                   unsigned int maxdepth)
{
  const BWTSeq *bwtseq = (const BWTSeq *) voidbwtseq;
  ArrayBoundsatdepth stack;
  Boundsatdepth parent, child;
  unsigned long idx;
  Seqpos *rangeOccs;
  AlphabetRangeSize rangesize;
  Pckbuckettable *pckbt;

  rangesize = MRAEncGetRangeSize(EISGetAlphabet(bwtseq->seqIdx),0);
  INITARRAY(&stack,Boundsatdepth);
  child.lowerbound = 0;
  child.upperbound = totallength+1;
  child.depth = 0;
  child.code = (Codetype) 0;
  STOREINARRAY(&stack,Boundsatdepth,128,child);
  rangeOccs = ma_malloc(sizeof(*rangeOccs) * MULT2(numofchars));
  pckbt = allocandinitpckbuckettable(numofchars,maxdepth);
  while (stack.nextfreeBoundsatdepth > 0)
  {
    parent = stack.spaceBoundsatdepth[--stack.nextfreeBoundsatdepth];
    assert(parent.lowerbound < parent.upperbound);
    BWTSeqPosPairRangeOcc(bwtseq,0,parent.lowerbound,parent.upperbound,
                          rangeOccs);
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
      assert(child.depth <= (unsigned long) maxdepth);
      child.code = parent.code * numofchars + idx;
      /*
      printf("depth=%lu code=%lu: %lu %lu\n",
             child.depth,child.code,(unsigned long) child.lowerbound,
                                    (unsigned long) child.upperbound);
      */
      storeBoundsatdepth(pckbt,&child);
      if (child.depth < (unsigned long) maxdepth &&
          child.lowerbound + 1 < child.upperbound)
      {
        STOREINARRAY(&stack,Boundsatdepth,128,child);
      }
    }
  }
  FREEARRAY(&stack,Boundsatdepth);
  ma_free(rangeOccs);
  return pckbt;
}
