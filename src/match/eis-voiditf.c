/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "procmatch.h"
#include "pckbucket.h"

/* XXX make use of the types declared by the EIS-tools, similar to
  the module bitpack-itf.h */

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

  bspi = gt_malloc(sizeof (*bspi));
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
    GtUchar cc;

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
  gt_free(*bspi);
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

  bsci = gt_malloc(sizeof (*bsci));
  initExtBitsRetrieval(&bsci->extBits);
  bsci->bwtseq = (const BWTSeq *) voidbwtseq;
  bsci->bound = bound;
  return bsci;
}

GtUchar bwtseqgetsymbol(Seqpos bound,const void *voidbwtseq)
{
  if (bound != BWTSeqTerminatorPos(voidbwtseq))
  {
    return BWTSeqGetSym(voidbwtseq, bound);
  }
  return SEPARATOR;
}

GtUchar nextBwtseqcontextiterator(Seqpos *bound,Bwtseqcontextiterator *bsci)
{
  GtUchar cc;

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
  gt_free(*bsci);
  *bsci = NULL;
}

void *loadvoidBWTSeqForSA(const GtStr *indexname,
                          const Suffixarray *suffixarray,
                          Seqpos totallength,
                          bool withpckbt,
                          GtError *err)
{
  BWTSeq *bwtseq;
  bool haserr = false;

  bwtseq = loadBWTSeqForSA(indexname,
                           BWT_ON_BLOCK_ENC,
                           BWTDEFOPT_MULTI_QUERY,
                           suffixarray,
                           totallength+1, err);
  if (bwtseq == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (withpckbt && pckbuckettableexists(indexname))
    {
      unsigned int numofchars
        = getencseqAlphabetnumofchars(suffixarray->encseq);
      bwtseq->pckbuckettable = mappckbuckettable(indexname,numofchars,err);
      if (bwtseq->pckbuckettable == NULL)
      {
        haserr = true;
      }
    } else
    {
      bwtseq->pckbuckettable = NULL;
    }
  }
  if (haserr && bwtseq != NULL)
  {
    deletevoidBWTSeq(bwtseq);
    bwtseq = NULL;
  }
  return haserr ? NULL : bwtseq;
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

const Mbtab **bwtseq2mbtab(const void *voidbwtseq)
{
  if (((const BWTSeq *) voidbwtseq)->pckbuckettable == NULL)
  {
    return NULL;
  }
  return (const Mbtab **)
         pcktb2mbtab(((const BWTSeq *) voidbwtseq)->pckbuckettable);
}

unsigned int bwtseq2maxdepth(const void *voidbwtseq)
{
  if (((const BWTSeq *) voidbwtseq)->pckbuckettable == NULL)
  {
    return 0;
  }
  return pcktb2maxdepth(((const BWTSeq *) voidbwtseq)->pckbuckettable);
}

unsigned long bwtrangesplitallwithoutspecial(Mbtab *mbtab,
                                             Seqpos *rangeOccs,
                                             const void *voidBwtSeq,
                                             Seqpos lbound,
                                             Seqpos ubound)
{
  unsigned long idx;
  const BWTSeq *bwtseq = (const BWTSeq *) voidBwtSeq;
  AlphabetRangeSize rangesize
    = MRAEncGetRangeSize(EISGetAlphabet(bwtseq->seqIdx),0);

  BWTSeqPosPairRangeOcc(bwtseq, 0, lbound, ubound,rangeOccs);
  for (idx = 0; idx < rangesize; idx++)
  {
    if (rangeOccs[idx] < rangeOccs[rangesize+idx])
    {
      mbtab[idx].lowerbound = bwtseq->count[idx] + rangeOccs[idx];
      mbtab[idx].upperbound = bwtseq->count[idx] + rangeOccs[rangesize+idx];
    } else
    {
      mbtab[idx].lowerbound = mbtab[idx].upperbound = 0;
    }
  }
  return rangesize;
}

/*
void bwtrangewithspecial(GT_UNUSED GtArrayBoundswithchar *bwci,
                         Seqpos *rangeOccs,
                         GT_UNUSED unsigned long numofchars,
                         const void *voidBwtSeq,
                         const Lcpinterval *parent)
{
  const BWTSeq *bwtseq = (const BWTSeq *) voidBwtSeq;

  AlphabetRangeSize rangesize
    = MRAEncGetRangeSize(EISGetAlphabet(bwtseq->seqIdx),1);
  gt_assert(rangesize < (AlphabetRangeSize) 4);
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

void deletevoidBWTSeq(void *voidbwtseq)
{
  BWTSeq *bwtseq = (BWTSeq *) voidbwtseq;

  if (bwtseq->pckbuckettable != NULL)
  {
    pckbuckettable_free(bwtseq->pckbuckettable);
    bwtseq->pckbuckettable = NULL;
  }
  deleteBWTSeq(bwtseq);
}

unsigned long voidpackedindexuniqueforward(const void *voidbwtseq,
                                           GT_UNUSED unsigned long offset,
                                           GT_UNUSED Seqpos left,
                                           GT_UNUSED Seqpos right,
                                           GT_UNUSED Seqpos *witnessposition,
                                           const GtUchar *qstart,
                                           const GtUchar *qend)
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
  gt_assert((bwtseq->seqIdx->seqLen-1) >= (startpos + matchlength));
  return (bwtseq->seqIdx->seqLen - 1) - (startpos + matchlength);
}

unsigned long voidpackedindexmstatsforward(const void *voidbwtseq,
                                           GT_UNUSED unsigned long offset,
                                           GT_UNUSED Seqpos left,
                                           GT_UNUSED Seqpos right,
                                           Seqpos *witnessposition,
                                           const GtUchar *qstart,
                                           const GtUchar *qend)
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

bool pck_exactpatternmatching(const void *voidbwtseq,
                              const GtUchar *pattern,
                              unsigned long patternlength,
                              Seqpos totallength,
                              const GtUchar *dbsubstring,
                              Processmatch processmatch,
                              void *processmatchinfo)
{
  BWTSeqExactMatchesIterator *bsemi;
  Seqpos dbstartpos, numofmatches;
  GtMatch match;

  bsemi = newEMIterator((const BWTSeq *) voidbwtseq,
                        pattern,(size_t) patternlength, true);
  gt_assert(bsemi != NULL);
  numofmatches = EMINumMatchesTotal(bsemi);
  match.dbabsolute = true;
  match.dblen = patternlength;
  match.dbsubstring = dbsubstring;
  match.querystartpos = 0;
  match.querylen = patternlength;
  match.distance = 0;
  match.alignment = NULL;
  while (EMIGetNextMatch(bsemi,&dbstartpos,(const BWTSeq *) voidbwtseq))
  {
    gt_assert(totallength >= (dbstartpos + patternlength));
    match.dbstartpos = totallength - (dbstartpos + patternlength);
    processmatch(processmatchinfo,&match);
  }
  if (bsemi != NULL)
  {
    deleteEMIterator(bsemi);
    bsemi = NULL;
  }
  return numofmatches > 0 ? true : false;
}
