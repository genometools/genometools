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

#include "core/divmodmul.h"
#include "eis-bwtseq.h"
#include "eis-bwtseq-construct.h"
#include "eis-voiditf.h"
#include "splititv.h"
#include "procmatch.h"
#include "pckbucket.h"

/* XXX make use of the types declared by the EIS-tools, similar to
  the module core/bitpackarray.h */

unsigned long gt_bwtseqfirstmatch(const FMindex *voidbwtseq,
                                  unsigned long bound)
{
  struct extBitsRetrieval extBits;
  unsigned long pos;

  initExtBitsRetrieval(&extBits);
  pos = gt_BWTSeqLocateMatch((const BWTSeq *) voidbwtseq,bound,&extBits);
  destructExtBitsRetrieval(&extBits);
  return pos;
}

struct Bwtseqpositioniterator
{
  struct extBitsRetrieval extBits;
  const BWTSeq *bwtseq;
  unsigned long currentbound, upperbound;
};

Bwtseqpositioniterator *gt_newBwtseqpositioniterator(const FMindex *voidbwtseq,
                                                     unsigned long lowerbound,
                                                     unsigned long upperbound)
{
  Bwtseqpositioniterator *bspi;

  bspi = gt_malloc(sizeof (*bspi));
  initExtBitsRetrieval(&bspi->extBits);
  bspi->bwtseq = (const BWTSeq *) voidbwtseq;
  bspi->currentbound = lowerbound;
  bspi->upperbound = upperbound;
  return bspi;
}

bool gt_nextBwtseqpositioniterator(unsigned long *pos,
                                   Bwtseqpositioniterator *bspi)
{
  if (bspi->currentbound < bspi->upperbound)
  {
    *pos = gt_BWTSeqLocateMatch(bspi->bwtseq,bspi->currentbound,&bspi->extBits);
    bspi->currentbound++;
    return true;
  }
  return false;
}

bool gt_nextBwtseqpositionwithoutSEPiterator(unsigned long *pos,
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
      *pos = gt_BWTSeqLocateMatch(bspi->bwtseq,
                                  bspi->currentbound,
                                  &bspi->extBits);
      bspi->currentbound++;
      return true;
    }
    bspi->currentbound++;
  }
  return false;
}

void gt_freeBwtseqpositioniterator(Bwtseqpositioniterator **bspi)
{
  destructExtBitsRetrieval(&(*bspi)->extBits);
  gt_free(*bspi);
  *bspi = NULL;
}

struct Bwtseqcontextiterator
{
  struct extBitsRetrieval extBits;
  const BWTSeq *bwtseq;
  unsigned long bound;
};

Bwtseqcontextiterator *gt_newBwtseqcontextiterator(const FMindex *voidbwtseq,
                                                   unsigned long bound)
{
  Bwtseqcontextiterator *bsci;

  bsci = gt_malloc(sizeof (*bsci));
  initExtBitsRetrieval(&bsci->extBits);
  bsci->bwtseq = (const BWTSeq *) voidbwtseq;
  bsci->bound = bound;
  return bsci;
}

GtUchar gt_bwtseqgetsymbol(unsigned long bound,const FMindex *voidbwtseq)
{
  if (bound != BWTSeqTerminatorPos((const BWTSeq *) voidbwtseq))
  {
    return BWTSeqGetSym((const BWTSeq *) voidbwtseq, bound);
  }
  return SEPARATOR;
}

GtUchar gt_nextBwtseqcontextiterator(unsigned long *bound,
                                  Bwtseqcontextiterator *bsci)
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

void gt_freeBwtseqcontextiterator(Bwtseqcontextiterator **bsci)
{
  destructExtBitsRetrieval(&(*bsci)->extBits);
  gt_free(*bsci);
  *bsci = NULL;
}

FMindex *gt_loadvoidBWTSeqForSA(const char *indexname,
                                const Suffixarray *suffixarray,
                                unsigned long totallength,
                                bool withpckbt,
                                GtError *err)
{
  BWTSeq *bwtseq;
  bool haserr = false;

  bwtseq = gt_loadBWTSeqForSA(indexname,
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
    if (withpckbt && gt_pckbuckettableexists(indexname))
    {
      unsigned int numofchars
        = gt_alphabet_num_of_chars(
                           gt_encseq_alphabet(suffixarray->encseq));
      bwtseq->pckbuckettable = gt_mappckbuckettable(indexname,numofchars,err);
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
    gt_deletevoidBWTSeq((FMindex *) bwtseq);
    bwtseq = NULL;
  }
  return haserr ? NULL : (FMindex *) bwtseq;
}

void gt_bwtrangesplitwithoutspecial(GtArrayBoundswithchar *bwci,
                                    unsigned long *rangeOccs,
                                    const FMindex *voidBwtSeq,
                                    unsigned long lbound,
                                    unsigned long ubound)
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

const Mbtab **gt_bwtseq2mbtab(const FMindex *voidbwtseq)
{
  if (((const BWTSeq *) voidbwtseq)->pckbuckettable == NULL)
  {
    return NULL;
  }
  return (const Mbtab **)
         gt_pcktb2mbtab(((const BWTSeq *) voidbwtseq)->pckbuckettable);
}

unsigned int gt_bwtseq2maxdepth(const FMindex *voidbwtseq)
{
  if (((const BWTSeq *) voidbwtseq)->pckbuckettable == NULL)
  {
    return 0;
  }
  return gt_pcktb2maxdepth(((const BWTSeq *) voidbwtseq)->pckbuckettable);
}

unsigned long gt_bwtrangesplitallwithoutspecial(Mbtab *mbtab,
                                             unsigned long *rangeOccs,
                                             const FMindex *voidBwtSeq,
                                             unsigned long lbound,
                                             unsigned long ubound)
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
                         unsigned long *rangeOccs,
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
      pos = gt_BWTSeqLocateMatch((const BWTSeq *) voidbwtseq,bound,&extBits);
    }
  }
}
*/

void gt_deletevoidBWTSeq(FMindex *voidbwtseq)
{
  BWTSeq *bwtseq = (BWTSeq *) voidbwtseq;

  if (bwtseq->pckbuckettable != NULL)
  {
    gt_pckbuckettable_free(bwtseq->pckbuckettable);
    bwtseq->pckbuckettable = NULL;
  }
  gt_deleteBWTSeq(bwtseq);
}

unsigned long gt_voidpackedindexuniqueforward(const void *voidbwtseq,
                                              GT_UNUSED unsigned long offset,
                                              GT_UNUSED unsigned long left,
                                              GT_UNUSED unsigned long right,
                                              GT_UNUSED unsigned long
                                                        *witnessposition,
                                              const GtUchar *qstart,
                                              const GtUchar *qend)
{
  return gt_packedindexuniqueforward((const BWTSeq *) voidbwtseq,qstart,
                                     qend);
}

unsigned long gt_voidpackedfindfirstmatchconvert(const FMindex *voidbwtseq,
                                                 unsigned long witnessbound,
                                                 unsigned long matchlength)
{
  const BWTSeq *bwtseq = (const BWTSeq *) voidbwtseq;
  unsigned long startpos;

  startpos = gt_bwtseqfirstmatch(voidbwtseq,witnessbound);
  gt_assert((bwtseq->seqIdx->seqLen-1) >= (startpos + matchlength));
  return (bwtseq->seqIdx->seqLen - 1) - (startpos + matchlength);
}

unsigned long gt_voidpackedindexmstatsforward(const void *voidbwtseq,
                                              GT_UNUSED unsigned long offset,
                                              GT_UNUSED unsigned long left,
                                              GT_UNUSED unsigned long right,
                                              unsigned long *witnessposition,
                                              const GtUchar *qstart,
                                              const GtUchar *qend)
{
  const BWTSeq *bwtseq = (const BWTSeq *) voidbwtseq;
  unsigned long matchlength;

  matchlength = gt_packedindexmstatsforward(bwtseq,witnessposition,qstart,qend);
  if (matchlength > 0 && witnessposition != NULL)
  {
    *witnessposition = gt_voidpackedfindfirstmatchconvert(voidbwtseq,
                                                       *witnessposition,
                                                       matchlength);
  }
  return matchlength;
}

bool gt_pck_exactpatternmatching(const FMindex *voidbwtseq,
                              const GtUchar *pattern,
                              unsigned long patternlength,
                              unsigned long totallength,
                              const GtUchar *dbsubstring,
                              Processmatch processmatch,
                              void *processmatchinfo)
{
  BWTSeqExactMatchesIterator *bsemi;
  unsigned long dbstartpos, numofmatches;
  GtMatch match;

  bsemi = gt_newEMIterator((const BWTSeq *) voidbwtseq,
                        pattern,(size_t) patternlength, true);
  gt_assert(bsemi != NULL);
  numofmatches = gt_EMINumMatchesTotal(bsemi);
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
    gt_deleteEMIterator(bsemi);
    bsemi = NULL;
  }
  return numofmatches > 0 ? true : false;
}
