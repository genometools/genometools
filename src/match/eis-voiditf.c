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
#include "core/encseq_metadata.h"
#include "core/log_api.h"
#include "eis-bwtseq-construct.h"
#include "eis-bwtseq-priv.h"
#include "eis-bwtseq.h"
#include "eis-voiditf.h"
#include "procmatch.h"

/* XXX make use of the types declared by the EIS-tools, similar to
  the module core/bitpackarray.h */

unsigned long gt_bwtseqfirstmatch(const FMindex *fmindex,
                                  unsigned long bound)
{
  struct extBitsRetrieval extBits;
  unsigned long pos;

  initExtBitsRetrieval(&extBits);
  pos = gt_BWTSeqLocateMatch((const BWTSeq *) fmindex,bound,&extBits);
  destructExtBitsRetrieval(&extBits);
  return pos;
}

struct Bwtseqpositioniterator
{
  struct extBitsRetrieval extBits;
  const BWTSeq *bwtseq;
  unsigned long currentbound, upperbound;
};

Bwtseqpositioniterator *gt_Bwtseqpositioniterator_new(const FMindex *fmindex,
                                                      unsigned long lowerbound,
                                                      unsigned long upperbound)
{
  Bwtseqpositioniterator *bspi;

  bspi = gt_malloc(sizeof (*bspi));
  initExtBitsRetrieval(&bspi->extBits);
  bspi->bwtseq = (const BWTSeq *) fmindex;
  bspi->currentbound = lowerbound;
  bspi->upperbound = upperbound;
  return bspi;
}

bool gt_Bwtseqpositioniterator_next(unsigned long *pos,
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

bool gt_BwtseqpositionwithoutSEPiterator_next(unsigned long *pos,
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

void gt_Bwtseqpositioniterator_delete(Bwtseqpositioniterator *bspi)
{
  destructExtBitsRetrieval(&bspi->extBits);
  gt_free(bspi);
}

struct BwtSeqpositionextractor
{
  struct extBitsRetrieval extBits;
  const BWTSeq *bwtseq;
  unsigned long upperbound;
};

BwtSeqpositionextractor *gt_newBwtSeqpositionextractor(
                                                     const FMindex *voidbwtseq,
                                                     unsigned long upperbound)
{
  BwtSeqpositionextractor *bspex;

  bspex = gt_malloc(sizeof (*bspex));
  initExtBitsRetrieval(&bspex->extBits);
  bspex->bwtseq = (const BWTSeq *) voidbwtseq;
  bspex->upperbound = upperbound;
  return bspex;
}

unsigned long gt_BwtSeqpositionextractor_extract(BwtSeqpositionextractor *bspex,
                                                 unsigned long lowerbound)
{
  unsigned long pos;

  gt_assert(lowerbound < bspex->upperbound);
  pos = gt_BWTSeqLocateMatch(bspex->bwtseq,lowerbound,&bspex->extBits);
  return pos;
}

void gt_freeBwtSeqpositionextractor(BwtSeqpositionextractor *bspex)
{
  destructExtBitsRetrieval(&(bspex->extBits));
  gt_free(bspex);
}

struct Bwtseqcontextiterator
{
  struct extBitsRetrieval extBits;
  const BWTSeq *bwtseq;
  unsigned long bound;
};

Bwtseqcontextiterator *gt_Bwtseqcontextiterator_new(const FMindex *fmindex,
                                                    unsigned long bound)
{
  Bwtseqcontextiterator *bsci;

  bsci = gt_malloc(sizeof (*bsci));
  initExtBitsRetrieval(&bsci->extBits);
  bsci->bwtseq = (const BWTSeq *) fmindex;
  bsci->bound = bound;
  return bsci;
}

GtUchar gt_bwtseqgetsymbol(unsigned long bound,const FMindex *fmindex)
{
  if (bound != BWTSeqTerminatorPos((const BWTSeq *) fmindex))
  {
    return BWTSeqGetSym((const BWTSeq *) fmindex, bound);
  }
  return SEPARATOR;
}

GtUchar gt_Bwtseqcontextiterator_next(unsigned long *bound,
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

void gt_Bwtseqcontextiterator_delete(Bwtseqcontextiterator *bsci)
{
  if (bsci != NULL)
  {
    destructExtBitsRetrieval(&bsci->extBits);
    gt_free(bsci);
  }
}

FMindex *gt_loadvoidBWTSeqForSA(const char *indexname,
                                bool withpckbt,
                                GtError *err)
{
  BWTSeq *bwtseq = NULL;
  bool haserr = false;
  GtEncseqMetadata *emd;
  GtAlphabet *alphabet;

  emd = gt_encseq_metadata_new(indexname, err);
  if (emd == NULL) {
    gt_assert(gt_error_is_set(err));
    haserr = true;
  }
  if (!haserr) {
    alphabet = gt_encseq_metadata_alphabet(emd);
    if (alphabet == NULL)
    {
      gt_assert(gt_error_is_set(err));
      haserr = true;
    }
  }
  if (!haserr)
  {
    bwtseq = gt_loadBWTSeqForSA(indexname,
                                BWT_ON_BLOCK_ENC,
                                BWTDEFOPT_MULTI_QUERY,
                                alphabet,
                                err);
    if (bwtseq == NULL)
    {
      gt_assert(gt_error_is_set(err));
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (withpckbt && gt_pckbuckettable_exists(indexname))
    {
      unsigned int numofchars = gt_alphabet_num_of_chars(alphabet);
      bwtseq->pckbuckettable = gt_pckbuckettable_map(indexname,numofchars,err);
      if (bwtseq->pckbuckettable == NULL)
      {
        gt_assert(gt_error_is_set(err));
        haserr = true;
      }
    } else
    {
      bwtseq->pckbuckettable = NULL;
    }
  }
  gt_encseq_metadata_delete(emd);
  if (haserr && bwtseq != NULL)
  {
    gt_deletevoidBWTSeq((FMindex *) bwtseq);
    bwtseq = NULL;
  }
  return haserr ? NULL : (FMindex *) bwtseq;
}

const Mbtab **gt_bwtseq2mbtab(const FMindex *fmindex)
{
  gt_assert(fmindex != NULL);
  if (((const BWTSeq *) fmindex)->pckbuckettable == NULL)
  {
    return NULL;
  }
  return (const Mbtab **)
         gt_pckbuckettable_mbtab_get(((const BWTSeq *) fmindex)
                                     ->pckbuckettable);
}

unsigned int gt_bwtseq2maxdepth(const FMindex *fmindex)
{
  gt_assert(fmindex != NULL);
  if (((const BWTSeq *) fmindex)->pckbuckettable == NULL)
  {
    return 0;
  }
  return gt_pckbuckettable_maxdepth_get(((const BWTSeq *) fmindex)
                                        ->pckbuckettable);
}

unsigned int gt_bwtseq2numofchars(const FMindex *fmindex)
{
  gt_assert(fmindex != NULL);
  if (((const BWTSeq *) fmindex)->pckbuckettable == NULL)
  {
    return 0;
  }
  return gt_pckbuckettable_numofchars_get(((const BWTSeq *) fmindex)
                                          ->pckbuckettable);
}

void gt_bwtrangesplitwithoutspecial(GtArrayBoundswithchar *bwci,
                                    unsigned long *rangeOccs,
                                    const FMindex *fmindex,
                                    unsigned long lbound,
                                    unsigned long ubound)
{
  const BWTSeq *bwtseq = (const BWTSeq *) fmindex;
  AlphabetRangeSize idx, rangesize;

  rangesize = MRAEncGetRangeSize(EISGetAlphabet(bwtseq->seqIdx),0);
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

unsigned long gt_bwtrangesplitallwithoutspecial(Mbtab *mbtab,
                                                unsigned long *rangeOccs,
                                                const FMindex *fmindex,
                                                unsigned long lbound,
                                                unsigned long ubound)
{
  const BWTSeq *bwtseq = (const BWTSeq *) fmindex;
  AlphabetRangeSize idx, rangesize
    = MRAEncGetRangeSize(EISGetAlphabet(bwtseq->seqIdx),0);

  gt_assert(sizeof (AlphabetRangeSize) <= sizeof (unsigned long));
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
  return (unsigned long) rangesize;
}

unsigned long gt_bwtrangesplitallwithspecial(Mbtab *mbtab,
                                             unsigned long *rangeOccs,
                                             const FMindex *voidBwtSeq,
                                             unsigned long lbound,
                                             unsigned long ubound)
{
  unsigned long char_idx, range_idx, rangebase;
  const BWTSeq *bwtseq = (const BWTSeq *) voidBwtSeq;
  const MRAEnc *alphabet = BWTSeqGetAlphabet(bwtseq);
  AlphabetRangeID numofranges = MRAEncGetNumRanges(alphabet);
  AlphabetRangeSize rangesize = 0, totalrange = 0;

  for (range_idx = 0; range_idx < (unsigned long) numofranges; range_idx++)
  {
    unsigned long rangeOcc_idx = 0;
    rangesize = MRAEncGetRangeSize(alphabet, range_idx);
    totalrange += rangesize;
    BWTSeqPosPairRangeOcc(bwtseq, range_idx, lbound, ubound,rangeOccs);
    rangebase = (unsigned long) MRAEncGetRangeBase(alphabet, range_idx);
    for (char_idx = rangebase;
         char_idx < rangebase + rangesize; char_idx++)
    {
      if (rangeOccs[rangeOcc_idx] < rangeOccs[rangesize+rangeOcc_idx])
      {
        mbtab[char_idx].lowerbound = bwtseq->count[char_idx] +
                                     rangeOccs[rangeOcc_idx];
        mbtab[char_idx].upperbound = bwtseq->count[char_idx] +
                                     rangeOccs[rangesize+rangeOcc_idx];
      } else
      {
        mbtab[char_idx].lowerbound = mbtab[char_idx].upperbound = 0;
      }
      rangeOcc_idx++;
    }
  }
  return totalrange;
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
      pos = gt_BWTSeqLocateMatch((const BWTSeq *) fmindex,bound,&extBits);
    }
  }
}
*/

void gt_deletevoidBWTSeq(FMindex *fmindex)
{
  BWTSeq *bwtseq = (BWTSeq *) fmindex;

  if (bwtseq->pckbuckettable != NULL)
  {
    gt_pckbuckettable_delete(bwtseq->pckbuckettable);
    bwtseq->pckbuckettable = NULL;
  }
  gt_deleteBWTSeq(bwtseq);
}

unsigned long gt_voidpackedindexuniqueforward(const void *fmindex,
                                              GT_UNUSED unsigned long offset,
                                              GT_UNUSED unsigned long left,
                                              GT_UNUSED unsigned long right,
                                              GT_UNUSED unsigned long
                                                        *witnessposition,
                                              const GtUchar *qstart,
                                              const GtUchar *qend)
{
  return gt_packedindexuniqueforward((const BWTSeq *) fmindex,qstart,qend);
}

unsigned long gt_voidpackedfindfirstmatchconvert(const FMindex *fmindex,
                                                 unsigned long witnessbound,
                                                 unsigned long matchlength)
{
  const BWTSeq *bwtseq = (const BWTSeq *) fmindex;
  unsigned long startpos;

  startpos = gt_bwtseqfirstmatch(fmindex,witnessbound);
  gt_assert((bwtseq->seqIdx->seqLen-1) >= (startpos + matchlength));
  return (bwtseq->seqIdx->seqLen - 1) - (startpos + matchlength);
}

unsigned long gt_voidpackedindexmstatsforward(const void *fmindex,
                                              GT_UNUSED unsigned long offset,
                                              GT_UNUSED unsigned long left,
                                              GT_UNUSED unsigned long right,
                                              unsigned long *witnessposition,
                                              const GtUchar *qstart,
                                              const GtUchar *qend)
{
  const BWTSeq *bwtseq = (const BWTSeq *) fmindex;
  unsigned long matchlength;

  matchlength = gt_packedindexmstatsforward(bwtseq,witnessposition,qstart,qend);
  if (matchlength > 0 && witnessposition != NULL)
  {
    *witnessposition = gt_voidpackedfindfirstmatchconvert(fmindex,
                                                          *witnessposition,
                                                          matchlength);
  }
  return matchlength;
}

bool gt_pck_exactpatternmatching(const FMindex *fmindex,
                                 const GtUchar *pattern,
                                 unsigned long patternlength,
                                 unsigned long totallength,
                                 const GtUchar *dbsubstring,
                                 ProcessIdxMatch processmatch,
                                 void *processmatchinfo)
{
  BWTSeqExactMatchesIterator *bsemi;
  unsigned long dbstartpos, numofmatches;
  GtIdxMatch match;

  bsemi = gt_newEMIterator((const BWTSeq *) fmindex,
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
  while (EMIGetNextMatch(bsemi,&dbstartpos,(const BWTSeq *) fmindex))
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

unsigned long gt_voidpackedindex_totallength_get(const FMindex *fmindex)
{
  unsigned long bwtlen = BWTSeqLength((const BWTSeq *) fmindex);

  gt_assert(bwtlen > 0);
  return bwtlen - 1;
}

unsigned long gt_pck_getShuStringLength(const FMindex *bwtSubject,
                                       const GtUchar *suffix,
                                       unsigned long suffixLength)
{
  GtUlongPair occPair;
  Symbol curChar;
  const GtUchar *qptr, *qend;
  const MRAEnc *alphabet;
  unsigned long start, end;

  gt_assert(bwtSubject && suffix);
  alphabet = BWTSeqGetAlphabet((const BWTSeq *) bwtSubject);

  qptr = suffix;
  qend = suffix + suffixLength;

  curChar = MRAEncMapSymbol(alphabet, *qptr);

  qptr++;
  start = ((const BWTSeq *) bwtSubject)->count[curChar];
  end = ((const BWTSeq *) bwtSubject)->count[curChar + 1];
  for (/* Nothing */; start < end && qptr < qend; qptr++)
  {
    curChar = MRAEncMapSymbol(alphabet, *qptr);
    occPair = BWTSeqTransformedPosPairOcc((const BWTSeq *) bwtSubject,
                                          curChar,
                                          start,
                                          end);
    start = ((const BWTSeq *) bwtSubject)->count[curChar] + occPair.a;
    end = ((const BWTSeq *) bwtSubject)->count[curChar] + occPair.b;
  }
  if (qptr == qend && start < end)
    return suffixLength + 1;
  else
    return qptr - suffix;
}

double gt_pck_getGCcontent(const FMindex *bwtSubject,
                           const GtAlphabet *alphabet)
{
  unsigned long  c, length;
  double gc;
  const MRAEnc *FM_alphabet;
  GtUchar c_sym;

  FM_alphabet = BWTSeqGetAlphabet((const BWTSeq *) bwtSubject);

  c_sym = MRAEncMapSymbol(FM_alphabet,
                          gt_alphabet_encode(alphabet, 'c'));
  length = ((const BWTSeq *) bwtSubject)->seqIdx->seqLen;

  c = ((const BWTSeq *) bwtSubject)->count[c_sym+1] -
    ((const BWTSeq *) bwtSubject)->count[c_sym];

  gc = c * 2 / (double) (length - 2);
  return gc;
}

unsigned long gt_pck_get_nonspecial_count(const FMindex *index)
{
  const BWTSeq *bwtseq = (const BWTSeq *) index;

  return BWTSeqAggTransformedCount(bwtseq,
                                   MRAEncGetRangeSize(BWTSeqGetAlphabet(bwtseq),
                                     0));
}

unsigned long gt_pck_special_occ_in_nonspecial_intervals(const FMindex *index)
{
  const BWTSeq *bwtseq = (const BWTSeq *) index;
  unsigned long count = 0,
                *rangeOccs, first_special_row;
  unsigned short idx, rangesize;

  rangesize = (unsigned short) MRAEncGetRangeSize(BWTSeqGetAlphabet(bwtseq), 1);
  rangeOccs = gt_calloc(rangesize, sizeof (unsigned long));
  first_special_row = gt_pck_get_nonspecial_count(index);
  BWTSeqRangeOcc(bwtseq, 1, first_special_row, rangeOccs);
  for (idx = 0; idx < rangesize; idx++)
    count += rangeOccs[idx];
  gt_free(rangeOccs);
  return count;
}

unsigned long gt_pck_exact_pattern_count(const FMindex *index,
                                         const GtUchar *pattern,
                                         unsigned long patternlength) {
  return gt_BWTSeqMatchCount((const BWTSeq *) index,
                              pattern,
                              (size_t) patternlength,
                              true); /* XXX check if this is right! */
}
