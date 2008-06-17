/*
  Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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
#include <stdlib.h>
#include <string.h>

#include "libgtcore/yarandom.h"

#include "libgtmatch/seqpos-def.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/esa-map.pr"

#include "libgtmatch/eis-encidxseq.h"
#include "libgtmatch/eis-encidxseq-priv.h"

enum {
  AVG_RANGE_RANK_INTERVAL = 128, /**< do range rank queries
                                  * approximately at this interval  */
};

void
deleteEncIdxSeq(EISeq *seq)
{
  seq->classInfo->delete(seq);
}

const char *EISIntegrityCheckResultStrings[] =
{
  "no error in bwt comparison.",
  "reading from the index produced an incorrect symbol.",
  "reading from the BWT reference file failed or delivered a symbol"
  " not in the alphabet (bwt file corrupt?) ",
  "loading/mapping of the suffix array project failed. \n"
  "(did you generate the BWT?)",
  "the rank operation delivered a wrong count"
};

#define verifyIntegrityErrRet(retcode)                                  \
  {                                                                     \
    switch (retcode) {                                                  \
    case EIS_INTEGRITY_CHECK_INVALID_SYMBOL:                            \
      fprintf(stderr, "Comparision failed at position "FormatSeqpos     \
              ", reference symbol: %u, symbol read: %u\n",              \
              pos, symOrig, symEnc);                                    \
      error_set(err, "Invalid symbol encountered.");                    \
      break;                                                            \
    case EIS_INTEGRITY_CHECK_BWT_READ_ERROR:                            \
      fprintf(stderr, "Read of symbol failed at position "              \
              FormatSeqpos"\n", pos);                                   \
      error_set(err, "Failed reading reference BWT source.");           \
      break;                                                            \
    case EIS_INTEGRITY_CHECK_RANK_FAILED:                               \
      fprintf(stderr, "At position "FormatSeqpos                        \
              ", rank operation yielded  wrong count: "FormatSeqpos     \
              ", expected "FormatSeqpos" for symbol %d\n",              \
              pos, rankQueryResult, rankExpect, rankCmpSym);            \
      error_set(err, "Invalid rank result.");                           \
      break;                                                            \
    }                                                                   \
    EISPrintDiagsForPos(seqIdx, pos, stderr, hint);                     \
    retval = retcode;                                                   \
    break;                                                              \
  } do {} while (0)

/**
 * @param tickPrint if not zero, print a . every tickPrint symbols to
 * fp
 * @param fp descriptor to write progress marks to.
 * @return -1 on error, 0 on identity, >0 on inconsistency
 */
extern enum EISIntegrityCheckResults
EISVerifyIntegrity(EISeq *seqIdx, const Str *projectName, Seqpos skip,
                   unsigned long tickPrint, FILE *fp, int chkFlags,
                   Verboseinfo *verbosity, Error *err)
{
  FILE *bwtFP;
  Seqpos pos = 0;
  Suffixarray suffixArray;
  Symbol symOrig;
  unsigned symEnc;
  EISHint hint;
  Seqpos seqLastPos, rankQueryResult, rankExpect;
  const MRAEnc *alphabet;
  AlphabetRangeSize alphabetSize;
  AlphabetRangeID numRanges;
  enum EISIntegrityCheckResults retval = EIS_INTEGRITY_CHECK_NO_ERROR;
  /* two part process: enumerate all positions of original sequence
   * and verify that the query functions return correct values */
  if (streamsuffixarray(&suffixArray, &seqLastPos,
                        SARR_BWTTAB, projectName, verbosity, err))
  {
    error_set(err, "Cannot load suffix array project with"
                  " demand for BWT file\n");
    return EIS_INTEGRITY_CHECK_SA_LOAD_ERROR;
  }
  bwtFP = suffixArray.bwttabstream.fp;
  hint = newEISHint(seqIdx);
  alphabet = EISGetAlphabet(seqIdx);
  alphabetSize = MRAEncGetSize(alphabet);
  numRanges = MRAEncGetNumRanges(alphabet);
  do
  {
    Seqpos rankTable[alphabetSize], rangeRanks[2][alphabetSize],
      pairRangeRanks[2* alphabetSize], lastRangeRankPos = 0;
    int symRead, rt = 0;
    AlphabetRangeID lastRangeID = 0;
    unsigned rankCmpSym;
    memset(rankTable, 0, sizeof (rankTable));
    memset(rangeRanks, 0, sizeof (rangeRanks));
    if (skip > 0)
    {
      Seqpos len = EISLength(seqIdx);
      unsigned sym;
      if (skip >= len)
      {
        showverbose(verbosity, "Invalid skip request: %lld,"
                    " too large for sequence length: "FormatSeqpos,
                    (long long)skip, len);
        return -1;
      }
      fseeko(bwtFP, skip, SEEK_SET);
      for (sym = 0; sym <= UCHAR_MAX; ++sym)
        if (MRAEncSymbolHasValidMapping(alphabet, sym))
          rankTable[sym] = EISRank(seqIdx, sym, skip, hint);
      pos = skip;
    }
    while (retval == EIS_INTEGRITY_CHECK_NO_ERROR
           && (symRead = getc(bwtFP)) != EOF)
    {
      symOrig = MRAEncMapSymbol(alphabet, symRead);
      symEnc = EISGetTransformedSym(seqIdx, pos, hint);
      if (symEnc >= alphabetSize)
        verifyIntegrityErrRet(EIS_INTEGRITY_CHECK_INVALID_SYMBOL);
      if (symOrig >= alphabetSize)
        verifyIntegrityErrRet(EIS_INTEGRITY_CHECK_BWT_READ_ERROR);
      if (symEnc != symOrig)
        verifyIntegrityErrRet(EIS_INTEGRITY_CHECK_INVALID_SYMBOL);
      ++rankTable[symOrig];
      if (chkFlags & EIS_VERIFY_EXT_RANK)
      {
        for (rankCmpSym = 0; rankCmpSym < alphabetSize; ++rankCmpSym)
          if ((rankExpect = rankTable[rankCmpSym])
              != (rankQueryResult
                  = EISSymTransformedRank(seqIdx, rankCmpSym, pos + 1, hint)))
            verifyIntegrityErrRet(EIS_INTEGRITY_CHECK_RANK_FAILED);
      }
      else
      {
        rankCmpSym = symEnc;
        if ((rankExpect = rankTable[symEnc])
            != (rankQueryResult = EISSymTransformedRank(seqIdx, symEnc,
                                                        pos + 1, hint)))
          verifyIntegrityErrRet(EIS_INTEGRITY_CHECK_RANK_FAILED);
      }
      /* do rank for full range on some occasions */
      if (!(random() % AVG_RANGE_RANK_INTERVAL))
      {
        unsigned i;
        AlphabetRangeSize rangeSize;
        AlphabetRangeID range = random() % numRanges;
        Symbol rangeBase = MRAEncGetRangeBase(alphabet, range);
        rangeSize = MRAEncGetRangeSize(alphabet, range);
        EISRangeRank(seqIdx, range, pos + 1, rangeRanks[rt], hint);
        for (i = 0; i < rangeSize; ++i)
        {
          rankCmpSym = rangeBase + i;
          if ((rankQueryResult = rangeRanks[rt][i])
              != (rankExpect = rankTable[rankCmpSym]))
            verifyIntegrityErrRet(EIS_INTEGRITY_CHECK_RANK_FAILED);
        }
        if (retval)
          break;
        if (range == lastRangeID)
        {
          EISPosPairRangeRank(seqIdx, range, lastRangeRankPos, pos + 1,
                              pairRangeRanks, hint);
          for (i = 0; i < rangeSize; ++i)
          {
            rankCmpSym = rangeBase + i;
            if ((rankQueryResult = pairRangeRanks[rangeSize + i])
                != (rankExpect = rankTable[rankCmpSym]))
              verifyIntegrityErrRet(EIS_INTEGRITY_CHECK_RANK_FAILED);
            if ((rankQueryResult = pairRangeRanks[i])
                != (rankExpect = rangeRanks[rt ^ 1][rankCmpSym]))
              verifyIntegrityErrRet(EIS_INTEGRITY_CHECK_RANK_FAILED);
          }
          if (retval)
            break;
        }
        lastRangeID =
        lastRangeRankPos = pos + 1;
        rt ^= 1;
      }
      ++pos;
      if (tickPrint && !(pos % tickPrint))
        putc('.', fp);
    }
    if (retval)
      break;
    if (tickPrint)
      putc('\n', fp);
    if (ferror(bwtFP))
      verifyIntegrityErrRet(EIS_INTEGRITY_CHECK_BWT_READ_ERROR);
  } while (0);
  deleteEISHint(seqIdx, hint);
  freesuffixarray(&suffixArray);
  return retval;
}

extern unsigned
estimateSegmentSize(const struct seqBaseParam *params)
{
  unsigned segmentLen = 0;
  switch (params->encType)
  {
  case BWT_ON_BLOCK_ENC:
    segmentLen =
      blockEncIdxSeqSegmentLen(&params->encParams.blockEnc);
    break;
  default:
    fputs("Illegal/unknown/unimplemented encoding requested!", stderr);
    abort();
    break;
  }
  return segmentLen;
}
