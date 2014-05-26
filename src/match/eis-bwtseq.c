/*
  Copyright (c) 2007,2008 Thomas Jahns <Thomas.Jahns@gmx.net>

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

#include <string.h>
#include "core/assert_api.h"
#include "match/dataalign.h"
#include "core/error.h"
#include "core/log.h"
#include "core/str.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/undef_api.h"
#include "core/yarandom.h"
#include "match/eis-bitpackseqpos.h"
#include "match/eis-bwtseq.h"
#include "match/eis-bwtseq-extinfo.h"
#include "match/eis-bwtseq-param.h"
#include "match/eis-bwtseq-priv.h"
#include "match/eis-bwtseq-context.h"
#include "match/eis-encidxseq.h"
#include "match/eis-mrangealphabet.h"
#include "match/eis-suffixerator-interface.h"
#include "match/eis-suffixarray-interface.h"
#include "match/esa-map.h"
#include "match/sarr-def.h"

/**
 * @param alphabet ownership of alphabet is with the newly produced
 * sequence object if return value is not 0
 */
static int
initBWTSeqFromEncSeqIdx(BWTSeq *bwtSeq, struct encIdxSeq *seqIdx,
                        MRAEnc *alphabet, GtUword *counts,
                        enum rangeSortMode *rangeSort,
                        const enum rangeSortMode *defaultRangeSort)
{
  size_t alphabetSize;
  Symbol bwtTerminatorFlat;
  EISHint hint;
  gt_assert(bwtSeq && seqIdx);
  bwtSeq->alphabet = alphabet;
  alphabetSize = gt_MRAEncGetSize(alphabet);
  if (!alphabetSize)
    /* weird error, shouldn't happen, but I prefer error return to
     * segfault in case someone tampered with the input */
    return 0;
  /* FIXME: this should probably be handled in chardef.h to have a
   * unique mapping */
  /* FIXME: this assumes there is exactly two ranges */
  gt_MRAEncAddSymbolToRange(alphabet, bwtTerminatorSym, 1);
  gt_assert(gt_MRAEncGetSize(alphabet) ==  alphabetSize + 1);
  alphabetSize = gt_MRAEncGetSize(alphabet);
  bwtSeq->bwtTerminatorFallback = bwtTerminatorFlat
    = MRAEncMapSymbol(alphabet, UNDEFBWTCHAR);
  bwtSeq->bwtTerminatorFallbackRange = 1;
  bwtSeq->count = counts;
  bwtSeq->rangeSort = rangeSort;
  bwtSeq->seqIdx = seqIdx;
  bwtSeq->alphabetSize = alphabetSize;
  bwtSeq->hint = hint = newEISHint(seqIdx);
  {
    Symbol i;
    GtUword len = EISLength(seqIdx), *count = bwtSeq->count;
    count[0] = 0;
    for (i = 0; i < bwtTerminatorFlat; ++i)
      count[i + 1] = count[i]
        + EISSymTransformedRank(seqIdx, i, len, hint);
    /* handle character which the terminator has been mapped to specially */
    count[i + 1] = count[i]
      + EISSymTransformedRank(seqIdx, i, len, hint) - 1;
    gt_assert(count[i + 1] >= count[i]);
    /* now we can finish the rest of the symbols */
    for (i += 2; i < alphabetSize; ++i)
      count[i] = count[i - 1]
        + EISSymTransformedRank(seqIdx, i - 1, len, hint);
    /* and finally place the 1-count for the terminator */
    count[i] = count[i - 1] + 1;
#ifdef EIS_DEBUG
    gt_log_log("count[alphabetSize]="GT_WU", len="GT_WU"",count[alphabetSize],
               len);
    for (i = 0; i <= alphabetSize; ++i)
      gt_log_log("count[%u]="GT_WU"", (unsigned)i, count[i]);
#endif
    gt_assert(count[alphabetSize] == len);
  }
  gt_BWTSeqInitLocateHandling(bwtSeq, defaultRangeSort);
  return 1;
}

/**
 * @param alphabet ownership of alphabet is with the newly produced
 * sequence object if return value is non-NULL
 */
BWTSeq *
gt_newBWTSeq(EISeq *seqIdx, MRAEnc *alphabet,
          const enum rangeSortMode *defaultRangeSort)
{
  BWTSeq *bwtSeq;
  GtUword *counts;
  size_t countsOffset, rangeSortOffset, totalSize;
  enum rangeSortMode *rangeSort;
  unsigned alphabetSize;
  gt_assert(seqIdx);
  /* alphabetSize is increased by one to handle the flattened
   * terminator symbol correctly */
  alphabetSize = gt_MRAEncGetSize(alphabet) + 1;
  countsOffset = offsetAlign(sizeof (struct BWTSeq), sizeof (GtUword));
  rangeSortOffset = offsetAlign(countsOffset
                                + sizeof (GtUword) * (alphabetSize + 1),
                                sizeof (enum rangeSortMode));
  totalSize = rangeSortOffset + sizeof (enum rangeSortMode)
    * MRAEncGetNumRanges(alphabet);
  bwtSeq = gt_malloc(totalSize);
  bwtSeq->pckbuckettable = NULL;
  counts = (GtUword *)((char  *)bwtSeq + countsOffset);
  rangeSort = (enum rangeSortMode *)((char *)bwtSeq + rangeSortOffset);
  if (!initBWTSeqFromEncSeqIdx(bwtSeq, seqIdx, alphabet, counts, rangeSort,
                               defaultRangeSort))
  {
    gt_free(bwtSeq);
    bwtSeq = NULL;
  }
  return bwtSeq;
}

void
gt_deleteBWTSeq(BWTSeq *bwtSeq)
{
  gt_MRAEncDelete(bwtSeq->alphabet);
  deleteEISHint(bwtSeq->seqIdx, bwtSeq->hint);
  gt_deleteEncIdxSeq(bwtSeq->seqIdx);
  gt_free(bwtSeq);
}

typedef struct
{
  const Mbtab **mbtab;
  unsigned int numofchars, maxdepth, depth;
  GtCodetype code;
} GtPrebwtstate;

static const Mbtab *gt_prebwt_next(GtPrebwtstate *prebwt,unsigned int cc)
{
  prebwt->code = prebwt->code * prebwt->numofchars + cc;
  prebwt->depth++;
  return prebwt->mbtab[prebwt->depth] + prebwt->code;
}

static inline void
getMatchBound(const BWTSeq *bwtSeq, const Symbol *query, size_t queryLen,
              struct matchBound *match, bool forward)
{
  const Symbol *qptr, *qend;
  unsigned int cc;
  const Mbtab *mbptr;
  GtPrebwtstate prebwt;

  gt_assert(bwtSeq && query);
  if (forward)
  {
    qptr = query;
    qend = query + queryLen;
  } else
  {
    qptr = query + queryLen - 1;
    qend = query - 1;
  }
  gt_assert(ISNOTSPECIAL(*qptr));
  cc = (unsigned int) *qptr;
  prebwt.mbtab = gt_bwtseq2mbtab((const FMindex *) bwtSeq);
  if (prebwt.mbtab != NULL)
  {
    prebwt.numofchars = gt_bwtseq2numofchars((const FMindex *) bwtSeq);
    prebwt.maxdepth = gt_bwtseq2maxdepth((const FMindex *) bwtSeq);
    prebwt.code = 0;
    prebwt.depth = 0;
    mbptr = gt_prebwt_next(&prebwt,cc);
    match->start = mbptr->lowerbound;
    match->end = mbptr->upperbound;
  } else
  {
    prebwt.numofchars = GT_UNDEF_UINT;
    prebwt.maxdepth = GT_UNDEF_UINT;
    prebwt.code = 0;
    prebwt.depth = GT_UNDEF_UINT;
    match->start = bwtSeq->count[cc];
    match->end   = bwtSeq->count[cc + 1];
  }
  qptr = forward ? (qptr+1) : (qptr-1);
  while (match->start < match->end && qptr != qend)
  {
    GtUlongPair occPair;

    gt_assert(ISNOTSPECIAL(*qptr));
    cc = (unsigned int) *qptr;
    if (prebwt.mbtab != NULL && prebwt.depth < prebwt.maxdepth)
    {
      mbptr = gt_prebwt_next(&prebwt,cc);
      match->start = mbptr->lowerbound;
      match->end = mbptr->upperbound;
    } else
    {
      occPair = BWTSeqTransformedPosPairOcc(bwtSeq, (Symbol) cc, match->start,
                                            match->end);
      match->start = bwtSeq->count[cc] + occPair.a;
      match->end   = bwtSeq->count[cc] + occPair.b;
    }
    qptr = forward ? (qptr+1) : (qptr-1);
  }
}

GtUword gt_packedindexuniqueforward(const BWTSeq *bwtSeq,
                                          const GtUchar *qstart,
                                          const GtUchar *qend)
{
  GtUchar cc;
  const GtUchar *qptr;
  struct matchBound bwtbound;
  GtUlongPair seqpospair;
  Symbol curSym;
  const MRAEnc *alphabet;

  gt_assert(bwtSeq && qstart);
  alphabet = BWTSeqGetAlphabet(bwtSeq);
  qptr = qstart;
  cc = *qptr++;
#undef SKDEBUG
#ifdef SKDEBUG
  printf("# start cc=%u\n",cc);
#endif
  if (ISSPECIAL(cc))
  {
    return 0;
  }
  curSym = MRAEncMapSymbol(alphabet, cc);
  bwtbound.start = bwtSeq->count[curSym];
  bwtbound.end = bwtSeq->count[curSym+1];
#ifdef SKDEBUG
  printf("# bounds="GT_WU","GT_WU" = "GT_WU""
          "occurrences\n",
         bwtbound.start,
         bwtbound.end,
         bwtbound.end - bwtbound.start);
#endif
  while (qptr < qend && bwtbound.start + 1 < bwtbound.end)
  {
    cc = *qptr;
#ifdef SKDEBUG
    printf("# cc=%u\n",cc);
#endif
    if (ISSPECIAL (cc))
    {
      return 0;
    }
    curSym = MRAEncMapSymbol(alphabet, cc);
    seqpospair = BWTSeqTransformedPosPairOcc(bwtSeq, curSym,
                                             bwtbound.start,bwtbound.end);
    bwtbound.start = bwtSeq->count[curSym] + seqpospair.a;
    bwtbound.end = bwtSeq->count[curSym] + seqpospair.b;
#ifdef SKDEBUG
    printf("# bounds="GT_WU","GT_WU" = "GT_WU""
            "occurrences\n",
           bwtbound.start,
           bwtbound.end,
           bwtbound.end - bwtbound.start);
#endif
    qptr++;
  }
  if (bwtbound.start + 1 == bwtbound.end)
  {
    return (GtUword) (qptr - qstart);
  }
  return 0;
}

GtUword gt_packedindexmstatsforward(const BWTSeq *bwtSeq,
                                       GtUword *witnessleftbound,
                                       const GtUchar *qstart,
                                       const GtUchar *qend)
{
  GtUchar cc;
  const GtUchar *qptr;
  GtUword prevlbound;
  struct matchBound bwtbound;
  GtUlongPair seqpospair;
  Symbol curSym;
  GtUword matchlength;
  const MRAEnc *alphabet;

  gt_assert(bwtSeq && qstart && qstart < qend);
  alphabet = BWTSeqGetAlphabet(bwtSeq);
  qptr = qstart;
  cc = *qptr;
#undef SKDEBUG
#ifdef SKDEBUG
  printf("# start cc=%u\n",cc);
#endif
  if (ISSPECIAL(cc))
  {
    return 0;
  }
  curSym = MRAEncMapSymbol(alphabet, cc);
  bwtbound.start = bwtSeq->count[curSym];
  bwtbound.end = bwtSeq->count[curSym+1];
  if (bwtbound.start >= bwtbound.end)
  {
    return 0;
  }
#ifdef SKDEBUG
  printf("# bounds="GT_WU","GT_WU" = "GT_WU""
          "occurrences\n",
         bwtbound.start,
         bwtbound.end,
         bwtbound.end - bwtbound.start);
#endif
  prevlbound = bwtbound.start;
  for (qptr++; qptr < qend; qptr++)
  {
    cc = *qptr;
#ifdef SKDEBUG
    printf("# cc=%u\n",cc);
#endif
    if (ISSPECIAL (cc))
    {
      break;
    }
    curSym = MRAEncMapSymbol(alphabet, cc);
    seqpospair = BWTSeqTransformedPosPairOcc(bwtSeq, curSym,
                                             bwtbound.start,bwtbound.end);
    bwtbound.start = bwtSeq->count[curSym] + seqpospair.a;
    bwtbound.end = bwtSeq->count[curSym] + seqpospair.b;
#ifdef SKDEBUG
    printf("# bounds="GT_WU","GT_WU" = "GT_WU""
            "occurrences\n",
           bwtbound.start,
           bwtbound.end,
           bwtbound.end - bwtbound.start);
#endif
    if (bwtbound.start >= bwtbound.end)
    {
      break;
    }
    prevlbound = bwtbound.start;
  }
  matchlength = (GtUword) (qptr - qstart);
  if (witnessleftbound != NULL)
  {
    *witnessleftbound = prevlbound;
  }
  return matchlength;
}

GtUword
gt_BWTSeqMatchCount(const BWTSeq *bwtSeq, const Symbol *query, size_t queryLen,
                 bool forward)
{
  struct matchBound match;
  gt_assert(bwtSeq && query);
  getMatchBound(bwtSeq, query, queryLen, &match, forward);
  if (match.end < match.start)
    return 0;
  else
    return match.end - match.start;
}

bool
gt_initEMIterator(BWTSeqExactMatchesIterator *iter, const BWTSeq *bwtSeq,
               const Symbol *query, size_t queryLen, bool forward)
{
  gt_assert(iter && bwtSeq && query);
  if (!bwtSeq->locateSampleInterval)
  {
    fputs("Index does not contain locate information.\n"
          "Localization of matches impossible!", stderr);
    return false;
  }
  getMatchBound(bwtSeq, query, queryLen, &iter->bounds, forward);
  iter->nextMatchBWTPos = iter->bounds.start;
  initExtBitsRetrieval(&iter->extBits);
  return true;
}

bool
gt_initEmptyEMIterator(BWTSeqExactMatchesIterator *iter, const BWTSeq *bwtSeq)
{
  gt_assert(iter && bwtSeq);
  if (!bwtSeq->locateSampleInterval)
  {
    fputs("Index does not contain locate information.\n"
          "Localization of matches impossible!", stderr);
    return false;
  }
  iter->bounds.start = iter->bounds.end = iter->nextMatchBWTPos = 0;
  initExtBitsRetrieval(&iter->extBits);
  return true;
}

struct BWTSeqExactMatchesIterator *
gt_newEMIterator(const BWTSeq *bwtSeq, const Symbol *query, size_t queryLen,
              bool forward)
{
  struct BWTSeqExactMatchesIterator *iter;
  gt_assert(bwtSeq && query);
  iter = gt_malloc(sizeof (*iter));
  if (gt_initEMIterator(iter, bwtSeq, query, queryLen,forward))
  {
    return iter;
  }
  gt_free(iter);
  return NULL;
}

bool
gt_reinitEMIterator(BWTSeqExactMatchesIterator *iter, const BWTSeq *bwtSeq,
                 const Symbol *query, size_t queryLen, bool forward)
{
  getMatchBound(bwtSeq, query, queryLen, &iter->bounds, forward);
  iter->nextMatchBWTPos = iter->bounds.start;
  return true;
}

void
gt_destructEMIterator(struct BWTSeqExactMatchesIterator *iter)
{
  destructExtBitsRetrieval(&iter->extBits);
}

void
gt_deleteEMIterator(struct BWTSeqExactMatchesIterator *iter)
{
  gt_free(iter);
}

GtUword
gt_EMINumMatchesTotal(const struct BWTSeqExactMatchesIterator *iter)
{
  gt_assert(iter);
  if (iter->bounds.start > iter->bounds.end)
    return 0;
  else
    return iter->bounds.end - iter->bounds.start;
}

GtUword
gt_EMINumMatchesLeft(const struct BWTSeqExactMatchesIterator *iter)
{
  gt_assert(iter);
  if (iter->nextMatchBWTPos > iter->bounds.end)
    return 0;
  else
    return iter->bounds.end - iter->bounds.start;
}

enum
{
  MAX_CONTEXT_LEN = 1000,
  MIN_CONTEXT_LEN = 1,
  CONTEXT_FRACTION = 128,
  MAX_NUM_CONTEXT_CHECKS = 1000,
  CONTEXT_INTERVAL = 128,
};

enum verifyBWTSeqErrCode
gt_BWTSeqVerifyIntegrity(BWTSeq *bwtSeq, const char *projectName,
                      int checkFlags,
                      GtUword tickPrint, FILE *fp,
                      GtLogger *verbosity, GtError *err)
{
  Suffixarray suffixArray;
  struct extBitsRetrieval extBits;
  bool suffixArrayIsInitialized = false, extBitsAreInitialized = false;
  enum verifyBWTSeqErrCode retval = VERIFY_BWTSEQ_NO_ERROR;
  do
  {
    GtUword seqLen;
    gt_assert(bwtSeq && projectName && err);
    gt_error_check(err);

    initExtBitsRetrieval(&extBits);
    extBitsAreInitialized = true;

    if (gt_mapsuffixarray(&suffixArray,
                       SARR_SUFTAB | SARR_ESQTAB, projectName, verbosity, err))
    {
      gt_error_set(err, "Cannot load reference suffix array project with"
                    " demand for suffix table file and encoded sequence"
                    " for project: %s", projectName);
      retval = VERIFY_BWTSEQ_REFLOAD_ERROR;
      break;
    }
    suffixArrayIsInitialized = true;
    seqLen = gt_encseq_total_length(suffixArray.encseq) + 1;
    if (BWTSeqLength(bwtSeq) != seqLen)
    {
      gt_error_set(err, "length mismatch for suffix array project %s and "
                "bwt sequence index", projectName);
      retval = VERIFY_BWTSEQ_LENCOMPARE_ERROR;
      break;
    }

    if (checkFlags & VERIFY_BWTSEQ_SUFVAL
        && BWTSeqHasLocateInformation(bwtSeq))
    {
      GtUword i;
      for (i = 0; i < seqLen && retval == VERIFY_BWTSEQ_NO_ERROR; ++i)
      {
        if (gt_BWTSeqPosHasLocateInfo(bwtSeq, i, &extBits))
        {
          GtUword sfxArrayValue = gt_BWTSeqLocateMatch(bwtSeq, i,
                                                             &extBits);
          if (sfxArrayValue != ESASUFFIXPTRGET(suffixArray.suftab,i))
          {
            gt_error_set(err, "Failed suffix array value comparison"
                          " at position "GT_WU": "GT_WU" != "GT_WU"",
                          i, sfxArrayValue,
                          ESASUFFIXPTRGET(suffixArray.suftab,i));
            retval = VERIFY_BWTSEQ_SUFVAL_ERROR;
            break;
          }
        }
        if (tickPrint && !((i + 1) % tickPrint))
          putc('.', fp);
      }
      if (tickPrint)
        putc('\n', fp);
      if (retval != VERIFY_BWTSEQ_NO_ERROR)
        break;
    }
    else if (checkFlags & VERIFY_BWTSEQ_SUFVAL)
    {
      gt_error_set(err, "check of suffix array values was requested,"
                " but index contains no  locate information!");
      retval = VERIFY_BWTSEQ_SUFVAL_ERROR;
      break;
    }
    else if (!(checkFlags & VERIFY_BWTSEQ_SUFVAL)
             && BWTSeqHasLocateInformation(bwtSeq))
    {
      fputs("Not checking suftab values.\n", stderr);
    }
    if (BWTSeqHasLocateInformation(bwtSeq))
    {
      GtUword nextLocate = BWTSeqTerminatorPos(bwtSeq);
      if (suffixArray.longest.defined &&
          suffixArray.longest.valueunsignedlong != nextLocate)
      {
        gt_error_set(err, "terminator/0-rotation position mismatch "GT_WU""
                  " vs. "GT_WU"", suffixArray.longest.valueunsignedlong,
                  nextLocate);
        retval = VERIFY_BWTSEQ_TERMPOS_ERROR;
        break;
      }
      if ((checkFlags & VERIFY_BWTSEQ_LFMAPWALK)
          && (bwtSeq->featureToggles & BWTReversiblySorted))
      {
        GtUword i = seqLen;
        /* handle first symbol specially because the encseq
         * will not return the terminator symbol */
        {
          Symbol sym = BWTSeqGetSym(bwtSeq, nextLocate);
          if (sym != UNDEFBWTCHAR)
          {
            gt_error_set(err, "symbol mismatch at position "GT_WU": "
                      "%d vs. reference symbol %d", i - 1, (int)sym,
                      (int)UNDEFBWTCHAR);
            retval = VERIFY_BWTSEQ_LFMAPWALK_ERROR;
            break;
          }
          --i;
          nextLocate = BWTSeqLFMap(bwtSeq, nextLocate, &extBits);
        }
        while (i > 0)
        {
          Symbol symRef =
                         gt_encseq_get_encoded_char(suffixArray.encseq,
                                                          --i,
                                                          suffixArray.readmode);
          Symbol symCmp = BWTSeqGetSym(bwtSeq, nextLocate);
          if (symCmp != symRef)
          {
            gt_error_set(err, "symbol mismatch at position "GT_WU": "
                      "%d vs. reference symbol %d", i, symCmp, symRef);
            retval = VERIFY_BWTSEQ_LFMAPWALK_ERROR;
            break;
          }
          nextLocate = BWTSeqLFMap(bwtSeq, nextLocate, &extBits);
        }
        if (retval != VERIFY_BWTSEQ_NO_ERROR)
          break;
      }
      else if ((checkFlags & VERIFY_BWTSEQ_LFMAPWALK)
               && !(bwtSeq->featureToggles & BWTReversiblySorted))
      {
        gt_error_set(err, "requested complete backwards regeneration in index"
                  " without regeneration capability");
        retval = VERIFY_BWTSEQ_LFMAPWALK_IMP_ERROR;
        break;
      }
    }
    if (checkFlags & VERIFY_BWTSEQ_CONTEXT)
    {
      BWTSeqContextRetriever *bwtSeqCR =
        gt_BWTSeqCRLoad(bwtSeq, projectName, CTX_MAP_ILOG_AUTOSIZE);
      if (!bwtSeqCR)
      {
        gt_error_set(err, "cannot load BWT sequence context access table"
                  " for project %s", projectName);
        retval = VERIFY_BWTSEQ_CONTEXT_LOADFAIL;
        break;
      }
      fputs("Checking context regeneration.\n", stderr);
      {
        GtUword i, start, subSeqLen,
          maxSubSeqLen = MIN(MAX(MIN_CONTEXT_LEN, seqLen/CONTEXT_FRACTION),
                             MAX_CONTEXT_LEN),
          numTries = MIN(MAX_NUM_CONTEXT_CHECKS,
                         MAX(2, seqLen/CONTEXT_INTERVAL));
        Symbol *contextBuf = gt_malloc(sizeof (Symbol) * MAX_CONTEXT_LEN);
        GtEncseqReader *esr =
           gt_encseq_create_reader_with_readmode(suffixArray.encseq,
                                                 suffixArray.readmode,
                                                 0);
        for (i = 0; i < numTries && retval == VERIFY_BWTSEQ_NO_ERROR; ++i)
        {
          GtUword j, end, inSubSeqLen;
          subSeqLen = random()%maxSubSeqLen + 1;
          start = random()%(seqLen - subSeqLen + 1);
          end = start + subSeqLen;
          inSubSeqLen = subSeqLen - ((end==seqLen)?1:0);
          gt_BWTSeqCRAccessSubseq(bwtSeqCR, start, subSeqLen, contextBuf);
          gt_encseq_reader_reinit_with_readmode(esr, suffixArray.encseq,
                                                suffixArray.readmode, start);
          for (j = 0; j < inSubSeqLen; ++j)
          {
            Symbol symRef = gt_encseq_reader_next_encoded_char(esr);
            Symbol symCmp = contextBuf[j];
            if (symCmp != symRef)
            {
              gt_error_set(err, "symbol mismatch at position "GT_WU": "
                        "%d vs. reference symbol %d", start + j, (int)symCmp,
                        (int)symRef);
              retval = VERIFY_BWTSEQ_CONTEXT_SYMFAIL;
              break;
            }
          }
          while (j < subSeqLen)
          {
            Symbol symRef = UNDEFBWTCHAR;
            Symbol symCmp = contextBuf[j];
            if (symCmp != symRef)
            {
              gt_error_set(err, "symbol mismatch at position "GT_WU": "
                        "%d vs. reference symbol %d", start + j, (int)symCmp,
                        (int)symRef);
              retval = VERIFY_BWTSEQ_CONTEXT_SYMFAIL;
              break;
            }
            ++j;
          }
        }
        if (retval == VERIFY_BWTSEQ_NO_ERROR)
          fputs("Context regeneration completed successfully.\n", stderr);
        gt_encseq_reader_delete(esr);
        gt_free(contextBuf);
      }
      gt_deleteBWTSeqCR(bwtSeqCR);
    }
  } while (0);
  if (suffixArrayIsInitialized) gt_freesuffixarray(&suffixArray);
  if (extBitsAreInitialized) destructExtBitsRetrieval(&extBits);
  return retval;
}
