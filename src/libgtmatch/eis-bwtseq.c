/*
  Copyright (c) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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

#include <assert.h>
#include <string.h>

#include "libgtcore/dataalign.h"
#include "libgtcore/error.h"
#include "libgtcore/str.h"
#include "libgtcore/unused.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/esa-map.pr"
#include "libgtmatch/seqpos-def.h"

#include "libgtmatch/eis-bitpackseqpos.h"
#include "libgtmatch/eis-bwtseq.h"
#include "libgtmatch/eis-bwtseqconstruct.h"
#include "libgtmatch/eis-bwtconstruct_params.h"
#include "libgtmatch/eis-bwtseqpriv.h"
#include "libgtmatch/eis-bwtseqcreate.h"
#include "libgtmatch/eis-encidxseq.h"
#include "libgtmatch/eis-encidxseqconstruct.h"
#include "libgtmatch/eis-mrangealphabet.h"
#include "libgtmatch/eis-suffixerator-interface.h"
#include "libgtmatch/eis-suffixarray-interface.h"

static int
initBWTSeqFromEncSeqIdx(struct BWTSeq *bwtSeq, struct encIdxSeq *baseSeqIdx,
                        MRAEnc *alphabet, Seqpos *counts);

static BWTSeq *
newBWTSeq(struct encIdxSeq *seqIdx, MRAEnc *alphabet);

extern BWTSeq *
availBWTSeq(const struct bwtParam *params, Error *err)
{
  struct BWTSeq *bwtSeq = NULL;
  Suffixarray suffixArray;
  Seqpos len;
  Verboseinfo *verbosity;
  assert(params && err);
  error_check(err);
  /* FIXME: handle verbosity in a more sane fashion */
  verbosity = newverboseinfo(false);
  if (streamsuffixarray(&suffixArray, &len, SARR_SUFTAB | SARR_BWTTAB,
                        params->projectName, verbosity, err))
  {
    freeverboseinfo(&verbosity);
    return NULL;
  }
  ++len;
  bwtSeq = availBWTSeqFromSA(params, &suffixArray, len, err);
  freesuffixarray(&suffixArray);
  freeverboseinfo(&verbosity);
  return bwtSeq;
}

extern BWTSeq *
trSuftab2BWTSeq(const struct bwtParam *params, Error *err)
{
  struct BWTSeq *bwtSeq = NULL;
  Suffixarray suffixArray;
  Seqpos len;
  Verboseinfo *verbosity;
  assert(params && err);
  error_check(err);
  /* FIXME: handle verbosity in a more sane fashion */
  verbosity = newverboseinfo(false);
  if (streamsuffixarray(&suffixArray, &len,
                        SARR_SUFTAB | SARR_BWTTAB | SARR_ESQTAB,
                        params->projectName, verbosity, err)
      && streamsuffixarray(&suffixArray, &len, SARR_SUFTAB | SARR_ESQTAB,
                           params->projectName, verbosity, err))
  {
    freeverboseinfo(&verbosity);
    return NULL;
  }
  ++len;
  bwtSeq = createBWTSeqFromSA(params, &suffixArray, len, err);
  freesuffixarray(&suffixArray);
  freeverboseinfo(&verbosity);
  return bwtSeq;
}

extern BWTSeq *
availBWTSeqFromSA(const struct bwtParam *params, Suffixarray *sa,
                  Seqpos totalLen, Error *err)
{
  BWTSeq *bwtSeq;
  assert(sa && params && err);
  error_check(err);
  /* try loading index */
  bwtSeq = loadBWTSeqForSA(params->projectName, params->baseType,
                           params->seqParams.blockEnc.EISFeatureSet,
                           sa, totalLen, err);
  /* if loading didn't work try on-demand creation */
  if (!bwtSeq)
  {
    error_unset(err);
    bwtSeq = createBWTSeqFromSA(params, sa, totalLen, err);
  }
  return bwtSeq;
}

static int GTAlphabetRangeHandling[] = { NORMAL_RANGE, SPECIAL_RANGE };

extern BWTSeq *
loadBWTSeq(const Str *projectName, int BWTOptFlags, Error *err)
{
  struct BWTSeq *bwtSeq = NULL;
  Suffixarray suffixArray;
  Seqpos len;
  Verboseinfo *verbosity;
  assert(projectName && err);
  error_check(err);
  /* FIXME: handle verbosity in a more sane fashion */
  verbosity = newverboseinfo(false);
  if (mapsuffixarray(&suffixArray, &len, 0, projectName, verbosity, err))
  {
    freeverboseinfo(&verbosity);
    return NULL;
  }
  ++len;
  bwtSeq = loadBWTSeqForSA(projectName, BWT_ON_BLOCK_ENC, BWTOptFlags,
                           &suffixArray, len, err);
  freesuffixarray(&suffixArray);
  freeverboseinfo(&verbosity);
  return bwtSeq;
}

extern BWTSeq *
loadBWTSeqForSA(const Str *projectName, enum seqBaseEncoding baseType,
                int BWTOptFlags, const Suffixarray *sa,
                Seqpos totalLen, Error *err)
{
  struct BWTSeq *bwtSeq = NULL;
  EISeq *seqIdx = NULL;
  MRAEnc *alphabet = NULL;
  assert(projectName && sa && err);
  alphabet = newMRAEncFromSA(sa);
  switch (baseType)
  {
  case BWT_ON_BLOCK_ENC:
    if ((seqIdx = loadBlockEncIdxSeqForSA(
           sa, totalLen, projectName,
           convertBWTOptFlags2EISFeatures(BWTOptFlags), err)))
    {
      if (!(bwtSeq = newBWTSeq(seqIdx, alphabet)))
        break;
      fputs("Using pre-computed sequence index.\n", stderr);
    }
    break;
  default:
    error_set(err, "Illegal/unknown/unimplemented encoding requested!");
    break;
  }
  if (!bwtSeq && seqIdx)
    deleteEncIdxSeq(seqIdx);
  return bwtSeq;
}

extern BWTSeq *
createBWTSeqFromSA(const struct bwtParam *params, Suffixarray *sa,
                   Seqpos totalLen, Error *err)
{
  BWTSeq *bwtSeq = NULL;
  if (!sa->longest.defined)
  {
    fprintf(stderr,
            "error: position of null-rotation/longest suffix not available"
            " for project %s\n", str_get(params->projectName));
  }
  else
  {
    struct suffixarrayFileInterface sai;
    initSuffixarrayFileInterface(&sai, sa);
    bwtSeq = createBWTSeqFromSAI(params, &sai, totalLen, err);
    destructSuffixarrayFileInterface(&sai);
  }
  return bwtSeq;
}

extern BWTSeq *
createBWTSeqFromSAI(const struct bwtParam *params,
                    struct suffixarrayFileInterface *sai,
                    Seqpos totalLen, Error *err)
{
  BWTSeq *bwtSeq = NULL;
  MRAEnc *alphabet = NULL;
  SeqDataReader readSfxIdx = { NULL, NULL };
  assert(sai && err && params);
  alphabet = newMRAEncFromSAI(sai);
  if (params->locateInterval
      && !SDRIsValid(readSfxIdx
                               = SAIMakeReader(sai, SFX_REQUEST_SUFTAB)))
  {
    fprintf(stderr, "error: locate sampling requested but not available"
            " for project %s\n", str_get(params->projectName));
  }
  else
  {
    EISeq *seqIdx = NULL;
    switch (params->baseType)
    {
    case BWT_ON_BLOCK_ENC:
      seqIdx =
        createBWTSeqGeneric(
          params, (indexCreateFunc)newBlockEncIdxSeqFromSAI, sai, totalLen,
          alphabet, GTAlphabetRangeHandling, SAIGetOrigSeqSym, sai, readSfxIdx,
          reportSAILongest, sai, err);
      break;
    default:
      error_set(err, "Illegal/unknown/unimplemented encoding requested!");
      break;
    }
    if (seqIdx)
      bwtSeq = newBWTSeq(seqIdx, alphabet);
    if (!bwtSeq)
    {
      if (seqIdx)
        deleteEncIdxSeq(seqIdx);
      if (alphabet)
        MRAEncDelete(alphabet);
    }
  }
  return bwtSeq;
}

extern BWTSeq *
createBWTSeqFromSfxI(const struct bwtParam *params, sfxInterface *sfxi,
                     Seqpos totalLen, Error *err)
{
  EISeq *seqIdx = NULL;
  BWTSeq *bwtSeq = NULL;
  MRAEnc *alphabet = NULL;
  SeqDataReader readSfxIdx = { NULL, NULL };
  assert(sfxi && params && err);
  if (params->locateInterval)
  {
    if (!SDRIsValid(readSfxIdx
                              = SfxIRegisterReader(sfxi,
                                                   SFX_REQUEST_SUFTAB)))
      return NULL;
  }
  alphabet = newMRAEncFromSfxI(sfxi);
  seqIdx= createBWTSeqGeneric(
    params, (indexCreateFunc)newBlockEncIdxSeqFromSfxI, sfxi, totalLen,
    alphabet, GTAlphabetRangeHandling,
    SfxIGetOrigSeq, sfxi, readSfxIdx,
    (reportLongest)getSfxILongestPos, sfxi, err);
  if (seqIdx)
  {
    bwtSeq = newBWTSeq(seqIdx, alphabet);
  }
  if (!bwtSeq && seqIdx)
    deleteEncIdxSeq(seqIdx);
  return bwtSeq;
}

/**
 * @param alphabet ownership of alphabet is with the newly produced
 * sequence object if return value is not 0
 */
static int
initBWTSeqFromEncSeqIdx(BWTSeq *bwtSeq, struct encIdxSeq *seqIdx,
                        MRAEnc *alphabet, Seqpos *counts)
{
  size_t alphabetSize;
  Symbol bwtTerminatorFlat;
  EISHint hint;
  assert(bwtSeq && seqIdx);
  bwtSeq->alphabet = alphabet;
  alphabetSize = MRAEncGetSize(alphabet);
  if (!alphabetSize)
    /* weird error, shouldn't happen, but I prefer error return to
     * segfault in case someone tampered with the input */
    return 0;
  /* FIXME: this should probably be handled in chardef.h to have a
   * unique mapping */
  /* FIXME: this assumes there is exactly two ranges */
  MRAEncAddSymbolToRange(alphabet, bwtTerminatorSym, 1);
  assert(MRAEncGetSize(alphabet) ==  alphabetSize + 1);
  alphabetSize = MRAEncGetSize(alphabet);
  bwtSeq->bwtTerminatorFallback = bwtTerminatorFlat =
    MRAEncMapSymbol(alphabet, UNDEFBWTCHAR);
  bwtSeq->bwtTerminatorFallbackRange = 1;
  bwtSeq->count = counts;
  bwtSeq->seqIdx = seqIdx;
  bwtSeq->alphabetSize = alphabetSize;
  {
    struct locateHeader header;
    if (!readLocateInfoHeader(seqIdx, &header)
        || !header.locateInterval)
    {
      fputs("Index does not contain locate information.\n"
            "Localization of matches will not be supported!\n", stderr);
      bwtSeq->locateSampleInterval = 0;
    }
    else
    {
      bwtSeq->locateSampleInterval = header.locateInterval;
      bwtSeq->longest = header.longest;
      /* FIXME: this really deserves its own header */
      bwtSeq->featureToggles = header.featureToggles;
    }
  }
  bwtSeq->hint = hint = newEISHint(seqIdx);
  {
    Symbol i;
    Seqpos len = EISLength(seqIdx), *count = bwtSeq->count;
    count[0] = 0;
    for (i = 0; i < bwtTerminatorFlat; ++i)
      count[i + 1] = count[i]
        + EISSymTransformedRank(seqIdx, i, len, hint);
    /* handle character which the terminator has been mapped to specially */
    count[i + 1] = count[i]
      + EISSymTransformedRank(seqIdx, i, len, hint) - 1;
    assert(count[i + 1] >= count[i]);
    /* now we can finish the rest of the symbols */
    for (i += 2; i < alphabetSize; ++i)
      count[i] = count[i - 1]
        + EISSymTransformedRank(seqIdx, i - 1, len, hint);
    /* and finally place the 1-count for the terminator */
    count[i] = count[i - 1] + 1;
#ifdef DEBUG
    fprintf(stderr, "count[alphabetSize]="FormatSeqpos
            ", len="FormatSeqpos"\n", count[alphabetSize], len);
    for (i = 0; i <= alphabetSize; ++i)
      fprintf(stderr, "count[%u]="FormatSeqpos"\n", (unsigned)i, count[i]);
#endif
    assert(count[alphabetSize] == len);
  }
  return 1;
}

/**
 * @param alphabet ownership of alphabet is with the newly produced
 * sequence object if return value is non-NULL
 */
static BWTSeq *
newBWTSeq(EISeq *seqIdx, MRAEnc *alphabet)
{
  BWTSeq *bwtSeq;
  Seqpos *counts;
  unsigned alphabetSize;
  assert(seqIdx);
  /* alphabetSize is increased by one to handle the flattened
   * terminator symbol correctly */
  alphabetSize = MRAEncGetSize(alphabet) + 1;
  bwtSeq = ma_malloc(offsetAlign(sizeof (struct BWTSeq), sizeof (Seqpos))
                     + sizeof (Seqpos) * (alphabetSize + 1));
  counts = (Seqpos *)((char  *)bwtSeq
                      + offsetAlign(sizeof (struct BWTSeq),
                                    sizeof (Seqpos)));
  if (!initBWTSeqFromEncSeqIdx(bwtSeq, seqIdx, alphabet, counts))
  {
    ma_free(bwtSeq);
    bwtSeq = NULL;
  }
  return bwtSeq;
}

void
deleteBWTSeq(BWTSeq *bwtSeq)
{
  MRAEncDelete(bwtSeq->alphabet);
  deleteEISHint(bwtSeq->seqIdx, bwtSeq->hint);
  deleteEncIdxSeq(bwtSeq->seqIdx);
  ma_free(bwtSeq);
}

static inline void
getMatchBound(const BWTSeq *bwtSeq, const Symbol *query, size_t queryLen,
              struct matchBound *match)
{
  size_t i = queryLen;
  const Seqpos *count;
  Symbol curSym;
  const MRAEnc *alphabet;
  assert(bwtSeq && query);
  count = bwtSeq->count;
  alphabet = BWTSeqGetAlphabet(bwtSeq);
  curSym = MRAEncMapSymbol(alphabet, query[--i]);
  match->start = count[curSym];
  match->end   = count[curSym + 1];
  while ((match->start <= match->end) && (i > 0))
  {
    struct SeqposPair occPair;
    curSym = MRAEncMapSymbol(alphabet, query[--i]);
    occPair = BWTSeqTransformedPosPairOcc(bwtSeq, curSym, match->start,
                                          match->end);
    match->start = count[curSym] + occPair.a;
    match->end   = count[curSym] + occPair.b;
/*     match->start = count[curSym] + BWTSeqTransformedOcc(bwtSeq, curSym,
 *     match->start); */
/*     match->end = count[curSym] + BWTSeqTransformedOcc(bwtSeq, curSym,
 *     match->end); */
  }
}

unsigned long packedindexuniqueforward(const void *genericindex,
                                       UNUSED unsigned long offset,
                                       UNUSED Seqpos left,
                                       UNUSED Seqpos right,
                                       UNUSED Seqpos *witnessposition,
                                       const Uchar *qstart,
                                       const Uchar *qend)
{
  Uchar cc;
  const Uchar *qptr;
  struct matchBound bwtbound;
  const BWTSeq *bwtSeq = (BWTSeq *) genericindex;
  Symbol curSym;
  const MRAEnc *alphabet;

  assert(bwtSeq && qstart);
  alphabet = BWTSeqGetAlphabet(bwtSeq);
  qptr = qstart;
  cc = *qptr++;
#undef mydebug
#ifdef mydebug
  printf("# start cc=%u\n",cc);
#endif
  if (ISSPECIAL(cc))
  {
    return 0;
  }
  curSym = MRAEncMapSymbol(alphabet, cc);
  bwtbound.end = bwtSeq->count[curSym];
  bwtbound.start = bwtSeq->count[curSym+1];
#ifdef mydebug
  printf("# bounds=" FormatSeqpos "," FormatSeqpos " = " FormatSeqos
          "occurrences\n",
         PRINTSeqposcast(bwtbound.end),
         PRINTSeqposcast(bwtbound.start),
         PRINTSeqposcast(bwtbound.start - bwtbound.end));
#endif
  while (qptr < qend && bwtbound.end + 1 < bwtbound.start)
  {
    cc = *qptr;
#ifdef mydebug
    printf("# cc=%u\n",cc);
#endif
    if (ISSPECIAL (cc))
    {
      return 0;
    }
    curSym = MRAEncMapSymbol(alphabet, cc);

    bwtbound.end = bwtSeq->count[curSym] +
                     BWTSeqOcc(bwtSeq, curSym, bwtbound.end);
    bwtbound.start = bwtSeq->count[curSym] +
                     BWTSeqOcc(bwtSeq, curSym, bwtbound.start);
    /*
      aber vorher keine Transformation;
      BWTSeqPosPairOcc(const BWTSeq *bwtSeq, Symbol sym,
                       Seqpos posA, Seqpos posB);
    */
#ifdef mydebug
    printf("# bounds=" FormatSeqpos "," FormatSeqpos " = " FormatSeqos
            "occurrences\n",
           PRINTSeqposcast(bwtbound.end),
           PRINTSeqposcast(bwtbound.start),
           PRINTSeqposcast(bwtbound.start - bwtbound.end));
#endif
    qptr++;
  }
  if (bwtbound.end + 1 == bwtbound.start)
  {
    return (unsigned long) (qptr - qstart);
  }
  return 0;
}

unsigned long packedindexmstatsforward(const void *genericindex,
                                       UNUSED unsigned long offset,
                                       UNUSED Seqpos left,
                                       UNUSED Seqpos right,
                                       Seqpos *witnessposition,
                                       const Uchar *qstart,
                                       const Uchar *qend)
{
  Uchar cc;
  const Uchar *qptr;
  Seqpos prevlbound;
  struct matchBound bwtbound;
  const BWTSeq *bwtSeq = (BWTSeq *) genericindex;
  Symbol curSym;
  unsigned long matchlength;
  const MRAEnc *alphabet;

  assert(bwtSeq && qstart && qstart < qend);
  alphabet = BWTSeqGetAlphabet(bwtSeq);
  qptr = qstart;
  cc = *qptr;
#undef mydebug
#ifdef mydebug
  printf("# start cc=%u\n",cc);
#endif
  if (ISSPECIAL(cc))
  {
    return 0;
  }
  curSym = MRAEncMapSymbol(alphabet, cc);
  bwtbound.end = bwtSeq->count[curSym];
  bwtbound.start = bwtSeq->count[curSym+1];
  if (bwtbound.end >= bwtbound.start)
  {
    return 0;
  }
#ifdef mydebug
  printf("# bounds=" FormatSeqpos "," FormatSeqpos " = " FormatSeqos
          "occurrences\n",
         PRINTSeqposcast(bwtbound.end),
         PRINTSeqposcast(bwtbound.start),
         PRINTSeqposcast(bwtbound.start - bwtbound.end));
#endif
  prevlbound = bwtbound.end;
  for (qptr++; qptr < qend; qptr++)
  {
    cc = *qptr;
#ifdef mydebug
    printf("# cc=%u\n",cc);
#endif
    if (ISSPECIAL (cc))
    {
      break;
    }
    curSym = MRAEncMapSymbol(alphabet, cc);
    bwtbound.end = bwtSeq->count[curSym] +
                     BWTSeqOcc(bwtSeq, curSym, bwtbound.end);
    bwtbound.start = bwtSeq->count[curSym] +
                     BWTSeqOcc(bwtSeq, curSym, bwtbound.start);
#ifdef mydebug
    printf("# bounds=" FormatSeqpos "," FormatSeqpos " = " FormatSeqos
            "occurrences\n",
           PRINTSeqposcast(bwtbound.end),
           PRINTSeqposcast(bwtbound.start),
           PRINTSeqposcast(bwtbound.start - bwtbound.end));
#endif
    if (bwtbound.end >= bwtbound.start)
    {
      break;
    }
    prevlbound = bwtbound.end;
  }
  matchlength = (unsigned long) (qptr - qstart);
  if (witnessposition != NULL)
  {
    Seqpos startpos = pckfindfirstmatch(bwtSeq,prevlbound);
    assert((bwtSeq->seqIdx->seqLen-1) >= (startpos + matchlength));
    *witnessposition = (bwtSeq->seqIdx->seqLen - 1) - (startpos + matchlength);
  }
  return matchlength;
}

extern Seqpos
BWTSeqMatchCount(const BWTSeq *bwtSeq, const Symbol *query, size_t queryLen)
{
  struct matchBound match;
  assert(bwtSeq && query);
  getMatchBound(bwtSeq, query, queryLen, &match);
  if (match.end < match.start)
    return 0;
  else
    return match.end - match.start;
}

extern bool
initEMIterator(BWTSeqExactMatchesIterator *iter, const BWTSeq *bwtSeq,
               const Symbol *query, size_t queryLen)
{
  assert(iter && bwtSeq && query);
  if (!bwtSeq->locateSampleInterval)
  {
    fputs("Index does not contain locate information.\n"
          "Localization of matches impossible!", stderr);
    return false;
  }
  getMatchBound(bwtSeq, query, queryLen, &iter->bounds);
  iter->nextMatchBWTPos = iter->bounds.start;
  initExtBitsRetrieval(&iter->extBits);
  return true;
}

extern bool
initEmptyEMIterator(BWTSeqExactMatchesIterator *iter, const BWTSeq *bwtSeq)
{
  assert(iter && bwtSeq);
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
newEMIterator(const BWTSeq *bwtSeq, const Symbol *query, size_t queryLen)
{
  struct BWTSeqExactMatchesIterator *iter;
  assert(bwtSeq && query);
  iter = ma_malloc(sizeof (*iter));
  if (initEMIterator(iter, bwtSeq, query, queryLen))
    return iter;
  else
  {
    ma_free(iter);
    return NULL;
  }
}

extern bool
reinitEMIterator(BWTSeqExactMatchesIterator *iter, const BWTSeq *bwtSeq,
                 const Symbol *query, size_t queryLen)
{
  getMatchBound(bwtSeq, query, queryLen, &iter->bounds);
  iter->nextMatchBWTPos = iter->bounds.start;
  return true;
}

extern void
destructEMIterator(struct BWTSeqExactMatchesIterator *iter)
{
  destructExtBitsRetrieval(&iter->extBits);
}

void
deleteEMIterator(struct BWTSeqExactMatchesIterator *iter)
{
  ma_free(iter);
}

Seqpos
EMINumMatchesTotal(const struct BWTSeqExactMatchesIterator *iter)
{
  assert(iter);
  if (iter->bounds.start > iter->bounds.end)
    return 0;
  else
    return iter->bounds.end - iter->bounds.start;
}

extern Seqpos
EMINumMatchesLeft(const struct BWTSeqExactMatchesIterator *iter)
{
  assert(iter);
  if (iter->nextMatchBWTPos > iter->bounds.end)
    return 0;
  else
    return iter->bounds.end - iter->bounds.start;
}

extern int
BWTSeqVerifyIntegrity(BWTSeq *bwtSeq, const Str *projectName,
                      unsigned long tickPrint, FILE *fp,
                      Error *err)
{
  Suffixarray suffixArray;
  struct extBitsRetrieval extBits;
  bool suffixArrayIsInitialized = false, extBitsAreInitialized = false;
  Verboseinfo *verbosity = NULL;
  enum verifyBWTSeqErrCode retval = VERIFY_BWTSEQ_NO_ERROR;
  do
  {
    Seqpos len;
    assert(bwtSeq && projectName && err);
    error_check(err);

    verbosity = newverboseinfo(true);
    initExtBitsRetrieval(&extBits);
    if (mapsuffixarray(&suffixArray, &len,
                       SARR_SUFTAB | SARR_ESQTAB, projectName, verbosity, err))
    {
      error_set(err, "Cannot load reference suffix array project with"
                    " demand for suffix table file and encoded sequence"
                    " for project: %s", str_get(projectName));
      freeverboseinfo(&verbosity);
      retval = VERIFY_BWTSEQ_REFLOAD_ERROR;
      break;
    }
    suffixArrayIsInitialized = true;
    ++len;
    if (BWTSeqLength(bwtSeq) != len)
    {
      error_set(err, "length mismatch for suffix array project %s and "
                "bwt sequence index", str_get(projectName));
      retval = VERIFY_BWTSEQ_LENCOMPARE_ERROR;
      break;
    }

    if (BWTSeqHasLocateInformation(bwtSeq))
    {
      Seqpos i;
      for (i = 0; i < len && retval == VERIFY_BWTSEQ_NO_ERROR; ++i)
      {
        if (BWTSeqPosHasLocateInfo(bwtSeq, i, &extBits))
        {
          Seqpos sfxArrayValue = BWTSeqLocateMatch(bwtSeq, i, &extBits);
          if (sfxArrayValue != suffixArray.suftab[i])
          {
            error_set(err, "Failed suffixarray value comparison"
                          " at position "FormatSeqpos": "FormatSeqpos" != "
                          FormatSeqpos,
                          i, sfxArrayValue, suffixArray.suftab[i]);
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
    else
    {
      fputs("Not checking suftab values (no locate information present)!\n",
            stderr);
    }
    if ((bwtSeq->featureToggles & BWTProperlySorted)
        && suffixArray.longest.defined && len)
    {
      Seqpos nextLocate = suffixArray.longest.valueseqpos,
        i = len;
      /* handle first symbol specially because the encodedsequence
       * will not return the terminator symbol */
      {
        Symbol sym = EISGetSym(bwtSeq->seqIdx, nextLocate, bwtSeq->hint);
        if (sym != UNDEFBWTCHAR)
        {
          error_set(err, "symbol mismatch at position "FormatSeqpos": "
                        "%d vs. reference symbol %d", i - 1, sym,
                        UNDEFBWTCHAR);
          retval = VERIFY_BWTSEQ_LFMAPWALK_ERROR;
          break;
        }
        --i;
        nextLocate = BWTSeqLFMap(bwtSeq, nextLocate);
      }
      while (i > 0)
      {
        Symbol symRef = getencodedchar(suffixArray.encseq,
                                       --i, suffixArray.readmode);
        Symbol symCmp = EISGetSym(bwtSeq->seqIdx, nextLocate, bwtSeq->hint);
        if (symCmp != symRef)
        {
          error_set(err, "symbol mismatch at position "FormatSeqpos": "
                        "%d vs. reference symbol %d", i, symCmp, symRef);
          retval = VERIFY_BWTSEQ_LFMAPWALK_ERROR;
          break;
        }
        nextLocate = BWTSeqLFMap(bwtSeq, nextLocate);
      }
    }
  } while (0);
  if (suffixArrayIsInitialized)
    freesuffixarray(&suffixArray);
  if (verbosity)
    freeverboseinfo(&verbosity);
  if (extBitsAreInitialized)
    destructExtBitsRetrieval(&extBits);
  return retval;
}
