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
#include "libgtcore/env.h"
#include "libgtcore/str.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/esa-map.pr"
#include "libgtmatch/seqpos-def.h"

#include "libgtmatch/eis-bitpackseqpos.h"
#include "libgtmatch/eis-bwtseq.h"
#include "libgtmatch/eis-bwtseqconstruct.h"
#include "libgtmatch/eis-bwtseqpriv.h"
#include "libgtmatch/eis-bwtseqcreate.h"
#include "libgtmatch/eis-encidxseq.h"
#include "libgtmatch/eis-encidxseqconstruct.h"
#include "libgtmatch/eis-mrangealphabet.h"
#include "libgtmatch/eis-suffixerator-interface.h"
#include "libgtmatch/eis-suffixarray-interface.h"

static int
initBWTSeqFromEncSeqIdx(struct BWTSeq *bwtSeq, struct encIdxSeq *baseSeqIdx,
                        MRAEnc *alphabet, Seqpos longest, Seqpos *counts,
                        Env *env);

static BWTSeq *
newBWTSeq(struct encIdxSeq *seqIdx, MRAEnc *alphabet, Seqpos longest, Env *env);

extern BWTSeq *
availBWTSeq(const struct bwtParam *params, Env *env)
{
  struct BWTSeq *bwtSeq = NULL;
  Suffixarray suffixArray;
  Seqpos len;
  Verboseinfo *verbosity;
  assert(params && env);
  env_error_check(env);
  /* FIXME: handle verbosity in a more sane fashion */
  verbosity = newverboseinfo(false, env);
  if (streamsuffixarray(&suffixArray, &len, SARR_SUFTAB | SARR_BWTTAB,
                        params->projectName, verbosity, env))
  {
    freeverboseinfo(&verbosity, env);
    return NULL;
  }
  ++len;
  bwtSeq = availBWTSeqFromSA(params, &suffixArray, len, env);
  freesuffixarray(&suffixArray, env);
  freeverboseinfo(&verbosity, env);
  return bwtSeq;
}

extern BWTSeq *
availBWTSeqFromSA(const struct bwtParam *params, Suffixarray *sa,
                  Seqpos totalLen, Env *env)
{
  BWTSeq *bwtSeq;
  assert(sa && params && env);
  env_error_check(env);
  /* try loading index */
  bwtSeq = loadBWTSeqForSA(params, sa, totalLen, env);
  /* if loading didn't work try on-demand creation */
  if (!bwtSeq)
  {
    env_error_unset(env);
    bwtSeq = createBWTSeqFromSA(params, sa, totalLen, env);
  }
  return bwtSeq;
}

static int GTAlphabetRangeHandling[] = { NORMAL_RANGE, SPECIAL_RANGE };

extern BWTSeq *
loadBWTSeqForSA(const struct bwtParam *params, Suffixarray *sa,
                Seqpos totalLen, Env *env)
{
  struct BWTSeq *bwtSeq = NULL;
  EISeq *seqIdx = NULL;
  MRAEnc *alphabet = NULL;
  assert(params && sa && env);
  if (!sa->longest.defined)
    return NULL;
  alphabet = newMRAEncFromSA(sa, env);
  switch (params->baseType)
  {
  case BWT_ON_BLOCK_ENC:
    if ((seqIdx = loadBlockEncIdxSeqForSA(
           sa, totalLen, params->projectName,
           params->seqParams.blockEnc.EISFeatureSet, env)))
    {
      if (!(bwtSeq = newBWTSeq(seqIdx, alphabet, sa->longest.valueseqpos, env)))
        break;
      fputs("Using pre-computed sequence index.\n", stderr);
    }
    break;
  default:
    env_error_set(env, "Illegal/unknown/unimplemented encoding requested!");
    break;
  }
  if (!bwtSeq && seqIdx)
    deleteEncIdxSeq(seqIdx, env);
  return bwtSeq;
}

extern BWTSeq *
createBWTSeqFromSA(const struct bwtParam *params, Suffixarray *sa,
                   Seqpos totalLen, Env *env)
{
  BWTSeq *bwtSeq = NULL;
  MRAEnc *alphabet = NULL;
  if (params->locateInterval &&
      (!sa->suftabstream.fp || !sa->longest.defined))
  {
    fprintf(stderr, "error: locate sampling requested but not available"
            " for project %s\n", str_get(params->projectName));
  }
  else
  {
    EISeq *seqIdx = NULL;
    alphabet = newMRAEncFromSA(sa, env);
    switch (params->baseType)
    {
    case BWT_ON_BLOCK_ENC:
      seqIdx =
        createBWTSeqGeneric(
          params, (indexCreateFunc)newBlockEncIdxSeqFromSA, sa, totalLen,
          alphabet, GTAlphabetRangeHandling, saGetOrigSeqSym, sa, saReadSeqpos,
          sa, reportSALongest, sa, env);
      break;
    default:
      env_error_set(env, "Illegal/unknown/unimplemented encoding requested!");
      break;
    }
    if (seqIdx)
      bwtSeq = newBWTSeq(seqIdx, alphabet, sa->longest.valueseqpos, env);
    if (!bwtSeq)
    {
      if (seqIdx)
        deleteEncIdxSeq(seqIdx, env);
      if (alphabet)
        MRAEncDelete(alphabet, env);
    }
  }
  return bwtSeq;
}

struct sfxIReadInfo
{
  sfxInterface *si;
  listenerID id;
};

static int
sfxIReadSeqpos(void *src, Seqpos *dest, size_t len, Env *env)
{
  return readSfxISufTabRange(((struct sfxIReadInfo *)src)->si,
                             ((struct sfxIReadInfo *)src)->id,
                             len, dest, env) == len;
}

#if 0
static int
sfxIReadBWTSym(void *src, Symbol *dest, size_t len, Env *env)
{
  return readSfxIBWTRange(((struct sfxIReadInfo *)src)->si,
                          ((struct sfxIReadInfo *)src)->id,
                          len, dest, env) == len;
}
#endif

extern BWTSeq *
createBWTSeqFromSfxI(const struct bwtParam *params, sfxInterface *si,
                     Seqpos totalLen, Env *env)
{
  struct sfxIReadInfo siriSeqpos;
  EISeq *seqIdx = NULL;
  BWTSeq *bwtSeq = NULL;
  MRAEnc *alphabet = NULL;
  assert(si && params && env);

  if (params->locateInterval)
  {
    siriSeqpos.si = si;
    if (!SfxIRegisterReader(si, &siriSeqpos.id, SFX_REQUEST_SUFTAB, env))
      return NULL;
  }
  alphabet = newMRAEncFromSfxI(si, env);
  seqIdx= createBWTSeqGeneric(
    params, (indexCreateFunc)newBlockEncIdxSeqFromSfxI, si, totalLen,
    alphabet, GTAlphabetRangeHandling,
    SfxIGetOrigSeq, si, sfxIReadSeqpos, &siriSeqpos,
    (reportLongest)getSfxILongestPos, si, env);
  if (seqIdx)
  {
    DefinedSeqpos longest = getSfxILongestPos(si);
    if (longest.defined)
      bwtSeq = newBWTSeq(seqIdx, alphabet, longest.valueseqpos, env);
    else
      env_error_set(env, "Position of terminator in BWT not found!");
  }
  if (!bwtSeq && seqIdx)
    deleteEncIdxSeq(seqIdx, env);
  return bwtSeq;
}

/**
 * @param alphabet ownership of alphabet is with the newly produced
 * sequence object if return value is not 0
 */
static int
initBWTSeqFromEncSeqIdx(BWTSeq *bwtSeq, struct encIdxSeq *seqIdx,
                        MRAEnc *alphabet, Seqpos longest, Seqpos *counts,
                        Env *env)
{
  size_t alphabetSize;
  Symbol bwtTerminatorFlat;
  EISHint hint;
  assert(bwtSeq && seqIdx && env);
  bwtSeq->alphabet = alphabet;
  alphabetSize = MRAEncGetSize(alphabet);
  if (!alphabetSize)
    /* weird error, shouldn't happen, but I prefer error return to
     * segfault in case someone tampered with the input */
    return 0;
  /* FIXME: this should probably be handled in chardef.h to have a
   * unique mapping */
  MRAEncAddSymbolToRange(alphabet, SEPARATOR - 3, 1);
  assert(MRAEncGetSize(alphabet) ==  alphabetSize + 1);
  alphabetSize = MRAEncGetSize(alphabet);
  bwtSeq->bwtTerminatorFallback = bwtTerminatorFlat =
    MRAEncMapSymbol(alphabet, UNDEFBWTCHAR);

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
  bwtSeq->hint = hint = newEISHint(seqIdx, env);
  {
    Symbol i;
    Seqpos len = EISLength(seqIdx), *count = bwtSeq->count;
    count[0] = 0;
    for (i = 0; i < bwtTerminatorFlat; ++i)
      count[i + 1] = count[i]
        + EISSymTransformedRank(seqIdx, i, len, hint, env);
    /* handle character which the terminator has been mapped to specially */
    count[i + 1] = count[i]
      + EISSymTransformedRank(seqIdx, i, len, hint, env) - 1;
    assert(count[i + 1] >= count[i]);
    /* now we can finish the rest of the symbols */
    for (i += 2; i < alphabetSize; ++i)
      count[i] = count[i - 1]
        + EISSymTransformedRank(seqIdx, i - 1, len, hint, env);
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
newBWTSeq(EISeq *seqIdx, MRAEnc *alphabet, Seqpos longest, Env *env)
{
  BWTSeq *bwtSeq;
  Seqpos *counts;
  unsigned alphabetSize;
  assert(seqIdx && env);
  env_error_check(env);
  /* alphabetSize is increased by one to handle the flattened
   * terminator symbol correctly */
  alphabetSize = MRAEncGetSize(alphabet) + 1;
  bwtSeq = ma_malloc(offsetAlign(sizeof (struct BWTSeq), sizeof (Seqpos))
                     + sizeof (Seqpos) * (alphabetSize + 1));
  counts = (Seqpos *)((char  *)bwtSeq
                      + offsetAlign(sizeof (struct BWTSeq),
                                    sizeof (Seqpos)));
  if (!initBWTSeqFromEncSeqIdx(bwtSeq, seqIdx, alphabet, longest, counts, env))
  {
    ma_free(bwtSeq);
    bwtSeq = NULL;
  }
  return bwtSeq;
}

void
deleteBWTSeq(BWTSeq *bwtSeq, Env *env)
{
  MRAEncDelete(bwtSeq->alphabet, env);
  deleteEISHint(bwtSeq->seqIdx, bwtSeq->hint, env);
  deleteEncIdxSeq(bwtSeq->seqIdx, env);
  ma_free(bwtSeq);
}

static inline void
getMatchBound(const BWTSeq *bwtSeq, const Symbol *query, size_t queryLen,
              struct matchBound *match, Env *env)
{
  size_t i = queryLen;
  const Seqpos *count;
  Symbol curSym;
  const MRAEnc *alphabet;
  assert(bwtSeq && query && env);
  env_error_check(env);
  count = bwtSeq->count;
  alphabet = BWTSeqGetAlphabet(bwtSeq);
  curSym = MRAEncMapSymbol(alphabet, query[--i]);
  match->upper = count[curSym];
  match->lower = count[curSym + 1];
  while ((match->upper <= match->lower) && (i > 0))
  {
    curSym = MRAEncMapSymbol(alphabet, query[--i]);
    match->upper = count[curSym]
      + BWTSeqOcc(bwtSeq, curSym, match->upper, env);
    match->lower = count[curSym]
      + BWTSeqOcc(bwtSeq, curSym, match->lower, env);
  }
}

extern Seqpos
BWTSeqMatchCount(const BWTSeq *bwtSeq, const Symbol *query, size_t queryLen,
                 Env *env)
{
  struct matchBound match;
  assert(bwtSeq && query && env);
  env_error_check(env);
  getMatchBound(bwtSeq, query, queryLen, &match, env);
  if (match.lower < match.upper)
    return 0;
  else
    return match.lower - match.upper;
}

struct BWTSeqExactMatchesIterator *
newEMIterator(const BWTSeq *bwtSeq, const Symbol *query, size_t queryLen,
              Env *env)
{
  struct BWTSeqExactMatchesIterator *newIter;
  assert(bwtSeq && query && env);
  env_error_check(env);
  if (!bwtSeq->locateSampleInterval)
  {
    fputs("Index does not contain locate information.\n"
          "Localization of matches impossible!", stderr);
    return NULL;
  }
  newIter = ma_malloc(sizeof (*newIter));
  getMatchBound(bwtSeq, query, queryLen, &newIter->bounds, env);
  newIter->nextMatchBWTPos = newIter->bounds.upper;
  initExtBitsRetrieval(&newIter->extBits, env);
  return newIter;
}

void
deleteEMIterator(struct BWTSeqExactMatchesIterator *iter, Env *env)
{
  destructExtBitsRetrieval(&iter->extBits, env);
  ma_free(iter);
}

Seqpos
EMINumMatchesTotal(const struct BWTSeqExactMatchesIterator *iter)
{
  assert(iter);
  if (iter->bounds.upper > iter->bounds.lower)
    return 0;
  else
    return iter->bounds.lower - iter->bounds.upper;
}

extern Seqpos
EMINumMatchesLeft(const struct BWTSeqExactMatchesIterator *iter)
{
  assert(iter);
  if (iter->nextMatchBWTPos > iter->bounds.lower)
    return 0;
  else
    return iter->bounds.lower - iter->bounds.upper;
}

extern int
BWTSeqVerifyIntegrity(BWTSeq *bwtSeq, const Str *projectName,
                      unsigned long tickPrint, FILE *fp, Env *env)
{
  Suffixarray suffixArray;
  struct extBitsRetrieval extBits;
  bool suffixArrayIsInitialized = false, extBitsAreInitialized = false;
  Verboseinfo *verbosity = NULL;
  enum verifyBWTSeqErrCode retval = VERIFY_BWTSEQ_NO_ERROR;
  do
  {
    Seqpos len;
    assert(bwtSeq && projectName && env);
    env_error_check(env);

    verbosity = newverboseinfo(true, env);
    initExtBitsRetrieval(&extBits, env);
    if (mapsuffixarray(&suffixArray, &len,
                       SARR_SUFTAB | SARR_ESQTAB, projectName, verbosity, env))
    {
      env_error_set(env, "Cannot load reference suffix array project with"
                    " demand for suffix table file and encoded sequence"
                    " for project: %s", str_get(projectName));
      freeverboseinfo(&verbosity, env);
      retval = VERIFY_BWTSEQ_REFLOAD_ERROR;
      break;
    }
    suffixArrayIsInitialized = true;
    ++len;
    if (BWTSeqLength(bwtSeq) != len)
    {
      env_error_set(env, "length mismatch for suffix array project %s and "
                    "bwt sequence index", str_get(projectName));
      retval = VERIFY_BWTSEQ_LENCOMPARE_ERROR;
      break;
    }

    if (BWTSeqHasLocateInformation(bwtSeq))
    {
      Seqpos i;
      for (i = 0; i < len && retval == VERIFY_BWTSEQ_NO_ERROR; ++i)
        if (BWTSeqPosHasLocateInfo(bwtSeq, i, &extBits, env))
        {
          Seqpos sfxArrayValue = BWTSeqLocateMatch(bwtSeq, i, &extBits, env);
          if (sfxArrayValue != suffixArray.suftab[i])
          {
            env_error_set(env, "Failed suffixarray value comparison"
                          " at position "FormatSeqpos": "FormatSeqpos" != "
                          FormatSeqpos,
                          i, sfxArrayValue, suffixArray.suftab[i]);
            retval = VERIFY_BWTSEQ_SUFVAL_ERROR;
            break;
          }
        }
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
        Symbol sym = EISGetSym(bwtSeq->seqIdx, nextLocate, bwtSeq->hint, env);
        if (sym != UNDEFBWTCHAR)
        {
          env_error_set(env, "symbol mismatch at position "FormatSeqpos": "
                        "%d vs. reference symbol %d", i - 1, sym,
                        UNDEFBWTCHAR);
          retval = VERIFY_BWTSEQ_LFMAPWALK_ERROR;
          break;
        }
        --i;
        nextLocate = BWTSeqLFMap(bwtSeq, nextLocate, env);
      }
      while (i > 0)
      {
        Symbol symRef = getencodedchar(suffixArray.encseq,
                                       --i, suffixArray.readmode);
        Symbol symCmp = EISGetSym(bwtSeq->seqIdx, nextLocate, bwtSeq->hint,
                                  env);
        if (symCmp != symRef)
        {
          env_error_set(env, "symbol mismatch at position "FormatSeqpos": "
                        "%d vs. reference symbol %d", i, symCmp, symRef);
          retval = VERIFY_BWTSEQ_LFMAPWALK_ERROR;
          break;
        }
        nextLocate = BWTSeqLFMap(bwtSeq, nextLocate, env);
      }
    }
  } while (0);
  if (suffixArrayIsInitialized)
    freesuffixarray(&suffixArray, env);
  if (verbosity)
    freeverboseinfo(&verbosity, env);
  if (extBitsAreInitialized)
    destructExtBitsRetrieval(&extBits, env);
  return retval;
}
