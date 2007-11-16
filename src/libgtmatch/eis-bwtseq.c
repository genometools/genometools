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

#include "libgtmatch/eis-bwtseq.h"
#include "libgtmatch/eis-bwtseqconstruct.h"
#include "libgtmatch/eis-bwtseqpriv.h"
#include "libgtmatch/eis-bwtseqcreate.h"
#include "libgtmatch/eis-encidxseq.h"
#include "libgtmatch/eis-encidxseqconstruct.h"
#include "libgtmatch/eis-encidxseqconstruct.h"

static int
initBWTSeqFromEncSeqIdx(struct BWTSeq *bwtSeq, struct encIdxSeq *baseSeqIdx,
                        Seqpos *counts, Env *env);

static BWTSeq *
newBWTSeq(struct encIdxSeq *seqIdx, Env *env);

DECLAREREADFUNCTION(Seqpos)

static int
streamReadSeqpos(Seqpos *dest, void *src, Env *env)
{
  Suffixarray *suffixArray = src;
  assert(suffixArray);
  return readnextSeqposfromstream(dest, &suffixArray->suftabstream, env);
}

extern BWTSeq *
availBWTSeq(const struct bwtParam *createParams, Env *env)
{
  struct BWTSeq *bwtSeq = NULL;
  Suffixarray suffixArray;
  Seqpos len;
  Verboseinfo *verbosity;
  /* FIXME: handle verbosity in a more sane fashion */
  verbosity = newverboseinfo(false, env);
  if (streamsuffixarray(&suffixArray, &len, SARR_SUFTAB | SARR_BWTTAB,
                        createParams->projectName, verbosity, env))
  {
    freeverboseinfo(&verbosity, env);
    return NULL;
  }
  ++len;
  bwtSeq = availBWTSeqFromSA(createParams, &suffixArray, len, env);
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
  /* try loading index */
  bwtSeq = loadBWTSeqForSA(params, sa, totalLen, env);
  /* if loading didn't work try on-demand creation */
  if (!bwtSeq)
  {
    bwtSeq = createBWTSeqFromSA(params, sa, totalLen, env);
  }
  return bwtSeq;
}

extern BWTSeq *
loadBWTSeqForSA(const struct bwtParam *params, Suffixarray *sa,
                Seqpos totalLen, Env *env)
{
  struct BWTSeq *bwtSeq = NULL;
  EISeq *seqIdx = NULL;
  switch (params->baseType)
  {
  case BWT_ON_BLOCK_ENC:
    if ((seqIdx = loadBlockEncIdxSeqForSA(
           sa, totalLen, params->projectName,
           params->seqParams.blockEnc.EISFeatureSet, env)))
    {
      if (!(bwtSeq = newBWTSeq(seqIdx, env)))
        break;
      fputs("Using pre-computed sequence index.\n", stderr);
    }
    break;
  default:
    fprintf(stderr, "Illegal/unknown/unimplemented encoding requested!"
            " (%s:%d)\n", __FILE__, __LINE__);
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
  if (params->locateInterval &&
      (!sa->suftabstream.fp || !sa->longest.defined))
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
          params, sa, (indexCreateFunc)newBlockEncIdxSeqFromSA,
          params->seqParams.blockEnc.blockSize *
          params->seqParams.blockEnc.bucketBlocks, streamReadSeqpos, totalLen,
          env);
      break;
    default:
      fprintf(stderr, "Illegal/unknown/unimplemented encoding requested!"
              " (%s:%d)\n", __FILE__, __LINE__);
      break;
    }
    if (seqIdx)
      bwtSeq = newBWTSeq(seqIdx, env);
    if (!bwtSeq && seqIdx)
      deleteEncIdxSeq(seqIdx, env);
  }
  return bwtSeq;
}

struct sfxIReadInfo
{
  sfxInterface *si;
  listenerID id;
};

static EISeq *
createBlockEncIdxSeqFromSfxIRI(void *src, Seqpos totalLen,
                               const Str *projectName,
                               const union bwtSeqParam *params,
                               size_t numExtHeaders, uint16_t *headerIDs,
                               uint32_t *extHeaderSizes,
                               headerWriteFunc *extHeaderCallbacks,
                               void **headerCBData,
                               bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                               BitOffset maxVarExtBitsPerPos, void *cbState,
                               Env *env)
{
  assert(src);
  return newBlockEncIdxSeqFromSfxI(((struct sfxIReadInfo *)src)->si, totalLen,
                                   projectName, &params->blockEnc,
                                   numExtHeaders,
                                   headerIDs, extHeaderSizes,
                                   extHeaderCallbacks, headerCBData, biFunc,
                                   cwExtBitsPerPos, maxVarExtBitsPerPos,
                                   cbState, env);
}

static int
sfxIReadSeqpos(Seqpos *dest, void *src, Env *env)
{
  return readSfxISufTabRange(((struct sfxIReadInfo *)src)->si,
                             ((struct sfxIReadInfo *)src)->id,
                             1, dest, env) == 1;
}

extern BWTSeq *
createBWTSeqFromSfxI(const struct bwtParam *params, sfxInterface *si,
                     Seqpos totalLen, Env *env)
{
  struct sfxIReadInfo siri;
  EISeq *seqIdx = NULL;
  BWTSeq *bwtSeq = NULL;
  assert(si && params && env);
  siri.si = si;
  if (params->locateInterval)
  {
    if (!SfxIRegisterReader(si, &siri.id, SFX_REQUEST_SUFTAB, env))
      return NULL;
  }
  seqIdx= createBWTSeqGeneric(params, &siri, createBlockEncIdxSeqFromSfxIRI,
                              params->seqParams.blockEnc.blockSize *
                              params->seqParams.blockEnc.bucketBlocks,
                              sfxIReadSeqpos, totalLen, env);
  if (seqIdx)
    bwtSeq = newBWTSeq(seqIdx, env);
  if (!bwtSeq && seqIdx)
    deleteEncIdxSeq(seqIdx, env);
  return bwtSeq;
}

static int
initBWTSeqFromEncSeqIdx(BWTSeq *bwtSeq, struct encIdxSeq *seqIdx,
                        Seqpos *counts, Env *env)
{
  const MRAEnc *alphabet;
  Symbol alphabetSize;
  EISHint hint;
  assert(bwtSeq && seqIdx);
  alphabet = EISGetAlphabet(seqIdx);
  alphabetSize = MRAEncGetSize(alphabet);
  if (!alphabetSize)
    /* weird error, shouldn't happen, but I prefer error return to
     * segfault in case someone tampered with the input */
    return 0;
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
    }
  }
  bwtSeq->hint = hint = newEISHint(seqIdx, env);
  {
    Symbol i;
    Seqpos len = EISLength(seqIdx), *count = bwtSeq->count;
    count[0] = 0;
    for (i = 0; i < alphabetSize; ++i)
      count[i + 1] = count[i]
        + EISSymTransformedRank(seqIdx, i, len, hint, env);
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

static BWTSeq *
newBWTSeq(EISeq *seqIdx, Env *env)
{
  BWTSeq *bwtSeq;
  Seqpos *counts;
  unsigned alphabetSize;
  assert(seqIdx && env);
  env_error_check(env);
  alphabetSize = MRAEncGetSize(EISGetAlphabet(seqIdx));
  bwtSeq = env_ma_malloc(env, offsetAlign(sizeof (struct BWTSeq),
                                          sizeof (Seqpos))
                         + sizeof (Seqpos) * (alphabetSize + 1));
  counts = (Seqpos *)((char  *)bwtSeq
                      + offsetAlign(sizeof (struct BWTSeq),
                                    sizeof (Seqpos)));
  if (!initBWTSeqFromEncSeqIdx(bwtSeq, seqIdx, counts, env))
  {
    env_ma_free(bwtSeq, env);
    bwtSeq = NULL;
  }
  return bwtSeq;
}

void
deleteBWTSeq(BWTSeq *bwtSeq, Env *env)
{
  deleteEISHint(bwtSeq->seqIdx, bwtSeq->hint, env);
  deleteEncIdxSeq(bwtSeq->seqIdx, env);
  env_ma_free(bwtSeq, env);
}

int
BWTSeqHasLocateInformation(const BWTSeq *bwtSeq)
{
  return bwtSeq->locateSampleInterval;
}

static inline int
BWTSeqPosHasLocateInfo(const BWTSeq *bwtSeq, Seqpos pos,
                       struct extBitsRetrieval *extBits, Env *env)
{
  EISRetrieveExtraBits(bwtSeq->seqIdx, pos, EBRF_RETRIEVE_CWBITS, extBits,
                       bwtSeq->hint, env);
  return bsGetBit(extBits->cwPart, extBits->cwOffset + pos - extBits->start);
}

static inline void
getMatchBound(const BWTSeq *bwtSeq, const Symbol *query, size_t queryLen,
              struct matchBound *match, Env *env)
{
  size_t i = queryLen;
  const Seqpos *count;
  Symbol curSym;
  const MRAEnc *alphabet;
  assert(bwtSeq && query);
  count = bwtSeq->count;
  alphabet = EISGetAlphabet(bwtSeq->seqIdx);
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
  getMatchBound(bwtSeq, query, queryLen, &match, env);
  if (match.lower < match.upper)
    return 0;
  else
    return match.lower - match.upper;
}

struct BWTSeqExactMatchesIterator
{
  struct matchBound bounds;
  Seqpos nextMatchBWTPos;
  struct MatchData nextMatch;
  struct extBitsRetrieval extBits;
};

struct BWTSeqExactMatchesIterator *
newEMIterator(const BWTSeq *bwtSeq, const Symbol *query, size_t queryLen,
              Env *env)
{
  struct BWTSeqExactMatchesIterator *newIter;
  assert(bwtSeq && query && env);
  if (!bwtSeq->locateSampleInterval)
  {
    fputs("Index does not contain locate information.\n"
          "Localization of matches impossible!", stderr);
    return NULL;
  }
  newIter = env_ma_malloc(env, sizeof (*newIter));
  getMatchBound(bwtSeq, query, queryLen, &newIter->bounds, env);
  newIter->nextMatchBWTPos = newIter->bounds.upper;
  initExtBitsRetrieval(&newIter->extBits, env);
  return newIter;
}

void
deleteEMIterator(struct BWTSeqExactMatchesIterator *iter, Env *env)
{
  env_ma_free(iter, env);
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

struct MatchData *
EMIGetNextMatch(struct BWTSeqExactMatchesIterator *iter, const BWTSeq *bwtSeq,
                Env *env)
{
  if (iter->nextMatchBWTPos < iter->bounds.lower)
  {
    Seqpos nextLocate = iter->nextMatchBWTPos, locateOffset = 0;
    while (!BWTSeqPosHasLocateInfo(bwtSeq, nextLocate, &iter->extBits, env))
      nextLocate = BWTSeqLFMap(bwtSeq, nextLocate, env), ++locateOffset;
    EISRetrieveExtraBits(bwtSeq->seqIdx, nextLocate,
                         EBRF_RETRIEVE_CWBITS | EBRF_RETRIEVE_VARBITS,
                         &iter->extBits, bwtSeq->hint, env);
    {
      unsigned bitsPerBWTPos = requiredUInt64Bits(iter->extBits.len - 1),
        bitsPerOrigPos = requiredUInt64Bits(
          (EISLength(bwtSeq->seqIdx) - 1)/bwtSeq->locateSampleInterval);
      BitOffset locateRecordIndex =
        bs1BitsCount(iter->extBits.cwPart, iter->extBits.cwOffset,
                     nextLocate - iter->extBits.start),
        locateRecordOffset = (bitsPerBWTPos + bitsPerOrigPos)
        * locateRecordIndex;
      iter->nextMatch.sfxArrayValue =
        bsGetUInt64(iter->extBits.varPart, iter->extBits.varOffset
                    + locateRecordOffset + bitsPerBWTPos, bitsPerOrigPos)
        * bwtSeq->locateSampleInterval + locateOffset;
      assert(bsGetUInt64(iter->extBits.varPart,
                         iter->extBits.varOffset + locateRecordOffset,
                         bitsPerBWTPos)
             == nextLocate - iter->extBits.start);
    }
    {
      /* FIXME: map position back to original encoded sequence */
      iter->nextMatch.dbFile = 0;
    }
    ++iter->nextMatchBWTPos;
    return &iter->nextMatch;
  }
  else
    return NULL;
}
