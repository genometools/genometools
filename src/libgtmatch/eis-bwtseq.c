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
#include "libgtmatch/eis-bwtseqpriv.h"
#include "libgtmatch/eis-encidxseq.h"

typedef int (*srcReadFunc)(Seqpos *dest, void *src, Env *env);
typedef struct encIdxSeq *(*indexCreateFunc)
  (void *src, Seqpos totalLen, const Str *projectName, unsigned blockSize,
   unsigned bucketBlocks, int features, size_t numExtHeaders,
   uint16_t *headerIDs, uint32_t *extHeaderSizes,
   headerWriteFunc *extHeaderCallbacks, void **headerCBData,
   bitInsertFunc biFunc, BitOffset cwBitsPerPos, BitOffset maxBitsPerPos,
   void *cbState, Env *env);
typedef struct encIdxSeq *(*indexLoadFunc)(void *src, Seqpos totalLen,
                                           const Str *projectName,
                                           int features, Env *env);

static BWTSeq *
newBWTSeqGen(const struct bwtParam *params, void *baseSrc,
             indexCreateFunc createIndex, indexLoadFunc loadIndex,
             srcReadFunc readCallback, Seqpos totalLen, Env *env);

enum {
  LOCATE_INFO_IN_INDEX_HEADERID = 1111,
};

struct locateHeader
{
  unsigned locateInterval;
};

enum {
  LOCATE_HEADER_SIZE = sizeof (uint32_t),
};

static int
writeLocateInfoHeader(FILE *fp, void *cbData)
{
  const struct locateHeader *headerData = cbData;
  uint32_t locateInterval;
  assert(cbData);
  locateInterval = headerData->locateInterval;
  return fwrite(&locateInterval, sizeof (locateInterval), 1, fp);
}

/**
 * @return 0 on error
 */
static inline int
readLocateInfoHeader(FILE *fp, struct locateHeader *headerData)
{
  uint32_t locateInterval;
  assert(headerData);
  if (fread(&locateInterval, sizeof (locateInterval), 1, fp) != 1)
    return 0;
  headerData->locateInterval = locateInterval;
  return 1;
}

struct seqRevMapEntry
{
  Seqpos bwtPos, origPos;
};

struct addLocateInfoState
{
  void *src;
  srcReadFunc readSeqpos;
  unsigned locateInterval, bitsPerOrigPos, bitsPerSeqpos;
  int locateFeatures;
  size_t revMapCacheSize;
  struct seqRevMapEntry *revMapCache;
};

static inline unsigned
bitsPerPosUpperBound(struct addLocateInfoState *state)
{
  return state->bitsPerOrigPos + state->bitsPerSeqpos;
}

static void
initAddLocateInfoState(struct addLocateInfoState *state,
                       void *src, srcReadFunc readFunc,
                       Seqpos srcLen, const Str *projectName,
                       unsigned locateInterval,
                       unsigned aggregationExpVal, Env *env)
{
  Seqpos lastPos;
  assert(src && projectName);
  state->src = src;
  state->readSeqpos = readFunc;
  lastPos = srcLen - 1;
  state->locateInterval = locateInterval;
  state->bitsPerSeqpos = requiredUInt64Bits(lastPos);
  if (locateInterval)
  {
    state->bitsPerOrigPos = requiredUInt64Bits(lastPos/locateInterval);
    state->revMapCacheSize = aggregationExpVal/locateInterval + 1;
    state->revMapCache = env_ma_malloc(env, sizeof (state->revMapCache[0])
                                       * state->revMapCacheSize);
  }
  else
  {
    state->bitsPerOrigPos = 0;
    state->revMapCacheSize = 0;
    state->revMapCache = NULL;
  }
}

static void
destructAddLocateInfoState(struct addLocateInfoState *state, Env *env)
{
  env_ma_free(state->revMapCache, env);
}

#if 1
DECLAREREADFUNCTION(Seqpos)

static int
streamReadSeqpos(Seqpos *dest, void *src, Env *env)
{
  Suffixarray *suffixArray = src;
  assert(suffixArray);
  return readnextSeqposfromstream(dest, &suffixArray->suftabstream, env);
}

static BitOffset
addLocateInfo(BitString cwDest, BitOffset cwOffset,
              BitString varDest, BitOffset varOffset,
              Seqpos start, Seqpos len, void *cbState,
              Env *env)
{
  BitOffset bitsWritten = 0;
  struct addLocateInfoState *state = cbState;
  unsigned bitsPerBWTPos, bitsPerOrigPos;
  int retcode;
  assert(varDest && cbState);
  /* 0. resize caches if necessary */
  if (len > state->revMapCacheSize)
  {
    state->revMapCache = env_ma_realloc(env, state->revMapCache,
                                        len * sizeof (state->revMapCache[0]));
    state->revMapCacheSize = len;
  }
  bitsPerBWTPos = requiredUInt64Bits(len - 1);
  bitsPerOrigPos = state->bitsPerOrigPos;
  /* 1. read rangeLen suffix array indices from suftab */
  {
    Seqpos i, mapVal;
    size_t revMapCacheLen = 0;
    for (i = 0; i < len; ++i)
    {
      /* 1.a read array index*/
      if ((retcode = state->readSeqpos(&mapVal, state->src, env)) != 1)
        return (BitOffset)-1;
      /* 1.b check wether the index into the original sequence is an
       * even multiple of the sampling interval */
      if (!(mapVal % state->locateInterval))
      {
        /* 1.b.1 enter index into cache */
        state->revMapCache[revMapCacheLen].bwtPos = i;
        state->revMapCache[revMapCacheLen].origPos
          = mapVal / state->locateInterval;
        ++revMapCacheLen;
        /* 1.b.2 mark position in bwt sequence */
        bsSetBit(cwDest, cwOffset + i);
      }
      else
        bsClearBit(cwDest, cwOffset + i);
    }
    /* 2. copy revMapCache into output */
    for (i = 0; i < revMapCacheLen; ++i)
    {
      bsStoreUInt64(varDest, varOffset + bitsWritten, bitsPerBWTPos,
                    state->revMapCache[i].bwtPos);
      bitsWritten += bitsPerBWTPos;
      bsStoreUInt64(varDest, varOffset + bitsWritten, bitsPerOrigPos,
                    state->revMapCache[i].origPos);
      bitsWritten += bitsPerOrigPos;
    }
  }
  return bitsWritten;
}
#endif

extern BWTSeq *
newBWTSeq(const struct bwtParam *params, Env *env)
{
  struct BWTSeq *bwtSeq = NULL;
  Suffixarray suffixArray;
  Seqpos len;
  Verboseinfo *verbosity;
  /* FIXME: handle verbosity in a more sane fashion */
  verbosity = newverboseinfo(false, env);
  if (streamsuffixarray(&suffixArray, &len, SARR_SUFTAB | SARR_BWTTAB,
                        params->projectName, verbosity, env))
  {
    freeverboseinfo(&verbosity, env);
    return NULL;
  }
  ++len;
  bwtSeq = newBWTSeqFromSA(params, &suffixArray, len, env);
  freesuffixarray(&suffixArray, env);
  freeverboseinfo(&verbosity, env);
  return bwtSeq;
}

BWTSeq *
newBWTSeqFromSA(const struct bwtParam *params, Suffixarray *sa,
                Seqpos totalLen, Env *env)
{
  assert(sa && params && env);
  if (params->locateInterval &&
      (!sa->suftabstream.fp || !sa->longest.defined))
  {
    fprintf(stderr, "error: locate sampling requested but not available"
            " for project %s\n", str_get(params->projectName));
    return NULL;
  }
  return newBWTSeqGen(params, sa,
                      (indexCreateFunc)newBlockEncIdxSeqFromSA,
                      (indexLoadFunc)loadBlockEncIdxSeqForSA,
                      streamReadSeqpos, totalLen, env);
}

static EISeq *
loadBlockEncIdxSeqForSfxI(void *srcNotUsed, Seqpos totalLen,
                          const Str *projectName, int features, Env *env)
{
  return loadBlockEncIdxSeq(projectName, features, env);
}

struct sfxIReadInfo
{
  sfxInterface *si;
  listenerID id;
};

static EISeq *
newBlockEncIdxSeqFromSfxIRI(void *src, Seqpos totalLen,
                            const Str *projectName, unsigned blockSize,
                            unsigned bucketBlocks, int features,
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
                                   projectName, blockSize, bucketBlocks,
                                   features, numExtHeaders, headerIDs,
                                   extHeaderSizes, extHeaderCallbacks,
                                   headerCBData, biFunc, cwExtBitsPerPos,
                                   maxVarExtBitsPerPos, cbState, env);
}

static int
sfxIReadSeqpos(Seqpos *dest, void *src, Env *env)
{
  return readSfxISufTabRange(((struct sfxIReadInfo *)src)->si,
                             ((struct sfxIReadInfo *)src)->id,
                             1, dest, env) == 1;
}

BWTSeq *
newBWTSeqFromSfxI(const struct bwtParam *params, sfxInterface *si,
                  Seqpos totalLen, Env *env)
{
  struct sfxIReadInfo siri;
  assert(si && params && env);
  siri.si = si;
  if (params->locateInterval)
  {
    if (!SfxIRegisterReader(si, &siri.id, SFX_REQUEST_SUFTAB, env))
      return NULL;
  }
  return newBWTSeqGen(params, &siri, newBlockEncIdxSeqFromSfxIRI,
                      (indexLoadFunc)loadBlockEncIdxSeqForSfxI,
                      sfxIReadSeqpos, totalLen, env);
}

#define newBWTSeqErrRet()                                \
  do {                                                   \
    if (varState.src)                                    \
      destructAddLocateInfoState(&varState, env);        \
    if (baseSeqIdx) deleteEncIdxSeq(baseSeqIdx, env);    \
    return NULL;                                         \
  } while (0)

/* FIXME: make flag bits optional */
static BWTSeq *
newBWTSeqGen(const struct bwtParam *params, void *baseSrc,
             indexCreateFunc createIndex, indexLoadFunc loadIndex,
             srcReadFunc readCallback, Seqpos totalLen, Env *env)
{
  struct BWTSeq *bwtSeq = NULL;
  struct encIdxSeq *baseSeqIdx = NULL;
  const MRAEnc *alphabet;
  Symbol alphabetSize;
  EISHint hint;
  struct addLocateInfoState varState;
  unsigned locateInterval;
  assert(baseSrc && params && env);
  locateInterval = params->locateInterval;
  env_error_check(env);
  varState.src = NULL;
  switch (params->baseType)
  {
  case BWT_ON_BLOCK_ENC:
    if (!(baseSeqIdx = loadIndex(baseSrc, totalLen, params->projectName,
                                 params->seqParams.blockEnc.EISFeatureSet,
                                 env)))
    {
      struct locateHeader headerData = { locateInterval };
      void *p[] = { &headerData };
      uint16_t headerIDs[] = { LOCATE_INFO_IN_INDEX_HEADERID };
      uint32_t headerSizes[] = { LOCATE_HEADER_SIZE };
      headerWriteFunc headerFuncs[] = { writeLocateInfoHeader };
      env_error_unset(env);
      initAddLocateInfoState(&varState, baseSrc, readCallback, totalLen,
                             params->projectName, locateInterval,
                             params->seqParams.blockEnc.blockSize
                             * params->seqParams.blockEnc.blockSize, env);
      if (locateInterval)
      {
        if (!(baseSeqIdx
             = createIndex(baseSrc, totalLen, params->projectName,
                           params->seqParams.blockEnc.blockSize,
                           params->seqParams.blockEnc.bucketBlocks,
                           params->seqParams.blockEnc.EISFeatureSet,
                           sizeof (p)/sizeof (p[0]), headerIDs,
                           headerSizes,
                           headerFuncs, p,
                           addLocateInfo, 1 /* one bit for flag */,
                           bitsPerPosUpperBound(&varState),
                           &varState, env)))
          newBWTSeqErrRet();
      }
      else
      {
        if (!(baseSeqIdx
             = createIndex(baseSrc, totalLen, params->projectName,
                           params->seqParams.blockEnc.blockSize,
                           params->seqParams.blockEnc.bucketBlocks,
                           params->seqParams.blockEnc.EISFeatureSet,
                           0, NULL, NULL, NULL, NULL, NULL, 0, 0,
                           &varState, env)))
          newBWTSeqErrRet();
      }
      destructAddLocateInfoState(&varState, env);
    }
    else
    {
      FILE *fp = EISSeekToHeader(baseSeqIdx, LOCATE_INFO_IN_INDEX_HEADERID,
                                 NULL);
      fputs("Using pre-computed sequence index.\n", stderr);
      if (!fp)
      {
        fputs("Specified index does not contain locate information.\n"
              "Localization of matches will not be supported!\n", stderr);
        locateInterval = 0;
      }
      else
      {
        struct locateHeader header;
        if (!readLocateInfoHeader(fp, &header))
          newBWTSeqErrRet();
        locateInterval = header.locateInterval;
      }
    }
    alphabet = EISGetAlphabet(baseSeqIdx);
    alphabetSize = MRAEncGetSize(alphabet);
    if (!alphabetSize)
      /* weird error, shouldn't happen, but I prefer error return to
       * segfault in case someone tampered with the input */
      newBWTSeqErrRet();
    bwtSeq = env_ma_malloc(env, offsetAlign(sizeof (struct BWTSeq),
                                            sizeof (Seqpos))
                           + sizeof (Seqpos) * (alphabetSize + 1));
    bwtSeq->count = (Seqpos *)((char  *)bwtSeq
                               + offsetAlign(sizeof (struct BWTSeq),
                                             sizeof (Seqpos)));
    bwtSeq->type = params->baseType;
    bwtSeq->seqIdx = baseSeqIdx;
    bwtSeq->alphabetSize = alphabetSize;
    bwtSeq->locateSampleInterval = locateInterval;
    bwtSeq->hint = hint = newEISHint(baseSeqIdx, env);
    {
      Symbol i;
      Seqpos len = EISLength(baseSeqIdx), *count = bwtSeq->count;
      count[0] = 0;
      for (i = 0; i < alphabetSize; ++i)
        count[i + 1] = count[i]
          + EISSymTransformedRank(baseSeqIdx, i, len, hint, env);
#ifdef DEBUG
      fprintf(stderr, "count[alphabetSize]="FormatSeqpos
              ", len="FormatSeqpos"\n", count[alphabetSize], len);
      for (i = 0; i <= alphabetSize; ++i)
        fprintf(stderr, "count[%u]="FormatSeqpos"\n", (unsigned)i, count[i]);
#endif
      assert(count[alphabetSize] == len);
    }
    break;
    /* FIXME: implement missing types. */
  case BWT_ON_RLE:
  case BWT_ON_WAVELET_TREE_ENC:
  default:
    fprintf(stderr, "Illegal/Unknown encoding requested! (%s:%d)\n",
            __FILE__, __LINE__);
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
