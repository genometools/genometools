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
#include "libgtmatch/eis-bwtseq.h"
#include "libgtmatch/eis-bwtseqcreate.h"
#include "libgtmatch/eis-bwtseqpriv.h"
/*************************************************************
 * generic method for bwt index creation
 *************************************************************/

enum {
  LOCATE_INFO_IN_INDEX_HEADERID = 1111,
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
extern int
readLocateInfoHeader(EISeq *seqIdx, struct locateHeader *headerData)
{
  FILE *fp;
  assert(seqIdx && headerData);
  if (!(fp = EISSeekToHeader(seqIdx, LOCATE_INFO_IN_INDEX_HEADERID, NULL)))
    return 0;
  if (fread(&headerData->locateInterval, sizeof (headerData->locateInterval),
            1, fp) != 1)
    return 0;
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

extern EISeq *
createBWTSeqGeneric(const struct bwtParam *params, void *baseSrc,
                    indexCreateFunc createIndex, unsigned aggInterval,
                    srcReadFunc readCallback, Seqpos totalLen, Env *env)
{
  struct encIdxSeq *baseSeqIdx = NULL;
  struct addLocateInfoState varState;
  unsigned locateInterval;
  assert(baseSrc && params && env);
  env_error_check(env);
  locateInterval = params->locateInterval;
  varState.src = NULL;
  do
  {
    struct locateHeader headerData = { locateInterval };
    void *p[] = { &headerData };
    uint16_t headerIDs[] = { LOCATE_INFO_IN_INDEX_HEADERID };
    uint32_t headerSizes[] = { LOCATE_HEADER_SIZE };
    headerWriteFunc headerFuncs[] = { writeLocateInfoHeader };
    env_error_unset(env);
    initAddLocateInfoState(&varState, baseSrc, readCallback, totalLen,
                           params->projectName, locateInterval,
                           aggInterval, env);
    if (locateInterval)
    {
      if (!(baseSeqIdx
            = createIndex(baseSrc, totalLen, params->projectName,
                          &params->seqParams, sizeof (p)/sizeof (p[0]),
                          headerIDs, headerSizes, headerFuncs, p,
                          addLocateInfo, 1 /* one bit for flag */,
                          bitsPerPosUpperBound(&varState),
                          &varState, env)))
        break;
    }
    else
    {
      if (!(baseSeqIdx
            = createIndex(baseSrc, totalLen, params->projectName,
                          &params->seqParams, 0, NULL, NULL, NULL,
                          NULL, NULL, 0, 0, &varState, env)))
        break;
    }
  } while (0);
  if (varState.src)
    destructAddLocateInfoState(&varState, env);
  return baseSeqIdx;
}
