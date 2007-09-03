/*
** Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
**  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**  
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**  
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
**  
*/

#include <assert.h>
#include <string.h>

#include <libgtcore/env.h>
#include <libgtcore/str.h>
#include <libgtmatch/sarr-def.h>
#include <libgtmatch/sfx-map.pr>

#include "biofmi2misc.h"
#include "bwtseq.h"
#include "encidxseq.h"

struct BWTSeq
{
  enum seqBaseEncoding type;
  struct encIdxSeq *seqIdx;
  size_t alphabetSize;
  EISHint hint;
  Seqpos *count;
};

#if 0
static BitOffset
doNothing(BitString dest, BitOffset sOffset,
          Seqpos start, Seqpos end, void *cbState,
          Env *env)
{
  return 0;
}
#endif

struct seqRevMapEntry
{
  Seqpos bwtPos, origPos;
};

struct addLocateInfoState
{
  Suffixarray suffixArray;
  unsigned sampleInterval, bitsPerSeqpos;
  size_t revMapCacheSize;
  struct seqRevMapEntry *revMapCache;
};

static inline unsigned
bitsPerPosUpperBound(struct addLocateInfoState *state)
{
  return 1                      /* one bit for flag in preceding bitmap */
    + state->bitsPerSeqpos * 2;
}

static void
initAddLocateInfoState(struct addLocateInfoState *state, const Str *projectName,
                       unsigned sampleInterval, unsigned aggregationExpVal,
                       Env *env)
{
  Seqpos len;
  streamsuffixarray(&state->suffixArray, &len, SARR_SUFTAB | SARR_BWTTAB,
                    projectName, true, env);
  state->sampleInterval = sampleInterval;
  state->bitsPerSeqpos = requiredUInt64Bits(len - 1);
  state->revMapCacheSize = aggregationExpVal/sampleInterval + 1;
  state->revMapCache = env_ma_malloc(env, sizeof(state->revMapCache[0])
                                     * state->revMapCacheSize);
}

static void
destructAddLocateInfoState(struct addLocateInfoState *state, Env *env)
{
  env_ma_free(state->revMapCache, env);
  freesuffixarray(&state->suffixArray, env);
}

DECLAREREADFUNCTION(Seqpos)

static BitOffset
addLocateInfo(BitString dest, BitOffset sOffset,
              Seqpos start, Seqpos len, void *cbState,
              Env *env)
{
  BitOffset bitsWritten = len; /* 1 bit for each position in the bit mask */
  struct addLocateInfoState *state = cbState;
  unsigned bitsPerBWTPos, bitsPerOrigPos;
  int retcode;
  assert(dest && cbState);
  /* 0. resize caches if necessary */
  if(len > state->revMapCacheSize)
  {
    state->revMapCache = env_ma_realloc(env, state->revMapCache,
                                        len * sizeof(state->revMapCache[0]));
    state->revMapCacheSize = len;
  }
  bitsPerBWTPos = requiredUInt64Bits(len - 1);
  bitsPerOrigPos = state->bitsPerSeqpos;
  /* 1. read rangeLen suffix array indices from suftab */
  {
    Seqpos i, mapVal;
    size_t revMapCachePos = 0;
    for(i = 0; i < len; ++i)
    {
      /* 1.a read array index*/
      if((retcode = 
          readnextSeqposfromstream(&mapVal, &state->suffixArray.suftabstream,
                                   env)) != 1)
        return (BitOffset)-1;
      /* 1.b check wether the index into the original sequence is an
       * even multiple of the sampling interval */
      if(!(mapVal % state->sampleInterval))
      {
        /* 1.b.1 enter index into cache */
        state->revMapCache[revMapCachePos].bwtPos = start + i;
        state->revMapCache[revMapCachePos].origPos = mapVal;
        ++revMapCachePos;
        /* 1.b.2 mark position in bwt sequence */
        bsSetBit(dest, sOffset + i);
      }
      else
        bsClearBit(dest, sOffset + i);
    }
    /* 2. copy revMapCache into output */
    for(i = 0; i < revMapCachePos; ++i)
    {
      bsStoreUInt64(dest, sOffset + bitsWritten, bitsPerBWTPos,
                    state->revMapCache[i].bwtPos);
      bitsWritten += bitsPerBWTPos;
      bsStoreUInt64(dest, sOffset + bitsWritten, bitsPerOrigPos,
                    state->revMapCache[i].origPos);
      bitsWritten += bitsPerOrigPos;
    }
  }
  return bitsWritten;
}

BWTSeq *
newBWTSeq(enum seqBaseEncoding baseType, union bwtSeqParam *extraParams,
          const Str *projectName, Env *env)
{
  struct BWTSeq *bwtSeq = NULL;
  struct encIdxSeq *baseSeqIdx;
  const MRAEnc *alphabet;
  Symbol alphabetSize;
  EISHint hint;
  struct addLocateInfoState varState;
  assert(projectName && env);
  env_error_check(env);
  initAddLocateInfoState(&varState, projectName, 16,
                         extraParams->blockEncParams.blockSize
                         * extraParams->blockEncParams.blockSize, env);
  switch(baseType)
  {
  case BWT_ON_BLOCK_ENC:
    if(!(baseSeqIdx = loadBlockEncIdxSeq(projectName, env)))
    {
      env_error_unset(env);
      if(!(baseSeqIdx
           = newBlockEncIdxSeq(projectName,
                               extraParams->blockEncParams.blockSize,
                               addLocateInfo, bitsPerPosUpperBound(&varState),
                               &varState, env)))
      {
        destructAddLocateInfoState(&varState, env);
        return NULL;
      }
    }
    destructAddLocateInfoState(&varState, env);
    alphabet = EISGetAlphabet(baseSeqIdx);
    alphabetSize = MRAEncGetSize(alphabet);
    if(!alphabetSize)
    {
      /* weird error, shouldn't happen, but I prefer error return to
       * segfault in case someone tampered with the input */
      deleteEncIdxSeq(baseSeqIdx, env);
      return NULL;
    }
    bwtSeq = env_ma_malloc(env, offsetAlign(sizeof(struct BWTSeq),
                                            sizeof(Seqpos))
                           + sizeof(Seqpos) * (alphabetSize + 1));
    bwtSeq->count = (Seqpos *)((char  *)bwtSeq + sizeof(struct BWTSeq));
    bwtSeq->type = baseType;
    bwtSeq->seqIdx = baseSeqIdx;
    bwtSeq->alphabetSize = alphabetSize;
    bwtSeq->hint = hint = newEISHint(baseSeqIdx, env);
    {
      Symbol i;
      Seqpos lastPos = EISLength(baseSeqIdx) - 1, *count = bwtSeq->count;
      count[0] = 0;
      for(i = 0; i < alphabetSize; ++i)
        count[i + 1] = count[i]
          + EISSymTransformedRank(baseSeqIdx, i, lastPos, hint, env);
      fprintf(stderr, "count[alphabetSize]="FormatSeqpos
              ", lastPos="FormatSeqpos"\n", count[alphabetSize], lastPos);
      for(i = 0; i <= alphabetSize; ++i)
        fprintf(stderr, "count[%u]="FormatSeqpos"\n", i, count[i]);
      assert(count[alphabetSize] == lastPos + 1);
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

static inline Seqpos
BWTOcc(const BWTSeq *bwtSeq, Symbol tsym, Seqpos pos, Env *env)
{
  assert(bwtSeq && env);
  return EISSymTransformedRank(bwtSeq->seqIdx, tsym, pos, bwtSeq->hint, env);
}

extern Seqpos
matchCount(const BWTSeq *bwtSeq, Symbol *query, size_t queryLen, Env *env)
{
  Seqpos i = queryLen - 1, first, last, *count;
  Symbol curSym;
  const MRAEnc *alphabet;
  assert(bwtSeq && query);
  count = bwtSeq->count;
  alphabet = EISGetAlphabet(bwtSeq->seqIdx);
  curSym = MRAEncMapSymbol(alphabet, query[0]);
  first = count[curSym];
  last = count[curSym] - 1;
  while((first <= last) && (i > 0))
  {
    curSym = MRAEncMapSymbol(alphabet, query[i]);
    /* FIXME: what value should first get when first == 0? */
    first = count[curSym] + BWTOcc(bwtSeq, curSym, first, env) + 1;
    last = count[curSym] + BWTOcc(bwtSeq, curSym, last - 1, env);
    --i;
  }
  if(last < first)
    return -1;
  else
    return last - first;
}
