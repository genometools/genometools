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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if (!defined DEBUG) && (!defined NDEBUG)
#define NDEBUG
#endif /* DEBUG */

#include <assert.h>
#include <string.h>

#include <libgtcore/env.h>
#include <libgtcore/str.h>
#include <libgtmatch/sarr-def.h>
#include <libgtmatch/esa-map.pr>

#include "biofmi2misc.h"
#include "bwtseq.h"
#include "bwtseqpriv.h"
#include "encidxseq.h"

enum {
  LOCATE_INFO_IN_INDEX_HEADERID = 1111,
};

struct locateHeader
{
  unsigned locateInterval;
};

enum {
  LOCATE_HEADER_SIZE = sizeof(uint32_t), 
};

static int
writeLocateInfoHeader(FILE *fp, void *cbData)
{
  const struct locateHeader *headerData = cbData;
  uint32_t locateInterval;
  assert(cbData);
  locateInterval = headerData->locateInterval;
  return fwrite(&locateInterval, sizeof(locateInterval), 1, fp);
}

/**
 * @return 0 on error
 */
static inline int
readLocateInfoHeader(FILE *fp, struct locateHeader *headerData)
{
  uint32_t locateInterval;
  assert(headerData);
  if(fread(&locateInterval, sizeof(locateInterval), 1, fp) != 1)
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
  Suffixarray *suffixArray;
  unsigned locateInterval, bitsPerOrigPos, bitsPerSeqpos;
  size_t revMapCacheSize;
  struct seqRevMapEntry *revMapCache;
};

static inline unsigned
bitsPerPosUpperBound(struct addLocateInfoState *state)
{
  return state->bitsPerOrigPos + state->bitsPerSeqpos;
}

static int
initAddLocateInfoState(struct addLocateInfoState *state,
                       Suffixarray *sa, const Str *projectName,
                       unsigned locateInterval, unsigned aggregationExpVal,
                       Env *env)
{
  Seqpos lastPos;
  if(locateInterval && (!sa->suftabstream.fp || !sa->longest.defined))
  {
    fprintf(stderr, "error: locate sampling requested but not available"
            " for project %s\n", str_get(projectName));
    return 0;
  }
  state->suffixArray = sa;
  lastPos = sa->longest.valueseqpos;
  state->locateInterval = locateInterval;
  state->bitsPerSeqpos = requiredUInt64Bits(lastPos);
  if(locateInterval)
  {
    state->bitsPerOrigPos = requiredUInt64Bits(lastPos/locateInterval);
    state->revMapCacheSize = aggregationExpVal/locateInterval + 1;
    state->revMapCache = env_ma_malloc(env, sizeof(state->revMapCache[0])
                                       * state->revMapCacheSize);
  }
  else
  {
    state->bitsPerOrigPos = 0;
    state->revMapCacheSize = 0;
    state->revMapCache = NULL;
  }
  return 1;
}

static void
destructAddLocateInfoState(struct addLocateInfoState *state, Env *env)
{
  env_ma_free(state->revMapCache, env);
  freesuffixarray(state->suffixArray, env);
}

#if 1
DECLAREREADFUNCTION(Seqpos)

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
  if(len > state->revMapCacheSize)
  {
    state->revMapCache = env_ma_realloc(env, state->revMapCache,
                                        len * sizeof(state->revMapCache[0]));
    state->revMapCacheSize = len;
  }
  bitsPerBWTPos = requiredUInt64Bits(len - 1);
  bitsPerOrigPos = state->bitsPerOrigPos;
  /* 1. read rangeLen suffix array indices from suftab */
  {
    Seqpos i, mapVal;
    size_t revMapCacheLen = 0;
    for(i = 0; i < len; ++i)
    {
      /* 1.a read array index*/
      if((retcode =
          readnextSeqposfromstream(&mapVal, &state->suffixArray->suftabstream,
                                   env)) != 1)
        return (BitOffset)-1;
      /* 1.b check wether the index into the original sequence is an
       * even multiple of the sampling interval */
      if(!(mapVal % state->locateInterval))
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
    for(i = 0; i < revMapCacheLen; ++i)
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

BWTSeq *
newBWTSeq(enum seqBaseEncoding baseType, unsigned locateInterval,
          union bwtSeqParam *extraParams, const Str *projectName, Env *env)
{
  struct BWTSeq *bwtSeq = NULL;
  Suffixarray suffixArray;
  Seqpos len;
  if(streamsuffixarray(&suffixArray, &len, SARR_SUFTAB | SARR_BWTTAB,
                       projectName, NULL, env))
    return NULL;
  ++len;
  bwtSeq = newBWTSeqFromSA(baseType, locateInterval, extraParams,
                           &suffixArray, len, projectName, env);
  freesuffixarray(&suffixArray, env);
  return bwtSeq;
}


#define newBWTSeqErrRet()                               \
  do {                                                  \
    if(varState.suffixArray)                            \
      destructAddLocateInfoState(&varState, env);       \
    if(baseSeqIdx) deleteEncIdxSeq(baseSeqIdx, env);    \
    return NULL;                                        \
  } while(0)


BWTSeq *
newBWTSeqFromSA(enum seqBaseEncoding baseType, unsigned locateInterval,
                union bwtSeqParam *extraParams, Suffixarray *sa,
                Seqpos totalLen, const Str *projectName, Env *env)
{
  struct BWTSeq *bwtSeq = NULL;
  struct encIdxSeq *baseSeqIdx = NULL;
  const MRAEnc *alphabet;
  Symbol alphabetSize;
  EISHint hint;
  struct addLocateInfoState varState;
  assert(sa && projectName && env);
  env_error_check(env);
  varState.suffixArray = sa;
  switch(baseType)
  {
  case BWT_ON_BLOCK_ENC:
    if(!(baseSeqIdx = loadBlockEncIdxSeqForSA(sa, totalLen, projectName, env)))
    {
      struct locateHeader headerData = { locateInterval };
      void *p[] = { &headerData };
      uint16_t headerIDs[] = { LOCATE_INFO_IN_INDEX_HEADERID };
      uint32_t headerSizes[] = { LOCATE_HEADER_SIZE };
      headerWriteFunc headerFuncs[] = { writeLocateInfoHeader };
      if(!initAddLocateInfoState(&varState, sa, projectName, locateInterval,
                                 extraParams->blockEncParams.blockSize
                                 * extraParams->blockEncParams.blockSize, env))
        newBWTSeqErrRet();        
      env_error_unset(env);
      if(locateInterval)
      {
        if(!(baseSeqIdx
             = newBlockEncIdxSeqFromSA(sa, totalLen, projectName,
                                       extraParams->blockEncParams.blockSize,
                                       extraParams->blockEncParams.bucketBlocks,
                                       sizeof(p)/sizeof(p[0]), headerIDs,
                                       headerSizes,
                                       headerFuncs, p,
/*                                0, NULL, NULL, NULL, NULL, */
                                       addLocateInfo, 1 /* one bit for flag */,
                                       bitsPerPosUpperBound(&varState),
                                       &varState, env)))
          newBWTSeqErrRet();
      }
      else
      {
        if(!(baseSeqIdx
             = newBlockEncIdxSeqFromSA(sa, totalLen, projectName,
                                       extraParams->blockEncParams.blockSize,
                                       extraParams->blockEncParams.bucketBlocks,
                                       0, NULL, NULL, NULL, NULL,
                                       NULL, 0, 0, &varState, env)))
          newBWTSeqErrRet();
      }
      destructAddLocateInfoState(&varState, env);
    }
    else
    {
      /* FIXME: read sample interval here, to make sure, the read
       * index has locate information */
      FILE *fp = EISSeekToHeader(baseSeqIdx, LOCATE_INFO_IN_INDEX_HEADERID,
                                 NULL);
      fputs("Using pre-computed sequence index.\n", stderr);
      if(!fp)
      {
        fputs("Specified index does not contain locate information.\n"
              "Localization of matches will not be supported!\n", stderr);
        locateInterval = 0;
      }
      else
      {
        struct locateHeader header;
        if(!readLocateInfoHeader(fp, &header))
          newBWTSeqErrRet();
        locateInterval = header.locateInterval;
      }
    }
    alphabet = EISGetAlphabet(baseSeqIdx);
    alphabetSize = MRAEncGetSize(alphabet);
    if(!alphabetSize)
      /* weird error, shouldn't happen, but I prefer error return to
       * segfault in case someone tampered with the input */
      newBWTSeqErrRet();
    bwtSeq = env_ma_malloc(env, offsetAlign(sizeof(struct BWTSeq),
                                            sizeof(Seqpos))
                           + sizeof(Seqpos) * (alphabetSize + 1));
    bwtSeq->count = (Seqpos *)((char  *)bwtSeq
                               + offsetAlign(sizeof(struct BWTSeq),
                                             sizeof(Seqpos)));
    bwtSeq->type = baseType;
    bwtSeq->seqIdx = baseSeqIdx;
    bwtSeq->alphabetSize = alphabetSize;
    bwtSeq->locateSampleInterval = locateInterval;
    bwtSeq->hint = hint = newEISHint(baseSeqIdx, env);
    {
      Symbol i;
      Seqpos len = EISLength(baseSeqIdx), *count = bwtSeq->count;
      count[0] = 0;
      for(i = 0; i < alphabetSize; ++i)
        count[i + 1] = count[i]
          + EISSymTransformedRank(baseSeqIdx, i, len, hint, env);
      fprintf(stderr, "count[alphabetSize]="FormatSeqpos
              ", len="FormatSeqpos"\n", count[alphabetSize], len);
      for(i = 0; i <= alphabetSize; ++i)
        fprintf(stderr, "count[%u]="FormatSeqpos"\n", i, count[i]);
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
getMatchBound(const BWTSeq *bwtSeq, Symbol *query, size_t queryLen,
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
  while((match->upper <= match->lower) && (i > 0))
  {
    curSym = MRAEncMapSymbol(alphabet, query[--i]);
    match->upper = count[curSym]
      + BWTSeqOcc(bwtSeq, curSym, match->upper, env);
    match->lower = count[curSym]
      + BWTSeqOcc(bwtSeq, curSym, match->lower, env);
  }
}

extern Seqpos
BWTSeqMatchCount(const BWTSeq *bwtSeq, Symbol *query, size_t queryLen,
                 Env *env)
{
  struct matchBound match;
  assert(bwtSeq && query && env);
  getMatchBound(bwtSeq, query, queryLen, &match, env);
  if(match.lower < match.upper)
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
newEMIterator(const BWTSeq *bwtSeq, Symbol *query, size_t queryLen, Env *env)
{
  struct BWTSeqExactMatchesIterator *newIter;
  assert(bwtSeq && query && env);
  if(!bwtSeq->locateSampleInterval)
  {
    fputs("Index does not contain locate information.\n"
          "Localization of matches impossible!", stderr);
    return NULL;
  }
  newIter = env_ma_malloc(env, sizeof(*newIter));
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
  if(iter->bounds.upper > iter->bounds.lower)
    return 0;
  else
    return iter->bounds.lower - iter->bounds.upper;
}

extern Seqpos
EMINumMatchesLeft(const struct BWTSeqExactMatchesIterator *iter)
{
  assert(iter);
  if(iter->nextMatchBWTPos > iter->bounds.lower)
    return 0;
  else
    return iter->bounds.lower - iter->bounds.upper;
}


struct MatchData *
EMIGetNextMatch(struct BWTSeqExactMatchesIterator *iter, const BWTSeq *bwtSeq,
                Env *env)
{
  if(iter->nextMatchBWTPos < iter->bounds.lower)
  {
    Seqpos nextLocate = iter->nextMatchBWTPos, locateOffset = 0;
    while(!BWTSeqPosHasLocateInfo(bwtSeq, nextLocate, &iter->extBits, env))
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


