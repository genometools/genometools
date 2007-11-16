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
/**
 * \file blockcomp.c
 * \brief Methods to build block-compressed representation of indexed
 * sequence and answer queries on said representation.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */
/*
  TODO:
  - normalize use  of  seqIdx variable naming (seq, bseq etc.)
  - split init/new functionality cleanly for all structs
 */
#include <assert.h>
#include <stddef.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>

#include "libgtcore/bitpackstring.h"
#include "libgtcore/dataalign.h"
#include "libgtcore/env.h"
#include "libgtcore/minmax.h"
#include "libgtcore/str.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/seqpos-def.h"
#include "libgtmatch/esa-map.pr"

#include "libgtmatch/eis-bitpackseqpos.h"
#include "libgtmatch/eis-encidxseq.h"
#include "libgtmatch/eis-encidxseqconstruct.h"
#include "libgtmatch/eis-encidxseqpriv.h"
#include "libgtmatch/eis-seqranges.h"
#include "libgtmatch/eis-seqblocktranslate.h"
#include "libgtmatch/eis-suffixerator-interface.h"

/**
 * generic method to aquire next readLen symbols of BWT string
 */
typedef int (*bwtReadFunc)(void *state, size_t readLen, Symbol *dest,
                           Env *env);

struct encIdxSeq *
newBlockEncIdxSeq(const Str *projectName, const struct blockEncParams *params,
                  size_t numExtHeaders, uint16_t *headerIDs,
                  uint32_t *extHeaderSizes,
                  headerWriteFunc *extHeaderCallbacks, void **headerCBData,
                  bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                  BitOffset maxVarExtBitsPerPos, void *cbState, Env *env)
{
  Suffixarray suffixArray;
  struct encIdxSeq *newSeqIdx;
  Seqpos length;
  Verboseinfo *verbosity;
  assert(projectName);
  /* map and interpret index project file */
  /* FIXME: handle verbosity in a more sane fashion */
  verbosity = newverboseinfo(false, env);
  if (streamsuffixarray(&suffixArray, &length,
                       SARR_SUFTAB | SARR_BWTTAB, projectName, verbosity, env))
  {
    freeverboseinfo(&verbosity, env);
    return NULL;
  }
  ++length;
  newSeqIdx = newBlockEncIdxSeqFromSA(&suffixArray, length,
                                      projectName, params,
                                      numExtHeaders, headerIDs,
                                      extHeaderSizes, extHeaderCallbacks,
                                      headerCBData, biFunc, cwExtBitsPerPos,
                                      maxVarExtBitsPerPos, cbState, env);
  freesuffixarray(&suffixArray, env);
  freeverboseinfo(&verbosity, env);
  return newSeqIdx;
}

static EISeq *
newGenBlockEncIdxSeq(Seqpos totalLen, const Str *projectName,
                     MRAEnc *alphabet, const struct seqStats *stats,
                     bwtReadFunc getNextBlock, void *bwtReadFuncState,
                     const struct blockEncParams *params,
                     size_t numExtHeaders, uint16_t *headerIDs,
                     uint32_t *extHeaderSizes,
                     headerWriteFunc *extHeaderCallbacks,
                     void **headerCBData,
                     bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                     BitOffset maxVarExtBitsPerPos, void *cbState, Env *env);

struct fileReadState
{
  FILE *fp;
  MRAEnc *alphabet;
};

static int
saBWTReadNext(void *state, size_t readLen, Symbol *dest, Env *env)
{
  struct fileReadState *FRState = state;
  assert(state);
  return MRAEncReadAndTransform(FRState->alphabet, FRState->fp,
                                readLen, dest);
}

extern EISeq *
newBlockEncIdxSeqFromSA(Suffixarray *sa, Seqpos totalLen,
                        const Str *projectName,
                        const struct blockEncParams *params,
                        size_t numExtHeaders, uint16_t *headerIDs,
                        uint32_t *extHeaderSizes,
                        headerWriteFunc *extHeaderCallbacks,
                        void **headerCBData,
                        bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                        BitOffset maxVarExtBitsPerPos, void *cbState, Env *env)
{
  struct encIdxSeq *newSeqIdx;
  struct fileReadState state;
  assert(sa && projectName && env);
  if (!sa->bwttabstream.fp)
  {
    fprintf(stderr, "error: bwt data not available for given project: %s",
            str_get(projectName));
    return NULL;
  }
  /* 1. get bwttab file pointer */
  state.fp = sa->bwttabstream.fp;
  /* convert alphabet */
  state.alphabet = MRAEncGTAlphaNew(sa->alpha, env);
  MRAEncAddSymbolToRange(state.alphabet, SEPARATOR, 1);
  newSeqIdx = newGenBlockEncIdxSeq(totalLen, projectName, state.alphabet,
                                   NULL, saBWTReadNext, &state, params,
                                   numExtHeaders, headerIDs, extHeaderSizes,
                                   extHeaderCallbacks, headerCBData, biFunc,
                                   cwExtBitsPerPos, maxVarExtBitsPerPos,
                                   cbState, env);
  if (!newSeqIdx)
    MRAEncDelete(state.alphabet, env);
  return newSeqIdx;
}

struct sfxInterfaceReadState
{
  sfxInterface *si;
  listenerID bwtReadId;
  MRAEnc *alphabet;
};

static int
sfxIBWTReadNext(void *state, size_t readLen, Symbol *dest, Env *env)
{
  struct sfxInterfaceReadState *SIState = state;
  size_t readCount;
  assert(state && SIState->si);
  readCount = readSfxIBWTRangeSym(SIState->si, SIState->bwtReadId, readLen,
                                  dest, env);
  MRAEncSymbolsTransform(SIState->alphabet, dest, readCount);
  return (readCount==readLen?1:-1);
}

extern EISeq *
newBlockEncIdxSeqFromSfxI(sfxInterface *si, Seqpos totalLen,
                          const Str *projectName,
                          const struct blockEncParams *params,
                          size_t numExtHeaders, uint16_t *headerIDs,
                          uint32_t *extHeaderSizes,
                          headerWriteFunc *extHeaderCallbacks,
                          void **headerCBData, bitInsertFunc biFunc,
                          BitOffset cwExtBitsPerPos,
                          BitOffset maxVarExtBitsPerPos, void *cbState,
                          Env *env)
{
  struct encIdxSeq *newSeqIdx;
  struct sfxInterfaceReadState state;
  assert(si && projectName && env);
  /* register bwttab reader */
  if (!SfxIRegisterReader(si, &state.bwtReadId, SFX_REQUEST_BWTTAB, env))
    return NULL;
  state.si = si;
  /* convert alphabet */
  state.alphabet = MRAEncGTAlphaNew(getSfxIAlphabet(si), env);
  MRAEncAddSymbolToRange(state.alphabet, SEPARATOR, 1);
  newSeqIdx = newGenBlockEncIdxSeq(totalLen, projectName,
                                   state.alphabet, getSfxISeqStats(si),
                                   sfxIBWTReadNext, &state, params,
                                   numExtHeaders, headerIDs, extHeaderSizes,
                                   extHeaderCallbacks, headerCBData, biFunc,
                                   cwExtBitsPerPos, maxVarExtBitsPerPos,
                                   cbState, env);
  if (!newSeqIdx)
    MRAEncDelete(state.alphabet, env);
  return newSeqIdx;
}

struct onDiskBlockCompIdx
{
  FILE *idxFP;                 /**< used for writing the index and
                                * if an mmap of the whole index
                                * isn't possible also for reading */
  char *idxMMap;               /**< if mmapping of the constant- and
                                * variable-width strings succeeded,
                                * stores the base of the
                                * corresponding memory mapped file
                                * portion */
  off_t cwDataPos,             /**< constant width data of the index */
    varDataPos,                /**< variable width part */
    rangeEncPos;               /**< in-file position of special
                                *  symbol representation */
};

/**
 * mode in which to encode a range of the alphabet
 */
enum rangeStoreMode {
  DIRECT_SYM_ENCODE,          /**< direct sequence encoding (not implemented) */
  BLOCK_COMPOSITION_INCLUDE,  /**< encode as two bitstrings, one
                               * containing a sequence of constant
                               * width composition indices, one
                               * containing a sequence of variable
                               * width permutation indices */
  REGIONS_LIST,               /**< encode as linear list of sequence regions */
};

struct blockCompositionSeq
{
  struct encIdxSeq baseClass;
  struct onDiskBlockCompIdx externalData;
  struct compList compositionTable;
  struct seqRangeList *rangeEncs;
  struct extHeaderPos *extHeaderPos;
  size_t numExtHeaders;
  BitOffset maxVarExtBitsPerBucket, cwExtBitsPerBucket;
  /* maxVarBitsPerBucket, cwBits */
  MRAEnc *blockMapAlphabet, *rangeMapAlphabet;
  int *modes;
  unsigned bucketBlocks, blockSize, callBackDataOffsetBits, bitsPerSeqpos,
    bitsPerVarDiskOffset;
  Symbol blockEncNumSyms, blockEncFallback, rangeEncFallback;
  int numModes;
  unsigned *partialSymSumBits, *partialSymSumBitsSums, symSumBits;
};

static inline size_t
blockEncIdxSeqHeaderLength(struct blockCompositionSeq *seqIdx,
                           size_t numExtHeaders, uint32_t *extHeaderSizes);

static int
openOnDiskData(const Str *projectName, struct onDiskBlockCompIdx *idx,
               char *mode, Env *env);

static void
initOnDiskBlockCompIdx(struct onDiskBlockCompIdx *idx,
                       size_t headerLen, off_t cwLen);

static void
destructOnDiskBlockCompIdx(struct onDiskBlockCompIdx *idx, Env *env);

static size_t
writeIdxHeader(struct blockCompositionSeq *seqIdx,
               size_t numExtHeaders, uint16_t *headerIDs,
               uint32_t *extHeaderSizes, headerWriteFunc *extHeaderCallbacks,
               void **headerCBData,
               Env *env);

static inline int
tryMMapOfIndex(struct onDiskBlockCompIdx *idxData);

static const struct encIdxSeqClass blockCompositionSeqClass;

static inline struct blockCompositionSeq *
encIdxSeq2blockCompositionSeq(struct encIdxSeq *seq)
{
  assert(seq && seq->classInfo == &blockCompositionSeqClass);
  return (struct blockCompositionSeq *)((char *)seq
    - offsetof(struct blockCompositionSeq, baseClass));
}

static inline struct encIdxSeq *
blockCompositionSeq2encIdxSeq(struct blockCompositionSeq *seq)
{
  assert(seq);
  return &(seq->baseClass);
}

static inline const struct blockCompositionSeq *
constEncIdxSeq2blockCompositionSeq(const struct encIdxSeq *seq)
{
  assert(seq && seq->classInfo == &blockCompositionSeqClass);
  return (struct blockCompositionSeq *)((char *)seq
    - offsetof(struct blockCompositionSeq, baseClass));
}

static inline const struct encIdxSeq *
constBlockCompositionSeq2encIdxSeq(const struct blockCompositionSeq *seq)
{
  assert(seq);
  return &(seq->baseClass);
}

struct appendState
{
  BitString compCache, permCache;
  BitOffset compCacheLen, permCacheLen, cwMemPos, varMemPos, varDiskOffset;
  off_t cwDiskOffset;
  unsigned varMemOldBits, cwMemOldBits;
};

static inline void
initAppendState(struct appendState *aState,
                const struct blockCompositionSeq *seqIdx, Env *env);
static void
destructAppendState(struct appendState *aState, Env *env);

static inline void
append2IdxOutput(struct blockCompositionSeq *newSeqIdx,
                 struct appendState *state,
                 PermCompIndex permCompIdx[2],
                 unsigned bitsOfCompositionIdx, unsigned bitsOfPermutationIdx);

static BitOffset
appendCallBackOutput(struct appendState *state,
                     const struct blockCompositionSeq *seqIdx,
                     bitInsertFunc biFunc, Seqpos start, Seqpos len,
                     unsigned callBackDataOffsetBits, void *cbState, Env *env);

typedef Seqpos *partialSymSums;

static inline partialSymSums
newPartialSymSums(unsigned alphabetSize, Env *env);

static inline void
deletePartialSymSums(partialSymSums sums, Env *env);

static inline void
addBlock2PartialSymSums(partialSymSums sums, Symbol *block, unsigned blockSize);

static inline void
copyPartialSymSums(unsigned alphabetSize, partialSymSums dest,
                   const partialSymSums src);

static inline Seqpos
numBuckets(Seqpos seqLen, Seqpos bucketLen);

static inline off_t
cwSize(const struct blockCompositionSeq *seqIdx);

static inline BitOffset
vwBits(Seqpos seqLen, unsigned blockSize, unsigned bucketBlocks,
       unsigned maxPermIdxBits, BitOffset maxVarExtBitsPerBucket);

static void
addRangeEncodedSyms(struct seqRangeList *rangeList, Symbol *block,
                    unsigned blockSize, Seqpos blockNum, MRAEnc *alphabet,
                    int selection, int *rangeSel, Env *env);

static int
updateIdxOutput(struct blockCompositionSeq *seqIdx,
                struct appendState *aState,
                partialSymSums buck);

static int
finalizeIdxOutput(struct blockCompositionSeq *seqIdx,
                  struct appendState *state,
                  partialSymSums buck, Env *env);

static inline void
symSumBitsDefaultSetup(struct blockCompositionSeq *seqIdx);

#define newBlockEncIdxSeqErrRet()                                       \
  do {                                                                  \
    if (newSeqIdx->externalData.idxFP)                                  \
      destructOnDiskBlockCompIdx(&newSeqIdx->externalData, env);        \
    if (newSeqIdx->compositionTable.bitsPerCount)                       \
      destructCompositionList(&newSeqIdx->compositionTable, env);       \
    if (newSeqIdx->rangeEncs)                                           \
      deleteSeqRangeList(newSeqIdx->rangeEncs, env);                    \
    if (newSeqIdx->extHeaderPos)                                        \
      env_ma_free(newSeqIdx->extHeaderPos, env);                        \
    if (newSeqIdx->partialSymSumBits)                                   \
      env_ma_free(newSeqIdx->partialSymSumBits, env);                   \
    if (newSeqIdx->extHeaderPos)                                        \
      env_ma_free(newSeqIdx->extHeaderPos, env);                        \
    if (newSeqIdx) env_ma_free(newSeqIdx, env);                         \
    if (modesCopy)                                                      \
      env_ma_free(modesCopy, env);                                      \
    if (blockMapAlphabet) MRAEncDelete(blockMapAlphabet, env);          \
    if (rangeMapAlphabet) MRAEncDelete(rangeMapAlphabet, env);          \
    return NULL;                                                        \
  } while (0)

#define newBlockEncIdxSeqLoopErr()                      \
  destructAppendState(&aState, env);                    \
  deletePartialSymSums(buck, env);                      \
  deletePartialSymSums(buckLast, env);                  \
  env_ma_free(compositionPreAlloc, env);                \
  env_ma_free(permCompBSPreAlloc, env);                 \
  env_ma_free(block, env);                              \
  break

/**
 * @param alphabet ownership of alphabet is transferred to the sequence
 * index produced unless NULL is returned
 */
static EISeq *
newGenBlockEncIdxSeq(Seqpos totalLen, const Str *projectName,
                     MRAEnc *alphabet, const struct seqStats *stats,
                     bwtReadFunc getNextBlock, void *bwtReadFuncState,
                     const struct blockEncParams *params,
                     size_t numExtHeaders, uint16_t *headerIDs,
                     uint32_t *extHeaderSizes,
                     headerWriteFunc *extHeaderCallbacks,
                     void **headerCBData,
                     bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                     BitOffset maxVarExtBitsPerPos, void *cbState, Env *env)
{
  struct blockCompositionSeq *newSeqIdx = NULL;
  Symbol blockMapAlphabetSize, totalAlphabetSize;
  size_t regionsEstimate=totalLen / 100;
  MRAEnc *blockMapAlphabet = NULL, *rangeMapAlphabet = NULL;
  BitOffset bitsPerComposition, bitsPerPermutation;
  unsigned compositionIdxBits, callBackDataOffsetBits,
    blockSize = params->blockSize, bucketBlocks = params->bucketBlocks;
  size_t bucketLen = (size_t)bucketBlocks * blockSize;
  int *modesCopy = NULL;
  enum rangeStoreMode modes[] = { BLOCK_COMPOSITION_INCLUDE,
                                  REGIONS_LIST };
  assert(projectName);
  assert(blockSize > 0);
  env_error_check(env);

  newSeqIdx = env_ma_calloc(env, sizeof (struct blockCompositionSeq), 1);
  newSeqIdx->bucketBlocks = bucketBlocks;
  newSeqIdx->bitsPerSeqpos = requiredSeqposBits(newSeqIdx->baseClass.seqLen
                                                = totalLen);
  newSeqIdx->baseClass.alphabet = alphabet;
  /* TODO: improve guessing of number of necessary ranges */
  {
    size_t range, numAlphabetRanges = newSeqIdx->numModes =
      MRAEncGetNumRanges(alphabet);
    totalAlphabetSize = MRAEncGetSize(alphabet);
    blockMapAlphabetSize = 0;
    newSeqIdx->modes = modesCopy =
      env_ma_malloc(env, sizeof (int) * numAlphabetRanges);
    for (range = 0; range < numAlphabetRanges; ++range)
    {
      modesCopy[range] = modes[range];
      switch (modes[range])
      {
      case BLOCK_COMPOSITION_INCLUDE:
        blockMapAlphabetSize += MRAEncGetRangeSize(alphabet, range);
        break;
      case DIRECT_SYM_ENCODE:
        /*< OUTLOOK: insert proper code to process ranges
         * to improve orthogonality */
        break;
      case REGIONS_LIST:
        /*< OUTLOOK: insert proper code to process ranges
         * to improve orthogonality */
        break;
      default:
        /* TODO: improve diagnostics */
        fprintf(stderr, "Invalid encoding request.\n");
        newBlockEncIdxSeqErrRet();
        break;
      }
    }
    newSeqIdx->blockEncNumSyms = blockMapAlphabetSize;
    newSeqIdx->blockMapAlphabet = blockMapAlphabet =
      MRAEncSecondaryMapping(alphabet, BLOCK_COMPOSITION_INCLUDE, modesCopy,
                             newSeqIdx->blockEncFallback = 0, env);
    newSeqIdx->rangeMapAlphabet = rangeMapAlphabet =
      MRAEncSecondaryMapping(alphabet, REGIONS_LIST, modesCopy,
                             newSeqIdx->rangeEncFallback = 0, env);
    assert(MRAEncGetSize(blockMapAlphabet) == blockMapAlphabetSize);
    assert(MRAEncGetSize(rangeMapAlphabet)
           == totalAlphabetSize - blockMapAlphabetSize);
  }
  newSeqIdx->partialSymSumBits
    = env_ma_malloc(env, sizeof (newSeqIdx->partialSymSumBits[0])
                    * blockMapAlphabetSize * 2);
  newSeqIdx->partialSymSumBitsSums
    = newSeqIdx->partialSymSumBits + blockMapAlphabetSize;
  if (stats)
  {
    Seqpos *symCounts;
    symCounts = env_ma_malloc(env, sizeof (symCounts[0])
                              * blockMapAlphabetSize);
    switch (stats->sourceAlphaType)
    {
    case sourceUInt8:
      memset(symCounts, 0,
             sizeof (newSeqIdx->partialSymSumBits[0]) * blockMapAlphabetSize);
      {
        Symbol eSym, bSym;
        unsigned i;
        for (i = 0; i <= UINT8_MAX; ++i)
          if (MRAEncSymbolIsInSelectedRanges(
               alphabet, eSym = MRAEncMapSymbol(alphabet, i),
               BLOCK_COMPOSITION_INCLUDE, modesCopy)
             && ((bSym = MRAEncMapSymbol(blockMapAlphabet, eSym))
                 < blockMapAlphabetSize))
            symCounts[bSym]
              += stats->symbolDistributionTable[i];
#if DEBUG > 1
        for (i = 0; i < blockMapAlphabetSize; ++i)
        {
          fprintf(stderr, "symCount[%u]="FormatSeqpos"\n", (unsigned)i,
                  stats->symbolDistributionTable[i]);
        }
#endif /* DEBUG > 1 */
        if (blockMapAlphabetSize)
        {
          newSeqIdx->partialSymSumBitsSums[0] = 0;
          newSeqIdx->partialSymSumBits[0] = requiredSeqposBits(symCounts[0]);
          for (i = 1; i < blockMapAlphabetSize; ++i)
          {
            newSeqIdx->partialSymSumBitsSums[i]
              = newSeqIdx->partialSymSumBitsSums[i - 1]
              + newSeqIdx->partialSymSumBits[i - 1];
            newSeqIdx->partialSymSumBits[i]
              = requiredSeqposBits(symCounts[i]);
          }
#ifdef DEBUG
          for (i = 0; i < blockMapAlphabetSize; ++i)
          {
            fprintf(stderr, "bitsPerSymSum[%u]=%u\n", (unsigned)i,
                    newSeqIdx->partialSymSumBits[i]);
          }
#endif  /* DEBUG */
          newSeqIdx->symSumBits
            = newSeqIdx->partialSymSumBitsSums[blockMapAlphabetSize - 1]
            + newSeqIdx->partialSymSumBits[blockMapAlphabetSize - 1];
#ifdef DEBUG
          fprintf(stderr, "symSumBits total: %u\n", newSeqIdx->symSumBits);
#endif  /* DEBUG */
        }
      }
      /* count special characters to estimate number of regions required */
      {
        Symbol eSym, bSym;
        Seqpos regionSymCount = 0;
        unsigned i;
        for (i = 0; i <= UINT8_MAX; ++i)
          if (MRAEncSymbolIsInSelectedRanges(
               alphabet, eSym = MRAEncMapSymbol(alphabet, i),
               REGIONS_LIST, modesCopy)
             && ((bSym = MRAEncMapSymbol(blockMapAlphabet, eSym))
                 < blockMapAlphabetSize))
            regionSymCount += stats->symbolDistributionTable[i];
        regionsEstimate = regionSymCount/20;
#ifdef DEBUG
        fprintf(stderr, "Expected "FormatSeqpos
                " symbols to encode in regions.\n",
                regionSymCount);
#endif
      }
      break;
    default:
      symSumBitsDefaultSetup(newSeqIdx);
    }
    env_ma_free(symCounts, env);
  }
  else
  {
    symSumBitsDefaultSetup(newSeqIdx);
  }
  {
    int regionFeatures = SRL_NO_FEATURES;
    if (params->EISFeatureSet & EIS_FEATURE_REGION_SUMS)
      regionFeatures |= SRL_PARTIAL_SYMBOL_SUMS;
    newSeqIdx->rangeEncs = newSeqRangeList(regionsEstimate, rangeMapAlphabet,
                                           regionFeatures, env);
  }
  newSeqIdx->baseClass.classInfo = &blockCompositionSeqClass;
  newSeqIdx->blockSize = blockSize;
  newSeqIdx->cwExtBitsPerBucket = cwExtBitsPerPos * bucketLen;
  newSeqIdx->maxVarExtBitsPerBucket = maxVarExtBitsPerPos * bucketLen;
  if (!initCompositionList(&newSeqIdx->compositionTable, blockSize,
                           blockMapAlphabetSize, env))
    newBlockEncIdxSeqErrRet();
  bitsPerComposition = newSeqIdx->compositionTable.bitsPerCount
    * blockMapAlphabetSize;
  compositionIdxBits = newSeqIdx->compositionTable.compositionIdxBits;
  bitsPerPermutation = newSeqIdx->compositionTable.bitsPerSymbol * blockSize;
  if (biFunc)
    newSeqIdx->callBackDataOffsetBits = callBackDataOffsetBits
      = requiredUInt64Bits(newSeqIdx->compositionTable.maxPermIdxBits
                           * bucketBlocks);
  else
    newSeqIdx->callBackDataOffsetBits = callBackDataOffsetBits = 0;
  newSeqIdx->bitsPerVarDiskOffset =
    requiredUInt64Bits(vwBits(totalLen, blockSize, bucketBlocks,
                              newSeqIdx->compositionTable.maxPermIdxBits,
                              newSeqIdx->maxVarExtBitsPerBucket));
  {
    size_t headerLen = blockEncIdxSeqHeaderLength(newSeqIdx, numExtHeaders,
                                                  extHeaderSizes);
    if (!openOnDiskData(projectName, &newSeqIdx->externalData, "wb+", env))
      newBlockEncIdxSeqErrRet();
    initOnDiskBlockCompIdx(&newSeqIdx->externalData,
                           headerLen, cwSize(newSeqIdx));
  }
  /* At this point everything should be ready to receive the actual sequence
   * information, steps: */
  {
    int hadError = 0;
    do
    {
      {
        Symbol *block;
        unsigned *compositionPreAlloc;
        BitString permCompBSPreAlloc;
        partialSymSums buck, buckLast;
        block = env_ma_malloc(env, sizeof (Symbol) * blockSize);
        compositionPreAlloc = env_ma_malloc(env, sizeof (unsigned)
                                            * blockMapAlphabetSize);
        permCompBSPreAlloc =
          env_ma_malloc(env, bitElemsAllocSize(bitsPerComposition
                                               + bitsPerPermutation)
                        * sizeof (BitElem));
        buck = newPartialSymSums(totalAlphabetSize, env);
        buckLast = newPartialSymSums(totalAlphabetSize, env);
        /* 2. read block sized chunks from bwttab and suffix array */
        {
          Seqpos numFullBlocks = totalLen / blockSize, blockNum,
            lastUpdatePos = 0;
          /* pos == totalLen - symbolsLeft */
          struct appendState aState;
          initAppendState(&aState, newSeqIdx, env);
          blockNum = 0;
          while (blockNum < numFullBlocks)
          {
            PermCompIndex permCompIdx[2];
            unsigned significantPermIdxBits;
            int readResult;
            /* 3. for each chunk: */
            readResult = getNextBlock(bwtReadFuncState, blockSize, block, env);
            if (readResult < 0)
            {
              hadError = 1;
              perror("error condition while reading index data");
              break;
            }
            /* 3.a. update superbucket table */
            addBlock2PartialSymSums(buck, block, blockSize);
            /* 3.b. add ranges of differently encoded symbols to
             * corresponding representation */
            addRangeEncodedSyms(newSeqIdx->rangeEncs, block, blockSize,
                                blockNum, alphabet, REGIONS_LIST,
                                modesCopy, env);
            /* 3.c. add to table of composition/permutation indices */
            MRAEncSymbolsTransform(blockMapAlphabet, block, blockSize);
            /* FIXME control remapping */
            block2IndexPair(&newSeqIdx->compositionTable, blockSize,
                            blockMapAlphabetSize, block, permCompIdx,
                            &significantPermIdxBits, env,
                            permCompBSPreAlloc, compositionPreAlloc);
            append2IdxOutput(newSeqIdx, &aState,
                             permCompIdx, compositionIdxBits,
                             significantPermIdxBits);
            /* update on-disk structure */
            if (!((++blockNum) % bucketBlocks))
            {
              Seqpos pos = blockNum * blockSize;
              if (biFunc)
                if (appendCallBackOutput(&aState, newSeqIdx, biFunc,
                                        lastUpdatePos, bucketLen,
                                        callBackDataOffsetBits, cbState, env)
                   == (BitOffset)-1)
                {
                  hadError = 1;
                  perror("error condition while writing block-compressed"
                         " index data");
                  break;
                }
              if (!updateIdxOutput(newSeqIdx, &aState, buckLast))
              {
                hadError = 1;
                perror("error condition while writing block-compressed"
                       " index data");
                break;
              }
              /* update retained data */
              copyPartialSymSums(totalAlphabetSize, buckLast, buck);
              lastUpdatePos = pos;
            }
          }
          /* handle last chunk */
          if (!hadError)
          {
            PermCompIndex permCompIdx[2];
            unsigned significantPermIdxBits;
            int readResult;
            Seqpos symbolsLeft = totalLen % blockSize;
            readResult = getNextBlock(bwtReadFuncState, symbolsLeft, block,
                                      env);
            if (readResult < 0)
            {
              hadError = 1;
              perror("error condition while reading index data");
              newBlockEncIdxSeqLoopErr();
            }
            else
            {
              memset(block + symbolsLeft, 0,
                     sizeof (Symbol) * (blockSize - symbolsLeft));
              addBlock2PartialSymSums(buck, block, blockSize);
              addRangeEncodedSyms(newSeqIdx->rangeEncs, block, blockSize,
                                  blockNum, alphabet, REGIONS_LIST,
                                  modesCopy, env);
              MRAEncSymbolsTransform(blockMapAlphabet, block, blockSize);
              block2IndexPair(&newSeqIdx->compositionTable, blockSize,
                              blockMapAlphabetSize, block, permCompIdx,
                              &significantPermIdxBits, env,
                              permCompBSPreAlloc, compositionPreAlloc);
              append2IdxOutput(newSeqIdx, &aState, permCompIdx,
                               compositionIdxBits, significantPermIdxBits);
              if (biFunc)
                if (appendCallBackOutput(&aState, newSeqIdx, biFunc,
                                        lastUpdatePos, totalLen - lastUpdatePos,
                                        callBackDataOffsetBits, cbState, env)
                   == (BitOffset)-1)
                {
                  hadError = 1;
                  perror("error condition while writing block-compressed"
                         " index data");
                  break;
                }
              if (!updateIdxOutput(newSeqIdx, &aState, buckLast))
              {
                hadError = 1;
                perror("error condition while writing block-compressed"
                       " index data");
                newBlockEncIdxSeqLoopErr();
              }
              if (!finalizeIdxOutput(newSeqIdx, &aState, buck, env))
              {
                hadError = 1;
                perror("error condition while writing block-compressed"
                       " index data");
                newBlockEncIdxSeqLoopErr();
              }
              if (!writeIdxHeader(newSeqIdx, numExtHeaders, headerIDs,
                                  extHeaderSizes, extHeaderCallbacks,
                                  headerCBData, env))
              {
                hadError = 1;
                perror("error condition while writing block-compressed"
                       " index header");
                newBlockEncIdxSeqLoopErr();
              }
              if (fflush(newSeqIdx->externalData.idxFP))
              {
                hadError = 1;
                perror("error condition while writing block-compressed"
                       " index header");
                newBlockEncIdxSeqLoopErr();
              }
              tryMMapOfIndex(&newSeqIdx->externalData);
            }
          }
          if (hadError)
          {
            newBlockEncIdxSeqLoopErr();
          }
          /* 4. dealloc resources no longer required */
          destructAppendState(&aState, env);
          deletePartialSymSums(buck, env);
          deletePartialSymSums(buckLast, env);
        }
        /* 4. dealloc resources no longer required */
        env_ma_free(compositionPreAlloc, env);
        env_ma_free(permCompBSPreAlloc, env);
        env_ma_free(block, env);
      }
    } while (0);
    /* close bwttab and suffix array */
    if (hadError)
      newBlockEncIdxSeqErrRet();
  }
  return &(newSeqIdx->baseClass);
}

static void
deleteBlockEncIdxSeq(struct encIdxSeq *seq, Env *env)
{
  struct blockCompositionSeq *bseq;
  assert(seq && seq->classInfo == &blockCompositionSeqClass);
  bseq = encIdxSeq2blockCompositionSeq(seq);
  env_ma_free(bseq->extHeaderPos, env);
  env_ma_free(bseq->partialSymSumBits, env);
  destructOnDiskBlockCompIdx(&bseq->externalData, env);
  destructCompositionList(&bseq->compositionTable, env);
  MRAEncDelete(bseq->baseClass.alphabet, env);
  MRAEncDelete(bseq->rangeMapAlphabet, env);
  MRAEncDelete(bseq->blockMapAlphabet, env);
  deleteSeqRangeList(bseq->rangeEncs, env);
  env_ma_free(bseq->modes, env);
  env_ma_free(bseq, env);
}

#define USE_SBLOCK_CACHE

static inline int
seqIdxUsesMMap(const struct blockCompositionSeq *seqIdx)
{
  return seqIdx->externalData.idxMMap != NULL;
}

/**
 * Since the sequence represented in the form of blocks is incomplete
 * without the information in corresponding regions, a superBlock
 * struct stores only information for the strings of composition and
 * permutation indices. The cwData string is prepended by partial
 * symbol sums and offsets in the global variable width string. Also
 * appended to the indices is extra information inserted by the
 * callback functions on index creation.
 */
struct superBlock
{
  unsigned varDataMemBase, cwIdxMemBase;
  BitString varData, cwData;
};

static inline size_t
superBlockCWMaxReadSize(const struct blockCompositionSeq *seqIdx);

static inline size_t
superBlockVarMaxReadSize(const struct blockCompositionSeq *seqIdx);

static inline BitOffset
superBlockCWBits(const struct blockCompositionSeq *seqIdx);

static inline size_t
superBlockMemSize(const struct blockCompositionSeq *seqIdx)
{
  if (seqIdxUsesMMap(seqIdx))
  {
    return sizeof (struct superBlock);
  }
  else
  {
    size_t offset;
    offset = offsetAlign(sizeof (struct superBlock), sizeof (BitElem));
    offset = offsetAlign(offset + superBlockCWMaxReadSize(seqIdx),
                         sizeof (BitElem));
    offset += superBlockVarMaxReadSize(seqIdx);
    return offsetAlign(offset, MAX_ALIGN_REQUIREMENT);
  }
}

static inline void
initEmptySuperBlock(struct superBlock *sBlock,
                    const struct blockCompositionSeq *seqIdx)
{
  if (!seqIdxUsesMMap(seqIdx))
  {
    size_t offset;
    sBlock->cwData = (BitString)
      ((char *)sBlock + (offset = offsetAlign(sizeof (*sBlock),
                                              sizeof (BitElem))));
    sBlock->varData = (BitString)
      ((char *)sBlock + (offset = offsetAlign(offset
                                              + superBlockCWMaxReadSize(seqIdx),
                                              sizeof (BitElem))));
  }
}

static struct superBlock *
newEmptySuperBlock(const struct blockCompositionSeq *seqIdx, Env *env)
{
  struct superBlock *sBlock;
  sBlock = env_ma_malloc(env, superBlockMemSize(seqIdx));
  initEmptySuperBlock(sBlock, seqIdx);
  return sBlock;
}

static void
deleteSuperBlock(struct superBlock *sBlock, Env *env)
{
  env_ma_free(sBlock, env);
}

static inline void
symSumBitsDefaultSetup(struct blockCompositionSeq *seqIdx)
{
  unsigned i;
  unsigned blockMapAlphabetSize = seqIdx->blockEncNumSyms;
  seqIdx->partialSymSumBitsSums[0] = 0;
  seqIdx->partialSymSumBits[0] = seqIdx->bitsPerSeqpos;
  for (i = 1; i < blockMapAlphabetSize; ++i)
    seqIdx->partialSymSumBitsSums[i] = seqIdx->partialSymSumBitsSums[i - 1]
      + (seqIdx->partialSymSumBits[i] = seqIdx->bitsPerSeqpos);
  seqIdx->symSumBits = blockMapAlphabetSize * seqIdx->bitsPerSeqpos;
#ifdef DEBUG
  fprintf(stderr, "symSumBits=%u, blockEncNumSyms=%u\n",
          seqIdx->symSumBits, seqIdx->blockEncNumSyms);
#endif
  assert(seqIdx->partialSymSumBitsSums[i - 1] + seqIdx->bitsPerSeqpos
         == seqIdx->symSumBits);
}

static inline Seqpos
sBlockGetPartialSymSum(struct superBlock *sBlock, Symbol sym,
                       const struct blockCompositionSeq *seqIdx)
{
  return bsGetSeqpos(sBlock->cwData, seqIdx->partialSymSumBitsSums[sym]
                     + sBlock->cwIdxMemBase, seqIdx->partialSymSumBits[sym]);
}

static inline BitOffset
cwPreVarIdxBits(const struct blockCompositionSeq *seqIdx)
{
  return seqIdx->symSumBits;
}

static inline BitOffset
sBlockGetVarIdxOffset(const struct superBlock *sBlock,
                      const struct blockCompositionSeq *seqIdx)
{
  BitOffset offset = cwPreVarIdxBits(seqIdx) + sBlock->cwIdxMemBase;
  return bsGetUInt64(sBlock->cwData, offset, seqIdx->bitsPerVarDiskOffset);
}

static inline BitOffset
cwPreCBOffsetBits(const struct blockCompositionSeq *seqIdx)
{
  BitOffset offset = cwPreVarIdxBits(seqIdx) + seqIdx->bitsPerVarDiskOffset;
  return offset;
}

static inline BitOffset
cwPreCompIdxBits(const struct blockCompositionSeq *seqIdx)
{
  BitOffset offset = cwPreCBOffsetBits(seqIdx) + seqIdx->callBackDataOffsetBits;
  return offset;
}

static inline BitOffset
sBlockGetCompIdxOffset(const struct superBlock *sBlock,
                       const struct blockCompositionSeq *seqIdx,
                       unsigned compIdxNum)
{
  BitOffset offset = sBlock->cwIdxMemBase + cwPreCompIdxBits(seqIdx)
    + compIdxNum * seqIdx->compositionTable.compositionIdxBits;
  return offset;
}

static inline BitOffset
cwPreCWExtBits(const struct blockCompositionSeq *seqIdx)
{
  BitOffset offset = cwPreCompIdxBits(seqIdx)
    + seqIdx->bucketBlocks * seqIdx->compositionTable.compositionIdxBits;
  return offset;
}

static inline BitOffset
sBlockCWExtBitsOffset(const struct superBlock *sBlock,
                      const struct blockCompositionSeq *seqIdx)
{
  BitOffset offset = cwPreCWExtBits(seqIdx) + sBlock->cwIdxMemBase;
  return offset;
}

static inline BitOffset
sBlockGetcbOffsetOffset(const struct superBlock *sBlock,
                        const struct blockCompositionSeq *seqIdx)
{
  BitOffset offset = cwPreCBOffsetBits(seqIdx) + sBlock->cwIdxMemBase;
  return offset;
}

static inline BitOffset
sBlockGetcbOffset(const struct superBlock *sBlock,
                  const struct blockCompositionSeq *seqIdx)
{
  BitOffset offset = sBlockGetcbOffsetOffset(sBlock, seqIdx);
  return bsGetUInt64(sBlock->cwData, offset, seqIdx->callBackDataOffsetBits);
}

static inline PermCompIndex
sBlockGetCompIdx(const struct superBlock *sBlock, unsigned compIdxNum,
                 const struct blockCompositionSeq *seqIdx)
{
  BitOffset offset = sBlockGetCompIdxOffset(sBlock, seqIdx, compIdxNum);
  unsigned bitsPerCompositionIdx = seqIdx->compositionTable.compositionIdxBits;
  return bsGetPermCompIndex(sBlock->cwData, offset, bitsPerCompositionIdx);
}

static inline Seqpos
blockNumFromPos(const struct blockCompositionSeq *seqIdx, Seqpos pos)
{
  return pos / seqIdx->blockSize;
}

static inline Seqpos
bucketNumFromPos(const struct blockCompositionSeq *seqIdx, Seqpos pos)
{
  return pos / (seqIdx->bucketBlocks * seqIdx->blockSize);
}

static inline Seqpos
bucketBasePos(const struct blockCompositionSeq *seqIdx, Seqpos bucketNum)
{
  return bucketNum * seqIdx->bucketBlocks * seqIdx->blockSize;
}

static inline Seqpos
bucketNumFromBlockNum(struct blockCompositionSeq *seqIdx, Seqpos blockNum)
{
  return blockNum / seqIdx->bucketBlocks;
}

#define fetchSuperBlockErrRet()                         \
  do {                                                  \
    if (!sBlockPreAlloc) deleteSuperBlock(retval, env); \
    return NULL;                                        \
  } while (0)

static struct superBlock *
fetchSuperBlock(const struct blockCompositionSeq *seqIdx, Seqpos bucketNum,
                struct superBlock *sBlockPreAlloc, Env *env)
{
  struct superBlock *retval = NULL;
  assert(seqIdx);
  assert(bucketNum * seqIdx->bucketBlocks * seqIdx->blockSize
         < seqIdx->baseClass.seqLen);
  if (sBlockPreAlloc)
    retval = sBlockPreAlloc;
  else
    retval = newEmptySuperBlock(seqIdx, env);
  if (seqIdxUsesMMap(seqIdx))
  {
    BitOffset bucketOffset = bucketNum * superBlockCWBits(seqIdx);
    BitOffset varDataOffset;
    retval->cwData = (BitString)(seqIdx->externalData.idxMMap
      + bucketOffset / bitElemBits * sizeof (BitElem));
    retval->cwIdxMemBase = bucketOffset%bitElemBits;

    varDataOffset = sBlockGetVarIdxOffset(retval, seqIdx);
    retval->varData = (BitString)(seqIdx->externalData.idxMMap
      + seqIdx->externalData.varDataPos - seqIdx->externalData.cwDataPos
      + varDataOffset/bitElemBits * sizeof (BitElem));
    retval->varDataMemBase = varDataOffset%bitElemBits;
  }
  else
  {
    FILE *idxFP;
    size_t superBlockCWDiskSize = superBlockCWMaxReadSize(seqIdx);
    BitOffset bucketOffset = bucketNum * superBlockCWBits(seqIdx);
    BitOffset varDataOffset;
    idxFP = seqIdx->externalData.idxFP;
    if (fseeko(idxFP, seqIdx->externalData.cwDataPos
              + bucketOffset / bitElemBits * sizeof (BitElem), SEEK_SET))
      fetchSuperBlockErrRet();
    if (fread(retval->cwData, 1, superBlockCWDiskSize, idxFP)
       != superBlockCWDiskSize)
      fetchSuperBlockErrRet();
    retval->cwIdxMemBase = bucketOffset%bitElemBits;
    varDataOffset = sBlockGetVarIdxOffset(retval, seqIdx);
    if (fseeko(idxFP, seqIdx->externalData.varDataPos
              + varDataOffset/bitElemBits * sizeof (BitElem), SEEK_SET))
      fetchSuperBlockErrRet();
    retval->varDataMemBase = varDataOffset%bitElemBits;
    fread(retval->varData, sizeof (BitElem), superBlockVarMaxReadSize(seqIdx),
          idxFP);
    if (ferror(idxFP))
      fetchSuperBlockErrRet();
  }
  return retval;
}

static Seqpos
blockCompSeqSelect(struct encIdxSeq *seq, Symbol sym, Seqpos count,
                   union EISHint *hint, Env *env)
{
  /* FIXME: implementation pending */
  abort();
  return 0;
}

/*
 * routines for management of super-Block-Cache, this does currently
 * use a simple direct-mapped caching
 */

static inline size_t
pos2CachePos(struct seqCache *sCache, Seqpos pos)
{
  return pos % sCache->numEntries;
}

static void
initSuperBlockSeqCache(struct seqCache *sBlockCache,
                       struct blockCompositionSeq *seqIdx,
                       size_t numEntries, Env *env)
{
  unsigned bucketBlocks, bucketLen, blockSize;
  size_t superBlockSize;

  assert(seqIdx && sBlockCache);
  blockSize = seqIdx->blockSize;
  bucketLen = (bucketBlocks = seqIdx->bucketBlocks) * blockSize;
  superBlockSize = superBlockMemSize(seqIdx);
  sBlockCache->numEntries = numEntries;
  {
    void *temp = env_ma_malloc(env, (sizeof (Seqpos) + superBlockSize
                                     + sizeof (void *)) * numEntries);
    sBlockCache->cachedPos = temp;
    sBlockCache->entriesPtr = (void **)((char *)temp
                                        + sizeof (Seqpos) * numEntries);
    sBlockCache->entries = (char *)temp + (sizeof (Seqpos) + sizeof (void *))
      * numEntries;
  }
  {
    size_t i;
    for (i = 0; i < numEntries; ++i)
    {
      sBlockCache->entriesPtr[i] = (char *)sBlockCache->entries
        + superBlockSize * i;
      initEmptySuperBlock(sBlockCache->entriesPtr[i], seqIdx);
      sBlockCache->cachedPos[i] = -1;
    }
  }
}

static void
destructSuperBlockSeqCache(struct seqCache *sBlockCache, Env *env)
{
  env_ma_free(sBlockCache->cachedPos, env);
}

static inline int
inSeqCache(struct seqCache *sCache, Seqpos pos)
{
  return sCache->cachedPos[pos2CachePos(sCache, pos)] == pos;
}

#ifdef USE_SBLOCK_CACHE
static struct superBlock *
cacheFetchSuperBlock(const struct blockCompositionSeq *seqIdx,
                     Seqpos superBlockNum, struct seqCache *sBlockCache,
                     Env *env)
{
  size_t cachePos = pos2CachePos(sBlockCache, superBlockNum);
  if (inSeqCache(sBlockCache, superBlockNum))
    return (struct superBlock *)sBlockCache->entriesPtr[cachePos];
  else
  {
    struct superBlock *sb =
      (struct superBlock *)sBlockCache->entriesPtr[cachePos];
    return fetchSuperBlock(seqIdx, superBlockNum, sb, env);
  }
}
#endif /* USE_SBLOCK_CACHE */

#define walkCompIndices(seqIdx, sBlock, blockNum, cwOffset,             \
                        codeForCompIndex, varOffset)                    \
  do {                                                                  \
    unsigned blocksLeft = blockNum;                                     \
    unsigned bitsPerCompositionIdx =                                    \
      seqIdx->compositionTable.compositionIdxBits;                      \
    cwOffset = sBlockGetCompIdxOffset(sBlock, seqIdx, 0);               \
    varOffset = sBlock->varDataMemBase;                                 \
    while (blocksLeft)                                                  \
    {                                                                   \
      PermCompIndex compIndex;                                          \
      compIndex = bsGetPermCompIndex(sBlock->cwData, cwIdxMemOffset,    \
                                     bitsPerCompositionIdx);            \
      codeForCompIndex;                                                 \
      varOffset +=                                                      \
        seqIdx->compositionTable.permutations[compIndex].permIdxBits;   \
      cwOffset += bitsPerCompositionIdx;                                \
      --blocksLeft;                                                     \
    }                                                                   \
  } while (0)

static inline void
unpackBlock(const struct blockCompositionSeq *seqIdx,
            const struct superBlock *sBlock,
            BitOffset cwOffset, BitOffset varOffset, Symbol *block,
            unsigned sublen)
{
  unsigned varIdxBits, bitsPerCompositionIdx;
  PermCompIndex compIndex, permIndex;
  bitsPerCompositionIdx = seqIdx->compositionTable.compositionIdxBits;
  compIndex = bsGetPermCompIndex(sBlock->cwData, cwOffset,
                                 bitsPerCompositionIdx);
  varIdxBits = seqIdx->compositionTable.permutations[compIndex].permIdxBits;
  permIndex = bsGetPermCompIndex(sBlock->varData, varOffset, varIdxBits);
  indexPair2block(&seqIdx->compositionTable, seqIdx->blockSize,
                  compIndex, permIndex, block, sublen);
}

/*
 * regular, user-accessible query functions
 */
static Symbol *
blockCompSeqGetBlock(struct blockCompositionSeq *seqIdx, Seqpos blockNum,
                     struct blockEncIdxSeqHint *hint, int queryRangeEnc,
                     struct superBlock *sBlockPreFetch,
                     Symbol *blockPA, Env *env)
{
  struct superBlock *sBlock;
  BitOffset varDataMemOffset, cwIdxMemOffset;
  Seqpos relBlockNum;
  unsigned blockSize;
  Symbol *block;
  assert(seqIdx);
  blockSize = seqIdx->blockSize;
  if (blockNum * blockSize >= EISLength(&seqIdx->baseClass))
    return NULL;
  if (sBlockPreFetch)
    sBlock = sBlockPreFetch;
  else
  {
    Seqpos bucketNum = bucketNumFromBlockNum(seqIdx, blockNum);
#ifdef USE_SBLOCK_CACHE
    sBlock = cacheFetchSuperBlock(seqIdx, bucketNum,
                                  &hint->sBlockCache, env);
#else
    sBlock = fetchSuperBlock(seqIdx, bucketNum, NULL, env);
#endif
  }
  relBlockNum = blockNum % seqIdx->bucketBlocks;
  if (blockPA)
    block = blockPA;
  else
    block = env_ma_calloc(env, sizeof (Symbol), blockSize);
  walkCompIndices(seqIdx, sBlock, blockNum % seqIdx->bucketBlocks,
                  cwIdxMemOffset, , varDataMemOffset);
  unpackBlock(seqIdx, sBlock, cwIdxMemOffset, varDataMemOffset, block,
              blockSize);
  if (queryRangeEnc)
    SRLapplyRangesToSubString(seqIdx->rangeEncs, seqIdx->rangeMapAlphabet,
                              block, blockNum * blockSize, blockSize,
                              blockNum * blockSize, &hint->rangeHint);
#ifndef USE_SBLOCK_CACHE
  if (!sBlockPreFetch)
    deleteSuperBlock(sBlock, env);
#endif
  return block;
}

/* Note: pos is meant exclusively, i.e. returns 0
   for any query where pos==0 because that corresponds to the empty prefix */
static Seqpos
blockCompSeqRank(struct encIdxSeq *eSeqIdx, Symbol eSym, Seqpos pos,
                 union EISHint *hint, Env *env)
{
  struct blockCompositionSeq *seqIdx;
  Seqpos blockNum, rankCount;
  assert(eSeqIdx && env && eSeqIdx->classInfo == &blockCompositionSeqClass);
  seqIdx = encIdxSeq2blockCompositionSeq(eSeqIdx);
  assert(MRAEncSymbolIsInSelectedRanges(seqIdx->baseClass.alphabet,
                                        eSym, BLOCK_COMPOSITION_INCLUDE,
                                        seqIdx->modes) >= 0);
  if (MRAEncSymbolIsInSelectedRanges(seqIdx->baseClass.alphabet, eSym,
                                     BLOCK_COMPOSITION_INCLUDE, seqIdx->modes))
  {
    BitOffset varDataMemOffset, cwIdxMemOffset;
    struct superBlock *sBlock;
    Symbol bSym = MRAEncMapSymbol(seqIdx->blockMapAlphabet, eSym);
    Seqpos bucketNum;
    unsigned blockSize = seqIdx->blockSize, bitsPerCompositionIdx
      = seqIdx->compositionTable.compositionIdxBits;
    bucketNum = bucketNumFromPos(seqIdx, pos);
#ifdef USE_SBLOCK_CACHE
    sBlock = cacheFetchSuperBlock(seqIdx, bucketNum,
                                  &hint->bcHint.sBlockCache, env);
#else
    sBlock = fetchSuperBlock(seqIdx, bucketNum, NULL, env);
#endif
    rankCount = sBlockGetPartialSymSum(sBlock, bSym, seqIdx);
    blockNum = blockNumFromPos(seqIdx, pos);
    walkCompIndices(
      seqIdx, sBlock, blockNum % seqIdx->bucketBlocks, cwIdxMemOffset,
      rankCount += symCountFromComposition(
        &seqIdx->compositionTable, seqIdx->blockEncNumSyms, compIndex, bSym);,
      varDataMemOffset);
    {
      Seqpos inBlockPos;
      if ((inBlockPos = pos % blockSize)
          && symCountFromComposition(
            &seqIdx->compositionTable, seqIdx->blockEncNumSyms,
            bsGetPermCompIndex(sBlock->cwData, cwIdxMemOffset,
                               bitsPerCompositionIdx), bSym))
      {
        Symbol block[blockSize];
        unsigned i;
        unpackBlock(seqIdx, sBlock, cwIdxMemOffset, varDataMemOffset, block,
                    inBlockPos);
        for (i = 0; i < inBlockPos; ++i)
        {
          if (block[i] == eSym)
            ++rankCount;
        }
      }
      if (bSym == seqIdx->blockEncFallback)
      {
        Seqpos base = bucketBasePos(seqIdx, bucketNum);
        rankCount -= SRLAllSymbolsCountInSeqRegion(
          seqIdx->rangeEncs, base, pos, &hint->bcHint.rangeHint);
      }
    }
#ifndef USE_SBLOCK_CACHE
    deleteSuperBlock(sBlock, env);
#endif
  }
  else
  {
    rankCount = SRLSymbolCountInSeqRegion(seqIdx->rangeEncs, 0, pos, eSym,
                                          &hint->bcHint.rangeHint);
  }
  return rankCount;
}

static void
blockCompSeqExpose(struct encIdxSeq *eSeqIdx, Seqpos pos, int flags,
                   struct extBitsRetrieval *retval, union EISHint *hint,
                   Env *env)
{
  struct blockCompositionSeq *seqIdx;
  assert(eSeqIdx && env && eSeqIdx->classInfo == &blockCompositionSeqClass);
  assert(retval);
  seqIdx = encIdxSeq2blockCompositionSeq(eSeqIdx);
  {
    struct superBlock *sBlock;
    Seqpos bucketNum, end;
    bucketNum = bucketNumFromPos(seqIdx, pos);
    retval->start = bucketBasePos(seqIdx, bucketNum);
    end = bucketBasePos(seqIdx, bucketNum + 1);
    if (end >= seqIdx->baseClass.seqLen)
      end = seqIdx->baseClass.seqLen;
    retval->len = end - retval->start;
#ifdef USE_SBLOCK_CACHE
    sBlock = cacheFetchSuperBlock(seqIdx, bucketNum,
                                  &hint->bcHint.sBlockCache, env);
#else
    sBlock = fetchSuperBlock(seqIdx, bucketNum, NULL, env);
#endif
    retval->cwOffset = sBlockCWExtBitsOffset(sBlock, seqIdx);
    if (flags & EBRF_PERSISTENT_CWBITS)
    {
      if (!(retval->flags & EBRF_PERSISTENT_CWBITS))
      {
          retval->varPart
            = env_ma_malloc(env, superBlockCWMaxReadSize(seqIdx));
      }
      memcpy(retval->cwPart, sBlock->cwData,
             superBlockCWMaxReadSize(seqIdx));
    }
    else
    {
      if (retval->flags & EBRF_PERSISTENT_CWBITS)
        env_ma_free(retval->cwPart, env);
      retval->cwPart = sBlock->cwData;
    }
    if (EBRF_RETRIEVE_VARBITS)
    {
      retval->varOffset = sBlockGetcbOffset(sBlock, seqIdx)
        + sBlock->varDataMemBase;
      if (flags & EBRF_PERSISTENT_VARBITS)
      {
        if (!(retval->flags & EBRF_PERSISTENT_VARBITS))
        {
          retval->varPart
            = env_ma_malloc(env, superBlockVarMaxReadSize(seqIdx));
        }
        memcpy(retval->varPart, sBlock->varData,
               superBlockVarMaxReadSize(seqIdx));
      }
      else
      {
        if (retval->flags & EBRF_PERSISTENT_VARBITS)
          env_ma_free(retval->varPart, env);
        retval->varPart = sBlock->varData;
      }
    }
    else
    {
      if (retval->flags & EBRF_PERSISTENT_VARBITS)
        env_ma_free(retval->varPart, env);
      retval->varPart = NULL;
      retval->varOffset = 0;
    }
#ifndef USE_SBLOCK_CACHE
    deleteSuperBlock(sBlock, env);
#endif
    retval->flags = flags;
  }
}

static Symbol
blockCompSeqGet(struct encIdxSeq *seq, Seqpos pos, union EISHint *hint,
                Env *env)
{
  Symbol sym;
  struct blockCompositionSeq *seqIdx;
  unsigned blockSize;
  assert(seq && seq->classInfo == &blockCompositionSeqClass);
  if (pos >= seq->seqLen)
    return ~(Symbol)0;
  seqIdx = encIdxSeq2blockCompositionSeq(seq);
  blockSize = seqIdx->blockSize;
  {
    Symbol block[blockSize];
    blockCompSeqGetBlock(seqIdx, pos/blockSize, &(hint->bcHint),
                         1, NULL, block, env);
    sym = block[pos%blockSize];
  }
  return sym;
}

static inline partialSymSums
newPartialSymSums(unsigned alphabetSize, Env *env)
{
  return env_ma_calloc(env, alphabetSize, sizeof (Seqpos));
}

static inline void
deletePartialSymSums(partialSymSums sums, Env *env)
{
  env_ma_free(sums, env);
}

static inline void
addBlock2PartialSymSums(partialSymSums sums, Symbol *block, unsigned blockSize)
{
  unsigned i;
  for (i = 0; i < blockSize; ++i)
    ++(sums[block[i]]);
}

static inline void
copyPartialSymSums(unsigned alphabetSize, partialSymSums dest,
                   const partialSymSums src)
{
  memcpy(dest, src, sizeof (dest[0]) * alphabetSize);
}

static inline BitOffset
superBlockCWBits(const struct blockCompositionSeq *seqIdx)
{
  return seqIdx->symSumBits + seqIdx->bitsPerVarDiskOffset
    + seqIdx->callBackDataOffsetBits
    + seqIdx->compositionTable.compositionIdxBits * seqIdx->bucketBlocks
    + seqIdx->cwExtBitsPerBucket;
}

static inline size_t
superBlockCWMaxReadSize(const struct blockCompositionSeq *seqIdx)
{
  return bitElemsAllocSize(superBlockCWBits(seqIdx) + bitElemBits - 1)
    * sizeof (BitElem);
}

static inline size_t
superBlockVarMaxBits(const struct blockCompositionSeq *seqIdx)
{
  return seqIdx->compositionTable.maxPermIdxBits * seqIdx->bucketBlocks
    + seqIdx->maxVarExtBitsPerBucket;
}

static inline size_t
superBlockVarMaxReadSize(const struct blockCompositionSeq *seqIdx)
{
  return bitElemsAllocSize(superBlockVarMaxBits(seqIdx)
                           + bitElemBits - 1);
}

static inline Seqpos
numBuckets(Seqpos seqLen, Seqpos bucketLen)
{
  return seqLen / bucketLen + ((seqLen % bucketLen)?1:0);
}

static inline off_t
cwSize(const struct blockCompositionSeq *seqIdx)
{
  off_t cwLen;
  cwLen = bitElemsAllocSize(
    superBlockCWBits(seqIdx)
    * numBuckets(seqIdx->baseClass.seqLen,
                 seqIdx->bucketBlocks * seqIdx->blockSize)) * sizeof (BitElem);
  return cwLen;
}

static inline BitOffset
vwBits(Seqpos seqLen, unsigned blockSize, unsigned bucketBlocks,
       unsigned maxPermIdxBits, BitOffset maxVarExtBitsPerBucket)
{
  return numBuckets(seqLen, bucketBlocks * blockSize)
    * (maxPermIdxBits * bucketBlocks + maxVarExtBitsPerBucket);
}

/**
 * @return 0 on error, 1 otherwise
 */
static int
openOnDiskData(const Str *projectName, struct onDiskBlockCompIdx *idx,
               char *mode, Env *env)
{
  Str *bdxName = str_clone(projectName, env);
  str_append_cstr(bdxName, ".bdx", env);
  idx->idxFP = env_fa_fopen(env, str_get(bdxName), mode);
  str_delete(bdxName, env);
  if (!idx->idxFP)
    return 0;
  else
    return 1;
}

enum {
  EXT_HEADER_PREFIX_SIZE = 8,
  HEADER_PAGESIZE_ROUNDUP = 8192,
  HEADER_ID_BLOCK_LEN = 8,
};

static void
initOnDiskBlockCompIdx(struct onDiskBlockCompIdx *idx,
                       size_t headerLen, off_t cwLen)
{
  idx->cwDataPos = roundUp(headerLen, HEADER_PAGESIZE_ROUNDUP);
  idx->varDataPos = cwLen + idx->cwDataPos;
  idx->rangeEncPos = 0;
}

static void
destructOnDiskBlockCompIdx(struct onDiskBlockCompIdx *idx, Env *env)
{
  if (idx->idxMMap)
    munmap(idx->idxMMap, idx->rangeEncPos - idx->cwDataPos);
  if (idx->idxFP)
    env_fa_xfclose(idx->idxFP, env);
}

static inline void
initAppendState(struct appendState *aState,
                const struct blockCompositionSeq *seqIdx, Env *env)
{
  BitOffset compCacheLen = superBlockCWBits(seqIdx) + bitElemBits - 1,
    permCacheLen = superBlockVarMaxBits(seqIdx) + bitElemBits - 1;
  aState->compCacheLen = compCacheLen;
  aState->permCacheLen = permCacheLen;
  aState->compCache = env_ma_malloc(env, sizeof (BitElem) *
                                    bitElemsAllocSize(compCacheLen));
  aState->permCache = env_ma_malloc(env, sizeof (BitElem) *
                                    bitElemsAllocSize(permCacheLen));
  aState->cwMemPos = cwPreCompIdxBits(seqIdx);
  aState->cwDiskOffset = aState->varMemPos = aState->cwMemOldBits =
    aState->varDiskOffset = aState->varMemOldBits = 0;
}

static void
destructAppendState(struct appendState *aState, Env *env)
{
  env_ma_free(aState->permCache, env);
  env_ma_free(aState->compCache, env);
}

/**
 * @return 0 on successfull update, <0 on error
 */
static inline void
append2IdxOutput(struct blockCompositionSeq *newSeqIdx,
                 struct appendState *state,
                 PermCompIndex permCompIdx[2],
                 unsigned bitsOfCompositionIdx, unsigned bitsOfPermutationIdx)
{
  assert(state->cwMemPos + bitsOfCompositionIdx <= state->compCacheLen);
  bsStorePermCompIndex(state->compCache, state->cwMemPos, bitsOfCompositionIdx,
                       permCompIdx[0]);
  state->cwMemPos += bitsOfCompositionIdx;
  assert(state->varMemPos + bitsOfPermutationIdx <= state->permCacheLen);
  bsStorePermCompIndex(state->permCache, state->varMemPos, bitsOfPermutationIdx,
                       permCompIdx[1]);
  state->varMemPos += bitsOfPermutationIdx;
}

static BitOffset
appendCallBackOutput(struct appendState *state,
                     const struct blockCompositionSeq *seqIdx,
                     bitInsertFunc biFunc, Seqpos start, Seqpos len,
                     unsigned callBackDataOffsetBits, void *cbState, Env *env)
{
  BitOffset bitsWritten;
  assert(state && env);
  if (callBackDataOffsetBits)
  {
    BitOffset offset = cwPreCBOffsetBits(seqIdx);
    bsStoreUInt64(state->compCache, state->cwMemOldBits + offset,
                  callBackDataOffsetBits,
                  state->varMemPos - state->varMemOldBits);
  }
  bitsWritten = biFunc(state->compCache, state->cwMemOldBits
                       + cwPreCWExtBits(seqIdx),
                       state->permCache, state->varMemPos, start, len, cbState,
                       env);
  if (bitsWritten == (BitOffset)-1)
    return bitsWritten;
  state->cwMemPos = cwPreCWExtBits(seqIdx) + state->cwMemOldBits
    + seqIdx->cwExtBitsPerBucket;
  state->varMemPos += bitsWritten;
  return bitsWritten;
}

/**
 * @return >0 on success, 0 on I/O-error
 */
static int
updateIdxOutput(struct blockCompositionSeq *seqIdx,
                struct appendState *aState,
                partialSymSums buck)
{
  size_t recordsExpected, cwBitElems;
  unsigned blockAlphabetSize;
/*   BitOffset sBucketBits; */
  assert(seqIdx && aState && buck);
  /* seek2/write constant width indices */
  assert(seqIdx->externalData.cwDataPos + aState->cwDiskOffset
         < seqIdx->externalData.varDataPos
         + aState->varDiskOffset/bitElemBits * sizeof (BitElem));
  blockAlphabetSize = seqIdx->blockEncNumSyms;
/*   sBucketBits = seqIdx->bitsPerSeqpos * blockAlphabetSize; */
  /* put bucket data for count up to the beginning of the current
   * block into the cw BitString */
  assert(sizeof (Seqpos) * CHAR_BIT >= seqIdx->bitsPerSeqpos);
  {
    Symbol i;
    for (i = 0; i < blockAlphabetSize; ++i)
    {
      bsStoreSeqpos(aState->compCache,
                    aState->cwMemOldBits + seqIdx->partialSymSumBitsSums[i],
                    seqIdx->partialSymSumBits[i], buck[i]);
    }
  }
  /* append variable width offset position */
  bsStoreUInt64(aState->compCache, aState->cwMemOldBits
                + cwPreVarIdxBits(seqIdx), seqIdx->bitsPerVarDiskOffset,
                aState->varDiskOffset);
  if (fseeko(seqIdx->externalData.idxFP,
            seqIdx->externalData.cwDataPos + aState->cwDiskOffset,
            SEEK_SET))
    return 0;
  cwBitElems = recordsExpected = aState->cwMemPos / bitElemBits;
  if (recordsExpected != fwrite(aState->compCache, sizeof (BitElem),
                               recordsExpected, seqIdx->externalData.idxFP))
    return 0;
  if ((aState->cwMemOldBits = aState->cwMemPos % bitElemBits))
    aState->compCache[0] = aState->compCache[recordsExpected];
  /* seek2/write variable width indices */
  if (fseeko(seqIdx->externalData.idxFP, seqIdx->externalData.varDataPos
            + aState->varDiskOffset/bitElemBits * sizeof (BitElem), SEEK_SET))
    return 0;
  recordsExpected = aState->varMemPos/bitElemBits;
  if (recordsExpected != fwrite(aState->permCache, sizeof (BitElem),
                               recordsExpected, seqIdx->externalData.idxFP))
    return 0;
  /* move last elem of permCache with unwritten bits to front */
  if (aState->varMemPos % bitElemBits)
    aState->permCache[0] = aState->permCache[recordsExpected];
  /* FIXME: this assumes that the string of variable width date does
   * indeed occupy at least one full bitElem */
  aState->cwDiskOffset += cwBitElems;
  aState->cwMemPos = cwPreCompIdxBits(seqIdx) + aState->cwMemOldBits;
  aState->varDiskOffset += aState->varMemPos - aState->varMemOldBits;
  aState->varMemOldBits = (aState->varMemPos %= bitElemBits);
  return 1;
}

/* Caution: EH??-headers are reserved for extension headers */
enum bdxHeader {
  BKSZ_HEADER_FIELD = 0x424b535a,
  BBLK_HEADER_FIELD = 0x42424c4b,
  VOFF_HEADER_FIELD = 0x564f4646,
  ROFF_HEADER_FIELD = 0x524f4646,
  NMRN_HEADER_FIELD = 0x4e4d524e,
  CBMB_HEADER_FIELD = 0x43424d42,
  MEXB_HEADER_FIELD = 0x4d455842,
  CEXB_HEADER_FIELD = 0x43455842,
  SPBT_HEADER_FIELD = 0x53504254,
  SSBT_HEADER_FIELD = 0x53534254,
  BEFB_HEADER_FIELD = 0x42454642,
  REFB_HEADER_FIELD = 0x52454642,
  EH_HEADER_PREFIX = 0x45480000,
};

static const char bdxHeader[] = "BDX";

static inline off_t
extHeadersSizeAggregate(size_t numExtHeaders, uint32_t *extHeaderSizes)
{
  off_t len = 0, i;
  for (i = 0; i < numExtHeaders; ++i)
    len += extHeaderSizes[i] + EXT_HEADER_PREFIX_SIZE;
  return len;
}

static inline size_t
blockEncIdxSeqHeaderLength(struct blockCompositionSeq *seqIdx,
                           size_t numExtHeaders, uint32_t *extHeaderSizes)
{
  size_t headerSize =
    4                           /* BDX identifier */
    + 4                         /* length field */
    + 8                         /* block size */
    + 8                         /* blocks per bucket */
    + 12                        /* offset of variable length data */
    + 12                        /* offset of range encodings */
    + 4 + 4                     /* bits used per seqpos */
    + 4 + 4 + 4 * seqIdx->blockEncNumSyms /* bit counts for partial sums */
    + 4 + 4                     /* block encoding fallback symbol */
    + 4 + 4                     /* range encoding fallback symbol */
    + 4 + 4                     /* num modes */
    + 4 * seqIdx->numModes      /* one uint32_t for every mode */
    ;
  if (seqIdx->callBackDataOffsetBits)
    headerSize += 4 + 4         /* extra offset bits per constant block */
      + 4 + 8                   /* extension bits stored in constant block */
      + 4 + 8;                  /* variable area bits added per bucket max */
  headerSize += extHeadersSizeAggregate(numExtHeaders, extHeaderSizes);
  return headerSize;
}

struct extHeaderPos
{
  off_t pos;
  uint32_t headerID;
};

static inline void
appendExtHeaderPos(struct extHeaderPos **headerList, size_t numHeaders,
                   off_t pos, uint32_t headerID,
                   Env *env)
{
  *headerList = env_ma_realloc(env, *headerList,
                               sizeof (**headerList) * (numHeaders + 1));
  (*headerList)[numHeaders].pos = pos;
  (*headerList)[numHeaders].headerID = headerID;
};

static inline int
writeExtIdxHeader(FILE *fp, uint16_t headerID, size_t len,
                  headerWriteFunc cbFunc, void *cbData)
{
  uint32_t expHeader[2] =
    { EH_HEADER_PREFIX | headerID, len };
  if (fwrite(expHeader, sizeof (uint32_t), 2, fp)!= 2)
    return 0;
  return cbFunc(fp, cbData);
}

#define writeIdxHeaderErrRet(retval)            \
  do {                                          \
    env_ma_free(buf, env);                      \
    return 0;                                   \
  } while (0)
/**
 * FIXME: this doesn't work on platforms with sizeof (uint32_t) != 4
 * or sizeof (int) < 4
 * @return 0 on error, header length in bytes on success
 */
static size_t
writeIdxHeader(struct blockCompositionSeq *seqIdx,
               size_t numExtHeaders, uint16_t *headerIDs,
               uint32_t *extHeaderSizes, headerWriteFunc *extHeaderCallbacks,
               void **headerCBData, Env *env)
{
  FILE *fp;
  /* construct memory buffer with header data */
  size_t i, bufLen;
  off_t offset, len;
  char *buf;
  assert(seqIdx && env);
  fp = seqIdx->externalData.idxFP;
  bufLen = blockEncIdxSeqHeaderLength(seqIdx, 0, NULL);
  buf = env_ma_malloc(env, bufLen);
  /* account for ext headers */
  len = roundUp(bufLen, HEADER_PAGESIZE_ROUNDUP);
  /* 1. 4 identifier bytes at offset 0 */
  strcpy(buf, bdxHeader);
  /* 2. header length */
  *(uint32_t *)(buf + 4) = len;
  /* 3. block size */
  offset = HEADER_ID_BLOCK_LEN;
  *(uint32_t *)(buf + offset) = BKSZ_HEADER_FIELD;
  *(uint32_t *)(buf + offset + 4) = seqIdx->blockSize;
  offset += 8;
  *(uint32_t *)(buf + offset) = BBLK_HEADER_FIELD;
  *(uint32_t *)(buf + offset + 4) = seqIdx->bucketBlocks;
  offset += 8;
  *(uint32_t *)(buf + offset) = VOFF_HEADER_FIELD;
  *(uint64_t *)(buf + offset + 4) = seqIdx->externalData.varDataPos;
  offset += 12;
  *(uint32_t *)(buf + offset) = ROFF_HEADER_FIELD;
  *(uint64_t *)(buf + offset + 4) = seqIdx->externalData.rangeEncPos;
  offset += 12;
  *(uint32_t *)(buf + offset) = SPBT_HEADER_FIELD;
  *(uint32_t *)(buf + offset + 4) = seqIdx->bitsPerSeqpos;
  offset += 8;
  *(uint32_t *)(buf + offset) = SSBT_HEADER_FIELD;
  *(uint32_t *)(buf + offset + 4) = seqIdx->blockEncNumSyms;
  for (i = 0; i < seqIdx->blockEncNumSyms; ++i)
    *(uint32_t *)(buf + offset + 8 + 4*i) = seqIdx->partialSymSumBits[i];
  offset += 8 + 4 * seqIdx->blockEncNumSyms;
  *(uint32_t *)(buf + offset) = BEFB_HEADER_FIELD;
  *(uint32_t *)(buf + offset + 4) = seqIdx->blockEncFallback;
  offset += 8;
  *(uint32_t *)(buf + offset) = REFB_HEADER_FIELD;
  *(uint32_t *)(buf + offset + 4) = seqIdx->rangeEncFallback;
  offset += 8;
  *(uint32_t *)(buf + offset) = NMRN_HEADER_FIELD;
  *(uint32_t *)(buf + offset + 4) = seqIdx->numModes;
  offset += 8;
  for (i = 0; i < seqIdx->numModes; ++i)
  {
    *(uint32_t *)(buf + offset) = seqIdx->modes[i];
    offset += 4;
  }
  if (seqIdx->callBackDataOffsetBits)
  {
    *(uint32_t *)(buf + offset) = CBMB_HEADER_FIELD;
    *(uint32_t *)(buf + offset + 4) = seqIdx->callBackDataOffsetBits;
    offset += 8;
    *(uint32_t *)(buf + offset) = CEXB_HEADER_FIELD;
    *(uint64_t *)(buf + offset + 4) = seqIdx->cwExtBitsPerBucket;
    offset += 12;
    *(uint32_t *)(buf + offset) = MEXB_HEADER_FIELD;
    *(uint64_t *)(buf + offset + 4) = seqIdx->maxVarExtBitsPerBucket;
    offset += 12;
  }
  assert(offset == bufLen);
  if (fseeko(fp, 0, SEEK_SET))
    writeIdxHeaderErrRet(0);
  if (fwrite(buf, bufLen, 1, fp) != 1)
    writeIdxHeaderErrRet(0);
  if (numExtHeaders)
  {
    off_t offsetTargetValue = bufLen;
    if (!writeExtIdxHeader(fp, headerIDs[0], extHeaderSizes[0],
                          extHeaderCallbacks[0], headerCBData[0]))
      writeIdxHeaderErrRet(0);
    else
      appendExtHeaderPos(&seqIdx->extHeaderPos, seqIdx->numExtHeaders++,
                         bufLen + EXT_HEADER_PREFIX_SIZE,
                         EH_HEADER_PREFIX | headerIDs[0], env);
    for (i = 1; i < numExtHeaders; ++i)
    {
      offsetTargetValue += EXT_HEADER_PREFIX_SIZE + extHeaderSizes[i - 1];
      if (fseeko(fp, offsetTargetValue, SEEK_SET))
        writeIdxHeaderErrRet(0);
      if (!writeExtIdxHeader(fp, headerIDs[i], extHeaderSizes[i],
                            extHeaderCallbacks[i], headerCBData[i]))
        writeIdxHeaderErrRet(0);
      else
        appendExtHeaderPos(&seqIdx->extHeaderPos, seqIdx->numExtHeaders++,
                           offsetTargetValue + EXT_HEADER_PREFIX_SIZE,
                           EH_HEADER_PREFIX | headerIDs[0], env);
    }
    offset = offsetTargetValue + EXT_HEADER_PREFIX_SIZE
      + extHeaderSizes[numExtHeaders - 1];
    assert(offset >= ftello(fp));
  }
  assert(seqIdx->externalData.cwDataPos == len);
  env_ma_free(buf, env);
  return len;
}

struct encIdxSeq *
loadBlockEncIdxSeq(const Str *projectName, int features, Env *env)
{
  struct encIdxSeq *newSeqIdx;
  Suffixarray suffixArray;
  Seqpos len;
  Verboseinfo *verbosity;
  /* FIXME: handle verbosity in a saner fashion */
  verbosity = newverboseinfo(false, env);
  if (streamsuffixarray(&suffixArray, &len,
                        0, projectName, verbosity, env))
  {
    freeverboseinfo(&verbosity, env);
    return NULL;
  }
  ++len;
  newSeqIdx = loadBlockEncIdxSeqForSA(&suffixArray, len, projectName,
                                      features, env);
  freesuffixarray(&suffixArray, env);
  freeverboseinfo(&verbosity, env);
  return newSeqIdx;
}

#define loadBlockEncIdxSeqErrRet()                                      \
  do {                                                                  \
    if (newSeqIdx->externalData.idxFP)                                  \
      destructOnDiskBlockCompIdx(&newSeqIdx->externalData, env);        \
    if (newSeqIdx->compositionTable.bitsPerCount)                       \
      destructCompositionList(&newSeqIdx->compositionTable, env);       \
    if (newSeqIdx->rangeEncs)                                           \
      deleteSeqRangeList(newSeqIdx->rangeEncs, env);                    \
    if (newSeqIdx->extHeaderPos)                                        \
      env_ma_free(newSeqIdx->extHeaderPos, env);                        \
    if (buf) env_ma_free(buf, env);                                     \
    if (alphabet) MRAEncDelete(alphabet, env);                          \
    if (modesCopy)                                                      \
      env_ma_free(modesCopy, env);                                      \
    if (blockMapAlphabet) env_ma_free(blockMapAlphabet, env);           \
    if (rangeMapAlphabet) env_ma_free(rangeMapAlphabet, env);           \
    if (newSeqIdx) env_ma_free(newSeqIdx, env);                         \
    return NULL;                                                        \
  } while (0)

struct encIdxSeq *
loadBlockEncIdxSeqForSA(Suffixarray *sa, Seqpos totalLen,
                        const Str *projectName, int features, Env *env)
{
  struct blockCompositionSeq *newSeqIdx = NULL;
  Symbol blockMapAlphabetSize, totalAlphabetSize;
  MRAEnc *alphabet = NULL, *blockMapAlphabet = NULL, *rangeMapAlphabet = NULL;
  size_t headerLen;
  int *modesCopy = NULL;
  char *buf = NULL;
  assert(projectName && env);
  newSeqIdx = env_ma_calloc(env, sizeof (struct blockCompositionSeq), 1);
  newSeqIdx->baseClass.seqLen = totalLen;
  newSeqIdx->baseClass.alphabet = alphabet =
    MRAEncGTAlphaNew(sa->alpha, env);
  MRAEncAddSymbolToRange(alphabet, SEPARATOR, 1);
  newSeqIdx->baseClass.classInfo = &blockCompositionSeqClass;

  if (!openOnDiskData(projectName, &newSeqIdx->externalData, "rb", env))
    loadBlockEncIdxSeqErrRet();
  {
    size_t offset = HEADER_ID_BLOCK_LEN;
    buf = env_ma_malloc(env, HEADER_ID_BLOCK_LEN);
    if (fread(buf, HEADER_ID_BLOCK_LEN, 1, newSeqIdx->externalData.idxFP) != 1)
      loadBlockEncIdxSeqErrRet();
    if (strcmp(buf, bdxHeader)!= 0)
      loadBlockEncIdxSeqErrRet();
    newSeqIdx->externalData.cwDataPos = headerLen = *(uint32_t *)(buf + 4);
    buf = env_ma_realloc(env, buf, headerLen);
    if (fread(buf + HEADER_ID_BLOCK_LEN, headerLen - HEADER_ID_BLOCK_LEN,
             1, newSeqIdx->externalData.idxFP) != 1)
      loadBlockEncIdxSeqErrRet();
    while (offset < headerLen)
    {
      uint32_t currentHeader;
      switch (currentHeader = *(uint32_t *)(buf + offset))
      {
      case BKSZ_HEADER_FIELD:
        newSeqIdx->blockSize = *(uint32_t *)(buf + offset + 4);
        offset += 8;
        break;
      case BBLK_HEADER_FIELD:
        newSeqIdx->bucketBlocks = *(uint32_t *)(buf + offset + 4);
        offset += 8;
        break;
      case VOFF_HEADER_FIELD:
        newSeqIdx->externalData.varDataPos
          = *(uint64_t *)(buf + offset + 4);
        offset += 12;
        break;
      case ROFF_HEADER_FIELD:
        newSeqIdx->externalData.rangeEncPos = *(uint64_t *)(buf + offset + 4);
        offset += 12;
        break;
      case NMRN_HEADER_FIELD:
        {
          size_t numModes = newSeqIdx->numModes
            = *(uint32_t *)(buf + (offset += 4));
          offset += 4;
          modesCopy = newSeqIdx->modes
            = env_ma_malloc(env, sizeof (int) * numModes);
          size_t j;
          for (j = 0; j < numModes; ++j)
          {
            modesCopy[j] = *(uint32_t *)(buf + offset);
            offset += 4;
          }
        }
        break;
      case CBMB_HEADER_FIELD:
        newSeqIdx->callBackDataOffsetBits = *(uint32_t *)(buf + offset + 4);
        offset += 8;
        break;
      case SPBT_HEADER_FIELD:
        newSeqIdx->bitsPerSeqpos = *(uint32_t *)(buf + offset + 4);
        offset += 8;
        break;
      case SSBT_HEADER_FIELD:
        {
          size_t i;
          Symbol blockMapAlphabetSize = *(uint32_t *)(buf + offset + 4);
          newSeqIdx->partialSymSumBits
            = env_ma_malloc(env, sizeof (newSeqIdx->partialSymSumBits[0])
                            * blockMapAlphabetSize * 2);
          newSeqIdx->partialSymSumBitsSums
            = newSeqIdx->partialSymSumBits + blockMapAlphabetSize;
          if (blockMapAlphabetSize)
          {
            newSeqIdx->partialSymSumBits[0]= *(uint32_t *)(buf + offset + 8);
            newSeqIdx->partialSymSumBitsSums[0] = 0;
#ifdef DEBUG
            fprintf(stderr, "partialSymSumBits[0]=%u\n",
                    newSeqIdx->partialSymSumBits[0]);
#endif
            for (i = 1; i < blockMapAlphabetSize; ++i)
            {
              newSeqIdx->partialSymSumBits[i]
                = *(uint32_t *)(buf + offset + 8 + 4*i);
#ifdef DEBUG
              fprintf(stderr, "partialSymSumBits[%u]=%u\n",
                      (unsigned)i, newSeqIdx->partialSymSumBits[i]);
#endif
              newSeqIdx->partialSymSumBitsSums[i] =
                newSeqIdx->partialSymSumBitsSums[i - 1]
                + newSeqIdx->partialSymSumBits[i - 1];
            }
            newSeqIdx->symSumBits
              = newSeqIdx->partialSymSumBitsSums[blockMapAlphabetSize - 1]
              + newSeqIdx->partialSymSumBits[blockMapAlphabetSize - 1];
          }
          offset += 8 + 4 * blockMapAlphabetSize;
        }
      case BEFB_HEADER_FIELD:
        newSeqIdx->blockEncFallback = *(uint32_t *)(buf + offset + 4);
        offset += 8;
        break;
      case REFB_HEADER_FIELD:
        newSeqIdx->rangeEncFallback = *(uint32_t *)(buf + offset + 4);
        offset += 8;
        break;
      case MEXB_HEADER_FIELD:
        newSeqIdx->maxVarExtBitsPerBucket = *(uint64_t *)(buf + offset + 4);
        offset += 12;
        break;
      case CEXB_HEADER_FIELD:
        newSeqIdx->cwExtBitsPerBucket = *(uint64_t *)(buf + offset + 4);
        offset += 12;
        break;
      case 0:
        /* empty header skip to next portion */
        offset = headerLen;
        break;
      default:
        if ((currentHeader & 0xffff0000L) == EH_HEADER_PREFIX)
        {
          uint32_t extHeaderLen = *(uint32_t *)(buf + offset + 4);
          offset += EXT_HEADER_PREFIX_SIZE;
          appendExtHeaderPos(&newSeqIdx->extHeaderPos,
                             newSeqIdx->numExtHeaders++,
                             offset, currentHeader, env);
          offset += extHeaderLen;
        }
        else
        {
          fprintf(stderr, "Unknown header field: %4s\n", buf + offset);
          loadBlockEncIdxSeqErrRet();
        }
      }
    }
    assert(newSeqIdx->modes
           && newSeqIdx->bucketBlocks
           && newSeqIdx->numModes
           && offset == headerLen);
    /* compute values not read from header */
    if (!newSeqIdx->bitsPerSeqpos)
      newSeqIdx->bitsPerSeqpos =
        requiredSeqposBits(newSeqIdx->baseClass.seqLen - 1);
    else if (newSeqIdx->bitsPerSeqpos !=
            requiredSeqposBits(newSeqIdx->baseClass.seqLen - 1))
      loadBlockEncIdxSeqErrRet();
  }
  {
    size_t range, numAlphabetRanges = newSeqIdx->numModes =
      MRAEncGetNumRanges(alphabet);
    totalAlphabetSize = MRAEncGetSize(alphabet);
    blockMapAlphabetSize = 0;
    for (range = 0; range < numAlphabetRanges; ++range)
    {
      switch (modesCopy[range])
      {
      case BLOCK_COMPOSITION_INCLUDE:
        blockMapAlphabetSize += MRAEncGetRangeSize(alphabet, range);
        break;
      case DIRECT_SYM_ENCODE:
        /*< FIXME: insert proper code to process direct encoding regions */
        break;
      case REGIONS_LIST:
        /*< FIXME: insert proper code to process ranges */
        break;
      default:
        /* TODO: improve diagnostics */
        fprintf(stderr, "Invalid encoding request.\n");
        loadBlockEncIdxSeqErrRet();
        break;
      }
    }
    newSeqIdx->blockMapAlphabet = blockMapAlphabet =
      MRAEncSecondaryMapping(alphabet, BLOCK_COMPOSITION_INCLUDE, modesCopy,
                             newSeqIdx->blockEncFallback, env);
    newSeqIdx->rangeMapAlphabet = rangeMapAlphabet =
      MRAEncSecondaryMapping(alphabet, REGIONS_LIST, modesCopy,
                             newSeqIdx->rangeEncFallback, env);
    newSeqIdx->blockEncNumSyms = blockMapAlphabetSize;
    assert(MRAEncGetSize(blockMapAlphabet) == blockMapAlphabetSize);
  }
  if (!newSeqIdx->partialSymSumBits && blockMapAlphabetSize)
  {
    newSeqIdx->partialSymSumBits
      = env_ma_malloc(env, sizeof (newSeqIdx->partialSymSumBits[0])
                      * blockMapAlphabetSize * 2);
    newSeqIdx->partialSymSumBitsSums
      = newSeqIdx->partialSymSumBits + blockMapAlphabetSize;
    symSumBitsDefaultSetup(newSeqIdx);
  }
  if (!initCompositionList(&newSeqIdx->compositionTable, newSeqIdx->blockSize,
                           blockMapAlphabetSize, env))
    loadBlockEncIdxSeqErrRet();
  newSeqIdx->bitsPerVarDiskOffset =
    requiredUInt64Bits(vwBits(newSeqIdx->baseClass.seqLen, newSeqIdx->blockSize,
                              newSeqIdx->bucketBlocks,
                              newSeqIdx->compositionTable.maxPermIdxBits,
                              newSeqIdx->maxVarExtBitsPerBucket));

  if (fseeko(newSeqIdx->externalData.idxFP,
             newSeqIdx->externalData.rangeEncPos, SEEK_SET))
    loadBlockEncIdxSeqErrRet();
  {
    int regionFeatures = SRL_NO_FEATURES;
    if (features & EIS_FEATURE_REGION_SUMS)
      regionFeatures |= SRL_PARTIAL_SYMBOL_SUMS;
    if (!(newSeqIdx->rangeEncs =
          SRLReadFromStream(newSeqIdx->externalData.idxFP, rangeMapAlphabet,
                            regionFeatures, env)))
      loadBlockEncIdxSeqErrRet();
  }
  tryMMapOfIndex(&newSeqIdx->externalData);
  env_ma_free(buf, env);
  return &newSeqIdx->baseClass;
}

static inline int
tryMMapOfIndex(struct onDiskBlockCompIdx *idxData)
{
  void *indexMMap;
  assert(idxData && idxData->idxFP);
  if ((indexMMap = mmap((void *)0, idxData->rangeEncPos - idxData->cwDataPos,
                        PROT_READ, MAP_SHARED, fileno(idxData->idxFP),
                        idxData->cwDataPos)))
  {
    idxData->idxMMap = indexMMap;
  }
  return idxData->idxMMap != NULL;
}

static FILE *
seekToHeader(const struct encIdxSeq *seqIdx, uint16_t headerID,
             uint32_t *lenRet)
{
  const struct blockCompositionSeq *bSeqIdx;
  size_t i;
  uint32_t query = headerID | EH_HEADER_PREFIX;
  assert(seqIdx);
  bSeqIdx =
    constEncIdxSeq2blockCompositionSeq(seqIdx);
  for (i = 0; i < bSeqIdx->numExtHeaders; ++i)
  {
    if (bSeqIdx->extHeaderPos[i].headerID == query)
    {
      uint32_t headerPrefix[2];
      if (lenRet)
      {
        if (fseeko(bSeqIdx->externalData.idxFP,
                  bSeqIdx->extHeaderPos[i].pos - EXT_HEADER_PREFIX_SIZE,
                  SEEK_SET))
           return NULL;
        if (fread(headerPrefix, sizeof (headerPrefix[0]), 2,
                 bSeqIdx->externalData.idxFP) != 2)
          return NULL;
        *lenRet = headerPrefix[1];
      }
      else
        if (fseeko(bSeqIdx->externalData.idxFP, bSeqIdx->extHeaderPos[i].pos,
                  SEEK_SET))
          return NULL;
      return bSeqIdx->externalData.idxFP;
    }
  }
  return NULL;
}

/**
 * @return 0 on error, 1 on success
 */
static int
finalizeIdxOutput(struct blockCompositionSeq *seqIdx,
                  struct appendState *aState,
                  partialSymSums buck, Env *env)
{
  off_t rangeEncPos;
  size_t recordsExpected;
  assert(aState && seqIdx);
  assert(aState->cwMemOldBits < bitElemBits
         && aState->cwMemOldBits + cwPreCompIdxBits(seqIdx)
         == aState->cwMemPos);
  assert(aState->varMemOldBits < bitElemBits
         && aState->varMemOldBits == aState->varMemPos);
  if (aState->cwMemOldBits)
  {
    /* seek2/write variable width indices */
    if (fseeko(seqIdx->externalData.idxFP,
              seqIdx->externalData.cwDataPos
              + aState->cwDiskOffset * sizeof (BitElem), SEEK_SET))
      return 0;
    recordsExpected = 1;
    if (recordsExpected != fwrite(aState->compCache, sizeof (BitElem),
                                 recordsExpected,
                                 seqIdx->externalData.idxFP))
      return 0;
    ++(aState->cwDiskOffset);
  }
  if (aState->varMemOldBits)
  {
    /* seek2/write variable width indices */
    if (fseeko(seqIdx->externalData.idxFP,
              seqIdx->externalData.varDataPos
              + aState->varDiskOffset/bitElemBits * sizeof (BitElem), SEEK_SET))
      return 0;
    recordsExpected = 1;
    if (recordsExpected != fwrite(aState->permCache, sizeof (BitElem),
                                  recordsExpected,
                                  seqIdx->externalData.idxFP))
      return 0;
  }
  rangeEncPos = seqIdx->externalData.varDataPos
    + aState->varDiskOffset / bitElemBits * sizeof (BitElem)
    + ((aState->varDiskOffset%bitElemBits)?sizeof (BitElem):0);
  seqIdx->externalData.rangeEncPos = rangeEncPos;
  /* insert terminator so every search for a next range will find a
   * range just beyond the sequence end */
  SRLAppendNewRange(seqIdx->rangeEncs,
                    seqIdx->baseClass.seqLen + seqIdx->blockSize, 1, 0, env);
  SRLCompact(seqIdx->rangeEncs, env);
  if (fseeko(seqIdx->externalData.idxFP, rangeEncPos, SEEK_SET))
    return 0;
  if (!(SRLSaveToStream(seqIdx->rangeEncs, seqIdx->externalData.idxFP)))
     return 0;
  return 1;
}

static void
addRangeEncodedSyms(struct seqRangeList *rangeList, Symbol *block,
                    unsigned blockSize, Seqpos blockNum, MRAEnc *alphabet,
                    int selection, int *rangeSel, Env *env)
{
  unsigned i;
  for (i = 0; i < blockSize; ++i)
  {
    assert(MRAEncSymbolIsInSelectedRanges(alphabet, block[i],
                                          selection, rangeSel) >= 0);
    if (MRAEncSymbolIsInSelectedRanges(alphabet, block[i], selection, rangeSel))
      SRLAddPosition(rangeList, blockNum * blockSize + i, block[i], env);
  }
}

static union EISHint *
newBlockCompSeqHint(struct encIdxSeq *seq, Env *env)
{
  union EISHint *hintret;
  struct blockCompositionSeq *seqIdx;
  assert(seq && seq->classInfo == &blockCompositionSeqClass);
  seqIdx = encIdxSeq2blockCompositionSeq(seq);
  hintret = env_ma_malloc(env, sizeof (union EISHint));
  SRLInitListSearchHint(seqIdx->rangeEncs, &hintret->bcHint.rangeHint);
  /* FIXME: make cache size user-configurable */
  initSuperBlockSeqCache(&hintret->bcHint.sBlockCache, seqIdx, 32, env);
  return hintret;
}

static void
deleteBlockCompSeqHint(struct encIdxSeq *seq, union EISHint *hint, Env *env)
{
  struct blockCompositionSeq *seqIdx;
  assert(seq && seq->classInfo == &blockCompositionSeqClass);
  seqIdx = encIdxSeq2blockCompositionSeq(seq);
  destructSuperBlockSeqCache(&hint->bcHint.sBlockCache, env);
  env_ma_free(hint, env);
}

static int
printBlock(Symbol *block, unsigned blockSize, FILE *fp)
{
  unsigned i;
  int outCount = 0;
  for (i = 0; i < blockSize; ++i)
    outCount += fprintf(fp, " %d", (int)block[i]);
  return outCount;
}

static int
printBucket(const struct blockCompositionSeq *seqIdx, Seqpos bucketNum,
            FILE *fp, union EISHint *hint, Env *env)
{
  Seqpos lastBucket = blockNumFromPos(seqIdx, EISLength(&seqIdx->baseClass)),
    start, end;
  unsigned i, blockMapAlphabetSize = seqIdx->blockEncNumSyms;
  int outCount = 0;
  assert(seqIdx && fp && hint && env);
  if (bucketBasePos(seqIdx, bucketNum) >= EISLength(&seqIdx->baseClass))
  {
    fprintf(stderr, "warning: querying bucket "FormatSeqpos
            " beyond end of sequence!\n", bucketNum);
    bucketNum = lastBucket;
  }
  start = bucketBasePos(seqIdx, bucketNum);
  end = ((bucketNum < lastBucket)?
         bucketBasePos(seqIdx, bucketNum + 1):
         EISLength(&seqIdx->baseClass));
  outCount +=
    fprintf(fp, "# Inspecting bucket: "FormatSeqpos"\n"
            "# bucket position start="FormatSeqpos", end="FormatSeqpos"\n"
            "# partial symbol sums up to start:\n",
            bucketNum, start, end - 1);
  {
    struct superBlock *sBlock;
    BitOffset varDataMemOffset, cwIdxMemOffset;
    Symbol *block;
    unsigned blockSize = seqIdx->blockSize;
#ifdef USE_SBLOCK_CACHE
    sBlock = cacheFetchSuperBlock(seqIdx, bucketNum,
                                  &hint->bcHint.sBlockCache, env);
#else
    sBlock = fetchSuperBlock(seqIdx, bucketNum, NULL, env);
#endif
    block = env_ma_malloc(env, sizeof (block[0]) * seqIdx->blockSize);
    for (i = 0; i < blockMapAlphabetSize; ++i)
    {
      outCount +=
        fprintf(fp, "# partial sum[%u]="FormatSeqpos"\n", i,
                sBlockGetPartialSymSum(sBlock, i, seqIdx));
    }
    fprintf(fp, "# sBlockVarIdxOffset=%llu\n",
            sBlockGetVarIdxOffset(sBlock, seqIdx));
    if (seqIdx->callBackDataOffsetBits)
      fprintf(fp, "# sBlockGetcbOffset=%llu\n",
              sBlockGetcbOffset(sBlock, seqIdx));
    i = 0;
    walkCompIndices(
      seqIdx, sBlock, seqIdx->bucketBlocks, cwIdxMemOffset,
      outCount +=
      fprintf(
        fp, "# block %u: comp idx: %lu, permIdxBits=%u, perm idx: %lu =>",
        i, (unsigned long)compIndex,
        (unsigned)seqIdx->compositionTable.permutations[compIndex].permIdxBits,
        (unsigned long)bsGetPermCompIndex(
          sBlock->varData, varDataMemOffset,
          seqIdx->compositionTable.permutations[compIndex].permIdxBits));
      unpackBlock(seqIdx, sBlock, cwIdxMemOffset, varDataMemOffset, block,
                  blockSize);
      outCount += printBlock(block, blockSize, fp);
      fputs("\n", fp);
      i++; , varDataMemOffset);
    env_ma_free(block, env);
  }
  return outCount;
}

static int
printBlockEncPosDiags(const EISeq *seq, Seqpos pos, FILE *fp, EISHint hint,
                      Env *env)
{
  const struct blockCompositionSeq *seqIdx;
  int outCount = 0;
  Seqpos bucketNum;
  assert(seq && seq->classInfo == &blockCompositionSeqClass);
  seqIdx = constEncIdxSeq2blockCompositionSeq(seq);
  bucketNum = bucketNumFromPos(seqIdx, pos);
  fputs("##################################################\n"
        "# This bucket:\n"
        "##################################################\n", fp);
  outCount += printBucket(seqIdx, bucketNum, fp, hint, env);
  if (bucketNum)
  {
    fputs("##################################################\n"
          "# Previous bucket:\n"
          "##################################################\n", fp);
    outCount += printBucket(seqIdx, bucketNum - 1, fp, hint, env);
  }
  if (bucketNum < bucketNumFromPos(seqIdx, EISLength(seq)))
  {
    fputs("##################################################\n"
          "# Next bucket:\n"
          "##################################################\n", fp);
    outCount += printBucket(seqIdx, bucketNum + 1, fp, hint, env);
  }
  return outCount;
}

static const struct encIdxSeqClass blockCompositionSeqClass =
{
  .delete = deleteBlockEncIdxSeq,
  .rank = blockCompSeqRank,
  .select = blockCompSeqSelect,
  .get = blockCompSeqGet,
  .newHint = newBlockCompSeqHint,
  .deleteHint = deleteBlockCompSeqHint,
  .expose = blockCompSeqExpose,
  .seekToHeader = seekToHeader,
  .printPosDiags = printBlockEncPosDiags,
};
