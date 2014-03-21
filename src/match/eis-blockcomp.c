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
 * \file eis-blockcomp.c
 * \brief Methods to build block-compressed representation of indexed
 * sequence and answer queries on said representation.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */
/*
 * TODO:
 * - normalize use  of  seqIdx variable naming (seq, bseq etc.)
 */
#include <stddef.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#ifndef S_SPLINT_S
#include <unistd.h>
#endif
#include "core/assert_api.h"
#include "core/bitpackstring.h"
#include "match/dataalign.h"
#include "core/error.h"
#include "core/fa.h"
#include "core/log.h"
#include "core/minmax.h"
#include "core/str.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "eis-blockcomp-construct.h"
#include "match/eis-bitpackseqpos.h"
#include "match/eis-encidxseq.h"
#include "match/eis-encidxseq-priv.h"
#include "match/eis-seqranges.h"
#include "match/eis-seqblocktranslate.h"
#include "match/eis-seqdatasrc.h"

struct onDiskBlockCompIdx
{
  FILE *idxFP;                 /**< used for writing the index and
                                * if an mmap of the whole index
                                * isn't possible also for reading */
  GtStr *idxFN;                /**< stores the filename of idxFP */
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
  unsigned bucketBlocks, blockSize, callBackDataOffsetBits, bitsPerUlong,
    bitsPerVarDiskOffset;
  AlphabetRangeSize blockMapAlphabetSize;
  Symbol blockEncFallback, rangeEncFallback;
  int numModes;
  unsigned *partialSymSumBits, *partialSymSumBitsSums, symSumBits;
};

static inline size_t
blockEncIdxSeqHeaderLength(const struct blockCompositionSeq *seqIdx,
                           size_t numExtHeaders,
                           const uint32_t *extHeaderSizes);

static int
openOnDiskData(const char *projectName, struct onDiskBlockCompIdx *idx,
               char *mode, GtError *err);

static void
initOnDiskBlockCompIdx(struct onDiskBlockCompIdx *idx,
                       size_t headerLen, off_t cwLen);

static void
destructOnDiskBlockCompIdx(struct onDiskBlockCompIdx *idx);

static size_t
writeIdxHeader(struct blockCompositionSeq *seqIdx,
               size_t numExtHeaders, const uint16_t *headerIDs,
               const uint32_t *extHeaderSizes,
               headerWriteFunc *extHeaderCallbacks,
               void **headerCBData,
               GtError *err);

static inline int
tryMMapOfIndex(struct onDiskBlockCompIdx *idxData);

static const struct encIdxSeqClass blockCompositionSeqClass;

static inline struct blockCompositionSeq *
encIdxSeq2blockCompositionSeq(struct encIdxSeq *seq)
{
  gt_assert(seq && seq->classInfo == &blockCompositionSeqClass);
  return (struct blockCompositionSeq *)((char *)seq
    - offsetof(struct blockCompositionSeq, baseClass));
}

static inline const struct blockCompositionSeq *
constEncIdxSeq2blockCompositionSeq(const struct encIdxSeq *seq)
{
  gt_assert(seq && seq->classInfo == &blockCompositionSeqClass);
  return (struct blockCompositionSeq *)((char *)seq
    - offsetof(struct blockCompositionSeq, baseClass));
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
                const struct blockCompositionSeq *seqIdx);
static void
destructAppendState(struct appendState *aState);

static inline void
append2IdxOutput(struct appendState *state,
                 PermCompIndex permCompIdx[2],
                 unsigned bitsOfCompositionIdx, unsigned bitsOfPermutationIdx);

static BitOffset
appendCallBackOutput(struct appendState *state,
                     const struct blockCompositionSeq *seqIdx,
                     bitInsertFunc biFunc, GtUword start,
                     GtUword len, unsigned callBackDataOffsetBits,
                     void *cbState);

typedef GtUword partialSymSum;

static inline partialSymSum *
newPartialSymSums(AlphabetRangeSize alphabetSize);

static inline void
deletePartialSymSums(partialSymSum *sums);

static inline void
addBlock2PartialSymSums(partialSymSum *sums, const Symbol *block,
                        unsigned blockSize);

static inline void
copyPartialSymSums(AlphabetRangeSize alphabetSize, partialSymSum *dest,
                   const partialSymSum *src);

static inline GtUword
numBuckets(GtUword seqLen, size_t bucketLen);

static inline off_t
cwSize(const struct blockCompositionSeq *seqIdx);

static inline BitOffset
vwBits(GtUword seqLen, unsigned blockSize, unsigned bucketBlocks,
       unsigned maxPermIdxBits, varExtBitsEstimator biVarBits, void *cbState,
       struct varBitsEstimate *extVarBitsUpperBound);

static void
addRangeEncodedSyms(struct seqRangeList *rangeList, const Symbol *block,
                    unsigned blockSize, GtUword blockNum,
                    const MRAEnc *alphabet, int selection, const int *rangeSel);

static int
updateIdxOutput(struct blockCompositionSeq *seqIdx,
                struct appendState *aState,
                const partialSymSum *buck);

static int
finalizeIdxOutput(struct blockCompositionSeq *seqIdx,
                  struct appendState *state);

static inline void
symSumBitsDefaultSetup(struct blockCompositionSeq *seqIdx);

#define newBlockEncIdxSeqErrRet()                                       \
  do {                                                                  \
    if (newSeqIdx->externalData.idxFP)                                  \
      destructOnDiskBlockCompIdx(&newSeqIdx->externalData);             \
    if (newSeqIdx->compositionTable.bitsPerCount)                       \
      gt_destructCompositionList(&newSeqIdx->compositionTable);         \
    if (newSeqIdx->rangeEncs)                                           \
      gt_deleteSeqRangeList(newSeqIdx->rangeEncs);                      \
    if (newSeqIdx->extHeaderPos)                                        \
      gt_free(newSeqIdx->extHeaderPos);                                 \
    if (newSeqIdx->partialSymSumBits)                                   \
      gt_free(newSeqIdx->partialSymSumBits);                            \
    if (newSeqIdx->extHeaderPos)                                        \
      gt_free(newSeqIdx->extHeaderPos);                                 \
    if (newSeqIdx) gt_free(newSeqIdx);                                  \
    if (modesCopy)                                                      \
      gt_free(modesCopy);                                               \
    if (blockMapAlphabet) gt_MRAEncDelete(blockMapAlphabet);            \
    if (rangeMapAlphabet) gt_MRAEncDelete(rangeMapAlphabet);            \
    return NULL;                                                        \
  } while (0)

static void
addBlock2OutputBuffer(
  struct blockCompositionSeq *newSeqIdx,
  partialSymSum *buck, GtUword blockNum,
  Symbol *block, unsigned blockSize,
  const MRAEnc *alphabet, const int *modes,
  const MRAEnc *blockMapAlphabet, AlphabetRangeSize blockMapAlphabetSize,
  BitString permCompBSPreAlloc, unsigned *compositionPreAlloc,
  unsigned compositionIdxBits, struct appendState *aState)
{
  PermCompIndex permCompIdx[2];
  unsigned significantPermIdxBits;
  /* a. update superbucket table */
  addBlock2PartialSymSums(buck, block, blockSize);
  /* b. add ranges of differently encoded symbols to
   * corresponding representation */
  addRangeEncodedSyms(newSeqIdx->rangeEncs, block, blockSize,
                      blockNum, alphabet, REGIONS_LIST,
                      modes);
  /* c. add to table of composition/permutation indices */
  gt_MRAEncSymbolsTransform(blockMapAlphabet, block, blockSize);
  /* FIXME control remapping */
  /* currently invalid characters in input can seriously break this */
  gt_block2IndexPair(&newSeqIdx->compositionTable, blockSize,
                  blockMapAlphabetSize, block, permCompIdx,
                  &significantPermIdxBits,
                  permCompBSPreAlloc, compositionPreAlloc);
  append2IdxOutput(aState, permCompIdx, compositionIdxBits,
                   significantPermIdxBits);
}

static int
writeOutputBuffer(struct blockCompositionSeq *newSeqIdx,
                  struct appendState *aState, bitInsertFunc biFunc,
                  GtUword lastUpdatePos, size_t bucketLen,
                  unsigned callBackDataOffsetBits, void *cbState,
                  const partialSymSum *buckLast)
{
  if (biFunc)
    if (appendCallBackOutput(aState, newSeqIdx, biFunc,
                             lastUpdatePos, bucketLen,
                             callBackDataOffsetBits, cbState)
        == (BitOffset)-1)
    {
      perror("error condition while writing block-compressed"
             " index data");
      return -1;
    }
  if (!updateIdxOutput(newSeqIdx, aState, buckLast))
  {
    perror("error condition while writing block-compressed"
           " index data");
    return -1;
  }
  return 1;
}

#define newBlockEncIdxSeqLoopErr()                      \
  destructAppendState(&aState);                         \
  deletePartialSymSums(buck);                           \
  deletePartialSymSums(buckLast);                       \
  gt_free(compositionPreAlloc);                         \
  gt_free(permCompBSPreAlloc);                          \
  gt_free(block);                                       \
  break

EISeq *
gt_newGenBlockEncIdxSeq(GtUword totalLen, const char *projectName,
                        MRAEnc *alphabet, const struct seqStats *stats,
                        SeqDataReader BWTGenerator,
                        const struct seqBaseParam *params,
                        size_t numExtHeaders, const uint16_t *headerIDs,
                        const uint32_t *extHeaderSizes,
                        headerWriteFunc *extHeaderCallbacks,
                        void **headerCBData,
                        bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                        varExtBitsEstimator biVarBits, void *cbState,
                        GtError *err)
{
  struct blockCompositionSeq *newSeqIdx = NULL;
  AlphabetRangeSize blockMapAlphabetSize, totalAlphabetSize;
  size_t regionsEstimate=totalLen / 100;
  MRAEnc *blockMapAlphabet = NULL, *rangeMapAlphabet = NULL;
  BitOffset bitsPerComposition, bitsPerPermutation;
  unsigned compositionIdxBits, callBackDataOffsetBits,
    blockSize = params->encParams.blockEnc.blockSize,
    bucketBlocks = params->encParams.blockEnc.bucketBlocks;
  size_t bucketLen = (size_t)bucketBlocks * blockSize;
  int *modesCopy = NULL;
  struct varBitsEstimate biMaxExtSize;
  enum rangeStoreMode modes[] = { BLOCK_COMPOSITION_INCLUDE,
                                  REGIONS_LIST };
  gt_assert(projectName);
  gt_assert(params->encType == BWT_ON_BLOCK_ENC);
  gt_assert(blockSize > 0);
  gt_assert((biFunc && biVarBits) || (biFunc == NULL && biVarBits == NULL));
  gt_error_check(err);

  newSeqIdx = gt_calloc(sizeof (struct blockCompositionSeq), 1);
  newSeqIdx->bucketBlocks = bucketBlocks;
  newSeqIdx->bitsPerUlong = requiredUlongBits((newSeqIdx->baseClass.seqLen
                                                 = totalLen) - 1);
  newSeqIdx->baseClass.alphabet = alphabet;
  {
    size_t range, numAlphabetRanges = newSeqIdx->numModes =
      MRAEncGetNumRanges(alphabet);
    totalAlphabetSize = gt_MRAEncGetSize(alphabet);
    blockMapAlphabetSize = 0;
    newSeqIdx->modes = modesCopy = gt_malloc(sizeof (int) * numAlphabetRanges);
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
        gt_log_log("Invalid encoding request.\n");
        newBlockEncIdxSeqErrRet();
        break;
      }
    }
    newSeqIdx->blockMapAlphabetSize = blockMapAlphabetSize;
    newSeqIdx->blockMapAlphabet = blockMapAlphabet =
      gt_MRAEncSecondaryMapping(alphabet, BLOCK_COMPOSITION_INCLUDE, modesCopy,
                             newSeqIdx->blockEncFallback = 0);
    newSeqIdx->rangeMapAlphabet = rangeMapAlphabet =
      gt_MRAEncSecondaryMapping(alphabet, REGIONS_LIST, modesCopy,
                             newSeqIdx->rangeEncFallback = 0);
    gt_assert(gt_MRAEncGetSize(blockMapAlphabet) == blockMapAlphabetSize);
    gt_assert(gt_MRAEncGetSize(rangeMapAlphabet)
           == totalAlphabetSize - blockMapAlphabetSize);
  }
  newSeqIdx->partialSymSumBits
    = gt_malloc(sizeof (newSeqIdx->partialSymSumBits[0])
                * blockMapAlphabetSize * 2);
  newSeqIdx->partialSymSumBitsSums
    = newSeqIdx->partialSymSumBits + blockMapAlphabetSize;
  if (stats)
  {
    GtUword *symCounts;
    symCounts = gt_malloc(sizeof (symCounts[0]) * blockMapAlphabetSize);
    switch (stats->sourceAlphaType)
    {
    case sourceUInt8:
      memset(symCounts, 0,
             sizeof (symCounts[0]) * blockMapAlphabetSize);
      {
        Symbol eSym, bSym;
        unsigned i;
        for (i = 0; i <= UINT8_MAX; ++i)
          if (gt_MRAEncSymbolIsInSelectedRanges(
                alphabet, eSym = MRAEncMapSymbol(alphabet, i),
                BLOCK_COMPOSITION_INCLUDE, modesCopy)
              && ((bSym = MRAEncMapSymbol(blockMapAlphabet, eSym))
                  < blockMapAlphabetSize))
            symCounts[bSym]
              += stats->symbolDistributionTable[i];
#ifdef EIS_DEBUG
        for (i = 0; i < blockMapAlphabetSize; ++i)
        {
          gt_log_log("symCount[%"PRIuSymbol"]="GT_WU"\n", (Symbol) i,
                  stats->symbolDistributionTable[i]);
        }
#endif /* EIS_DEBUG */
        if (blockMapAlphabetSize)
        {
          newSeqIdx->partialSymSumBitsSums[0] = 0;
          newSeqIdx->partialSymSumBits[0] = requiredUlongBits(symCounts[0]);
          for (i = 1; i < blockMapAlphabetSize; ++i)
          {
            newSeqIdx->partialSymSumBitsSums[i]
              = newSeqIdx->partialSymSumBitsSums[i - 1]
              + newSeqIdx->partialSymSumBits[i - 1];
            newSeqIdx->partialSymSumBits[i]
              = requiredUlongBits(symCounts[i]);
          }
#ifdef EIS_DEBUG
          for (i = 0; i < blockMapAlphabetSize; ++i)
          {
            gt_log_log("bitsPerSymSum[%"PRIuSymbol"]=%u\n", (Symbol)i,
                    newSeqIdx->partialSymSumBits[i]);
          }
#endif  /* EIS_DEBUG */
          newSeqIdx->symSumBits
            = newSeqIdx->partialSymSumBitsSums[blockMapAlphabetSize - 1]
            + newSeqIdx->partialSymSumBits[blockMapAlphabetSize - 1];
#ifdef EIS_DEBUG
          gt_log_log("symSumBits total: %u\n", newSeqIdx->symSumBits);
#endif  /* EIS_DEBUG */
        }
      }
      /* count special characters to estimate number of regions required */
      {
        Symbol eSym, rSym;
        GtUword regionSymCount = 0;
        AlphabetRangeSize rangeMapAlphabetSize
          = gt_MRAEncGetSize(rangeMapAlphabet);
        unsigned i;
        for (i = 0; i <= UINT8_MAX; ++i)
          if (gt_MRAEncSymbolIsInSelectedRanges(
               alphabet, eSym = MRAEncMapSymbol(alphabet, i),
               REGIONS_LIST, modesCopy)
             && ((rSym = MRAEncMapSymbol(rangeMapAlphabet, eSym))
                 < rangeMapAlphabetSize))
            regionSymCount += stats->symbolDistributionTable[i];
        regionsEstimate = regionSymCount/20;
#ifdef EIS_DEBUG
        gt_log_log("Expected "GT_WU" symbols to encode in regions.\n",
                    regionSymCount);
#endif
      }
      break;
    default:
      symSumBitsDefaultSetup(newSeqIdx);
    }
    gt_free(symCounts);
  }
  else
  {
    symSumBitsDefaultSetup(newSeqIdx);
  }
  {
    int regionFeatures = SRL_NO_FEATURES;
    if (params->EISFeatureSet & EIS_FEATURE_REGION_SUMS)
      regionFeatures |= SRL_PARTIAL_SYMBOL_SUMS;
    newSeqIdx->rangeEncs = gt_newSeqRangeList(regionsEstimate, rangeMapAlphabet,
                                           regionFeatures);
  }
  newSeqIdx->baseClass.classInfo = &blockCompositionSeqClass;
  if (!gt_initCompositionList(&newSeqIdx->compositionTable, blockSize,
                           blockMapAlphabetSize))
  {
    gt_error_set(err, "Insufficient memory for selected block size %u and "
              "alphabet size %u, try smaller block size?\n", blockSize,
              blockMapAlphabetSize);
    newBlockEncIdxSeqErrRet();
  }
  bitsPerComposition = newSeqIdx->compositionTable.bitsPerCount
    * blockMapAlphabetSize;
  compositionIdxBits = newSeqIdx->compositionTable.compositionIdxBits;
  bitsPerPermutation = newSeqIdx->compositionTable.bitsPerSymbol * blockSize;
  newSeqIdx->blockSize = blockSize;
  newSeqIdx->cwExtBitsPerBucket = cwExtBitsPerPos * bucketLen;
  newSeqIdx->callBackDataOffsetBits = callBackDataOffsetBits
    = biFunc ? gt_requiredUInt64Bits(newSeqIdx->compositionTable.maxPermIdxBits
                                  * bucketBlocks) : 0;
  {
    BitOffset maxVarBitsTotal =
      vwBits(totalLen, blockSize, bucketBlocks,
             newSeqIdx->compositionTable.maxPermIdxBits,
             biVarBits, cbState, &biMaxExtSize);
    newSeqIdx->bitsPerVarDiskOffset = gt_requiredUInt64Bits(maxVarBitsTotal);
  }
  newSeqIdx->maxVarExtBitsPerBucket = biMaxExtSize.maxBitsPerBucket;
  {
    size_t headerLen = blockEncIdxSeqHeaderLength(newSeqIdx, numExtHeaders,
                                                  extHeaderSizes);
    if (!openOnDiskData(projectName, &newSeqIdx->externalData, "wb+",err))
      newBlockEncIdxSeqErrRet();
    initOnDiskBlockCompIdx(&newSeqIdx->externalData,
                           headerLen, cwSize(newSeqIdx));
  }
  /* At this point everything should be ready to receive the actual sequence
   * information, steps: */
  {
    int hadGtError = 0;
    do
    {
      {
        Symbol *block;
        unsigned *compositionPreAlloc;
        BitString permCompBSPreAlloc;
        partialSymSum *buck, *buckLast;
        block = gt_malloc(sizeof (Symbol) * blockSize);
        compositionPreAlloc = gt_malloc(sizeof (compositionPreAlloc[0])
                                        * blockMapAlphabetSize);
        permCompBSPreAlloc =
          gt_malloc(bitElemsAllocSize(bitsPerComposition + bitsPerPermutation)
                    * sizeof (BitElem));
        buck = newPartialSymSums(totalAlphabetSize);
        buckLast = newPartialSymSums(totalAlphabetSize);
        /* 2. read block sized chunks from bwttab and suffix array */
        {
          GtUword numFullBlocks = totalLen / blockSize, blockNum,
            lastUpdatePos = 0;
          /* pos == totalLen - symbolsLeft */
          struct appendState aState;
          initAppendState(&aState, newSeqIdx);
          blockNum = 0;
          while (blockNum < numFullBlocks)
          {
            size_t readResult;
            /* 3. for each chunk: */
            readResult = SDRRead(BWTGenerator, block, blockSize);
            if (readResult != blockSize)
            {
              hadGtError = 1;
              perror("error condition while reading index data");
              break;
            }
            gt_MRAEncSymbolsTransform(alphabet, block, blockSize);
            addBlock2OutputBuffer(newSeqIdx, buck, blockNum,
                                  block, blockSize,
                                  alphabet, modesCopy,
                                  blockMapAlphabet, blockMapAlphabetSize,
                                  permCompBSPreAlloc, compositionPreAlloc,
                                  compositionIdxBits, &aState);
            /* update on-disk structure */
            if (!((++blockNum) % bucketBlocks))
            {
              GtUword pos = blockNum * blockSize;
              if (writeOutputBuffer(newSeqIdx, &aState, biFunc, lastUpdatePos,
                                    bucketLen, callBackDataOffsetBits, cbState,
                                    buckLast) < 0)
              {
                hadGtError = 1;
                break;
              }
              /* update retained data */
              copyPartialSymSums(totalAlphabetSize, buckLast, buck);
              lastUpdatePos = pos;
            }
          }
          /* handle last chunk */
          if (!hadGtError)
          {
            GtUword symbolsLeft = totalLen % blockSize;
            if (symbolsLeft)
            {
              size_t readResult;
              readResult = SDRRead(BWTGenerator, block, symbolsLeft);
              if (readResult < symbolsLeft)
              {
                hadGtError = 1;
                perror("error condition while reading index data");
                newBlockEncIdxSeqLoopErr();
              }
              else
              {
                gt_MRAEncSymbolsTransform(alphabet, block, symbolsLeft);
                memset(block + symbolsLeft, 0,
                       sizeof (Symbol) * (blockSize - symbolsLeft));
                addBlock2OutputBuffer(newSeqIdx, buck, blockNum,
                                      block, blockSize,
                                      alphabet, modesCopy,
                                      blockMapAlphabet, blockMapAlphabetSize,
                                      permCompBSPreAlloc, compositionPreAlloc,
                                      compositionIdxBits, &aState);
              }
            }
            if (lastUpdatePos <= totalLen)
            {
              /* one bucket still unfinished */
              if (writeOutputBuffer(newSeqIdx, &aState, biFunc, lastUpdatePos,
                                    totalLen - lastUpdatePos,
                                    callBackDataOffsetBits, cbState,
                                    buckLast) < 0)
              {
                hadGtError = 1;
                break;
              }
            }
            if (!finalizeIdxOutput(newSeqIdx, &aState))
            {
              hadGtError = 1;
              perror("error condition while writing block-compressed"
                     " index data");
              newBlockEncIdxSeqLoopErr();
            }
            if (!writeIdxHeader(newSeqIdx, numExtHeaders, headerIDs,
                                extHeaderSizes, extHeaderCallbacks,
                                headerCBData, err))
            {
              hadGtError = 1;
              perror("error condition while writing block-compressed"
                     " index header");
              newBlockEncIdxSeqLoopErr();
            }
            if (fflush(newSeqIdx->externalData.idxFP))
            {
              hadGtError = 1;
              perror("error condition while writing block-compressed"
                     " index header");
              newBlockEncIdxSeqLoopErr();
            }
            tryMMapOfIndex(&newSeqIdx->externalData);
          }
          if (hadGtError)
          {
            newBlockEncIdxSeqLoopErr();
          }
          /* 4. dealloc resources no longer required */
          destructAppendState(&aState);
          deletePartialSymSums(buck);
          deletePartialSymSums(buckLast);
        }
        /* 4. dealloc resources no longer required */
        gt_free(compositionPreAlloc);
        gt_free(permCompBSPreAlloc);
        gt_free(block);
      }
    } while (0);
    /* close bwttab and suffix array */
    if (hadGtError)
      newBlockEncIdxSeqErrRet();
  }
  return &(newSeqIdx->baseClass);
}

static void
deleteBlockEncIdxSeq(struct encIdxSeq *seq)
{
  struct blockCompositionSeq *bseq;
  gt_assert(seq && seq->classInfo == &blockCompositionSeqClass);
  bseq = encIdxSeq2blockCompositionSeq(seq);
  gt_free(bseq->extHeaderPos);
  gt_free(bseq->partialSymSumBits);
  destructOnDiskBlockCompIdx(&bseq->externalData);
  gt_destructCompositionList(&bseq->compositionTable);
  gt_MRAEncDelete(bseq->baseClass.alphabet);
  gt_MRAEncDelete(bseq->rangeMapAlphabet);
  gt_MRAEncDelete(bseq->blockMapAlphabet);
  gt_deleteSeqRangeList(bseq->rangeEncs);
  gt_free(bseq->modes);
  gt_free(bseq);
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
newEmptySuperBlock(const struct blockCompositionSeq *seqIdx)
{
  struct superBlock *sBlock;
  sBlock = gt_malloc(superBlockMemSize(seqIdx));
  initEmptySuperBlock(sBlock, seqIdx);
  return sBlock;
}

static void
deleteSuperBlock(struct superBlock *sBlock)
{
  gt_free(sBlock);
}

static inline void
symSumBitsDefaultSetup(struct blockCompositionSeq *seqIdx)
{
  unsigned i;
  AlphabetRangeSize blockMapAlphabetSize = seqIdx->blockMapAlphabetSize;
  seqIdx->partialSymSumBitsSums[0] = 0;
  seqIdx->partialSymSumBits[0] = seqIdx->bitsPerUlong;
  for (i = 1; i < blockMapAlphabetSize; ++i)
    seqIdx->partialSymSumBitsSums[i] = seqIdx->partialSymSumBitsSums[i - 1]
      + (seqIdx->partialSymSumBits[i] = seqIdx->bitsPerUlong);
  seqIdx->symSumBits = blockMapAlphabetSize * seqIdx->bitsPerUlong;
#ifdef EIS_DEBUG
  gt_log_log("symSumBits=%u, blockMapAlphabetSize=%u\n",
          seqIdx->symSumBits, seqIdx->blockMapAlphabetSize);
#endif
  gt_assert(seqIdx->partialSymSumBitsSums[i - 1] + seqIdx->bitsPerUlong
         == seqIdx->symSumBits);
}

static inline GtUword
sBlockGetPartialSymSum(struct superBlock *sBlock, Symbol sym,
                       const struct blockCompositionSeq *seqIdx)
{
  return gt_bsGetUlong(sBlock->cwData, seqIdx->partialSymSumBitsSums[sym]
                     + sBlock->cwIdxMemBase, seqIdx->partialSymSumBits[sym]);
}

static inline void
sBlockGetPartialSymSums(struct superBlock *sBlock,
                        const struct blockCompositionSeq *seqIdx,
                        GtUword *sums)
{
  gt_bsGetNonUniformUlongArray(
    sBlock->cwData, sBlock->cwIdxMemBase, seqIdx->blockMapAlphabetSize,
    seqIdx->symSumBits, seqIdx->partialSymSumBits,
#if defined (_LP64) || defined (_WIN64)
                             (uint64_t*) sums);
#else
                             (uint32_t*) sums);
#endif
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
  return gt_bsGetUInt64(sBlock->cwData, offset, seqIdx->bitsPerVarDiskOffset);
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
  return gt_bsGetUInt64(sBlock->cwData, offset, seqIdx->callBackDataOffsetBits);
}

static inline GtUword
blockNumFromPos(const struct blockCompositionSeq *seqIdx, GtUword pos)
{
  return pos / seqIdx->blockSize;
}

static inline GtUword
bucketNumFromPos(const struct blockCompositionSeq *seqIdx, GtUword pos)
{
  return pos / (seqIdx->bucketBlocks * seqIdx->blockSize);
}

static inline GtUword
bucketBasePos(const struct blockCompositionSeq *seqIdx, GtUword bucketNum)
{
  return bucketNum * seqIdx->bucketBlocks * seqIdx->blockSize;
}

static inline GtUword
bucketNumFromBlockNum(struct blockCompositionSeq *seqIdx,
                      GtUword blockNum)
{
  return blockNum / seqIdx->bucketBlocks;
}

#define fetchSuperBlockErrRet()                    \
  do {                                             \
    if (!sBlockPreAlloc) deleteSuperBlock(retval); \
    return NULL;                                   \
  } while (0)

static struct superBlock *
fetchSuperBlock(const struct blockCompositionSeq *seqIdx,
                GtUword bucketNum, struct superBlock *sBlockPreAlloc)
{
  struct superBlock *retval = NULL;
  gt_assert(seqIdx);
  gt_assert(bucketNum * seqIdx->bucketBlocks * seqIdx->blockSize
         <= seqIdx->baseClass.seqLen);
  if (sBlockPreAlloc)
    retval = sBlockPreAlloc;
  else
    retval = newEmptySuperBlock(seqIdx);
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
    GT_UNUSED int ret;
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
    ret = fread(retval->varData, sizeof (BitElem),
                superBlockVarMaxReadSize(seqIdx),
                idxFP);
    if (ferror(idxFP))
      fetchSuperBlockErrRet();
  }
  return retval;
}

static GtUword
blockCompSeqSelect(GT_UNUSED struct encIdxSeq *seq, GT_UNUSED Symbol sym,
                   GT_UNUSED GtUword count, GT_UNUSED union EISHint *hint)
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
pos2CachePos(struct seqCache *sCache, GtUword pos)
{
  return pos % sCache->numEntries;
}

static void
initSuperBlockSeqCache(struct seqCache *sBlockCache,
                       const struct blockCompositionSeq *seqIdx,
                       size_t numEntries)
{
  unsigned bucketBlocks, GT_UNUSED bucketLen, blockSize;
  size_t superBlockSize;

  gt_assert(seqIdx && sBlockCache);
  blockSize = seqIdx->blockSize;
  bucketLen = (bucketBlocks = seqIdx->bucketBlocks) * blockSize;
  superBlockSize = superBlockMemSize(seqIdx);
  sBlockCache->numEntries = numEntries;
  {
    void *temp = gt_malloc((sizeof (GtUword) + superBlockSize
                           + sizeof (void *)) * numEntries);
    sBlockCache->cachedPos = temp;
    sBlockCache->entriesPtr = (void **)((char *)temp
                                        + sizeof (GtUword) * numEntries);
    sBlockCache->entries =
                      (char *)temp + (sizeof (GtUword) + sizeof (void *))
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
destructSuperBlockSeqCache(struct seqCache *sBlockCache)
{
  gt_free(sBlockCache->cachedPos);
}

static inline int
inSeqCache(struct seqCache *sCache, GtUword pos)
{
  return sCache->cachedPos[pos2CachePos(sCache, pos)] == pos;
}

#ifdef USE_SBLOCK_CACHE
static struct superBlock *
cacheFetchSuperBlock(const struct blockCompositionSeq *seqIdx,
                     GtUword superBlockNum, struct seqCache *sBlockCache)
{
  size_t cachePos = pos2CachePos(sBlockCache, superBlockNum);
  if (inSeqCache(sBlockCache, superBlockNum))
    return (struct superBlock *)sBlockCache->entriesPtr[cachePos];
  else
  {
    struct superBlock *sb =
      (struct superBlock *)sBlockCache->entriesPtr[cachePos];
    return fetchSuperBlock(seqIdx, superBlockNum, sb);
  }
}
#endif /* USE_SBLOCK_CACHE */

#define walkCompIndices(seqIdx, sBlock, numBlocks, cwOffset,              \
                        codeForCompIndex, varOffset)                      \
  do {                                                                    \
    unsigned blocksLeft = (numBlocks),                                    \
      bitsPerCompositionIdx =                                             \
      (seqIdx)->compositionTable.compositionIdxBits;                      \
    while (blocksLeft)                                                    \
    {                                                                     \
      PermCompIndex compIndex;                                            \
      compIndex = gt_bsGetPermCompIndex((sBlock)->cwData, cwIdxMemOffset, \
                                     bitsPerCompositionIdx);              \
      codeForCompIndex;                                                   \
      (varOffset) +=                                                      \
        (seqIdx)->compositionTable.permutations[compIndex].permIdxBits;   \
      (cwOffset) += bitsPerCompositionIdx;                                \
      --blocksLeft;                                                       \
    }                                                                     \
  } while (0)

#define walkCompIndicesPrefix(seqIdx, sBlock, blockNum, cwOffset,       \
                              codeForCompIndex, varOffset)              \
  do {                                                                  \
    cwOffset = sBlockGetCompIdxOffset(sBlock, seqIdx, 0);               \
    varOffset = sBlock->varDataMemBase;                                 \
    walkCompIndices(seqIdx, sBlock, blockNum, cwOffset,                 \
                    codeForCompIndex, varOffset);                       \
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
  compIndex = gt_bsGetPermCompIndex(sBlock->cwData, cwOffset,
                                 bitsPerCompositionIdx);
  varIdxBits = seqIdx->compositionTable.permutations[compIndex].permIdxBits;
  permIndex = gt_bsGetPermCompIndex(sBlock->varData, varOffset, varIdxBits);
  indexPair2block(&seqIdx->compositionTable, seqIdx->blockSize,
                  compIndex, permIndex, block, sublen);
}

/*
 * regular, user-accessible query functions
 */
static Symbol *
blockCompSeqGetBlock(struct blockCompositionSeq *seqIdx, GtUword blockNum,
                     struct blockEncIdxSeqHint *hint, int queryRangeEnc,
                     struct superBlock *sBlockPreFetch,
                     Symbol *blockPA)
{
  struct superBlock *sBlock;
  BitOffset varDataMemOffset, cwIdxMemOffset;
  GT_UNUSED GtUword relBlockNum;
  unsigned blockSize;
  Symbol *block;
  gt_assert(seqIdx);
  blockSize = seqIdx->blockSize;
  if (blockNum * blockSize >= EISLength(&seqIdx->baseClass))
    return NULL;
  if (sBlockPreFetch)
    sBlock = sBlockPreFetch;
  else
  {
    GtUword bucketNum = bucketNumFromBlockNum(seqIdx, blockNum);
#ifdef USE_SBLOCK_CACHE
    sBlock = cacheFetchSuperBlock(seqIdx, bucketNum,
                                  &hint->sBlockCache);
#else
    sBlock = fetchSuperBlock(seqIdx, bucketNum, NULL);
#endif
  }
  relBlockNum = blockNum % seqIdx->bucketBlocks;
  if (blockPA)
    block = blockPA;
  else
    block = gt_calloc(sizeof (Symbol), blockSize);
  walkCompIndicesPrefix(seqIdx, sBlock, blockNum % seqIdx->bucketBlocks,
                        cwIdxMemOffset, , varDataMemOffset);
  unpackBlock(seqIdx, sBlock, cwIdxMemOffset, varDataMemOffset, block,
              blockSize);
  if (queryRangeEnc)
    gt_SRLApplyRangesToSubString(seqIdx->rangeEncs,
                              block, blockNum * blockSize, blockSize,
                              blockNum * blockSize, &hint->rangeHint);
#ifndef USE_SBLOCK_CACHE
  if (!sBlockPreFetch)
    deleteSuperBlock(sBlock);
#endif
  return block;
}

static inline GtUword
adjustPosRankForBlock(struct blockCompositionSeq *seqIdx,
                      struct superBlock *sBlock, GtUword pos, Symbol bSym,
                      unsigned blockSize, GtUword preBlockRankCount,
                      BitOffset cwIdxMemOffset, BitOffset varDataMemOffset,
                      unsigned bitsPerCompositionIdx)
{
  GtUword rankCount = preBlockRankCount;
  unsigned inBlockPos;
  if ((inBlockPos = pos % blockSize)
      && symCountFromComposition(
        &seqIdx->compositionTable, seqIdx->blockMapAlphabetSize,
        gt_bsGetPermCompIndex(sBlock->cwData, cwIdxMemOffset,
                           bitsPerCompositionIdx), bSym))
  {
    Symbol block[blockSize];
    unsigned i;
    unpackBlock(seqIdx, sBlock, cwIdxMemOffset, varDataMemOffset, block,
                inBlockPos);
    for (i = 0; i < inBlockPos; ++i)
    {
      if (block[i] == bSym)
        ++rankCount;
    }
  }
  return rankCount;
}

/* Note: pos is meant exclusively, i.e. returns 0
   for any query where pos==0 because that corresponds to the empty prefix */
static GtUword
blockCompSeqRank(struct encIdxSeq *eSeqIdx, Symbol eSym, GtUword pos,
                 union EISHint *hint)
{
  struct blockCompositionSeq *seqIdx;
  GtUword rankCount;
  gt_assert(eSeqIdx && eSeqIdx->classInfo == &blockCompositionSeqClass);
  seqIdx = encIdxSeq2blockCompositionSeq(eSeqIdx);
  gt_assert(gt_MRAEncSymbolIsInSelectedRanges(seqIdx->baseClass.alphabet,
                                        eSym, BLOCK_COMPOSITION_INCLUDE,
                                        seqIdx->modes) >= 0);
  if (gt_MRAEncSymbolIsInSelectedRanges(seqIdx->baseClass.alphabet, eSym,
                                     BLOCK_COMPOSITION_INCLUDE, seqIdx->modes))
  {
    BitOffset varDataMemOffset, cwIdxMemOffset;
    struct superBlock *sBlock;
    Symbol bSym = MRAEncMapSymbol(seqIdx->blockMapAlphabet, eSym);
    GtUword blockNum, bucketNum;
    unsigned blockSize = seqIdx->blockSize, bitsPerCompositionIdx
      = seqIdx->compositionTable.compositionIdxBits;
    bucketNum = bucketNumFromPos(seqIdx, pos);
#ifdef USE_SBLOCK_CACHE
    sBlock = cacheFetchSuperBlock(seqIdx, bucketNum,
                                  &hint->bcHint.sBlockCache);
#else
    sBlock = fetchSuperBlock(seqIdx, bucketNum, NULL);
#endif
    rankCount = sBlockGetPartialSymSum(sBlock, bSym, seqIdx);
    blockNum = blockNumFromPos(seqIdx, pos);
    walkCompIndicesPrefix(
      seqIdx, sBlock, blockNum % seqIdx->bucketBlocks, cwIdxMemOffset,
      rankCount += symCountFromComposition(
        &seqIdx->compositionTable, seqIdx->blockMapAlphabetSize, compIndex,
        bSym);,
      varDataMemOffset);
    rankCount = adjustPosRankForBlock(
      seqIdx, sBlock, pos, bSym, blockSize, rankCount, cwIdxMemOffset,
      varDataMemOffset, bitsPerCompositionIdx);
    if (bSym == seqIdx->blockEncFallback)
    {
      GtUword base = bucketBasePos(seqIdx, bucketNum);
      rankCount -= gt_SRLAllSymbolsCountInSeqRegion(
        seqIdx->rangeEncs, base, pos, &hint->bcHint.rangeHint);
    }
#ifndef USE_SBLOCK_CACHE
    deleteSuperBlock(sBlock);
#endif
  }
  else
  {
    rankCount = gt_SRLSymbolCountInSeqRegion(seqIdx->rangeEncs, 0, pos, eSym,
                                          &hint->bcHint.rangeHint);
  }
  return rankCount;
}

/* Note: pos is meant exclusively, i.e. returns 0
   for any query where pos==0 because that corresponds to the empty prefix */
static GtUlongPair
blockCompSeqPosPairRank(struct encIdxSeq *eSeqIdx, Symbol eSym,
                        GtUword posA, GtUword posB,
                        union EISHint *hint)
{
  struct blockCompositionSeq *seqIdx;
  GtUlongPair rankCounts;
  GtUword bucketNum;
  gt_assert(eSeqIdx && eSeqIdx->classInfo == &blockCompositionSeqClass);
  gt_assert(posA <= posB);
  seqIdx = encIdxSeq2blockCompositionSeq(eSeqIdx);
  gt_assert(gt_MRAEncSymbolIsInSelectedRanges(seqIdx->baseClass.alphabet,
                                        eSym, BLOCK_COMPOSITION_INCLUDE,
                                        seqIdx->modes) >= 0);
  /* Only when both positions are in same bucket, special treatment
   * makes sense. */
  {
    GtUword bucketNumA, bucketNumB;
    bucketNumA = bucketNumFromPos(seqIdx, posA);
    bucketNumB = bucketNumFromPos(seqIdx, posB);
    if (bucketNumA != bucketNumB)
    {
      rankCounts.a = blockCompSeqRank(eSeqIdx, eSym, posA, hint);
      rankCounts.b = blockCompSeqRank(eSeqIdx, eSym, posB, hint);
      return rankCounts;
    }
    bucketNum = bucketNumA;
  }
  if (gt_MRAEncSymbolIsInSelectedRanges(seqIdx->baseClass.alphabet, eSym,
                                     BLOCK_COMPOSITION_INCLUDE, seqIdx->modes))
  {
    BitOffset varDataMemOffset, cwIdxMemOffset;
    struct superBlock *sBlock;
    GtUword blockNumA = blockNumFromPos(seqIdx, posA),
      blockNumB = blockNumFromPos(seqIdx, posB);
    Symbol bSym = MRAEncMapSymbol(seqIdx->blockMapAlphabet, eSym);
    unsigned blockSize = seqIdx->blockSize, bitsPerCompositionIdx
      = seqIdx->compositionTable.compositionIdxBits,
      inBucketBlockNumA = blockNumA % seqIdx->bucketBlocks,
      inBucketBlockNumB = blockNumB % seqIdx->bucketBlocks;
#ifdef USE_SBLOCK_CACHE
    sBlock = cacheFetchSuperBlock(seqIdx, bucketNum,
                                  &hint->bcHint.sBlockCache);
#else
    sBlock = fetchSuperBlock(seqIdx, bucketNum, NULL);
#endif
    rankCounts.a = sBlockGetPartialSymSum(sBlock, bSym, seqIdx);
    walkCompIndicesPrefix(
      seqIdx, sBlock, inBucketBlockNumA, cwIdxMemOffset,
      rankCounts.a += symCountFromComposition(
        &seqIdx->compositionTable, seqIdx->blockMapAlphabetSize, compIndex,
        bSym),
      varDataMemOffset);
    rankCounts.b = rankCounts.a;
    rankCounts.a =
      adjustPosRankForBlock(seqIdx, sBlock, posA, bSym, blockSize,
                            rankCounts.a, cwIdxMemOffset, varDataMemOffset,
                            bitsPerCompositionIdx);
    walkCompIndices(
      seqIdx, sBlock, inBucketBlockNumB - inBucketBlockNumA, cwIdxMemOffset,
      rankCounts.b += symCountFromComposition(
        &seqIdx->compositionTable, seqIdx->blockMapAlphabetSize, compIndex,
        bSym),
      varDataMemOffset);
    rankCounts.b =
      adjustPosRankForBlock(seqIdx, sBlock, posB, bSym, blockSize,
                            rankCounts.b, cwIdxMemOffset, varDataMemOffset,
                            bitsPerCompositionIdx);
    if (bSym == seqIdx->blockEncFallback)
    {
      GtUword base = bucketBasePos(seqIdx, bucketNum);
      rankCounts.a -= gt_SRLAllSymbolsCountInSeqRegion(
        seqIdx->rangeEncs, base, posA, &hint->bcHint.rangeHint);
      rankCounts.b -= gt_SRLAllSymbolsCountInSeqRegion(
        seqIdx->rangeEncs, base, posB, &hint->bcHint.rangeHint);
    }
#ifndef USE_SBLOCK_CACHE
    deleteSuperBlock(sBlock);
#endif
  }
  else
  {
    rankCounts.a = gt_SRLSymbolCountInSeqRegion(seqIdx->rangeEncs, 0, posA,
                                                eSym, &hint->bcHint.rangeHint);
    rankCounts.b = gt_SRLSymbolCountInSeqRegion(seqIdx->rangeEncs, 0, posA,
                                                eSym, &hint->bcHint.rangeHint);
  }
  return rankCounts;
}

static void
blockCompSeqExpose(struct encIdxSeq *eSeqIdx, GtUword pos, int flags,
                   struct extBitsRetrieval *retval, union EISHint *hint)
{
  struct blockCompositionSeq *seqIdx;
  gt_assert(eSeqIdx && eSeqIdx->classInfo == &blockCompositionSeqClass);
  gt_assert(retval);
  seqIdx = encIdxSeq2blockCompositionSeq(eSeqIdx);
  {
    struct superBlock *sBlock;
    GtUword bucketNum, end;
    bucketNum = bucketNumFromPos(seqIdx, pos);
    retval->start = bucketBasePos(seqIdx, bucketNum);
    end = bucketBasePos(seqIdx, bucketNum + 1);
    if (end >= seqIdx->baseClass.seqLen)
      end = seqIdx->baseClass.seqLen;
    retval->len = end - retval->start;
#ifdef USE_SBLOCK_CACHE
    sBlock = cacheFetchSuperBlock(seqIdx, bucketNum,
                                  &hint->bcHint.sBlockCache);
#else
    sBlock = fetchSuperBlock(seqIdx, bucketNum, NULL);
#endif
    retval->cwOffset = sBlockCWExtBitsOffset(sBlock, seqIdx);
    if (flags & EBRF_PERSISTENT_CWBITS)
    {
      if (!(retval->flags & EBRF_PERSISTENT_CWBITS))
      {
          retval->varPart = gt_malloc(superBlockCWMaxReadSize(seqIdx));
      }
      memcpy(retval->cwPart, sBlock->cwData,
             superBlockCWMaxReadSize(seqIdx));
    }
    else
    {
      if (retval->flags & EBRF_PERSISTENT_CWBITS)
        gt_free(retval->cwPart);
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
          retval->varPart = gt_malloc(superBlockVarMaxReadSize(seqIdx));
        }
        memcpy(retval->varPart, sBlock->varData,
               superBlockVarMaxReadSize(seqIdx));
      }
      else
      {
        if (retval->flags & EBRF_PERSISTENT_VARBITS)
          gt_free(retval->varPart);
        retval->varPart = sBlock->varData;
      }
    }
    else
    {
      if (retval->flags & EBRF_PERSISTENT_VARBITS)
        gt_free(retval->varPart);
      retval->varPart = NULL;
      retval->varOffset = 0;
    }
#ifndef USE_SBLOCK_CACHE
    deleteSuperBlock(sBlock);
#endif
    retval->flags = flags;
  }
}

static inline void
adjustRanksForBlock(struct blockCompositionSeq *seqIdx,
                    struct superBlock *sBlock, GtUword pos,
                    unsigned blockSize, GtUword rankCounts[],
                    BitOffset cwIdxMemOffset, BitOffset varDataMemOffset)
{
  unsigned inBlockPos = pos % blockSize;
  if (inBlockPos)
  {
    Symbol block[blockSize];
    unsigned i;
    unpackBlock(seqIdx, sBlock, cwIdxMemOffset, varDataMemOffset, block,
                inBlockPos);
    for (i = 0; i < inBlockPos; ++i)
      ++(rankCounts[block[i]]);
  }
}

static void
blockCompSeqRangeRank(struct encIdxSeq *eSeqIdx, AlphabetRangeID range,
                      GtUword pos, GtUword *rankCounts,
                      union EISHint *hint)
{
  struct blockCompositionSeq *seqIdx;
  gt_assert(eSeqIdx && eSeqIdx->classInfo == &blockCompositionSeqClass);
  seqIdx = encIdxSeq2blockCompositionSeq(eSeqIdx);
  gt_assert(range < MRAEncGetNumRanges(EISGetAlphabet(eSeqIdx)));
  switch (seqIdx->modes[range])
  {
  case BLOCK_COMPOSITION_INCLUDE:
    {
      BitOffset varDataMemOffset, cwIdxMemOffset;
      struct superBlock *sBlock;
      GtUword blockNum, bucketNum;
      unsigned blockSize = seqIdx->blockSize;
      bucketNum = bucketNumFromPos(seqIdx, pos);
#ifdef USE_SBLOCK_CACHE
      sBlock = cacheFetchSuperBlock(seqIdx, bucketNum,
                                    &hint->bcHint.sBlockCache);
#else
      sBlock = fetchSuperBlock(seqIdx, bucketNum, NULL);
#endif
      {
        sBlockGetPartialSymSums(sBlock, seqIdx, rankCounts);
        blockNum = blockNumFromPos(seqIdx, pos);
        walkCompIndicesPrefix(
          seqIdx, sBlock, blockNum % seqIdx->bucketBlocks, cwIdxMemOffset,
          addSymCountsFromComposition(&seqIdx->compositionTable,
                                      seqIdx->blockMapAlphabetSize, compIndex,
                                      rankCounts);,
          varDataMemOffset);
        adjustRanksForBlock(seqIdx, sBlock, pos, blockSize, rankCounts,
                            cwIdxMemOffset, varDataMemOffset);
        {
          GtUword base = bucketBasePos(seqIdx, bucketNum);
          rankCounts[seqIdx->blockEncFallback]
            -= gt_SRLAllSymbolsCountInSeqRegion(
              seqIdx->rangeEncs, base, pos, &hint->bcHint.rangeHint);
        }
      }
#ifndef USE_SBLOCK_CACHE
      deleteSuperBlock(sBlock);
#endif
    }
    break;
  case REGIONS_LIST:
    {
      AlphabetRangeSize sym,
        rangeEncNumSyms = gt_MRAEncGetSize(seqIdx->rangeMapAlphabet);
      for (sym = 0; sym < rangeEncNumSyms; ++sym)
        rankCounts[sym] = gt_SRLSymbolCountInSeqRegion(
          seqIdx->rangeEncs, 0, pos,
          MRAEncRevMapSymbol(seqIdx->rangeMapAlphabet, sym),
          &hint->bcHint.rangeHint);
    }
    break;
  }
}

static void
blockCompSeqPosPairRangeRank(struct encIdxSeq *eSeqIdx, AlphabetRangeID range,
                             GtUword posA, GtUword posB,
                             GtUword *rankCounts, union EISHint *hint)
{
  struct blockCompositionSeq *seqIdx;
  gt_assert(eSeqIdx && eSeqIdx->classInfo == &blockCompositionSeqClass);
  seqIdx = encIdxSeq2blockCompositionSeq(eSeqIdx);
  gt_assert(range < MRAEncGetNumRanges(EISGetAlphabet(eSeqIdx)));
  switch (seqIdx->modes[range])
  {
  case BLOCK_COMPOSITION_INCLUDE:
    {
      BitOffset varDataMemOffset, cwIdxMemOffset;
      struct superBlock *sBlock;
      AlphabetRangeSize rsize = MRAEncGetRangeSize(eSeqIdx->alphabet, range);
      /* Only when both positions are in same bucket, special treatment
       * makes sense. */
      GtUword bucketNum = bucketNumFromPos(seqIdx, posA);
      if (bucketNum != bucketNumFromPos(seqIdx, posB))
      {
        blockCompSeqRangeRank(eSeqIdx, range, posA, rankCounts, hint);
        blockCompSeqRangeRank(eSeqIdx, range, posB, rankCounts + rsize, hint);
        return;
      }
#ifdef USE_SBLOCK_CACHE
      sBlock = cacheFetchSuperBlock(seqIdx, bucketNum,
                                    &hint->bcHint.sBlockCache);
#else
      sBlock = fetchSuperBlock(seqIdx, bucketNum, NULL);
#endif
      {
        GtUword blockNumA = blockNumFromPos(seqIdx, posA),
          blockNumB = blockNumFromPos(seqIdx, posB);
        unsigned blockSize = seqIdx->blockSize,
          inBucketBlockNumA = blockNumA % seqIdx->bucketBlocks,
          inBucketBlockNumB = blockNumB % seqIdx->bucketBlocks;
        sBlockGetPartialSymSums(sBlock, seqIdx, rankCounts);
        walkCompIndicesPrefix(
          seqIdx, sBlock, inBucketBlockNumA, cwIdxMemOffset,
          addSymCountsFromComposition(&seqIdx->compositionTable,
                                      seqIdx->blockMapAlphabetSize, compIndex,
                                      rankCounts);,
          varDataMemOffset);
        memcpy(rankCounts + rsize, rankCounts, rsize * sizeof (rankCounts[0]));
        adjustRanksForBlock(seqIdx, sBlock, posA, blockSize, rankCounts,
                            cwIdxMemOffset, varDataMemOffset);
        walkCompIndices(
          seqIdx, sBlock, inBucketBlockNumB - inBucketBlockNumA, cwIdxMemOffset,
          addSymCountsFromComposition(&seqIdx->compositionTable,
                                      seqIdx->blockMapAlphabetSize, compIndex,
                                      rankCounts + rsize);,
          varDataMemOffset);
        adjustRanksForBlock(seqIdx, sBlock, posB, blockSize,
                            rankCounts + rsize,
                            cwIdxMemOffset, varDataMemOffset);
        {
          GtUword base = bucketBasePos(seqIdx, bucketNum);
          rankCounts[seqIdx->blockEncFallback]
            -= gt_SRLAllSymbolsCountInSeqRegion(
              seqIdx->rangeEncs, base, posA, &hint->bcHint.rangeHint);
          rankCounts[rsize + seqIdx->blockEncFallback]
            -= gt_SRLAllSymbolsCountInSeqRegion(
              seqIdx->rangeEncs, base, posB, &hint->bcHint.rangeHint);
        }
      }
#ifndef USE_SBLOCK_CACHE
      deleteSuperBlock(sBlock);
#endif
    }
    break;
  case REGIONS_LIST:
    {
      AlphabetRangeSize sym,
        rangeEncNumSyms = gt_MRAEncGetSize(seqIdx->rangeMapAlphabet);
      for (sym = 0; sym < rangeEncNumSyms; ++sym)
        rankCounts[sym] = gt_SRLSymbolCountInSeqRegion(
          seqIdx->rangeEncs, 0, posA,
          MRAEncRevMapSymbol(seqIdx->rangeMapAlphabet, sym),
          &hint->bcHint.rangeHint);
      for (sym = 0; sym < rangeEncNumSyms; ++sym)
        rankCounts[sym + rangeEncNumSyms]
          = gt_SRLSymbolCountInSeqRegion(
            seqIdx->rangeEncs, 0, posB,
            MRAEncRevMapSymbol(seqIdx->rangeMapAlphabet, sym),
            &hint->bcHint.rangeHint);
    }
    break;
  }
}

static Symbol
blockCompSeqGet(struct encIdxSeq *seq, GtUword pos, union EISHint *hint)
{
  Symbol sym;
  struct blockCompositionSeq *seqIdx;
  unsigned blockSize;
  gt_assert(seq && seq->classInfo == &blockCompositionSeqClass);
  if (pos >= seq->seqLen)
    return ~(Symbol)0;
  seqIdx = encIdxSeq2blockCompositionSeq(seq);
  blockSize = seqIdx->blockSize;
  {
    Symbol block[blockSize];
    blockCompSeqGetBlock(seqIdx, pos/blockSize, &(hint->bcHint), 1, NULL,
                         block);
    sym = block[pos%blockSize];
  }
  return sym;
}

static inline partialSymSum *
newPartialSymSums(AlphabetRangeSize alphabetSize)
{
  return gt_calloc(alphabetSize, sizeof (GtUword));
}

static inline void
deletePartialSymSums(partialSymSum *sums)
{
  gt_free(sums);
}

static inline void
addBlock2PartialSymSums(partialSymSum *sums, const Symbol *block,
                        unsigned blockSize)
{
  unsigned i;
  for (i = 0; i < blockSize; ++i)
    ++(sums[block[i]]);
}

static inline void
copyPartialSymSums(AlphabetRangeSize alphabetSize, partialSymSum *dest,
                   const partialSymSum *src)
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

static inline GtUword
numBuckets(GtUword seqLen, size_t bucketLen)
{
  /* seqLen + 1 because the partial sums for seqLen are used  */
  return (seqLen + 1) / bucketLen + (((seqLen + 1) % bucketLen)?1:0);
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
vwBitsSimple(GtUword seqLen, unsigned blockSize, unsigned bucketBlocks,
             unsigned maxPermIdxBits, BitOffset maxVarExtBitsPerBucket)
{
  return numBuckets(seqLen, bucketBlocks * blockSize)
    * (maxPermIdxBits * bucketBlocks + maxVarExtBitsPerBucket);
}

static inline BitOffset
vwBits(GtUword seqLen, unsigned blockSize, unsigned bucketBlocks,
       unsigned maxPermIdxBits, varExtBitsEstimator biVarBits, void *cbState,
       struct varBitsEstimate *extVarBitsUpperBound)
{
  size_t bucketLen = (size_t)bucketBlocks * blockSize;
  BitOffset maxVarBits = numBuckets(seqLen, bucketLen)
    * (maxPermIdxBits * bucketBlocks);
  if (biVarBits)
  {
    struct segmentDesc desc[2];
    struct varBitsEstimate *extVarBits, extVarBitsTemp;
    extVarBits = extVarBitsUpperBound?extVarBitsUpperBound:&extVarBitsTemp;
    desc[0].repeatCount = (seqLen + 1) / bucketLen;
    desc[0].len = bucketLen;
    desc[1].repeatCount = ((seqLen + 1) % bucketLen)?1:0;
    desc[1].len = seqLen % bucketLen;
    if (biVarBits(cbState, desc, sizeof (desc)/sizeof (desc[0]), extVarBits))
      maxVarBits += extVarBits->maxBitsTotal;
    else
      maxVarBits += numBuckets(seqLen, bucketLen)
        * extVarBits->maxBitsPerBucket;
  }
  else if (extVarBitsUpperBound)
  {
    extVarBitsUpperBound->maxBitsTotal
      = extVarBitsUpperBound->maxBitsPerBucket
      = extVarBitsUpperBound->maxBitsPerPos = 0;
  }
  return maxVarBits;
}

/**
 * @return 0 on error, 1 otherwise
 */
static int
openOnDiskData(const char *projectName, struct onDiskBlockCompIdx *idx,
               char *mode,GtError *err)
{
  GtStr *bdxName = gt_str_new_cstr(projectName);
  gt_str_append_cstr(bdxName, ".bdx");
  idx->idxFP = gt_fa_fopen(gt_str_get(bdxName), mode, err);
  idx->idxFN = gt_str_ref(bdxName);
  gt_str_delete(bdxName);
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
destructOnDiskBlockCompIdx(struct onDiskBlockCompIdx *idx)
{
  if (idx->idxMMap)
    gt_fa_xmunmap(idx->idxMMap);
  if (idx->idxFP)
    gt_fa_xfclose(idx->idxFP);
  if (idx->idxFN)
    gt_str_delete(idx->idxFN);
}

static inline void
initAppendState(struct appendState *aState,
                const struct blockCompositionSeq *seqIdx)
{
  BitOffset compCacheLen = superBlockCWBits(seqIdx) + bitElemBits - 1,
    permCacheLen = superBlockVarMaxBits(seqIdx) + bitElemBits - 1;
  aState->compCacheLen = compCacheLen;
  aState->permCacheLen = permCacheLen;
  aState->compCache = gt_calloc(sizeof (BitElem),
                                bitElemsAllocSize(compCacheLen));
  aState->permCache = gt_calloc(sizeof (BitElem),
                                bitElemsAllocSize(permCacheLen));
  aState->cwMemPos = cwPreCompIdxBits(seqIdx);
  aState->cwDiskOffset = aState->varMemPos = aState->cwMemOldBits =
    aState->varDiskOffset = aState->varMemOldBits = 0;
}

static void
destructAppendState(struct appendState *aState)
{
  gt_free(aState->permCache);
  gt_free(aState->compCache);
}

/**
 * @return 0 on successfull update, <0 on error
 */
static inline void
append2IdxOutput(struct appendState *state,
                 PermCompIndex permCompIdx[2],
                 unsigned bitsOfCompositionIdx, unsigned bitsOfPermutationIdx)
{
  gt_assert(state->cwMemPos + bitsOfCompositionIdx <= state->compCacheLen);
  gt_bsStorePermCompIndex(state->compCache, state->cwMemPos,
                          bitsOfCompositionIdx, permCompIdx[0]);
  state->cwMemPos += bitsOfCompositionIdx;
  gt_assert(state->varMemPos + bitsOfPermutationIdx <= state->permCacheLen);
  gt_bsStorePermCompIndex(state->permCache, state->varMemPos,
                          bitsOfPermutationIdx, permCompIdx[1]);
  state->varMemPos += bitsOfPermutationIdx;
}

static BitOffset
appendCallBackOutput(struct appendState *state,
                     const struct blockCompositionSeq *seqIdx,
                     bitInsertFunc biFunc, GtUword start,
                     GtUword len, unsigned callBackDataOffsetBits,
                     void *cbState)
{
  BitOffset bitsWritten;
  gt_assert(state);
  if (callBackDataOffsetBits)
  {
    BitOffset offset = cwPreCBOffsetBits(seqIdx);
    gt_bsStoreUInt64(state->compCache, state->cwMemOldBits + offset,
                  callBackDataOffsetBits,
                  state->varMemPos - state->varMemOldBits);
  }
  bitsWritten = biFunc(state->compCache, state->cwMemOldBits
                       + cwPreCWExtBits(seqIdx),
                       state->permCache, state->varMemPos, start, len, cbState);
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
                const partialSymSum *buck)
{
  size_t recordsExpected, cwBitElems;
  AlphabetRangeSize blockAlphabetSize;
  gt_assert(seqIdx && aState && buck);
  /* seek2/write constant width indices */
  if (!(seqIdx->externalData.cwDataPos + aState->cwDiskOffset
         < seqIdx->externalData.varDataPos
         + aState->varDiskOffset/bitElemBits * sizeof (BitElem)))
  {
    fprintf(stderr,"cwDatapos="GT_WU"\n",
           (GtUword) seqIdx->externalData.cwDataPos);
    fprintf(stderr,"cwDiskOffset="GT_WU"\n",(GtUword) aState->cwDiskOffset);
    fprintf(stderr,"varDataPos="GT_WU"\n",
                   (GtUword) seqIdx->externalData.varDataPos);
    fprintf(stderr,"bitElemBits="GT_WU"\n",(GtUword) bitElemBits);
    fprintf(stderr,"aState->varDiskOffset="GT_WU"\n",
                    (GtUword) aState->varDiskOffset);
    fprintf(stderr,"aState->varDiskOffset/bitElemBits="GT_WU"\n",
                    (GtUword) (aState->varDiskOffset/bitElemBits));
    fprintf(stderr,"sizeof (BitElem)="GT_WU"\n",(GtUword) sizeof (BitElem));
    fprintf(stderr,"aState->varDiskOffset/bitElemBits * sizeof (BitElem)="
            GT_WU"\n", (GtUword) (aState->varDiskOffset/bitElemBits * sizeof
                                  (BitElem)));
    fprintf(stderr,"seqIdx->externalData.varDataPos + "
            "aState->varDiskOffset/bitElemBits * sizeof (BitElem)=" GT_WU"\n",
            (GtUword) (seqIdx->externalData.varDataPos +
                       aState->varDiskOffset/bitElemBits * sizeof (BitElem)));
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_assert(seqIdx->externalData.cwDataPos + aState->cwDiskOffset
         < seqIdx->externalData.varDataPos
         + aState->varDiskOffset/bitElemBits * sizeof (BitElem));
  blockAlphabetSize = seqIdx->blockMapAlphabetSize;
  /* put bucket data for count up to the beginning of the current
   * block into the cw BitString */
  gt_assert(sizeof (GtUword) * CHAR_BIT >= seqIdx->bitsPerUlong);
  {
    Symbol i;
    for (i = 0; i < blockAlphabetSize; ++i)
    {
      gt_bsStoreUlong(aState->compCache,
                    aState->cwMemOldBits + seqIdx->partialSymSumBitsSums[i],
                    seqIdx->partialSymSumBits[i], buck[i]);
    }
  }
  /* append variable width offset position */
  gt_bsStoreUInt64(aState->compCache, aState->cwMemOldBits
                + cwPreVarIdxBits(seqIdx), seqIdx->bitsPerVarDiskOffset,
                aState->varDiskOffset);
  if (fseeko(seqIdx->externalData.idxFP,
            seqIdx->externalData.cwDataPos + aState->cwDiskOffset,
            SEEK_SET))
    return 0;
  cwBitElems = recordsExpected = aState->cwMemPos / bitElemBits;
  gt_xfwrite(aState->compCache, sizeof (BitElem), recordsExpected,
             seqIdx->externalData.idxFP);
  if ((aState->cwMemOldBits = aState->cwMemPos % bitElemBits))
    aState->compCache[0] = aState->compCache[recordsExpected];
  /* seek2/write variable width indices */
  if (fseeko(seqIdx->externalData.idxFP, seqIdx->externalData.varDataPos
            + aState->varDiskOffset/bitElemBits * sizeof (BitElem), SEEK_SET))
    return 0;
  recordsExpected = aState->varMemPos/bitElemBits;
  gt_xfwrite(aState->permCache, sizeof (BitElem), recordsExpected,
             seqIdx->externalData.idxFP);
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
  BKSZ_HEADER_FIELD = 0x424b535a, /* block size */
  BBLK_HEADER_FIELD = 0x42424c4b, /* blocks per bucket */
  VOFF_HEADER_FIELD = 0x564f4646, /* variable string offset */
  ROFF_HEADER_FIELD = 0x524f4646, /* range encoding offset */
  NMRN_HEADER_FIELD = 0x4e4d524e, /* number of ranges */
  CBMB_HEADER_FIELD = 0x43424d42, /* block internal offset for ext
                                   * bits provided by callback */
  MEXB_HEADER_FIELD = 0x4d455842, /* maxVarExtBitsPerBucket */
  CEXB_HEADER_FIELD = 0x43455842, /* cwExtBitsPerBucket */
  SPBT_HEADER_FIELD = 0x53504254, /* bits stored for GtUword values */
  SSBT_HEADER_FIELD = 0x53534254, /* block map alphabet size */
  BEFB_HEADER_FIELD = 0x42454642, /* block encoding fallback symbol */
  REFB_HEADER_FIELD = 0x52454642, /* range encoding fallback symbol */
  VDOB_HEADER_FIELD = 0x56444f42, /* bitsPerVarDiskOffset */
  SELE_HEADER_FIELD = 0x53454c45, /* sequence length */
  EH_HEADER_PREFIX = 0x45480000,  /* extension headers */
};

static const char bdxHeader[] = "BDX";

static inline off_t
extHeadersSizeAggregate(size_t numExtHeaders, const uint32_t *extHeaderSizes)
{
  off_t len = 0, i;
  for (i = 0; i < numExtHeaders; ++i)
    len += extHeaderSizes[i] + EXT_HEADER_PREFIX_SIZE;
  return len;
}

static inline size_t
blockEncIdxSeqHeaderLength(const struct blockCompositionSeq *seqIdx,
                           size_t numExtHeaders, const uint32_t *extHeaderSizes)
{
  size_t headerSize =
    4                           /* BDX identifier */
    + 4                         /* length field */
    + 8                         /* block size */
    + 8                         /* blocks per bucket */
    + 12                        /* offset of variable length data */
    + 12                        /* offset of range encodings */
    + 4 + 4                     /* bits used per seqpos */
    + 4 + 4                     /* bits used per variable bit offset */
    + 4 + 4 + 4 * seqIdx->blockMapAlphabetSize /* bit counts for partial sums */
    + 4 + 4                     /* block encoding fallback symbol */
    + 4 + 4                     /* range encoding fallback symbol */
    + 4 + 4                     /* num modes */
    + 4 + 8                     /* length of sequence */
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
                   off_t pos, uint32_t headerID)
{
  *headerList = gt_realloc(*headerList,
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
  gt_xfwrite(expHeader, sizeof (uint32_t), 2, fp);
  return cbFunc(fp, cbData);
}

#define writeIdxHeaderErrRet(retval)            \
  do {                                          \
    gt_free(buf);                               \
    return 0;                                   \
  } while (0)
/**
 * FIXME: this doesn't work on platforms with sizeof (uint32_t) != 4
 * or sizeof (int) < 4
 * @return 0 on error, header length in bytes on success
 */
static size_t
writeIdxHeader(struct blockCompositionSeq *seqIdx,
               size_t numExtHeaders, const uint16_t *headerIDs,
               const uint32_t *extHeaderSizes,
               headerWriteFunc *extHeaderCallbacks,
               void **headerCBData, GT_UNUSED GtError *err)
{
  FILE *fp;
  /* construct memory buffer with header data */
  size_t i, bufLen;
  off_t offset, len;
  char *buf;
  gt_assert(seqIdx && err);
  fp = seqIdx->externalData.idxFP;
  bufLen = blockEncIdxSeqHeaderLength(seqIdx, 0, NULL);
  buf = gt_malloc(bufLen);
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
  *(uint32_t *)(buf + offset) = SELE_HEADER_FIELD;
  *(uint64_t *)(buf + offset + 4) = seqIdx->baseClass.seqLen;
  offset += 12;
  *(uint32_t *)(buf + offset) = SPBT_HEADER_FIELD;
  *(uint32_t *)(buf + offset + 4) = seqIdx->bitsPerUlong;
  offset += 8;
  *(uint32_t *)(buf + offset) = VDOB_HEADER_FIELD;
  *(uint32_t *)(buf + offset + 4) = seqIdx->bitsPerVarDiskOffset;
  offset += 8;
  *(uint32_t *)(buf + offset) = SSBT_HEADER_FIELD;
  *(uint32_t *)(buf + offset + 4) = seqIdx->blockMapAlphabetSize;
  for (i = 0; i < seqIdx->blockMapAlphabetSize; ++i)
    *(uint32_t *)(buf + offset + 8 + 4*i) = seqIdx->partialSymSumBits[i];
  offset += 8 + 4 * seqIdx->blockMapAlphabetSize;
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
  gt_assert(offset == bufLen);
  if (fseeko(fp, 0, SEEK_SET))
    writeIdxHeaderErrRet(0);
  gt_xfwrite(buf, bufLen, 1, fp);
  if (numExtHeaders)
  {
    off_t offsetTargetValue = bufLen;
    if (!writeExtIdxHeader(fp, headerIDs[0], extHeaderSizes[0],
                          extHeaderCallbacks[0], headerCBData[0]))
      writeIdxHeaderErrRet(0);
    else
      appendExtHeaderPos(&seqIdx->extHeaderPos, seqIdx->numExtHeaders++,
                         bufLen + EXT_HEADER_PREFIX_SIZE,
                         EH_HEADER_PREFIX | headerIDs[0]);
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
                           EH_HEADER_PREFIX | headerIDs[0]);
    }
    offset = offsetTargetValue + EXT_HEADER_PREFIX_SIZE
      + extHeaderSizes[numExtHeaders - 1];
    gt_assert(offset >= ftello(fp));
  }
  gt_assert(seqIdx->externalData.cwDataPos == len);
  gt_free(buf);
  return len;
}

#define loadBlockEncIdxSeqErrRet()                                \
  do {                                                            \
    destructOnDiskBlockCompIdx(&newSeqIdx->externalData);         \
    if (newSeqIdx->compositionTable.bitsPerCount)                 \
      gt_destructCompositionList(&newSeqIdx->compositionTable);   \
    if (newSeqIdx->rangeEncs)                                     \
      gt_deleteSeqRangeList(newSeqIdx->rangeEncs);                \
    gt_free(newSeqIdx->extHeaderPos);                             \
    gt_free(newSeqIdx->partialSymSumBits);                        \
    gt_free(buf);                                                 \
    if (alphabet) gt_MRAEncDelete(alphabet);                      \
    gt_free(modesCopy);                                           \
    gt_free(blockMapAlphabet);                                    \
    gt_free(rangeMapAlphabet);                                    \
    gt_free(newSeqIdx);                                           \
    return NULL;                                                  \
  } while (0)

struct encIdxSeq *
gt_loadBlockEncIdxSeqGen(MRAEnc *alphabet,
                         const char *projectName, int features, GtError *err)
{
  struct blockCompositionSeq *newSeqIdx = NULL;
  Symbol blockMapAlphabetSize, GT_UNUSED totalAlphabetSize;
  MRAEnc *blockMapAlphabet = NULL, *rangeMapAlphabet = NULL;
  size_t headerLen;
  int *modesCopy = NULL;
  char *buf = NULL;
  gt_assert(projectName && err);
  newSeqIdx = gt_calloc(sizeof (struct blockCompositionSeq), 1);
  newSeqIdx->baseClass.alphabet = alphabet;
  newSeqIdx->baseClass.classInfo = &blockCompositionSeqClass;

  if (!openOnDiskData(projectName, &newSeqIdx->externalData, "rb",err)) {
    gt_error_set(err, "error reading data from disk");
    loadBlockEncIdxSeqErrRet();
  }
  {
    size_t offset = HEADER_ID_BLOCK_LEN;
    buf = gt_malloc(HEADER_ID_BLOCK_LEN);
    if (fread(buf, HEADER_ID_BLOCK_LEN, 1,
          newSeqIdx->externalData.idxFP) != 1) {
      gt_error_set(err, "error reading header id");
      loadBlockEncIdxSeqErrRet();
    }
    if (strcmp(buf, bdxHeader)!= 0) {
      gt_error_set(err, "header is not BDX found %s instead", buf);
      loadBlockEncIdxSeqErrRet();
    }
    newSeqIdx->externalData.cwDataPos = headerLen = *(uint32_t *)(buf + 4);
    buf = gt_realloc(buf, headerLen);
    if (fread(buf + HEADER_ID_BLOCK_LEN, headerLen - HEADER_ID_BLOCK_LEN,
             1, newSeqIdx->externalData.idxFP) != 1) {
      gt_error_set(err, "error reading header");
      loadBlockEncIdxSeqErrRet();
    }
    while (offset < headerLen)
    {
      uint32_t currentHeader;
      switch (currentHeader = *(uint32_t *)(buf + offset))
      {
      case SELE_HEADER_FIELD:
        newSeqIdx->baseClass.seqLen = *(uint64_t *)(buf + offset + 4);
        offset += 12;
        break;
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
          modesCopy = newSeqIdx->modes = gt_malloc(sizeof (int) * numModes);
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
        newSeqIdx->bitsPerUlong = *(uint32_t *)(buf + offset + 4);
        offset += 8;
        break;
      case VDOB_HEADER_FIELD:
        newSeqIdx->bitsPerVarDiskOffset = *(uint32_t *)(buf + offset + 4);
        offset += 8;
        break;
      case SSBT_HEADER_FIELD:
        {
          size_t i;
          Symbol blockMapAlphabetSize = *(uint32_t *)(buf + offset + 4);
          newSeqIdx->partialSymSumBits
            = gt_malloc(sizeof (newSeqIdx->partialSymSumBits[0])
                        * blockMapAlphabetSize * 2);
          newSeqIdx->partialSymSumBitsSums
            = newSeqIdx->partialSymSumBits + blockMapAlphabetSize;
          if (blockMapAlphabetSize)
          {
            newSeqIdx->partialSymSumBits[0]= *(uint32_t *)(buf + offset + 8);
            newSeqIdx->partialSymSumBitsSums[0] = 0;
#ifdef EIS_DEBUG
            gt_log_log("partialSymSumBits[0]=%u\n",
                    newSeqIdx->partialSymSumBits[0]);
#endif
            for (i = 1; i < blockMapAlphabetSize; ++i)
            {
              newSeqIdx->partialSymSumBits[i]
                = *(uint32_t *)(buf + offset + 8 + 4*i);
#ifdef EIS_DEBUG
              gt_log_log("partialSymSumBits[%"PRIuSymbol"]=%u\n",
                      (Symbol)i, newSeqIdx->partialSymSumBits[i]);
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
                             offset, currentHeader);
          offset += extHeaderLen;
        }
        else
        {
          gt_error_set(err, "Unknown header field: %4s", buf + offset);
          loadBlockEncIdxSeqErrRet();
        }
      }
    }
    gt_assert(newSeqIdx->modes
           && newSeqIdx->bucketBlocks
           && newSeqIdx->numModes
           && offset == headerLen);
    /* compute values not read from header */
    if (!newSeqIdx->bitsPerUlong)
      newSeqIdx->bitsPerUlong =
        requiredUlongBits(newSeqIdx->baseClass.seqLen - 1);
    else if (newSeqIdx->bitsPerUlong !=
            requiredUlongBits(newSeqIdx->baseClass.seqLen - 1)) {
      gt_error_set(err, "bits per ulong wrong in header found %u, expected %u",
                   newSeqIdx->bitsPerUlong,
                   requiredUlongBits(newSeqIdx->baseClass.seqLen - 1));
      loadBlockEncIdxSeqErrRet();
    }
  }
  {
    AlphabetRangeID range, numAlphabetRanges = newSeqIdx->numModes =
      MRAEncGetNumRanges(alphabet);
    totalAlphabetSize = gt_MRAEncGetSize(alphabet);
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
        gt_error_set(err, "Invalid encoding request.");
        loadBlockEncIdxSeqErrRet();
        break;
      }
    }
    newSeqIdx->blockMapAlphabet = blockMapAlphabet =
      gt_MRAEncSecondaryMapping(alphabet, BLOCK_COMPOSITION_INCLUDE, modesCopy,
                             newSeqIdx->blockEncFallback);
    newSeqIdx->rangeMapAlphabet = rangeMapAlphabet =
      gt_MRAEncSecondaryMapping(alphabet, REGIONS_LIST, modesCopy,
                             newSeqIdx->rangeEncFallback);
    newSeqIdx->blockMapAlphabetSize = blockMapAlphabetSize;
    gt_assert(gt_MRAEncGetSize(blockMapAlphabet) == blockMapAlphabetSize);
  }
  if (!newSeqIdx->partialSymSumBits && blockMapAlphabetSize)
  {
    newSeqIdx->partialSymSumBits
      = gt_malloc(sizeof (newSeqIdx->partialSymSumBits[0])
                  * blockMapAlphabetSize * 2);
    newSeqIdx->partialSymSumBitsSums
      = newSeqIdx->partialSymSumBits + blockMapAlphabetSize;
    symSumBitsDefaultSetup(newSeqIdx);
  }
  if (!gt_initCompositionList(&newSeqIdx->compositionTable,
                              newSeqIdx->blockSize,
                              blockMapAlphabetSize)) {
    gt_error_set(err, "error initialising composition list");
    loadBlockEncIdxSeqErrRet();
  }
  if (newSeqIdx->bitsPerVarDiskOffset == 0)
    newSeqIdx->bitsPerVarDiskOffset =
      gt_requiredUInt64Bits(
        vwBitsSimple(newSeqIdx->baseClass.seqLen, newSeqIdx->blockSize,
                     newSeqIdx->bucketBlocks,
                     newSeqIdx->compositionTable.maxPermIdxBits,
                     newSeqIdx->maxVarExtBitsPerBucket));

  if (fseeko(newSeqIdx->externalData.idxFP,
             newSeqIdx->externalData.rangeEncPos, SEEK_SET)) {
    gt_error_set(err, "error in fseeko");
    loadBlockEncIdxSeqErrRet();
  }
  {
    int regionFeatures = SRL_NO_FEATURES;
    if (features & EIS_FEATURE_REGION_SUMS)
      regionFeatures |= SRL_PARTIAL_SYMBOL_SUMS;
    if (!(newSeqIdx->rangeEncs =
          gt_SRLReadFromStream(newSeqIdx->externalData.idxFP, rangeMapAlphabet,
                            regionFeatures, err))) {
      gt_error_set(err, "error creating rangeEncs");
      loadBlockEncIdxSeqErrRet();
    }
  }
  tryMMapOfIndex(&newSeqIdx->externalData);
  gt_free(buf);
  return &newSeqIdx->baseClass;
}

static inline int
tryMMapOfIndex(struct onDiskBlockCompIdx *idxData)
{
  size_t len = idxData->rangeEncPos - idxData->cwDataPos;
  gt_assert(idxData && idxData->idxFP && idxData->idxMMap == NULL);
  idxData->idxMMap = gt_fa_mmap_generic_fd(fileno(idxData->idxFP),
                                           gt_str_get(idxData->idxFN),len,
                                           idxData->cwDataPos, false, false,
                                           NULL);
  return (idxData->idxMMap != NULL ? 0 : 1);
}

static FILE *
seekToHeader(const struct encIdxSeq *seqIdx, uint16_t headerID,
             uint32_t *lenRet)
{
  const struct blockCompositionSeq *bSeqIdx;
  size_t i;
  uint32_t query = headerID | EH_HEADER_PREFIX;
  gt_assert(seqIdx);
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
                  struct appendState *aState)
{
  off_t rangeEncPos;
  size_t recordsExpected;
  gt_assert(aState && seqIdx);
  gt_assert(aState->cwMemOldBits < bitElemBits
         && aState->cwMemOldBits + cwPreCompIdxBits(seqIdx)
         == aState->cwMemPos);
  gt_assert(aState->varMemOldBits < bitElemBits
         && aState->varMemOldBits == aState->varMemPos);
  if (aState->cwMemOldBits)
  {
    /* seek2/write variable width indices */
    if (fseeko(seqIdx->externalData.idxFP,
               seqIdx->externalData.cwDataPos
               + aState->cwDiskOffset * sizeof (BitElem), SEEK_SET))
      return 0;
    recordsExpected = 1;
    gt_xfwrite(aState->compCache, sizeof (BitElem), recordsExpected,
               seqIdx->externalData.idxFP);
    ++(aState->cwDiskOffset);
  }
  if (aState->varMemOldBits)
  {
    /* seek2/write variable width indices */
    if (fseeko(seqIdx->externalData.idxFP,
               seqIdx->externalData.varDataPos
               + aState->varDiskOffset/bitElemBits * sizeof (BitElem),
               SEEK_SET))
      return 0;
    recordsExpected = 1;
    gt_xfwrite(aState->permCache, sizeof (BitElem),
               recordsExpected,
               seqIdx->externalData.idxFP);
  }
  rangeEncPos = seqIdx->externalData.varDataPos
    + aState->varDiskOffset / bitElemBits * sizeof (BitElem)
    + ((aState->varDiskOffset%bitElemBits)?sizeof (BitElem):0);
  seqIdx->externalData.rangeEncPos = rangeEncPos;
  /* insert terminator so every search for a next range will find a
   * range just beyond the sequence end */
  gt_SRLAppendNewRange(seqIdx->rangeEncs,
                    seqIdx->baseClass.seqLen + seqIdx->blockSize, 1, 0);
  gt_SRLCompact(seqIdx->rangeEncs);
  if (fseeko(seqIdx->externalData.idxFP, rangeEncPos, SEEK_SET))
    return 0;
  if (!(gt_SRLSaveToStream(seqIdx->rangeEncs, seqIdx->externalData.idxFP)))
     return 0;
  return 1;
}

static void
addRangeEncodedSyms(struct seqRangeList *rangeList, const Symbol *block,
                    unsigned blockSize, GtUword blockNum,
                    const MRAEnc *alphabet, int selection, const int *rangeSel)
{
  unsigned i;
  for (i = 0; i < blockSize; ++i)
  {
    gt_assert(gt_MRAEncSymbolIsInSelectedRanges(alphabet, block[i],
                                          selection, rangeSel) >= 0);
    if (gt_MRAEncSymbolIsInSelectedRanges(alphabet, block[i], selection,
                                          rangeSel)) {
      gt_SRLAddPosition(rangeList, blockNum * blockSize + i, block[i]);
    }
  }
}

static union EISHint *
newBlockCompSeqHint(const struct encIdxSeq *seq)
{
  union EISHint *hintret;
  const struct blockCompositionSeq *seqIdx;
  gt_assert(seq && seq->classInfo == &blockCompositionSeqClass);
  seqIdx = constEncIdxSeq2blockCompositionSeq(seq);
  hintret = gt_malloc(sizeof (union EISHint));
  gt_SRLInitListSearchHint(seqIdx->rangeEncs, &hintret->bcHint.rangeHint);
  /* FIXME: make cache size user-configurable */
  initSuperBlockSeqCache(&hintret->bcHint.sBlockCache, seqIdx, 32);
  return hintret;
}

static void
deleteBlockCompSeqHint(struct encIdxSeq *seq, union EISHint *hint)
{
  GT_UNUSED struct blockCompositionSeq *seqIdx;
  gt_assert(seq && seq->classInfo == &blockCompositionSeqClass);
  seqIdx = encIdxSeq2blockCompositionSeq(seq);
  destructSuperBlockSeqCache(&hint->bcHint.sBlockCache);
  gt_free(hint);
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

enum {
  BUCKET_PLAIN                     = 0,
  BUCKET_PRINT_BITSTRING           = 1 <<  0,
  BUCKET_PRINT_BITSTRING_SEPARATOR = 1 <<  1,
  BUCKET_PRINT_BITSIZES            = 1 <<  2,
  BUCKET_PRINT_RANGES_OVERLAY      = 1 <<  3,
};

static int
printBucket(const struct blockCompositionSeq *seqIdx, GtUword bucketNum,
            int flags, FILE *fp, union EISHint *hint)
{
  GtUword lastBucket = bucketNumFromPos(seqIdx,
                                             EISLength(&seqIdx->baseClass)),
                start, end;
  AlphabetRangeSize i, blockMapAlphabetSize = seqIdx->blockMapAlphabetSize;
  int outCount = 0;
  gt_assert(seqIdx && fp && hint);
  if (bucketBasePos(seqIdx, bucketNum) >= EISLength(&seqIdx->baseClass))
  {
    gt_log_log("warning: querying bucket "GT_WU""
            " beyond end of sequence!\n", bucketNum);
    bucketNum = lastBucket;
  }
  start = bucketBasePos(seqIdx, bucketNum);
  end = ((bucketNum < lastBucket)?
         bucketBasePos(seqIdx, bucketNum + 1):
         EISLength(&seqIdx->baseClass));
  outCount +=
    fprintf(fp, "# Inspecting bucket: "GT_WU"\n"
            "# bucket position start="GT_WU", end="GT_WU"\n"
            "# partial symbol sums up to start:\n",
            bucketNum, start, end - 1);
  {
    struct superBlock *sBlock;
    BitOffset varDataMemOffset, cwIdxMemOffset, varIdxOffset;
    Symbol *block;
    unsigned blockSize = seqIdx->blockSize,
      idxBits = seqIdx->compositionTable.compositionIdxBits;
#ifdef USE_SBLOCK_CACHE
    sBlock = cacheFetchSuperBlock(seqIdx, bucketNum,
                                  &hint->bcHint.sBlockCache);
#else
    sBlock = fetchSuperBlock(seqIdx, bucketNum, NULL);
#endif
    block = gt_malloc(sizeof (block[0]) * seqIdx->blockSize);
    for (i = 0; i < blockMapAlphabetSize; ++i)
    {
      if (flags & BUCKET_PRINT_BITSIZES)
      {
        outCount += fprintf(fp, "# partial sum[%u] bits: %u\n", i,
                            seqIdx->partialSymSumBits[i]);
      }
      outCount +=
        fprintf(fp, "# partial sum[%u]="GT_WU"\n", i,
                sBlockGetPartialSymSum(sBlock, i, seqIdx));
    }
    if (flags & BUCKET_PRINT_BITSIZES)
    {
      fprintf(fp, "# bit size of sBlockVarIdxOffset: %u\n",
              seqIdx->bitsPerVarDiskOffset);
    }
    fprintf(fp, "# sBlockVarIdxOffset="GT_LLU"\n",
            (varIdxOffset = sBlockGetVarIdxOffset(sBlock, seqIdx)));
    if (seqIdx->callBackDataOffsetBits)
      fprintf(fp, "# sBlockGetcbOffset="GT_LLU"\n",
              sBlockGetcbOffset(sBlock, seqIdx));
    i = 0;
    if (flags & BUCKET_PRINT_BITSIZES)
    {
      fprintf(fp, "# bits of constant width extension data: "GT_LLU"\n",
              (GtUint64)seqIdx->cwExtBitsPerBucket);
    }
    walkCompIndicesPrefix(
      seqIdx, sBlock, seqIdx->bucketBlocks, cwIdxMemOffset,
      outCount +=
      fprintf(
         fp, "# block %u: comp idx: "GT_WU", permIdxBits=%u, perm idx: " GT_WU
         "=>", i, (GtUword)compIndex,
         (unsigned)seqIdx->compositionTable.permutations[compIndex].permIdxBits,
         (GtUword)gt_bsGetPermCompIndex(
                   sBlock->varData, varDataMemOffset,
                   seqIdx->compositionTable.permutations[compIndex].permIdxBits)
         );
      unpackBlock(seqIdx, sBlock, cwIdxMemOffset, varDataMemOffset, block,
                  blockSize);
      outCount += printBlock(block, blockSize, fp);
      fputs("\n", fp);
      i++; , varDataMemOffset);
    gt_free(block);
    if (flags & BUCKET_PRINT_BITSTRING)
    {
      for (i = 0; i < blockMapAlphabetSize; ++i)
      {
        if (gt_bsPrint(fp, sBlock->cwData, seqIdx->partialSymSumBitsSums[i]
                    + sBlock->cwIdxMemBase, seqIdx->partialSymSumBits[i]))
          outCount += seqIdx->partialSymSumBits[i];
        if (flags & BUCKET_PRINT_BITSTRING_SEPARATOR)
          outCount += fputs("&", fp);
      }
      gt_bsPrint(fp, sBlock->cwData,
              cwPreVarIdxBits(seqIdx) + sBlock->cwIdxMemBase,
              seqIdx->bitsPerVarDiskOffset);
      if (flags & BUCKET_PRINT_BITSTRING_SEPARATOR)
        outCount += fputs("&", fp);
      if (seqIdx->callBackDataOffsetBits)
      {
        gt_bsPrint(fp, sBlock->cwData, sBlockGetcbOffsetOffset(sBlock, seqIdx),
                seqIdx->callBackDataOffsetBits);
        if (flags & BUCKET_PRINT_BITSTRING_SEPARATOR)
          outCount += fputs("&", fp);
      }
      walkCompIndicesPrefix(
        seqIdx, sBlock, seqIdx->bucketBlocks, cwIdxMemOffset,
        gt_bsPrint(fp, sBlock->cwData, cwIdxMemOffset, idxBits);
        if (flags & BUCKET_PRINT_BITSTRING_SEPARATOR)
          outCount += fputs("&", fp);, varDataMemOffset);
      gt_bsPrint(fp, sBlock->cwData, sBlockCWExtBitsOffset(sBlock, seqIdx),
              seqIdx->cwExtBitsPerBucket);
      fputs("\n# variable width string: \n", fp);
      walkCompIndicesPrefix(
        seqIdx, sBlock, seqIdx->bucketBlocks, cwIdxMemOffset,
        idxBits =
        seqIdx->compositionTable.permutations[compIndex].permIdxBits;
        gt_bsPrint(fp, sBlock->varData, varDataMemOffset, idxBits);
        if (flags & BUCKET_PRINT_BITSTRING_SEPARATOR)
          outCount += fputs("&", fp);, varDataMemOffset);
      if (bucketNum != lastBucket)
      {
        struct superBlock *nextSBlock;
        BitOffset nextVarIdxOffset;
#ifdef USE_SBLOCK_CACHE
        nextSBlock = cacheFetchSuperBlock(seqIdx, bucketNum + 1,
                                          &hint->bcHint.sBlockCache);
#else
        nextSBlock = fetchSuperBlock(seqIdx, bucketNum + 1, NULL);
#endif
        nextVarIdxOffset = sBlockGetVarIdxOffset(nextSBlock, seqIdx);
        gt_bsPrint(fp, sBlock->varData, varDataMemOffset,
                nextVarIdxOffset - varIdxOffset);
        fprintf(fp, "\n# varIdxOffset for next block: "GT_LLU"",
                (GtUint64)nextVarIdxOffset);
#ifndef USE_SBLOCK_CACHE
        deleteSuperBlock(sBlock);
#endif
      }
      fputs("\n", fp);
    }
    if (flags & BUCKET_PRINT_RANGES_OVERLAY)
    {
      fputs("# overlapping symbol ranges:\n", fp);
      gt_SRLPrintRangesInfo(seqIdx->rangeEncs, fp, start, end - start,
                         &hint->bcHint.rangeHint);
    }
  }
  return outCount;
}

unsigned
gt_blockEncIdxSeqSegmentLen(const struct blockEncParams *params)
{
  return params->blockSize * params->bucketBlocks;
}

static int
printBlockEncPosDiags(const EISeq *seq, GtUword pos, FILE *fp,
                      EISHint hint)
{
  const struct blockCompositionSeq *seqIdx;
  int outCount = 0;
  GtUword bucketNum;
  gt_assert(seq && seq->classInfo == &blockCompositionSeqClass);
  seqIdx = constEncIdxSeq2blockCompositionSeq(seq);
  bucketNum = bucketNumFromPos(seqIdx, pos);
  fputs("##################################################\n"
        "# This bucket:\n"
        "##################################################\n", fp);
  outCount += printBucket(seqIdx, bucketNum, BUCKET_PLAIN, fp, hint);
  if (bucketNum)
  {
    fputs("##################################################\n"
          "# Previous bucket:\n"
          "##################################################\n", fp);
    outCount += printBucket(seqIdx, bucketNum - 1, BUCKET_PLAIN, fp, hint);
  }
  if (bucketNum < bucketNumFromPos(seqIdx, EISLength(seq)))
  {
    fputs("##################################################\n"
          "# Next bucket:\n"
          "##################################################\n", fp);
    outCount += printBucket(seqIdx, bucketNum + 1, BUCKET_PLAIN, fp, hint);
  }
  return outCount;
}

static int
displayBlockEncBlock(const EISeq *seq, GtUword pos, FILE *fp,
                     EISHint hint)
{
  const struct blockCompositionSeq *seqIdx;
  int outCount;
  GtUword bucketNum;
  gt_assert(seq && seq->classInfo == &blockCompositionSeqClass);
  seqIdx = constEncIdxSeq2blockCompositionSeq(seq);
  bucketNum = bucketNumFromPos(seqIdx, pos);
  outCount = printBucket(seqIdx, bucketNum, BUCKET_PRINT_BITSTRING_SEPARATOR
                         | BUCKET_PRINT_BITSTRING | BUCKET_PRINT_BITSIZES
                         | BUCKET_PRINT_RANGES_OVERLAY, fp, hint);
  return outCount;
}

static const struct encIdxSeqClass blockCompositionSeqClass =
{
  .delete = deleteBlockEncIdxSeq,
  .rank = blockCompSeqRank,
  .posPairRank = blockCompSeqPosPairRank,
  .rangeRank = blockCompSeqRangeRank,
  .posPairRangeRank = blockCompSeqPosPairRangeRank,
  .select = blockCompSeqSelect,
  .get = blockCompSeqGet,
  .newHint = newBlockCompSeqHint,
  .deleteHint = deleteBlockCompSeqHint,
  .expose = blockCompSeqExpose,
  .seekToHeader = seekToHeader,
  .printPosDiags = printBlockEncPosDiags,
  .printExtPosDiags = displayBlockEncBlock,
};
