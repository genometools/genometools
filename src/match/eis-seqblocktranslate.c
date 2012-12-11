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

#ifndef NDEBUG
#include <math.h>
#endif

#include "core/bitpackstring.h"
#include "core/combinatorics.h"
#include "core/mathsupport.h"
#include "core/dynalloc.h"
#include "core/log.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "match/eis-seqblocktranslate.h"

#if 1
/**
 * Setup composition array for lexically maximal composition.
 * @param maxSym last symbol to represent (lexical maximum)
 * @param blockSize sum of symbol counts in composition
 * (i.e. composition contains this many symbols)
 * @param composition vector of symbol counts
 * @param symRMNZ points to index of rightmost non-zero symbol count
 * in composition vector (state of generating procedure)
 */
static inline void
initComposition(Symbol maxSym, unsigned blockSize,
                unsigned composition[maxSym + 1], unsigned *symRMNZ)
{
  memset(composition, 0, sizeof (composition[0]) * (maxSym));
  composition[*symRMNZ = maxSym] = blockSize;
}
#else
/**
 * Setup composition array for lexically minimal composition.
 * @param maxSym last symbol to represent (lexical maximum)
 * @param blockSize sum of symbol counts in composition
 * (i.e. composition contains this many symbols)
 * @param composition vector of symbol counts
 * @param symRMNZ points to index of rightmost non-zero symbol count
 * in composition vector (state of generating procedure)
 */
static inline void
initComposition(Symbol maxSym, unsigned blockSize,
                unsigned composition[maxSym + 1], unsigned *symLMNZ)
{
  memset(composition + 1, 0, sizeof (composition[0]) * (maxSym));
  composition[*symLMNZ = 0] = blockSize;
}
#endif

#if 1
/**
 * Compute lexically previous composition given any composition in
 * \link composition \endlink.
 * @param composition array[0..maxSym] of occurrence counts for each symbol
 * @param maxSym last symbol index of composition
 * @param symMNZ rightmost non-zero entry of composition[]
 */
static inline void
nextComposition(unsigned composition[],
                Symbol maxSym,
                unsigned *symRMNZ)
{
  ++composition[*symRMNZ - 1];
  if (!(--composition[*symRMNZ]))
  {
    --*symRMNZ;
  }
  else if (!(composition[maxSym]))
  {
    composition[maxSym] = composition[*symRMNZ];
    composition[*symRMNZ] = 0;
    *symRMNZ = maxSym;
  }
}
#else
/**
 * Compute lexically next composition given any composition in
 * \link composition \endlink.
 * @param composition array[0..maxSym] of occurrence counts for each symbol
 * @param maxSym last symbol index of composition
 * @param symMNZ rightmost non-zero entry of composition[]
 */
static inline void
nextComposition(unsigned composition[],
                Symbol maxSym,
                unsigned *symLMNZ)
{
  ++composition[*symLMNZ + 1];
  if (!(--composition[*symLMNZ]))
  {
    ++*symLMNZ;
  }
  else if (!(composition[0]))
  {
    composition[0] = composition[*symLMNZ];
    composition[*symLMNZ] = 0;
    *symLMNZ = 0;
  }
}
#endif

static int
initPermutationsList(const unsigned *composition, struct permList *permutation,
                     BitString permStore, BitOffset permOffset,
                     unsigned blockSize, unsigned alphabetSize,
                     unsigned bitsPerSymbol);

static void
destructPermutationsList(struct permList *permutation);

#ifdef EIS_DEBUG
static void
printComposition(FILE *fp, const unsigned *composition,
                 unsigned alphabetSize, unsigned blockSize);
#endif /* EIS_DEBUG */

static inline BitOffset
compListPermStartOffset(struct compList *list, unsigned numSyms)
{
  return list->numCompositions * list->bitsPerCount * numSyms;
}

#define gt_initCompositionListErrRet()                                  \
  do {                                                                  \
    if (newList->permutations)                                          \
    {                                                                   \
      unsigned i;                                                       \
      for (i = 0; i < cmpIdx; ++i)                                      \
        destructPermutationsList(newList->permutations + i );           \
      gt_free(newList->permutations);                                   \
    }                                                                   \
    if (newList->catCompsPerms)                                         \
      gt_free(newList->catCompsPerms);                                  \
    if (composition) gt_free(composition);                              \
    return 0;                                                           \
  } while (0)

int
gt_initCompositionList(struct compList *newList, unsigned blockSize,
                    unsigned alphabetSize)
{
  unsigned *composition = NULL;
  Symbol maxSym = alphabetSize - 1;
  BitOffset bitsPerComp, bitsPerCount, bitsPerPerm;
  size_t numCompositions, cmpIdx = 0;
  size_t maxNumPermutations = 0, numTotalPermutations;
  gt_assert(newList);
  newList->permutations = NULL;
  newList->catCompsPerms = NULL;
  if (!(composition = gt_malloc(sizeof (composition[0]) * alphabetSize)))
    gt_initCompositionListErrRet();
  bitsPerComp = alphabetSize * (bitsPerCount = gt_requiredUIntBits(blockSize));
  newList->bitsPerCount = bitsPerCount;
  numCompositions = newList->numCompositions =
    gt_combinatorics_binomial_simple(blockSize + maxSym, maxSym);
  numTotalPermutations = gt_combinatorics_i_pow(alphabetSize, blockSize);
  newList->compositionIdxBits = gt_requiredUInt64Bits(numCompositions - 1);
  newList->bitsPerSymbol = gt_requiredUIntBits(maxSym);
  bitsPerPerm = newList->bitsPerSymbol * blockSize;
  {
    size_t size = bitElemsAllocSize(numCompositions * bitsPerComp
                                    + numTotalPermutations * bitsPerPerm)
      * sizeof (BitElem);
    if (size == SIZE_MAX)
      gt_initCompositionListErrRet();
    newList->catCompsPerms = gt_malloc(size);
  }
  newList->permutations = gt_calloc(sizeof (struct permList), numCompositions);
  {
    unsigned symRMNZ;
    BitOffset offset = 0, permOffset =
      compListPermStartOffset(newList, alphabetSize);
#ifndef NDEBUG
    size_t permSum = 0;
#endif
    initComposition(maxSym, blockSize, composition, &symRMNZ);
    do
    {
#ifdef EIS_DEBUG
      printComposition(stderr, composition, alphabetSize, blockSize);
#endif /* EIS_DEBUG */
      gt_bsStoreUniformUIntArray(newList->catCompsPerms, offset,
                              bitsPerCount, alphabetSize, composition);
      gt_assert(cmpIdx > 1?(gt_bsCompare(newList->catCompsPerms, offset,
                                      bitsPerComp,
                              newList->catCompsPerms, offset - bitsPerComp,
                              bitsPerComp)>0):1);
      if (initPermutationsList(composition, newList->permutations + cmpIdx,
                               newList->catCompsPerms, permOffset, blockSize,
                               alphabetSize, newList->bitsPerSymbol))
        gt_initCompositionListErrRet();
#if defined(EIS_DEBUG) && !defined(NDEBUG)
      gt_log_log("%lu",
              (unsigned long)newList->permutations[cmpIdx].numPermutations);
#endif
#ifndef NDEBUG
      permSum += newList->permutations[cmpIdx].numPermutations;
#endif
      permOffset += newList->permutations[cmpIdx].numPermutations * bitsPerPerm;
      if (newList->permutations[cmpIdx].numPermutations > maxNumPermutations)
        maxNumPermutations = newList->permutations[cmpIdx].numPermutations;
      if (++cmpIdx < numCompositions)
        ;
      else
        break;
      offset += bitsPerComp;
      nextComposition(composition, maxSym, &symRMNZ);
    } while (1);
    /* verify that the last composition is indeed the lexically maximally */
    gt_assert(composition[0] == blockSize);
#if defined(EIS_DEBUG) && !defined(NDEBUG)
    gt_log_log("permSum=%lu, alphabetSize=%lu, blockSize=%d, "
            "pow(alphabetSize, blockSize)=%f",
            (unsigned long)permSum, (unsigned long)alphabetSize, blockSize,
            gt_power_for_small_exponents(alphabetSize, blockSize));
#endif
    gt_assert(permSum == gt_power_for_small_exponents(alphabetSize, blockSize));
  }
  newList->maxPermIdxBits = gt_requiredUInt64Bits(maxNumPermutations - 1);
  gt_free(composition);
  return 1;
}

void
gt_destructCompositionList(struct compList *clist)
{
  {
    unsigned i;
    for (i = 0; i < clist->numCompositions; ++i)
      destructPermutationsList(clist->permutations + i);
  }
  gt_free(clist->permutations);
  gt_free(clist->catCompsPerms);
}

struct compList *
gt_newCompositionList(unsigned blockSize, unsigned alphabetSize)
{
  struct compList *newList = NULL;
  if (!(newList = gt_calloc(1, sizeof (struct compList))))
    return NULL;
  if (!gt_initCompositionList(newList, blockSize, alphabetSize))
  {
    gt_free(newList);
    return NULL;
  }
  return newList;
}

void
gt_deleteCompositionList(struct compList *clist)
{
  gt_destructCompositionList(clist);
  gt_free(clist);
}

#ifdef EIS_DEBUG
/**
 * Compute number of digits a value would require when displayed in
 * selected number system (i.e. 2=binary, 10=decimal).
 * @param d value to be displayed
 * @param base number of symbols in output alphabet
 * @return number of digits required for output
 */
static inline int
digitPlaces(long d, int base)
{
  int l=1;
  while (d/=base)
    ++l;
  return l;
}

static void
printComposition(FILE *fp, const unsigned *composition,
                 unsigned alphabetSize, unsigned blockSize)
{
  Symbol sym;
  unsigned width = digitPlaces(blockSize, 10);
  for (sym = 0; sym < alphabetSize; ++sym)
  {
    fprintf(fp, "%*d ", width, composition[sym]);
  }
  fputs("\n", fp);
}
#endif /* EIS_DEBUG */

static void
initPermutation(Symbol *permutation, const unsigned *composition,
                unsigned alphabetSize);
static inline void
nextPermutation(Symbol *permutation, unsigned blockSize);

#ifdef EIS_DEBUG
static void
printPermutation(FILE *fp, Symbol *permutation, unsigned blockSize);
#endif /* EIS_DEBUG */

static int
initPermutationsList(const unsigned *composition, struct permList *permutation,
                     BitString permStore, BitOffset permOffset,
                     unsigned blockSize, unsigned alphabetSize,
                     unsigned bitsPerSymbol)
{
  size_t numPermutations = permutation->numPermutations =
    gt_combinatorics_multinomial(blockSize, alphabetSize, composition);
  if (numPermutations > 1)
    permutation->permIdxBits = gt_requiredUInt64Bits(numPermutations - 1);
  else
    permutation->permIdxBits = 0;
  Symbol *currentPermutation;
  BitOffset bitsPerPermutation = bitsPerSymbol * blockSize;
  if (!(currentPermutation = gt_malloc(sizeof (Symbol) * blockSize)))
    return -1;
  permutation->catPermsOffset = permOffset;
  initPermutation(currentPermutation, composition, alphabetSize);
  {
    size_t i = 0;
    BitOffset offset = permOffset;
    do
    {
      gt_bsStoreUniformSymbolArray(permStore, offset,
                                bitsPerSymbol, blockSize, currentPermutation);
#ifdef EIS_DEBUG
      printPermutation(stderr, currentPermutation, blockSize);
#endif /* EIS_DEBUG */
      gt_assert(i > 0?(gt_bsCompare(permStore, offset, bitsPerPermutation,
                              permStore,
                              offset - bitsPerPermutation,
                              bitsPerPermutation)>0):1);
      if (++i < numPermutations)
        ;
      else
        break;
      offset += bitsPerPermutation;
      nextPermutation(currentPermutation, blockSize);
    } while (1);
  }
  gt_free(currentPermutation);
  return 0;
}

static void
destructPermutationsList(GT_UNUSED struct permList *permutation)
{
/*   gt_free(permutation->catPerms); */
}

static void
initPermutation(Symbol *permutation, const unsigned *composition,
                unsigned alphabetSize)
{
  Symbol sym, *p = permutation;
  unsigned j;
  for (sym = 0; sym < alphabetSize; ++sym)
    for (j = 0; j < composition[sym]; ++j)
      *(p++) = sym;
}

static void
nextPermutation(Symbol *permutation, unsigned blockSize)
{
  /*
   * Every permutation represents a leaf in the decision tree
   * generating all permutations, thus given one permutation, we
   * ascend in the tree until we can make a decision resulting in a
   * lexically larger value.
   */
  unsigned upLvl = blockSize - 1;
  {
    Symbol maxSymInSuf = permutation[blockSize - 1];
    do
    {
      if (upLvl == 0 || maxSymInSuf > permutation[--upLvl])
        break;
      if (permutation[upLvl] > maxSymInSuf)
        maxSymInSuf = permutation[upLvl];
    } while (1);
  }
  /* now that we have found the branch point, descend again,
   * selecting first the next largest symbol */
  {
    Symbol saved = permutation[upLvl];
    Symbol swap = permutation[upLvl + 1];
    unsigned swapIdx = upLvl + 1;
    unsigned i;
    for (i = swapIdx + 1; i < blockSize; ++i)
    {
      if (permutation[i] > saved && permutation[i] < swap)
        swap = permutation[i], swapIdx = i;
    }
    permutation[upLvl] = swap;
    permutation[swapIdx] = saved;
  }
  for (++upLvl;upLvl < blockSize - 1; ++upLvl)
  {
    /* find minimal remaining symbol */
    unsigned i, minIdx = upLvl;
    Symbol minSym = permutation[upLvl];
    for (i = upLvl + 1; i < blockSize; ++i)
      if (permutation[i] < minSym)
        minSym = permutation[i], minIdx = i;
    permutation[minIdx] = permutation[upLvl];
    permutation[upLvl] = minSym;
  }
}

#ifdef EIS_DEBUG
static void
printPermutation(FILE *fp, Symbol *permutation, unsigned blockSize)
{
  unsigned i;
  if (blockSize)
    fprintf(fp, "%lu", (unsigned long)permutation[0]);
  for (i = 1; i < blockSize; ++i)
  {
    fprintf(fp, "%lu", (unsigned long)permutation[i]);
  }
  fputs("\n", fp);
}
#endif /* EIS_DEBUG */

int
gt_block2IndexPair(const struct compList *compositionTable,
                unsigned blockSize, unsigned alphabetSize,
                const Symbol *block, PermCompIndex idxOutput[2],
                unsigned *bitsOfPermIdx,
                BitString permCompPA, unsigned *compPA)
{
  unsigned bitsPerCount;
  BitOffset bitsPerComposition, bitsPerPermutation;
  BitString permCompBitString;
  gt_assert(compositionTable && idxOutput && block);
  gt_assert(blockSize > 0);
  bitsPerComposition = (bitsPerCount = compositionTable->bitsPerCount)
    * alphabetSize;
  bitsPerPermutation = compositionTable->bitsPerSymbol * blockSize;
  if (permCompPA)
    permCompBitString = permCompPA;
  else
    permCompBitString =
      gt_malloc(bitElemsAllocSize(bitsPerComposition + bitsPerPermutation)
                * sizeof (BitElem));
  {
    /* first compute composition from block */
    size_t i;
    unsigned *composition;
    if (compPA)
    {
      composition = compPA;
      memset(composition, 0, sizeof (composition[0])*alphabetSize);
    }
    else
      composition = gt_calloc(sizeof (composition[0]), alphabetSize);
    for (i = 0; i < blockSize; ++i)
    {
      ++composition[block[i]];
    }
    gt_bsStoreUniformUIntArray(permCompBitString, 0, bitsPerCount,
                            alphabetSize, composition);
    if (!compPA)
      gt_free(composition);
  }
  {
    /* do binary search for composition (simplified because the list
     * contains every composition possible and thus cmpresult will
     * become 0 at some point) */
    size_t compIndex = compositionTable->numCompositions / 2;
    size_t divStep = compIndex;
    int cmpresult;
    while ((cmpresult = gt_bsCompare(permCompBitString, 0, bitsPerComposition,
                                 compositionTable->catCompsPerms,
                                 compIndex * bitsPerComposition,
                                 bitsPerComposition)))
    {
      if (divStep > 1)
        divStep >>= 1; /* divStep /= 2 */
      if (cmpresult > 0)
        compIndex += divStep;
      else /* cmpresult < 0 */
        compIndex -= divStep;
    }
    idxOutput[0] = compIndex;
  }
  {
    const struct permList *permutation = compositionTable->permutations
      + idxOutput[0];
    if (bitsOfPermIdx)
      *bitsOfPermIdx = permutation->permIdxBits;
    if (permutation->numPermutations > 1)
    {
      /* build permutation bitstring */
      gt_bsStoreUniformSymbolArray(permCompBitString, bitsPerComposition,
                                compositionTable->bitsPerSymbol, blockSize,
                                block);
      /* do binary search for permutation */
      {
        size_t permIndex = permutation->numPermutations / 2;
        size_t divStep = permIndex;
        int cmpresult;
        while ((cmpresult = gt_bsCompare(permCompBitString, bitsPerComposition,
                                     bitsPerPermutation,
                                     compositionTable->catCompsPerms,
                                     permutation->catPermsOffset
                                     + permIndex * bitsPerPermutation,
                                     bitsPerPermutation)))
        {
          if (divStep > 1)
            divStep >>= 1; /* divStep /= 2 */
          if (cmpresult > 0)
            permIndex += divStep;
          else /* cmpresult < 0 */
            permIndex -= divStep;
        }
        idxOutput[1] = permIndex;
      }
    }
    else
        idxOutput[1] = 0;
  }
  if (!permCompPA)
    gt_free(permCompBitString);
  return 0;
}
