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
#ifndef EIS_SEQBLOCKTRANSLATE_H
#define EIS_SEQBLOCKTRANSLATE_H

/**
 * \file eis-seqblocktranslate.h
 * Translate fixed length substrings (q-words) to a pair of indices: one to
 * designate the composition counts and one to define the actual
 * sequence
 */
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#include "libgtmatch/seqpos-def.h"
#include "libgtmatch/eis-bitpackseqpos.h"
#include "libgtmatch/eis-mrangealphabet.h"

/** store indices in this scalar */
typedef uint32_t PermCompIndex;
/** method to retrieve one index from a BitString */
#define bsGetPermCompIndex bsGetUInt32
/** method to store one index in a BitString */
#define bsStorePermCompIndex bsStoreUInt32

/**
 * For permutations common to one composition.
 */
struct permList
{
  size_t numPermutations;       /**< number of permutations this
                                 * composition produces  */
  unsigned permIdxBits;         /**< bits required to store values
                                 * from 0..numPermuations-1  */
  BitOffset catPermsOffset;     /**< offset into the table of all
                                 * permutations at which the
                                 * permutations corresponding to this
                                 * composition are stored */
};

/**
 * Holds the lists of compositions and permutations.
 */
struct compList
{
  size_t numCompositions;       /**< how many ways are there to choose
                                 *   q symbols from the given alphabet */
  struct permList *permutations;/**< reference of list of permutations */
  BitString catCompsPerms;      /**< bitstring with all compositions,
                                 *   followed by all permutations
                                 *   encoded in minimal space
                                 *   concatenated together */
  unsigned bitsPerCount,        /**< bits required to hold one value 0..q */
    bitsPerSymbol,              /**< bits for each symbol */
    compositionIdxBits,         /**< log_2 of numcompositions */
    maxPermIdxBits;             /**< maximum bit length of permutation
                                 *   indices */
};

/**
 * @brief construct object to translate q-words to
 * composition/permutation index tuple
 * @param compList
 * @param blockSize length of q-words to encode
 * @param numSyms size of alphabet
 * @return 0 on error, !0 otherwise
 */
extern int
initCompositionList(struct compList *compList, unsigned blockSize,
                    unsigned numSyms);
/**
 * Deallocate resources of composition list object, storage struct is not freed.
 * @param clist
 */
extern void
destructCompositionList(struct compList *clist);

/**
 * @brief create new object to map q-words to composition/permutation
 * index tuples
 * @param blockSize length of q-words to map
 * @param alphabetSize
 * @return NULL on error
 */
extern struct compList *
newCompositionList(unsigned blockSize, unsigned alphabetSize);

/**
 * Delete composition list object.
 * @param clist
 */
extern void
deleteCompositionList(struct compList *clist);

/**
 * \brief Transforms a block-sized sequence of symbols to corresponding
 * index-pair representation of block-composition index.
 * @param compositionTable tables to use for lookup
 * @param blockSize length of encoded sequence blocks
 * @param alphabetSize encoded symbols are from range 0..(alphabetSize - 1)
 * @param block symbol sequence holding at least as many symbols as
 * required by seqIdx
 * @param idxOutput on return idxOutput[0] holds the composition
 * index, idxOutput[1] the permutation index.
 * @param bitsOfPermIdx if non-NULL, the number of significant bits
 * for the permutation index is stored here.
 * @param permCompPA if not NULL must point to a memory region of
 * sufficient size to hold the concatenated bistring representations
 * of composition and permutation, composition at offset 0,
 * permutation at offset compositionTable->bitsPerCount * alphabetSize.
 * @param compPA if not NULL must reference a memory region of
 * sufficient size to hold the composition representation as alphabet
 * range sized sequence.
 * @return -1 in case of memory exhaustion, cannot happen if both
 * permCompPA and compPA are valid preallocated memory regions, 0
 * otherwise.
 */
extern int
block2IndexPair(const struct compList *compositionTable,
                unsigned blockSize, unsigned alphabetSize,
                const Symbol *block, PermCompIndex idxOutput[2],
                unsigned *bitsOfPermIdx,
                BitString permCompPA, unsigned *compPA);

/**
 * @brief Give q-word corresponding to pair of indices.
 * @param compositionTable
 * @param blockSize q-word length (must be same as used on construction)
 * @param compIdx composition index
 * @param permIdx permutation index
 * @param block write sequence of symbols here
 * @param subLen only unpack this many symbols
 */
static inline void
indexPair2block(const struct compList *compositionTable, unsigned blockSize,
                PermCompIndex compIdx, PermCompIndex permIdx,
                Symbol *block, unsigned subLen)
{
  unsigned bitsPerPermutation;
  struct permList *permutationList;
  assert(compositionTable && block);
  assert(subLen <= blockSize);
  bitsPerPermutation = compositionTable->bitsPerSymbol * blockSize;
  permutationList = compositionTable->permutations + compIdx;
  bsGetUniformSymbolArray(
    compositionTable->catCompsPerms,
    permutationList->catPermsOffset + bitsPerPermutation * permIdx,
    compositionTable->bitsPerSymbol, subLen, block);
}

/**
 * @brief Find how often a symbol occurs in a composition given by index
 * @param compositionTable
 * @param alphabetSize
 * @param compIndex composition index
 * @param sym
 * @return number of symbol occurrences
 */
static inline unsigned
symCountFromComposition(struct compList *compositionTable,
                        unsigned alphabetSize,
                        PermCompIndex compIndex, Symbol sym)
{
  BitOffset bitsPerComp, bitsPerCount;
  assert(compositionTable);
  bitsPerCount = compositionTable->bitsPerCount;
  bitsPerComp = bitsPerCount * alphabetSize;
  assert(compIndex < compositionTable->numCompositions);
  return bsGetUInt(compositionTable->catCompsPerms,
                   compIndex * bitsPerComp + sym * bitsPerCount,
                   bitsPerCount);
}

/**
 * @brief Add symbol occurrences for a composition given by index
 * @param compositionTable
 * @param alphabetSize
 * @param compIndex composition index
 * @param counts add occurrence counts to this array
 */
static inline void
addSymCountsFromComposition(struct compList *compositionTable,
                            unsigned alphabetSize,
                            PermCompIndex compIndex, Seqpos *counts)
{
  BitOffset bitsPerComp, bitsPerCount;
  assert(compositionTable);
  bitsPerCount = compositionTable->bitsPerCount;
  bitsPerComp = bitsPerCount * alphabetSize;
  assert(compIndex < compositionTable->numCompositions);
  bsGetUniformSeqposArrayAdd(compositionTable->catCompsPerms,
                             compIndex * bitsPerComp, bitsPerCount,
                             alphabetSize, counts);
}

#endif
