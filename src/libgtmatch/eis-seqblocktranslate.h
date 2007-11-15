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

#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#include "libgtmatch/eis-mrangealphabet.h"

typedef uint32_t PermCompIndex;
#define bsGetPermCompIndex bsGetUInt32
#define bsStorePermCompIndex bsStoreUInt32

struct permList
{
  size_t numPermutations;
  unsigned permIdxBits;
  BitOffset catPermsOffset;
};

struct compList
{
  size_t numCompositions;
  struct permList *permutations;
  BitString catCompsPerms;             /**< bitstring with all
                                        *   compositions, followed by
                                        *   all permutations encoded
                                        *   in minimal space
                                        *   concatenated together */
  unsigned bitsPerCount, bitsPerSymbol, compositionIdxBits, maxPermIdxBits;
  size_t maxPermCount;
};

/**
 * @return 0 on error, !0 otherwise
 */
extern int
initCompositionList(struct compList *compList, unsigned blockSize,
                    unsigned numSyms, Env *env);
extern void
destructCompositionList(struct compList *clist, Env *env);

/**
 * @return NULL on error
 */
extern struct compList *
newCompositionList(unsigned blockSize, unsigned alphabetSize, Env *env);

extern void
deleteCompositionList(struct compList *clist, Env *env);

static inline BitOffset
compListPermStartOffset(struct compList *list, unsigned numSyms)
{
  return list->numCompositions * list->bitsPerCount * numSyms;
}

/**
 * \brief Transforms a block-sized sequence of symbols to corresponding
 * index-pair representation of block-composition index.
 * @param tables to use for lookup
 * @param blockSize length of encoded sequence blocks
 * @param alphabetSize encoded symbols are from range 0..(alphabetSize - 1)
 * @param block symbol sequence holding at least as many symbols as
 * required by seqIdx
 * @param idxOutput on return idxOutput[0] holds the composition
 * index, idxOutput[1] the permutation index.
 * @param bitsOfPermIdx if non-NULL, the number of significant bits
 * for the permutation index is stored here.
 * @param env Environment to use for memory allocation etc.
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
                unsigned *bitsOfPermIdx, Env *env,
                BitString permCompPA, unsigned *compPA);

static inline void
indexPair2block(const struct compList *compositionTable, unsigned blockSize,
                PermCompIndex compIdx, PermCompIndex permIdx,
                Symbol *block, unsigned sublen)
{
  unsigned bitsPerPermutation;
  struct permList *permutationList;
  assert(compositionTable && block);
  bitsPerPermutation = compositionTable->bitsPerSymbol * blockSize;
  permutationList = compositionTable->permutations + compIdx;
  bsGetUniformSymbolArray(
    compositionTable->catCompsPerms,
    permutationList->catPermsOffset + bitsPerPermutation * permIdx,
    compositionTable->bitsPerSymbol, sublen, block);
}

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

#endif
