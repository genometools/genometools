/*
** Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
**
** See LICENSE file or http://genometools.org/license.html for license details.
**
*/
#ifndef BITPACKARRAY_H_INCLUDED
#define BITPACKARRAY_H_INCLUDED
/**
 * \file bitpackarray.h
 * \brief The class presented in this file encapsulates a bitstring in
 * a data structure to store/retrieve integers at fixed bitlength and
 * integrate well with genome tools library.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */

#include <assert.h>
#include <stdlib.h>

#include "bitpackstring.h"
#include "libgtcore/env.h"

struct BitPackArray
{
  BitString store;
  BitOffset numElems;
  unsigned bitsPerElem;
};

typedef struct BitPackArray BitPackArray;

/**
 * Create new BitPackArray structure.
 * @param bits number of bits to encode each value stored with
 * @param numValues number of values to store
 * @return pointer to new BitPackArray structure or NULL on failure
 */
static inline BitPackArray *
newBitPackArray(unsigned bits, BitOffset numValues, Env *env)
{
  BitPackArray *newBPA = env_ma_malloc(env, sizeof(*newBPA));
  if(newBPA)
  {
    if(!(newBPA->store = env_ma_malloc(env, bitElemsAllocSize(bits*numValues)
                                       * sizeof(BitElem))))
    {
      env_ma_free(newBPA, env);
      return NULL;
    }
    newBPA->bitsPerElem = bits;
    newBPA->numElems = numValues;
  }
  return newBPA;
}

static inline void
deleteBitPackArray(BitPackArray *bpa, Env *env)
{
  env_ma_free(bpa->store, env);
  env_ma_free(bpa, env);
}

/**
 * Stores unsigned integer in BitPackArray at given index.
 * @param bparray array to use
 * @param index index to store value at.
 * @param val value to store (only as many bits as specified on
 * BitPackArray construction are stored).
 */
static inline void
bpaStoreUInt32(BitPackArray *array, BitOffset index, uint32_t val)
{
  assert(array && index < array->numElems
         && array->bitsPerElem <= sizeof(val)*CHAR_BIT);
  bsStoreUInt32(array->store, array->bitsPerElem * index,
                array->bitsPerElem, val);
}

/**
 * Stores unsigned integer in BitPackArray at given index.
 * @param bparray array to use
 * @param index index to store value at.
 * @param val value to store (only as many bits as specified on
 * BitPackArray construction are stored).
 */
static inline uint32_t
bpaGetUInt32(const BitPackArray *array, BitOffset index)
{
  assert(array && index < array->numElems
         && array->bitsPerElem <= sizeof(uint32_t)*CHAR_BIT);
  return bsGetUInt32(array->store, array->bitsPerElem * index,
                     array->bitsPerElem);
}

/**
 * Stores unsigned integer in BitPackArray at given index.
 * @param bparray array to use
 * @param index index to store value at.
 * @param val value to store (only as many bits as specified on
 * BitPackArray construction are stored).
 */
static inline void
bpaStoreUInt64(BitPackArray *array, BitOffset index, uint64_t val)
{
  assert(array && index < array->numElems
         && array->bitsPerElem <= sizeof(val)*CHAR_BIT);
  bsStoreUInt64(array->store, array->bitsPerElem * index,
                array->bitsPerElem, val);
}

/**
 * Stores unsigned integer in BitPackArray at given index.
 * @param bparray array to use
 * @param index index to store value at.
 * @param val value to store (only as many bits as specified on
 * BitPackArray construction are stored).
 */
static inline uint64_t
bpaGetUInt64(const BitPackArray *array, BitOffset index)
{
  assert(array && index < array->numElems
         && array->bitsPerElem <= sizeof(uint64_t)*CHAR_BIT);
  return bsGetUInt64(array->store, array->bitsPerElem * index,
                     array->bitsPerElem);
}

/**
 * Unit test function for BitPackArray.
 * @return 0 on success, -1 on error.
 */
extern int
bitPackArray_unit_test(Env *env);

#endif /* BITPACKARRAY_H_INCLUDED */
