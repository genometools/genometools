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

struct bitPackArray
{
  bitString store;
  size_t numElems;
  int bitsPerElem;
};

/**
 * Create new bitPackArray structure.
 * @param bits number of bits to encode each value stored with
 * @param numValues number of values to store
 * @return pointer to new bitPackArray structure or NULL on failure
 */
static inline struct bitPackArray *
newBitPackArray(bitOffset bits, bitOffset numValues, Env *env)
{
  struct bitPackArray *newBPA = env_ma_malloc(env, sizeof(*newBPA));
  if(newBPA)
  {
    if(!(newBPA->store = env_ma_malloc(env, bitElemsAllocSize(bits*numValues)
                                       * sizeof(bitElem))))
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
deleteBitPackArray(struct bitPackArray *bpa, Env *env)
{
  env_ma_free(bpa->store, env);
  env_ma_free(bpa, env);
}

/**
 * Stores unsigned integer in bitPackArray at given index.
 * @param bparray array to use
 * @param index index to store value at.
 * @param val value to store (only as many bits as specified on
 * bitPackArray construction are stored).
 */
static inline void
bpaStoreUInt32(struct bitPackArray *array, bitOffset index, uint32_t val)
{
  assert(array && index < array->numElems
         && array->bitsPerElem <= sizeof(val)*CHAR_BIT);
  bsStoreUInt32(array->store, array->bitsPerElem * index,
                array->bitsPerElem, val);
}

/**
 * Stores unsigned integer in bitPackArray at given index.
 * @param bparray array to use
 * @param index index to store value at.
 * @param val value to store (only as many bits as specified on
 * bitPackArray construction are stored).
 */
static inline uint32_t
bpaGetUInt32(struct bitPackArray *array, bitOffset index)
{
  assert(array && index < array->numElems
         && array->bitsPerElem <= sizeof(uint32_t)*CHAR_BIT);
  return bsGetUInt32(array->store, array->bitsPerElem * index,
                     array->bitsPerElem);
}

/**
 * Stores unsigned integer in bitPackArray at given index.
 * @param bparray array to use
 * @param index index to store value at.
 * @param val value to store (only as many bits as specified on
 * bitPackArray construction are stored).
 */
static inline void
bpaStoreUInt64(struct bitPackArray *array, bitOffset index, uint64_t val)
{
  assert(array && index < array->numElems
         && array->bitsPerElem <= sizeof(val)*CHAR_BIT);
  bsStoreUInt64(array->store, array->bitsPerElem * index,
                array->bitsPerElem, val);
}

/**
 * Stores unsigned integer in bitPackArray at given index.
 * @param bparray array to use
 * @param index index to store value at.
 * @param val value to store (only as many bits as specified on
 * bitPackArray construction are stored).
 */
static inline uint64_t
bpaGetUInt64(struct bitPackArray *array, bitOffset index)
{
  assert(array && index < array->numElems
         && array->bitsPerElem <= sizeof(uint64_t)*CHAR_BIT);
  return bsGetUInt64(array->store, array->bitsPerElem * index,
                     array->bitsPerElem);
}

/**
 * Unit test function for bitPackArray.
 * @return 0 on success, -1 on error.
 */
extern int
bitPackArray_unit_test(Env *env);

#endif /* BITPACKARRAY_H_INCLUDED */
