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

#ifndef BITPACKARRAY_H
#define BITPACKARRAY_H

/**
 * \file bitpackarray.h
 * \brief The class presented in this file encapsulates a bitstring in
 * a data structure to store/retrieve integers at fixed bitlength and
 * integrate well with genome tools library.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */

#include <assert.h>
#include <stdlib.h>
#include "libgtcore/bitpackstring.h"
#include "libgtcore/error.h"
#include "libgtcore/ma.h"

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
newBitPackArray(unsigned bits, BitOffset numValues)
{
  BitPackArray *newBPA = ma_malloc(sizeof (*newBPA));
  if (newBPA)
  {
    if (!(newBPA->store = ma_malloc(bitElemsAllocSize(bits*numValues)
                                    * sizeof (BitElem))))
    {
      ma_free(newBPA);
      return NULL;
    }
    newBPA->bitsPerElem = bits;
    newBPA->numElems = numValues;
  }
  return newBPA;
}

static inline void
deleteBitPackArray(BitPackArray *bpa)
{
  ma_free(bpa->store);
  ma_free(bpa);
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
         && array->bitsPerElem <= sizeof (val)*CHAR_BIT);
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
         && array->bitsPerElem <= sizeof (uint32_t)*CHAR_BIT);
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
         && array->bitsPerElem <= sizeof (val)*CHAR_BIT);
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
         && array->bitsPerElem <= sizeof (uint64_t)*CHAR_BIT);
  return bsGetUInt64(array->store, array->bitsPerElem * index,
                     array->bitsPerElem);
}

/**
 * Unit test function for BitPackArray.
 * @return 0 on success, -1 on error.
 */
extern int
bitPackArray_unit_test(Error*);

#endif
