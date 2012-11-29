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

#ifndef S_SPLINT_S
/**
 * \file bitpackarray.h
 * \brief The class presented in this file encapsulates a bitstring in
 * a data structure to store/retrieve integers at fixed bitlength and
 * integrate well with genome tools library.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */

#include <stdlib.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/bitpackstring.h"
#include "core/error.h"
#include "core/ma.h"

struct BitPackArray
{
  BitString store;
  BitOffset numElems;
  unsigned bitsPerElem;
};

/* SK added the following macro, as we need the structure component,
   to store a pointer to the address, which is used when mapping the
   bitpackarray as part of a larger memory map */

#define BITPACKARRAYSTOREVAR(BPA) (BPA)->store

typedef struct BitPackArray BitPackArray;

/**
 * determine size of BitPackArray structure.
 * @param bits number of bits to encode each value stored with
 * @param numValues number of values to store
 * @return size
 */

static inline size_t sizeofbitarray(unsigned bits, BitOffset numValues)
{
  return bitElemsAllocSize(bits * numValues) * sizeof (BitElem);
}

/**
 * Create new BitPackArray structure.
 * @param bits number of bits to encode each value stored with
 * @param numValues number of values to store
 * @return pointer to new BitPackArray structure or NULL on failure
 */
static inline BitPackArray *
bitpackarray_new(unsigned bits, BitOffset numValues, bool withstorealloc)
{
  BitPackArray *newBPA = gt_malloc(sizeof (*newBPA));
  if (newBPA) {
    if (withstorealloc) {
      if (!(newBPA->store = gt_calloc(bitElemsAllocSize(bits*numValues),
                                      sizeof (BitElem)))) {
        gt_free(newBPA);
        return NULL;
      }
    } else {
      newBPA->store = NULL;
    }
    newBPA->bitsPerElem = bits;
    newBPA->numElems = numValues;
  }
  return newBPA;
}

static inline void
bitpackarray_delete(BitPackArray *bpa)
{
  if (!bpa) return;
  gt_free(bpa->store);
  gt_free(bpa);
}

/**
 * Stores unsigned integer in BitPackArray at given index.
 * @param bparray array to use
 * @param index index to store value at.
 * @param val value to store (only as many bits as specified on
 * BitPackArray construction are stored).
 */
static inline void
bitpackarray_store_uint32(BitPackArray *array, BitOffset index, uint32_t val)
{
  gt_assert(array && index < array->numElems
            && array->bitsPerElem <= sizeof (val)*CHAR_BIT);
  gt_bsStoreUInt32(array->store,
                   array->bitsPerElem * index,
                   array->bitsPerElem,
                   val);
}

/**
 * Stores unsigned integer in BitPackArray at given index.
 * @param bparray array to use
 * @param index index to store value at.
 * @param val value to store (only as many bits as specified on
 * BitPackArray construction are stored).
 */
static inline uint32_t
bitpackarray_get_uint32(const BitPackArray *array, BitOffset index)
{
  gt_assert(array && index < array->numElems
            && array->bitsPerElem <= sizeof (uint32_t)*CHAR_BIT);
  return gt_bsGetUInt32(array->store, array->bitsPerElem * index,
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
bitpackarray_store_uint64(BitPackArray *array, BitOffset index, uint64_t val)
{
  gt_assert(array && index < array->numElems
            && array->bitsPerElem <= sizeof (val)*CHAR_BIT);
  gt_bsStoreUInt64(array->store, array->bitsPerElem * index,
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
bitpackarray_get_uint64(const BitPackArray *array, BitOffset index)
{
  gt_assert(array && index < array->numElems
            && array->bitsPerElem <= sizeof (uint64_t)*CHAR_BIT);
  return gt_bsGetUInt64(array->store, array->bitsPerElem * index,
                        array->bitsPerElem);
}

/**
 * Unit test function for BitPackArray.
 * @return 0 on success, -1 on error.
 */
int gt_bitpackarray_unit_test(GtError*);

/*
static inline void showbitpackarray(const BitPackArray *bitpackarray)
{
  unsigned long numofunits, idx;

  gt_assert(bitpackarray != NULL);
  gt_assert(bitpackarray->store != NULL);
  numofunits = (unsigned long) sizeofbitarray(bitpackarray->bitsPerElem,
                                              bitpackarray->numElems);
  printf("numofunits=%lu\n",numofunits);
  for (idx=0; idx < numofunits; idx++)
  {
    printf("%lu: %u\n",idx,(unsigned int) bitpackarray->store[idx]);
    fflush(stdout);
  }
}
*/

#else

/* some minimal set of declaration to satisfy splint. We cannot use the
   original implementation, as this produces many errors when fed into
   splint
*/

typedef unsigned char BitElem;
typedef BitElem *BitString;

typedef struct
{
  BitString store;
  /* and others */
} BitPackArray;

typedef unsigned long long BitOffset;
size_t sizeofbitarray(unsigned bits, BitOffset numValues);
void bitpackarray_delete(BitPackArray *bpa);
BitPackArray *bitpackarray_new(unsigned bits, BitOffset numValues,
                               bool withstorealloc);
void bitpackarray_store_uint32(BitPackArray *array, BitOffset index,
                               uint32_t val);
void bitpackarray_store_uint64(BitPackArray *array, BitOffset index,
                               uint64_t val);
uint32_t bitpackarray_get_uint32(const BitPackArray *array, BitOffset index);
uint64_t bitpackarray_get_uint64(const BitPackArray *array, BitOffset index);
#endif

#endif
