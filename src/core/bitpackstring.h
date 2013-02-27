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

#ifndef BITPACKSTRING_H
#define BITPACKSTRING_H

/**
 * \file bitpackstring.h
 *
 * \brief This file introduces methods and types to handle a sequence
 * of bytes as concatenation of variable size integers, encoded at
 * user-requested number of bits.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */

#include <inttypes.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "core/error.h"

/** basic unit of addressing BitStrings */
/* Caution: sizeof (BitElem) <= sizeof (unsigned long) must be met */
typedef unsigned char BitElem;
/** Stores arbitrary bit sequences */
typedef BitElem *BitString;
/** Immutable BitStrings */
typedef const BitElem *constBitString;
/** Since even on 32-bit architectures 512MiB hold over 2 billion
 * bits, offsets into BitStrings must be stored as unsigned long long x*/
typedef unsigned long long BitOffset;

enum {
  /** bits held in one BitElem */
  bitElemBits = sizeof (BitElem)*CHAR_BIT,
};

/*
 * Uglyness ahead: (uint32_t *)(void *) double-casts are necessary to prevent
 * aliasing warnings on platforms that don't treat uint32_t as
 * compatible to unsigned even though both are internally the same
 * type (I'm looking at you Windows!)
 */

/** bits required to store an unsigned value ranging from 0..val */
#define gt_requiredUIntBits(val) gt_requiredUInt32Bits(val)
/** bits required to store int value val exactly */
#define gt_requiredIntBits(val) gt_requiredInt32Bits(val)
/** extract unsigned from BitString */
#define gt_bsGetUInt(str, offset, numBits) gt_bsGetUInt32(str, offset, numBits)
/** extract int from BitString */
#define gt_bsGetInt(str, offset, numBits) gt_bsGetInt32(str, offset, numBits)
/** store unsigned in BitString */
#define gt_bsStoreUInt(str, offset, numBits, val) \
  gt_bsStoreUInt32(str, offset, numBits, val)
/** store int in BitString */
#define gt_bsStoreInt(str, offset, numBits, val) \
  gt_bsStoreInt32(str, offset, numBits, val)
/** store array of unsigned ints in BitString */
#define gt_bsStoreUniformUIntArray(str, offset, numBits, numValues, val) \
  gt_bsStoreUniformUInt32Array(str, offset, numBits, numValues, \
                               (uint32_t *)(void *)val)
/** store array of unsigned ints in BitString */
#define gt_bsStoreNonUniformUIntArray(str, offset, numValues, bitsTotal, \
                                      numBitsList, val) \
  gt_bsStoreNonUniformUInt32Array(str, offset, numValues, bitsTotal, \
                                  numBitsList, (uint32_t *)(void *)val)
/** store array of unsigned ints in BitString */
#define gt_bsStoreNonUniformIntArray(str, offset, numValues, bitsTotal, \
                                     numBitsList, val) \
  gt_bsStoreNonUniformInt32Array(str, offset, numValues, bitsTotal, \
                                 numBitsList, (int32_t *)val)
/** store array of ints in BitString */
#define gt_bsStoreUniformIntArray(str, offset, numBits, numValues, val) \
  gt_bsStoreUniformInt32Array(str, offset, numBits, numValues, (int32_t *)val)
/** get array of unsigned ints from BitString */
#define gt_bsGetUniformUIntArray(str, offset, numBits, numValues, val) \
  gt_bsGetUniformUInt32Array(str, offset, numBits, numValues, \
                             (uint32_t *)(void *)val)
/** get array of unsigned ints from BitString */
#define gt_bsGetNonUniformUIntArray(str, offset, numValues, bitsTotal, \
                                    numBitsList, val) \
  gt_bsGetNonUniformUInt32Array(str, offset, numValues, bitsTotal, \
                                numBitsList, (uint32_t *)(void *)val)
/** get array of ints from BitString */
#define gt_bsGetUniformIntArray(str, offset, numBits, numValues, val) \
  gt_bsGetUniformInt32Array(str, offset, numBits, numValues, \
                            (int32_t *)(void *)val)
#define gt_bsGetNonUniformIntArray(str, offset, numValues, bitsTotal, \
                                   numBitsList, val) \
  gt_bsGetNonUniformInt32Array(str, offset, numValues, bitsTotal, \
                               numBitsList, (int32_t *)(void *)val)

/**
 * \brief Computes number of BitElem objects needed to store requested number
 * of consecutive bits.
 * @param numBits number of bits to be stored in total
 * @return number of BitElem objects needed to store numBits
 * consecutive bits or SIZE_MAX / sizeof (BitElem) if numBits exceeds
 * the addressable storage
 */
static inline size_t
bitElemsAllocSize(BitOffset numBits);
/**
 * \brief Computes \f$\gt_log_2v + 1\f$, where \f$v\f$ is of unsigned
 * 8-bit type.
 *
 * An unsigned integer \f$v\f$ can be represented in this many bits.
 * @param val value to compute log2 of.
 * @return \f$\lfloor{}gt_log_2v\rfloor + 1\f$
 */
static inline int
gt_requiredUInt8Bits(uint8_t val);
/**
 * \brief Computes \f$\gt_log_2v + 1\f$, where \f$v\f$ is of unsigned
 * 16-bit type.
 *
 * An unsigned integer \f$v\f$ can be represented in this many bits.
 * @param val value to compute log2 of.
 * @return \f$\lfloor{}gt_log_2v\rfloor + 1\f$
 */
static inline int
gt_requiredUInt16Bits(uint16_t val);
/**
 * \brief Computes \f$\gt_log_2v + 1\f$, where \f$v\f$ is of unsigned
 * 32-bit type.
 *
 * An unsigned integer \f$v\f$ can be represented in this many bits.
 * @param val value to compute log2 of.
 * @return \f$\lfloor{}gt_log_2v\rfloor + 1\f$
 */
int gt_requiredUInt32Bits(uint32_t val);
/**
 * \brief Computes \f$\gt_log_2v + 1\f$, where \f$v\f$ is of unsigned
 * 64-bit type.
 *
 * An unsigned integer \f$v\f$ can be represented in this many bits.
 * @param val value to compute log2 of.
 * @return \f$\lfloor{}gt_log_2v\rfloor + 1\f$
 */
int gt_requiredUInt64Bits(uint64_t val);
/**
 * \brief Computes \f$\gt_log_2v + 2\f$, where \f$v\f$ is of signed 8-bit
 * type.
 *
 * An integer \f$v\f$ can be represented in this many bits.
 * @param v value to compute log2 of.
 * @return \f$\lfloor{}|gt_log_2v|\rfloor + 2\f$. Special case: INT8_MIN
 * is representable in 8 bits although
 * \f[\left\lfloor{}gt_log_2\left|-(2^{7})\right|\right\rfloor + 2 = 9.\f]
 */
static inline int
gt_requiredInt8Bits(int8_t v);
/**
 * \brief Computes \f$\gt_log_2v + 2\f$, where \f$v\f$ is of signed
 * 16-bit type.
 *
 * An integer \f$v\f$ can be represented in this many bits.
 * @param v value to compute log2 of.
 * @return \f$\lfloor{}|gt_log_2v|\rfloor + 2\f$. Special case: INT16_MIN
 * is representable in 16 bits although
 * \f[\left\lfloor{}gt_log_2\left|-(2^{15})\right|\right\rfloor + 2 = 17.\f]
 */
static inline int
gt_requiredInt16Bits(int16_t v);
/**
 * \brief Computes \f$\gt_log_2v + 2\f$, where \f$v\f$ is of signed
 * 32-bit type.
 *
 * An integer \f$v\f$ can be represented in this many bits.
 * @param v value to compute log2 of.
 * @return \f$\lfloor{}|gt_log_2v|\rfloor + 2\f$. Special case: INT32_MIN
 * is representable in 32 bits although
 * \f[\left\lfloor{}gt_log_2\left|-(2^{31})\right|\right\rfloor + 2 = 33.\f]
 */
static inline int
gt_requiredInt32Bits(int32_t v);
/**
 * \brief Computes \f$\gt_log_2v + 2\f$, where \f$v\f$ is of signed
 * 64-bit type.
 *
 * An integer \f$v\f$ can be represented in this many bits.
 * @param v value to compute log2 of.
 * @return \f$\left\lfloor{}gt_log_2|v|\right\rfloor + 2\f$
 * (Special case: INT64_MIN is representable in 64 bits)
 */
static inline int
gt_requiredInt64Bits(int64_t v);
/**
 * \brief Retrieve unsigned integer of specified length from
 * bitstring at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
uint8_t gt_bsGetUInt8(constBitString str, BitOffset offset, unsigned numBits);
/**
 * \brief Retrieve unsigned integer of specified length from bitstring
 * at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
uint16_t gt_bsGetUInt16(constBitString str, BitOffset offset, unsigned numBits);
/**
 * \brief Retrieve unsigned integer of specified length from bitstring
 * at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
uint32_t gt_bsGetUInt32(constBitString str, BitOffset offset, unsigned numBits);
/**
 * \brief Retrieve unsigned integer of specified length from bitstring
 * at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
uint64_t gt_bsGetUInt64(constBitString str, BitOffset offset, unsigned numBits);
/**
 * \brief Store unsigned integer of specified length in bitstring at
 * given position.
 *
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 * @param val value to store
 */
void gt_bsStoreUInt8(BitString str, BitOffset offset, unsigned numBits,
                     uint8_t val);
/**
 * \brief Store unsigned integer of specified length in bitstring at
 * given position.
 *
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 * @param val value to store
 */
void gt_bsStoreUInt16(BitString str, BitOffset offset, unsigned numBits,
                      uint16_t val);
/**
 * \brief Store unsigned integer of specified length in bitstring at
 * given position.
 *
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 * @param val value to store
 */
void gt_bsStoreUInt32(BitString str, BitOffset offset, unsigned numBits,
                      uint32_t val);
/**
 * \brief Store unsigned integer of specified length in bitstring at
 * given position.
 *
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 * @param val value to store
 */
void gt_bsStoreUInt64(BitString str, BitOffset offset, unsigned numBits,
                      uint64_t val);
/**
 * \brief Retrieve integer of specified length from bitstring at given
 * position.
 *
 * @param str bitstring to read from
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 * @return value read
 */
static inline int8_t
gt_bsGetInt8(constBitString str, BitOffset offset, unsigned numBits);
/**
 * \brief Retrieve integer of specified length from bitstring at given
 * position.
 *
 * @param str bitstring to read from
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 * @return value read
 */
static inline int16_t
gt_bsGetInt16(constBitString str, BitOffset offset, unsigned numBits);
/**
 * \brief Retrieve integer of specified length from bitstring at given
 * position.
 *
 * @param str bitstring to read from
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 * @return value read
 */
static inline int32_t
gt_bsGetInt32(constBitString str, BitOffset offset, unsigned numBits);
/**
 * \brief Retrieve integer of specified length from bitstring at given
 * position.
 *
 * @param str bitstring to read from
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 * @return value read
 */
static inline int64_t
gt_bsGetInt64(constBitString str, BitOffset offset, unsigned numBits);
/**
 * \brief Store integer of specified length in bitstring at given
 * position.
 *
 * @param str bitstring to write to
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @param val value to store
 */
static inline void
gt_bsStoreInt8(BitString str, BitOffset offset, unsigned numBits, uint8_t val);
/**
 * \brief Store integer of specified length in bitstring at given
 * position.
 *
 * @param str bitstring to write to
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @param val value to store
 */
static inline void
gt_bsStoreInt16(BitString str, BitOffset offset, unsigned numBits,
                uint16_t val);
/**
 * \brief Store integer of specified length in bitstring at given
 * position.
 *
 * @param str bitstring to write to
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @param val value to store
 */
static inline void
gt_bsStoreInt32(BitString str, BitOffset offset, unsigned numBits,
                uint32_t val);
/**
 * \brief Store integer of specified length in bitstring at given
 * position.
 *
 * @param str bitstring to write to
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @param val value to store
 */
static inline void
gt_bsStoreInt64(BitString str, BitOffset offset, unsigned numBits,
                uint64_t val);
/* higher level functions */

/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param val read integers from this array
 */
void
gt_bsStoreUniformUInt8Array(BitString str, BitOffset offset, unsigned numBits,
                            size_t numValues, const uint8_t val[]);
/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param val read integers from this array
 */
void
gt_bsStoreUniformUInt16Array(BitString str, BitOffset offset, unsigned numBits,
                             size_t numValues, const uint16_t val[]);
/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param val read integers from this array
 */
void
gt_bsStoreUniformUInt32Array(BitString str, BitOffset offset, unsigned numBits,
                             size_t numValues, const uint32_t val[]);
/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param val read integers from this array
 */
void
gt_bsStoreUniformUInt64Array(BitString str, BitOffset offset, unsigned numBits,
                             size_t numValues, const uint64_t val[]);

/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList number of bits composing each integer, one bit
 * count per value
 * @param val read integers from this array
 */
void
gt_bsStoreNonUniformUInt8Array(BitString str, BitOffset offset,
                               size_t numValues, BitOffset totalBitsLeft,
                               unsigned *numBitsList, const uint8_t val[]);
/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList number of bits composing each integer, one bit
 * count per value
 * @param val read integers from this array
 */
void
gt_bsStoreNonUniformUInt16Array(BitString str, BitOffset offset,
                                size_t numValues, BitOffset totalBitsLeft,
                                unsigned *numBitsList, const uint16_t val[]);
/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList number of bits composing each integer, one bit
 * count per value
 * @param val read integers from this array
 */
void
gt_bsStoreNonUniformUInt32Array(BitString str, BitOffset offset,
                                size_t numValues, BitOffset totalBitsLeft,
                                unsigned *numBitsList, const uint32_t val[]);
/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList number of bits composing each integer, one bit
 * count per value
 * @param val read integers from this array
 */
void
gt_bsStoreNonUniformUInt64Array(BitString str, BitOffset offset,
                                size_t numValues, BitOffset totalBitsLeft,
                                unsigned *numBitsList, const uint64_t val[]);

/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param val read integers from this array
 */
static inline void
gt_bsStoreUniformInt8Array(BitString str, BitOffset offset,
                           unsigned numBits, size_t numValues,
                           const int8_t val[]);
/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param val read integers from this array
 */
static inline void
gt_bsStoreUniformInt16Array(BitString str, BitOffset offset,
                            unsigned numBits, size_t numValues,
                            const int16_t val[]);
/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param val read integers from this array
 */
static inline void
gt_bsStoreUniformInt32Array(BitString str, BitOffset offset,
                            unsigned numBits, size_t numValues,
                            const int32_t val[]);
/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param val read integers from this array
 */
static inline void
gt_bsStoreUniformInt64Array(BitString str, BitOffset offset,
                            unsigned numBits, size_t numValues,
                            const int64_t val[]);

/**
 * \brief Retrieve n unsigned integers of specified length from
 * bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param val store integers read in this array
 */
void
gt_bsGetUniformUInt8Array(constBitString str, BitOffset offset,
                          unsigned numBits, size_t numValues, uint8_t val[]);
/**
 * \brief Retrieve n unsigned integers of specified length from
 * bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param val store integers read in this array
 */
void
gt_bsGetUniformUInt16Array(constBitString str, BitOffset offset,
                           unsigned numBits, size_t numValues, uint16_t val[]);
/**
 * \brief Retrieve n unsigned integers of specified length from
 * bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param val store integers read in this array
 */
void
gt_bsGetUniformUInt32Array(constBitString str, BitOffset offset,
                           unsigned numBits, size_t numValues, uint32_t val[]);
/**
 * \brief Retrieve n unsigned integers of specified length from
 * bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param val store integers read in this array
 */
void
gt_bsGetUniformUInt64Array(constBitString str, BitOffset offset,
                           unsigned numBits, size_t numValues, uint64_t val[]);

/**
 * \brief Retrieve n integers of specified length from bitstring,
 * starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @param numValues number of integers to read
 * @param val store integers read in this array
 */
static inline void
gt_bsGetUniformInt8Array(constBitString str, BitOffset offset,
                         unsigned numBits, size_t numValues, int8_t val[]);
/**
 * \brief Retrieve n integers of specified length from bitstring,
 * starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @param numValues number of integers to read
 * @param val store integers read in this array
 */
static inline void
gt_bsGetUniformInt16Array(constBitString str, BitOffset offset,
                          unsigned numBits, size_t numValues, int16_t val[]);
/**
 * \brief Retrieve n integers of specified length from bitstring,
 * starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @param numValues number of integers to read
 * @param val store integers read in this array
 */
static inline void
gt_bsGetUniformInt32Array(constBitString str, BitOffset offset,
                          unsigned numBits, size_t numValues, int32_t val[]);
/**
 * \brief Retrieve n integers of specified length from bitstring,
 * starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @param numValues number of integers to read
 * @param val store integers read in this array
 */
static inline void
gt_bsGetUniformInt64Array(constBitString str, BitOffset offset,
                          unsigned numBits, size_t numValues, int64_t val[]);

/**
 * \brief Retrieve n unsigned integers of specified length from
 * bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList vector of number of bits composing each integer
 * @param val store integers read in this array
 */
void
gt_bsGetNonUniformUInt8Array(constBitString str, BitOffset offset,
                             size_t numValues, BitOffset numBitsTotal,
                             unsigned numBitsList[], uint8_t val[]);
/**
 * \brief Retrieve n unsigned integers of specified length from
 * bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList vector of number of bits composing each integer
 * @param val store integers read in this array
 */
void
gt_bsGetNonUniformUInt16Array(constBitString str, BitOffset offset,
                              size_t numValues, BitOffset numBitsTotal,
                              unsigned numBitsList[], uint16_t val[]);
/**
 * \brief Retrieve n unsigned integers of specified length from
 * bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList vector of number of bits composing each integer
 * @param val store integers read in this array
 */
void
gt_bsGetNonUniformUInt32Array(constBitString str, BitOffset offset,
                              size_t numValues, BitOffset numBitsTotal,
                              unsigned numBitsList[], uint32_t val[]);
/**
 * \brief Retrieve n unsigned integers of specified length from
 * bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList vector of number of bits composing each integer
 * @param val store integers read in this array
 */
void
gt_bsGetNonUniformUInt64Array(constBitString str, BitOffset offset,
                              size_t numValues, BitOffset numBitsTotal,
                              unsigned numBitsList[], uint64_t val[]);

/**
 * \brief Retrieve n unsigned integers of specified length from
 * bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList vector of number of bits composing each integer
 * @param val store integers read in this array
 */
void
gt_bsGetNonUniformInt8Array(constBitString str, BitOffset offset,
                            size_t numValues, BitOffset numBitsTotal,
                            unsigned numBitsList[], int8_t val[]);
/**
 * \brief Retrieve n unsigned integers of specified length from
 * bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList vector of number of bits composing each integer
 * @param val store integers read in this array
 */
void
gt_bsGetNonUniformInt16Array(constBitString str, BitOffset offset,
                             size_t numValues, BitOffset numBitsTotal,
                             unsigned numBitsList[], int16_t val[]);
/**
 * \brief Retrieve n unsigned integers of specified length from
 * bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList vector of number of bits composing each integer
 * @param val store integers read in this array
 */
void
gt_bsGetNonUniformInt32Array(constBitString str, BitOffset offset,
                             size_t numValues, BitOffset numBitsTotal,
                             unsigned numBitsList[], int32_t val[]);
/**
 * \brief Retrieve n unsigned integers of specified length from
 * bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList vector of number of bits composing each integer
 * @param val store integers read in this array
 */
void
gt_bsGetNonUniformInt64Array(constBitString str, BitOffset offset,
                             size_t numValues, BitOffset numBitsTotal,
                             unsigned numBitsList[], int64_t val[]);

/**
 * \brief Add n unsigned integers of specified length from
 * bitstring to the given vector of values, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param val sum of integer read and the value in val will be
 * returned in this array
 */
void
gt_bsGetUniformUInt8ArrayAdd(constBitString str, BitOffset offset,
                             unsigned numBits, size_t numValues,
                             uint8_t val[]);
/**
 * \brief Add n unsigned integers of specified length from
 * bitstring to the given vector of values, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param val sum of integer read and the value in val will be
 * returned in this array
 */
void
gt_bsGetUniformUInt16ArrayAdd(constBitString str, BitOffset offset,
                              unsigned numBits, size_t numValues,
                              uint16_t val[]);
/**
 * \brief Add n unsigned integers of specified length from
 * bitstring to the given vector of values, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param val sum of integer read and the value in val will be
 * returned in this array
 */
void
gt_bsGetUniformUInt32ArrayAdd(constBitString str, BitOffset offset,
                              unsigned numBits, size_t numValues,
                              uint32_t val[]);
/**
 * \brief Add n unsigned integers of specified length from
 * bitstring to the given vector of values, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param val sum of integer read and the value in val will be
 * returned in this array
 */
void
gt_bsGetUniformUInt64ArrayAdd(constBitString str, BitOffset offset,
                              unsigned numBits, size_t numValues,
                              uint64_t val[]);

/**
 * \brief Add n unsigned integers of specified length from
 * bitstring to the given vector of values, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList vector of number of bits composing each integer
 * @param val sum of integer read and the value in val will be
 * returned in this array
 */
void
gt_bsGetNonUniformUInt8ArrayAdd(constBitString str, BitOffset offset,
                                size_t numValues, BitOffset numBitsTotal,
                                unsigned numBitsList[], uint8_t val[]);
/**
 * \brief Add n unsigned integers of specified length from
 * bitstring to the given vector of values, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList vector of number of bits composing each integer
 * @param val sum of integer read and the value in val will be
 * returned in this array
 */
void
gt_bsGetNonUniformUInt16ArrayAdd(constBitString str, BitOffset offset,
                                 size_t numValues, BitOffset numBitsTotal,
                                 unsigned numBitsList[], uint16_t val[]);
/**
 * \brief Add n unsigned integers of specified length from
 * bitstring to the given vector of values, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList vector of number of bits composing each integer
 * @param val sum of integer read and the value in val will be
 * returned in this array
 */
void
gt_bsGetNonUniformUInt32ArrayAdd(constBitString str, BitOffset offset,
                                 size_t numValues, BitOffset numBitsTotal,
                                 unsigned numBitsList[], uint32_t val[]);
/**
 * \brief Add n unsigned integers of specified length from
 * bitstring to the given vector of values, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList vector of number of bits composing each integer
 * @param val sum of integer read and the value in val will be
 * returned in this array
 */
void
gt_bsGetNonUniformUInt64ArrayAdd(constBitString str, BitOffset offset,
                                 size_t numValues, BitOffset numBitsTotal,
                                 unsigned numBitsList[], uint64_t val[]);

/**
 * \brief Add n unsigned integers of specified length from
 * bitstring to the given vector of values, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList vector of number of bits composing each integer
 * @param val sum of integer read and the value in val will be
 * returned in this array
 */
void
gt_bsGetNonUniformInt8ArrayAdd(constBitString str, BitOffset offset,
                               size_t numValues, BitOffset numBitsTotal,
                               unsigned numBitsList[], int8_t val[]);
/**
 * \brief Add n unsigned integers of specified length from
 * bitstring to the given vector of values, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList vector of number of bits composing each integer
 * @param val sum of integer read and the value in val will be
 * returned in this array
 */
void
gt_bsGetNonUniformInt16ArrayAdd(constBitString str, BitOffset offset,
                                size_t numValues, BitOffset numBitsTotal,
                                unsigned numBitsList[], int16_t val[]);
/**
 * \brief Add n unsigned integers of specified length from
 * bitstring to the given vector of values, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList vector of number of bits composing each integer
 * @param val sum of integer read and the value in val will be
 * returned in this array
 */
void
gt_bsGetNonUniformInt32ArrayAdd(constBitString str, BitOffset offset,
                                size_t numValues, BitOffset numBitsTotal,
                                unsigned numBitsList[], int32_t val[]);
/**
 * \brief Add n unsigned integers of specified length from
 * bitstring to the given vector of values, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numValues number of integers to read
 * @param numBitsTotal number of bits to be read in total, must equal
 * \f$\sum^\mathrm{numValues}_i numBitsList[i]\f$
 * @param numBitsList vector of number of bits composing each integer
 * @param val sum of integer read and the value in val will be
 * returned in this array
 */
void
gt_bsGetNonUniformInt64ArrayAdd(constBitString str, BitOffset offset,
                                size_t numValues, BitOffset numBitsTotal,
                                unsigned numBitsList[], int64_t val[]);

/**
 * \brief Compares substrings of bitstrings, starting at respective offsets,
 * both of lenght numBits.
 *
 * The bitstrings are treated as MSB-first encodings (which they are
 * if produced by above functions) and thus after the first
 * most-significant bit set in only one of the strings the comparision
 * can terminate.
 *
 * This function treats two sub-bitstrings as
 * @param a first bitstring to compare
 * @param offsetA corresponding start position
 * @param numBitsA length of substring \f$a'\f$ in a to use for comparison.
 * @param b second bitstring
 * @param offsetB corresponding start position
 * @param numBitsB length of substring \f$b'\f$ in b to use for comparison.
 * @return 0 for equality, \f$-1\f$ if \f$a < b\f$, \f$1\f$ if \f$b > a\f$
 */
int gt_bsCompare(constBitString a, BitOffset offsetA, BitOffset numBitsA,
                 constBitString b, BitOffset offsetB, BitOffset numBitsB);

/**
 * \brief copy (sub-)bitstring to another position (in same or other
 * bitstring)
 * @param src bitsring to copy from
 * @param offsetSrc bit position in src to start copying from
 * @param dest bitstring to copy to
 * @param offsetDest bit position in dest to start writing at
 * @param numBits number of bits to copy
 */
void gt_bsCopy(constBitString src, BitOffset offsetSrc,
               BitString dest, BitOffset offsetDest, BitOffset numBits);

/**
 * \brief set (sub-)bitstring to all one or zero bits
 * @param str bitsring to reset (portion of)
 * @param offset bit position in str to start at
 * @param numBits number of bits to copy
 * @param bitVal set all numBits bits to 0 if 0, to 1 otherwise
 */
void gt_bsClear(BitString str, BitOffset offset, BitOffset numBits, int bitVal);

/**
 * \brief set singular bit in bitstring to 1
 * @param str bitstring to modify
 * @param pos selects bit to set
 */
static inline void
bsSetBit(BitString str, BitOffset pos);

/**
 * \brief clear, i.e. set to 0 singular bit in bitstring
 * @param str bitstring to modify
 * @param pos selects bit to clear
 */
static inline void
gt_bsClearBit(BitString str, BitOffset pos);

/**
 * \brief XOR singular bit in bitstring
 * @param str bitstring to modify
 * @param pos selects bit to invert
 */
static inline void
bsToggleBit(BitString str, BitOffset pos);

/**
 * \brief Query value of single bit in bitstring.
 * @param str bitstring to read from
 * @param pos selects bit to retrieve
 * @return 1 if selected bit is set, 0 if not set
 */
static inline int
gt_bsGetBit(constBitString str, BitOffset pos);
/**
 * \brief Compute Hamming weight of (sub-)bitstring
 * @return number of bits set in (sub-)bitstring
 */
BitOffset gt_bs1BitsCount(constBitString str, BitOffset offset,
                          BitOffset numBits);

/**
 * \brief Print sequence of 0s and 1s to stream to display BitString
 * contents.
 * @param fp
 * @param str
 * @param offset start reading at this position
 * @param numBits print 0,1 chars for this many bits
 * @return -1 in case of error, 0 otherwise
 */
int gt_bsPrint(FILE *fp, constBitString str, BitOffset offset,
               BitOffset numBits);

/**
 * \brief Meta-Unit test function for bitPackString, calls all functions
 * mentioned below.
 * @return 0 on success, -1 on error.
 */
int gt_bitPackString_unit_test(GtError*);
/**
 * \brief Unit test function for bitPackString, integer functions.
 * @return 0 on success, -1 on error.
 */
int gt_bitPackStringInt_unit_test(GtError*);
/**
 * \brief Unit test function for bitPackString, 8-bit functions.
 * @return 0 on success, -1 on error.
 */
int gt_bitPackStringInt8_unit_test(GtError*);
/**
 * \brief Unit test function for bitPackString, 16-bit functions.
 * @return 0 on success, -1 on error.
 */
int gt_bitPackStringInt16_unit_test(GtError*);
/**
 * \brief Unit test function for bitPackString, 32-bit functions.
 * @return 0 on success, -1 on error.
 */
int gt_bitPackStringInt32_unit_test(GtError*);
/**
 * \brief Unit test function for bitPackString, 64-bit functions.
 * @return 0 on success, -1 on error.
 */
int gt_bitPackStringInt64_unit_test(GtError*);

#include "core/bitpackstringsimpleop.h"

#endif
