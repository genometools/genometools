/*
** Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
**  
** See LICENSE file or http://genometools.org/license.html for license details.
** 
*/
#ifndef BITPACKSTRING_H_INCLUDED
#define BITPACKSTRING_H_INCLUDED
/**
 * \file bitpackstring.h
 *
 * \brief This file introduces methods and types to handle a sequence
 * of bytes as concatenation of variable size integers, encoded at
 * some number of bits.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */
#include <limits.h>
#include <inttypes.h>
#include <stdlib.h>

/* Caution: sizeof(bitElem) <= sizeof(uint_fast32_t) must be met */
typedef unsigned char bitElem;
typedef bitElem *bitString;
typedef size_t bitOffset;

enum {
  bitElemBits = sizeof(bitElem)*CHAR_BIT,
};

#define requiredUIntBits(val) requiredUInt32Bits(val)
#define requiredIntBits(val) requiredInt32Bits(val)
#define bsGetUInt(str, offset, numBits) bsGetUInt32(str, offset, numBits)
#define bsGetInt(str, offset, numBits) bsGetInt32(str, offset, numBits)
#define bsStoreUInt(str, offset, numBits, val) \
  bsStoreUInt32(str, offset, numBits, val)
#define bsStoreInt(str, offset, numBits, val) \
  bsStoreInt32(str, offset, numBits, val)
#define bsStoreUniformUIntArray(str, offset, numBits, numValues, val) \
  bsStoreUniformUInt32Array(str, offset, numBits, numValues, val)
#define bsStoreUniformIntArray(str, offset, numBits, numValues, val) \
  bsStoreUniformInt32Array(str, offset, numBits, numValues, val)
#define bsGetUniformUIntArray(str, offset, numBits, numValues, val) \
  bsGetUniformUInt32Array(str, offset, numBits, numValues, val)
#define bsGetUniformIntArray(str, offset, numBits, numValues, val) \
  bsGetUniformInt32Array(str, offset, numBits, numValues, val)

static inline size_t
bitElemsAllocSize(bitOffset numBits);
/**
 * Computes \f$\log_2v + 1\f$, where \f$v\f$ is of unsigned 8-bit type.
 * An unsigned integer \f$v\f$ can be represented in this many bits.
 * @param val value to find log2 of.
 * @return \f$\lfloor{}log_2v\rfloor + 1\f$
 */
static inline int
requiredUInt8Bits(uint8_t val);
/**
 * Computes \f$\log_2v + 1\f$, where \f$v\f$ is of unsigned 16-bit type.
 * An unsigned integer \f$v\f$ can be represented in this many bits.
 * @param val value to find log2 of.
 * @return \f$\lfloor{}log_2v\rfloor + 1\f$
 */
static inline int
requiredUInt16Bits(uint16_t val);
/**
 * Computes \f$\log_2v + 1\f$, where \f$v\f$ is of unsigned 32-bit type.
 * An unsigned integer \f$v\f$ can be represented in this many bits.
 * @param val value to find log2 of.
 * @return \f$\lfloor{}log_2v\rfloor + 1\f$
 */
extern int
requiredUInt32Bits(uint32_t val);
/**
 * Computes \f$\log_2v + 1\f$, where \f$v\f$ is of unsigned 64-bit type.
 * An unsigned integer \f$v\f$ can be represented in this many bits.
 * @param val value to find log2 of.
 * @return \f$\lfloor{}log_2v\rfloor + 1\f$
 */
extern int
requiredUInt64Bits(uint64_t val);

/**
 * Computes \f$\log_2v + 2\f$, where \f$v\f$ is of signed 8-bit type.
 * An integer \f$v\f$ can be represented in this many bits.
 * @param val value to find log2 of.
 * @return \f$\lfloor{}|log_2v|\rfloor + 2\f$. Special case: INT8_MIN
 * is representable in 8 bits although
 * \f[\left\lfloor{}log_2\left|-(2^{7})\right|\right\rfloor + 2 = 9.\f]
 */
static inline int
requiredInt8Bits(int8_t v);
/**
 * Computes \f$\log_2v + 2\f$, where \f$v\f$ is of signed 16-bit type.
 * An integer \f$v\f$ can be represented in this many bits.
 * @param val value to find log2 of.
 * @return \f$\lfloor{}|log_2v|\rfloor + 2\f$. Special case: INT16_MIN
 * is representable in 16 bits although
 * \f[\left\lfloor{}log_2\left|-(2^{15})\right|\right\rfloor + 2 = 17.\f]
 */
static inline int
requiredInt16Bits(int16_t v);
/**
 * Computes \f$\log_2v + 2\f$, where \f$v\f$ is of signed 32-bit type.
 * An integer \f$v\f$ can be represented in this many bits.
 * @param val value to find log2 of.
 * @return \f$\lfloor{}|log_2v|\rfloor + 2\f$. Special case: INT32_MIN
 * is representable in 32 bits although
 * \f[\left\lfloor{}log_2\left|-(2^{31})\right|\right\rfloor + 2 = 33.\f]
 */
static inline int
requiredInt32Bits(int32_t val);
/**
 * Computes \f$\log_2v + 2\f$, where \f$v\f$ is of signed 64-bit type.
 * An integer \f$v\f$ can be represented in this many bits.
 * @param val value to find log2 of.
 * @return \f$\left\lfloor{}log_2|v|\right\rfloor + 2\f$ (Special case: INT64_MIN
 * is representable in 64 bits)
 */
static inline int
requiredInt64Bits(int64_t val);
/**
 * Retrieve unsigned integer of specified length from bitstring at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
extern uint8_t
bsGetUInt8(const bitString str, bitOffset offset, unsigned numBits);
/**
 * Retrieve unsigned integer of specified length from bitstring at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
extern uint16_t
bsGetUInt16(const bitString str, bitOffset offset, unsigned numBits);
/**
 * Retrieve unsigned integer of specified length from bitstring at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
extern uint32_t
bsGetUInt32(const bitString str, bitOffset offset, unsigned numBits);
/**
 * Retrieve unsigned integer of specified length from bitstring at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
extern uint64_t
bsGetUInt64(const bitString str, bitOffset offset, unsigned numBits);
/**
 * Store unsigned integer of specified length in bitstring at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 */
extern void
bsStoreUInt8(bitString str, bitOffset offset, unsigned numBits, uint8_t val);
/**
 * Store unsigned integer of specified length in bitstring at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 */
extern void
bsStoreUInt16(bitString str, bitOffset offset, unsigned numBits, uint16_t val);
/**
 * Store unsigned integer of specified length in bitstring at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 */
extern void
bsStoreUInt32(bitString str, bitOffset offset, unsigned numBits, uint32_t val);
/**
 * Store unsigned integer of specified length in bitstring at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 */
extern void
bsStoreUInt64(bitString str, bitOffset offset, unsigned numBits, uint64_t val);
/**
 * Store integer of specified length in bitstring at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 */
static inline int8_t
bsGetInt8(bitString str, bitOffset offset, unsigned numBits);
/**
 * Store integer of specified length in bitstring at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 */
static inline int16_t
bsGetInt16(bitString str, bitOffset offset, unsigned numBits);
/**
 * Store integer of specified length in bitstring at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 */
static inline int32_t
bsGetInt32(bitString str, bitOffset offset, unsigned numBits);
/**
 * Store integer of specified length in bitstring at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 */
static inline int64_t
bsGetInt64(bitString str, bitOffset offset, unsigned numBits);
/**
 * Retrieve integer of specified length from bitstring at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
static inline void
bsStoreInt8(bitString str, bitOffset offset, unsigned numBits, uint8_t val);
/**
 * Retrieve integer of specified length from bitstring at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
static inline void
bsStoreInt16(bitString str, bitOffset offset, unsigned numBits, uint16_t val);
/**
 * Retrieve integer of specified length from bitstring at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
static inline void
bsStoreInt32(bitString str, bitOffset offset, unsigned numBits, uint32_t val);
/**
 * Retrieve integer of specified length from bitstring at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
static inline void
bsStoreInt64(bitString str, bitOffset offset, unsigned numBits, uint64_t val);
/* higher level functions */
/**
 * Store n unsigned integers of specified length from array, in bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param read integers from this array
 */
void
bsStoreUniformUInt8Array(bitString str, bitOffset offset, unsigned numBits,
                         size_t numValues, const uint8_t val[]);
/**
 * Store n unsigned integers of specified length from array, in bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param read integers from this array
 */
void
bsStoreUniformUInt16Array(bitString str, bitOffset offset, unsigned numBits,
                          size_t numValues, const uint16_t val[]);
/**
 * Store n unsigned integers of specified length from array, in bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param read integers from this array
 */
void
bsStoreUniformUInt32Array(bitString str, bitOffset offset, unsigned numBits,
                          size_t numValues, const uint32_t val[]);
/**
 * Store n unsigned integers of specified length from array, in bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param read integers from this array
 */
void
bsStoreUniformUInt64Array(bitString str, bitOffset offset, unsigned numBits,
                          size_t numValues, const uint64_t val[]);
/**
 * Store n unsigned integers of specified length from array, in bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param read integers from this array
 */
static inline void
bsStoreUniformInt8Array(bitString str, bitOffset offset, unsigned numBits,
                        size_t numValues, const int8_t val[]);
/**
 * Store n unsigned integers of specified length from array, in bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param read integers from this array
 */
static inline void
bsStoreUniformInt16Array(bitString str, bitOffset offset, unsigned numBits,
                         size_t numValues, const int16_t val[]);
/**
 * Store n unsigned integers of specified length from array, in bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param read integers from this array
 */
static inline void
bsStoreUniformInt32Array(bitString str, bitOffset offset, unsigned numBits,
                         size_t numValues, const int32_t val[]);
/**
 * Store n unsigned integers of specified length from array, in bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param read integers from this array
 */
static inline void
bsStoreUniformInt64Array(bitString str, bitOffset offset, unsigned numBits,
                         size_t numValues, const int64_t val[]);
/**
 * Retrieve n unsigned integers of specified length from bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param store integers read in this array
 */
void
bsGetUniformUInt8Array(const bitString str, bitOffset offset, unsigned numBits,
                       size_t numValues, uint8_t val[]);
/**
 * Retrieve n unsigned integers of specified length from bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param store integers read in this array
 */
void
bsGetUniformUInt16Array(const bitString str, bitOffset offset, unsigned numBits,
                        size_t numValues, uint16_t val[]);
/**
 * Retrieve n unsigned integers of specified length from bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param store integers read in this array
 */
void
bsGetUniformUInt32Array(const bitString str, bitOffset offset, unsigned numBits,
                        size_t numValues, uint32_t val[]);
/**
 * Retrieve n unsigned integers of specified length from bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param store integers read in this array
 */
void
bsGetUniformUInt64Array(const bitString str, bitOffset offset, unsigned numBits,
                        size_t numValues, uint64_t val[]);
/**
 * Retrieve n integers of specified length from bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @param numValues number of integers to read
 * @param store integers read in this array
 */
static inline void
bsGetUniformInt8Array(const bitString str, bitOffset offset, unsigned numBits,
                      size_t numValues, int8_t val[]);
/**
 * Retrieve n integers of specified length from bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @param numValues number of integers to read
 * @param store integers read in this array
 */
static inline void
bsGetUniformInt16Array(const bitString str, bitOffset offset, unsigned numBits,
                       size_t numValues, int16_t val[]);
/**
 * Retrieve n integers of specified length from bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @param numValues number of integers to read
 * @param store integers read in this array
 */
static inline void
bsGetUniformInt32Array(const bitString str, bitOffset offset, unsigned numBits,
                       size_t numValues, int32_t val[]);
/**
 * Retrieve n integers of specified length from bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @param numValues number of integers to read
 * @param store integers read in this array
 */
static inline void
bsGetUniformInt64Array(const bitString str, bitOffset offset, unsigned numBits,
                       size_t numValues, int64_t val[]);

extern int
bsCompare(const bitString a, bitOffset offsetA, bitOffset numBitsA,
          const bitString b, bitOffset offsetB, bitOffset numBitsB);

#include "bitpackstringsimpleop.h"

#endif /* BITPACKSTRING_H_INCLUDED */
