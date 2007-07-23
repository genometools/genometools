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
 * user-requested number of bits.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */
#include <limits.h>
#include <inttypes.h>
#include <stdlib.h>

#include <libgtcore/env.h>

/* Caution: sizeof (BitElem) <= sizeof (uint_fast32_t) must be met */
typedef unsigned char BitElem;
typedef BitElem *BitString;
typedef unsigned long long BitOffset;

enum {
  bitElemBits = sizeof (BitElem)*CHAR_BIT,
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

/**
 * \brief Computes number of BitElem objects needed to store requested number
 * of consecutive bits.
 * @param numBits number of bits to be stored in total
 * @return number of BitElem objects needed to store numBits
 * consecutive bits
 */
static inline size_t
bitElemsAllocSize(BitOffset numBits);
/**
 * \brief Computes \f$\log_2v + 1\f$, where \f$v\f$ is of unsigned
 * 8-bit type.
 *
 * An unsigned integer \f$v\f$ can be represented in this many bits.
 * @param val value to find log2 of.
 * @return \f$\lfloor{}log_2v\rfloor + 1\f$
 */
static inline int
requiredUInt8Bits(uint8_t val);
/**
 * \brief Computes \f$\log_2v + 1\f$, where \f$v\f$ is of unsigned
 * 16-bit type.
 *
 * An unsigned integer \f$v\f$ can be represented in this many bits.
 * @param val value to find log2 of.
 * @return \f$\lfloor{}log_2v\rfloor + 1\f$
 */
static inline int
requiredUInt16Bits(uint16_t val);
/**
 * \brief Computes \f$\log_2v + 1\f$, where \f$v\f$ is of unsigned
 * 32-bit type.
 *
 * An unsigned integer \f$v\f$ can be represented in this many bits.
 * @param val value to find log2 of.
 * @return \f$\lfloor{}log_2v\rfloor + 1\f$
 */
extern int
requiredUInt32Bits(uint32_t val);
/**
 * \brief Computes \f$\log_2v + 1\f$, where \f$v\f$ is of unsigned
 * 64-bit type.
 *
 * An unsigned integer \f$v\f$ can be represented in this many bits.
 * @param val value to find log2 of.
 * @return \f$\lfloor{}log_2v\rfloor + 1\f$
 */
extern int
requiredUInt64Bits(uint64_t val);

/**
 * \brief Computes \f$\log_2v + 2\f$, where \f$v\f$ is of signed 8-bit
 * type.
 *
 * An integer \f$v\f$ can be represented in this many bits.
 * @param val value to find log2 of.
 * @return \f$\lfloor{}|log_2v|\rfloor + 2\f$. Special case: INT8_MIN
 * is representable in 8 bits although
 * \f[\left\lfloor{}log_2\left|-(2^{7})\right|\right\rfloor + 2 = 9.\f]
 */
static inline int
requiredInt8Bits(int8_t v);
/**
 * \brief Computes \f$\log_2v + 2\f$, where \f$v\f$ is of signed
 * 16-bit type.
 *
 * An integer \f$v\f$ can be represented in this many bits.
 * @param val value to find log2 of.
 * @return \f$\lfloor{}|log_2v|\rfloor + 2\f$. Special case: INT16_MIN
 * is representable in 16 bits although
 * \f[\left\lfloor{}log_2\left|-(2^{15})\right|\right\rfloor + 2 = 17.\f]
 */
static inline int
requiredInt16Bits(int16_t v);
/**
 * \brief Computes \f$\log_2v + 2\f$, where \f$v\f$ is of signed
 * 32-bit type.
 *
 * An integer \f$v\f$ can be represented in this many bits.
 * @param val value to find log2 of.
 * @return \f$\lfloor{}|log_2v|\rfloor + 2\f$. Special case: INT32_MIN
 * is representable in 32 bits although
 * \f[\left\lfloor{}log_2\left|-(2^{31})\right|\right\rfloor + 2 = 33.\f]
 */
static inline int
requiredInt32Bits(int32_t val);
/**
 * \brief Computes \f$\log_2v + 2\f$, where \f$v\f$ is of signed
 * 64-bit type.
 *
 * An integer \f$v\f$ can be represented in this many bits.
 * @param val value to find log2 of.
 * @return \f$\left\lfloor{}log_2|v|\right\rfloor + 2\f$
 * (Special case: INT64_MIN is representable in 64 bits)
 */
static inline int
requiredInt64Bits(int64_t val);
/**
 * \brief Retrieve unsigned integer of specified length from
 * bitstring at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
extern uint8_t
bsGetUInt8(const BitString str, BitOffset offset, unsigned numBits);
/**
 * \brief Retrieve unsigned integer of specified length from bitstring
 * at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
extern uint16_t
bsGetUInt16(const BitString str, BitOffset offset, unsigned numBits);
/**
 * \brief Retrieve unsigned integer of specified length from bitstring
 * at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
extern uint32_t
bsGetUInt32(const BitString str, BitOffset offset, unsigned numBits);
/**
 * \brief Retrieve unsigned integer of specified length from bitstring
 * at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
extern uint64_t
bsGetUInt64(const BitString str, BitOffset offset, unsigned numBits);
/**
 * \brief Store unsigned integer of specified length in bitstring at
 * given position.
 *
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 */
extern void
bsStoreUInt8(BitString str, BitOffset offset, unsigned numBits, uint8_t val);
/**
 * \brief Store unsigned integer of specified length in bitstring at
 * given position.
 *
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 */
extern void
bsStoreUInt16(BitString str, BitOffset offset, unsigned numBits, uint16_t val);
/**
 * \brief Store unsigned integer of specified length in bitstring at
 * given position.
 *
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 */
extern void
bsStoreUInt32(BitString str, BitOffset offset, unsigned numBits, uint32_t val);
/**
 * \brief Store unsigned integer of specified length in bitstring at
 * given position.
 *
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 */
extern void
bsStoreUInt64(BitString str, BitOffset offset, unsigned numBits, uint64_t val);
/**
 * Store integer of specified length in bitstring at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 */
static inline int8_t
bsGetInt8(BitString str, BitOffset offset, unsigned numBits);
/**
 * \brief Store integer of specified length in bitstring at given
 * position.
 *
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 */
static inline int16_t
bsGetInt16(BitString str, BitOffset offset, unsigned numBits);
/**
 * \brief Store integer of specified length in bitstring at given
 * position.
 *
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 */
static inline int32_t
bsGetInt32(BitString str, BitOffset offset, unsigned numBits);
/**
 * \brief Store integer of specified length in bitstring at given
 * position.
 *
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing integer to be written
 */
static inline int64_t
bsGetInt64(BitString str, BitOffset offset, unsigned numBits);
/**
 * Retrieve integer of specified length from bitstring at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
static inline void
bsStoreInt8(BitString str, BitOffset offset, unsigned numBits, uint8_t val);
/**
 * \brief Retrieve integer of specified length from bitstring at given
 * position.
 *
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
static inline void
bsStoreInt16(BitString str, BitOffset offset, unsigned numBits, uint16_t val);
/**
 * \brief Retrieve integer of specified length from bitstring at given
 * position.
 *
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
static inline void
bsStoreInt32(BitString str, BitOffset offset, unsigned numBits, uint32_t val);
/**
 * \brief Retrieve integer of specified length from bitstring at given
 * position.
 *
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @return value read
 */
static inline void
bsStoreInt64(BitString str, BitOffset offset, unsigned numBits, uint64_t val);
/* higher level functions */
/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param read integers from this array
 */
void
bsStoreUniformUInt8Array(BitString str, BitOffset offset, unsigned numBits,
                         size_t numValues, const uint8_t val[]);
/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param read integers from this array
 */
void
bsStoreUniformUInt16Array(BitString str, BitOffset offset, unsigned numBits,
                          size_t numValues, const uint16_t val[]);
/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param read integers from this array
 */
void
bsStoreUniformUInt32Array(BitString str, BitOffset offset, unsigned numBits,
                          size_t numValues, const uint32_t val[]);
/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param read integers from this array
 */
void
bsStoreUniformUInt64Array(BitString str, BitOffset offset, unsigned numBits,
                          size_t numValues, const uint64_t val[]);
/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param read integers from this array
 */
static inline void
bsStoreUniformInt8Array(BitString str, BitOffset offset, unsigned numBits,
                        size_t numValues, const int8_t val[]);
/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param read integers from this array
 */
static inline void
bsStoreUniformInt16Array(BitString str, BitOffset offset, unsigned numBits,
                         size_t numValues, const int16_t val[]);
/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param read integers from this array
 */
static inline void
bsStoreUniformInt32Array(BitString str, BitOffset offset, unsigned numBits,
                         size_t numValues, const int32_t val[]);
/**
 * \brief Store n unsigned integers of specified length from array, in
 * bitstring, starting at given position.
 * @param str bitstring to write to
 * @param offset position to start writing at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param read integers from this array
 */
static inline void
bsStoreUniformInt64Array(BitString str, BitOffset offset, unsigned numBits,
                         size_t numValues, const int64_t val[]);
/**
 * \brief Retrieve n unsigned integers of specified length from
 * bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param store integers read in this array
 */
void
bsGetUniformUInt8Array(const BitString str, BitOffset offset, unsigned numBits,
                       size_t numValues, uint8_t val[]);
/**
 * \brief Retrieve n unsigned integers of specified length from
 * bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param store integers read in this array
 */
void
bsGetUniformUInt16Array(const BitString str, BitOffset offset, unsigned numBits,
                        size_t numValues, uint16_t val[]);
/**
 * \brief Retrieve n unsigned integers of specified length from
 * bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param store integers read in this array
 */
void
bsGetUniformUInt32Array(const BitString str, BitOffset offset, unsigned numBits,
                        size_t numValues, uint32_t val[]);
/**
 * \brief Retrieve n unsigned integers of specified length from
 * bitstring, starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing each integer
 * @param numValues number of integers to read
 * @param store integers read in this array
 */
void
bsGetUniformUInt64Array(const BitString str, BitOffset offset, unsigned numBits,
                        size_t numValues, uint64_t val[]);
/**
 * \brief Retrieve n integers of specified length from bitstring,
 * starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @param numValues number of integers to read
 * @param store integers read in this array
 */
static inline void
bsGetUniformInt8Array(const BitString str, BitOffset offset, unsigned numBits,
                      size_t numValues, int8_t val[]);
/**
 * \brief Retrieve n integers of specified length from bitstring,
 * starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @param numValues number of integers to read
 * @param store integers read in this array
 */
static inline void
bsGetUniformInt16Array(const BitString str, BitOffset offset, unsigned numBits,
                       size_t numValues, int16_t val[]);
/**
 * \brief Retrieve n integers of specified length from bitstring,
 * starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @param numValues number of integers to read
 * @param store integers read in this array
 */
static inline void
bsGetUniformInt32Array(const BitString str, BitOffset offset, unsigned numBits,
                       size_t numValues, int32_t val[]);
/**
 * \brief Retrieve n integers of specified length from bitstring,
 * starting at given position.
 * @param str bitstring to read from
 * @param offset position to start reading at (bit exact)
 * @param numBits number of bits composing integer to be read
 * @param numValues number of integers to read
 * @param store integers read in this array
 */
static inline void
bsGetUniformInt64Array(const BitString str, BitOffset offset, unsigned numBits,
                       size_t numValues, int64_t val[]);

extern int
bsCompare(const BitString a, BitOffset offsetA, BitOffset numBitsA,
          const BitString b, BitOffset offsetB, BitOffset numBitsB);

/**
 * \brief copy (sub-)bitstring to another position (in same or other
 * bitstring)
 * @param src bitsring to copy from
 * @param offsetSrc bit position in src to start copying from
 * @param dest bitstring to copy to
 * @param offsetDest bit position in dest to start writing at
 * @param numBits number of bits to copy
 */
extern void
bsCopy(const BitString src, BitOffset offsetSrc, 
       const BitString dest, BitOffset offsetDest, BitOffset numBits);
/**
 * \brief Meta-Unit test function for bitPackString, calls all functions
 * mentioned below.
 * @return 0 on success, -1 on error.
 */
extern int
bitPackString_unit_test(Env *env);
/**
 * \brief Unit test function for bitPackString, integer functions.
 * @return 0 on success, -1 on error.
 */
extern int
bitPackStringInt_unit_test(Env *env);
/**
 * \brief Unit test function for bitPackString, unsigned functions.
 * @return 0 on success, -1 on error.
 */
extern int
bitPackStringUInt_unit_test(Env *env);
/**
 * \brief Unit test function for bitPackString, 8-bit functions.
 * @return 0 on success, -1 on error.
 */
extern int
bitPackStringInt8_unit_test(Env *env);
/**
 * \brief Unit test function for bitPackString, 16-bit functions.
 * @return 0 on success, -1 on error.
 */
extern int
bitPackStringInt16_unit_test(Env *env);
/**
 * \brief Unit test function for bitPackString, 32-bit functions.
 * @return 0 on success, -1 on error.
 */
extern int
bitPackStringInt32_unit_test(Env *env);
/**
 * \brief Unit test function for bitPackString, 64-bit functions.
 * @return 0 on success, -1 on error.
 */
extern int
bitPackStringInt64_unit_test(Env *env);

#include "bitpackstringsimpleop.h"

#endif /* BITPACKSTRING_H_INCLUDED */
