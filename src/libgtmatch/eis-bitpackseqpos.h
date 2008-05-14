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
#ifndef EIS_BITPACKSEQPOS_H
#define EIS_BITPACKSEQPOS_H

/**
 * @file eis-bitpackseqpos.h
 * @brief Defines appropriate width routines for storing and retrieving
 * Seqpos values in/from BitString.
 */

#include "libgtmatch/seqpos-def.h"
#include "libgtcore/bitpackstring.h"

#ifdef Seqposequalsunsignedint
/** retrieve Seqpos from BitString */
#define bsGetSeqpos bsGetUInt32
/** store Seqpos in BitString */
#define bsStoreSeqpos bsStoreUInt32
/** read back array of Seqpos values in BitString */
#define bsGetUniformSeqposArray bsGetUniformUInt32Array
/** read back array of Seqpos values in BitString */
#define bsGetUniformSeqposArrayAdd bsGetUniformUInt32ArrayAdd
/** read back array of Seqpos values in BitString */
#define bsGetNonUniformSeqposArray bsGetNonUniformUInt32Array
/** read back array of Seqpos values in BitString */
#define bsGetNonUniformSeqposArrayAdd bsGetNonUniformUInt32ArrayAdd
/** store array of Seqpos values in BitString */
#define bsStoreUniformSeqposArray bsStoreUniformUInt32Array
/** how many bits are required to store given Seqpos value */
#define requiredSeqposBits requiredUInt32Bits
#else
/** retrieve Seqpos from BitString */
#define bsGetSeqpos bsGetUInt64
/** store Seqpos in BitString */
#define bsStoreSeqpos bsStoreUInt64
/** read back array of Seqpos values in BitString */
#define bsGetUniformSeqposArray bsGetUniformUInt64Array
/** read back array of Seqpos values in BitString */
#define bsGetUniformSeqposArrayAdd bsGetUniformUInt64ArrayAdd
/** read back array of Seqpos values in BitString */
#define bsGetNonUniformSeqposArray bsGetNonUniformUInt64Array
/** read back array of Seqpos values in BitString */
#define bsGetNonUniformSeqposArrayAdd bsGetNonUniformUInt64ArrayAdd
/** store array of Seqpos values in BitString  */
#define bsStoreUniformSeqposArray bsStoreUniformUInt64Array
/** how many bits are required to store given Seqpos value */
#define requiredSeqposBits requiredUInt64Bits
#endif

#endif
