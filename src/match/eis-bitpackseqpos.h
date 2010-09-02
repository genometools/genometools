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
 * unsigned long values in/from BitString.
 */

#include "core/bitpackstring.h"

#ifndef _LP64
/** retrieve unsigned long from BitString */
#define gt_bsGetUlong gt_bsGetUInt32
/** store unsigned long in BitString */
#define gt_bsStoreUlong gt_bsStoreUInt32
/** read back array of unsigned long values in BitString */
#define gt_bsGetUniformUlongArray gt_bsGetUniformUInt32Array
/** read back array of unsigned long values in BitString */
#define gt_bsGetUniformUlongArrayAdd gt_bsGetUniformUInt32ArrayAdd
/** read back array of unsigned long values in BitString */
#define gt_bsGetNonUniformUlongArray gt_bsGetNonUniformUInt32Array
/** read back array of unsigned long values in BitString */
#define gt_bsGetNonUniformUlongArrayAdd gt_bsGetNonUniformUInt32ArrayAdd
/** store array of unsigned long values in BitString */
#define gt_bsStoreUniformUlongArray gt_bsStoreUniformUInt32Array
/** how many bits are required to store given unsigned long value */
#define requiredUlongBits gt_requiredUInt32Bits
#else
/** retrieve unsigned long from BitString */
#define gt_bsGetUlong gt_bsGetUInt64
/** store unsigned long in BitString */
#define gt_bsStoreUlong gt_bsStoreUInt64
/** read back array of unsigned long values in BitString */
#define gt_bsGetUniformUlongArray gt_bsGetUniformUInt64Array
/** read back array of unsigned long values in BitString */
#define gt_bsGetUniformUlongArrayAdd gt_bsGetUniformUInt64ArrayAdd
/** read back array of unsigned long values in BitString */
#define gt_bsGetNonUniformUlongArray gt_bsGetNonUniformUInt64Array
/** read back array of unsigned long values in BitString */
#define gt_bsGetNonUniformUlongArrayAdd gt_bsGetNonUniformUInt64ArrayAdd
/** store array of unsigned long values in BitString  */
#define gt_bsStoreUniformUlongArray gt_bsStoreUniformUInt64Array
/** how many bits are required to store given unsigned long value */
#define requiredUlongBits gt_requiredUInt64Bits
#endif

#endif
