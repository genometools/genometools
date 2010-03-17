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

#include "core/seqpos.h"
#include "core/bitpackstring.h"

#ifdef Seqposequalsunsignedint
/** retrieve Seqpos from BitString */
#define gt_bsGetSeqpos gt_bsGetUInt32
/** store Seqpos in BitString */
#define gt_bsStoreSeqpos gt_bsStoreUInt32
/** read back array of Seqpos values in BitString */
#define gt_bsGetUniformSeqposArray gt_bsGetUniformUInt32Array
/** read back array of Seqpos values in BitString */
#define gt_bsGetUniformSeqposArrayAdd gt_bsGetUniformUInt32ArrayAdd
/** read back array of Seqpos values in BitString */
#define gt_bsGetNonUniformSeqposArray gt_bsGetNonUniformUInt32Array
/** read back array of Seqpos values in BitString */
#define gt_bsGetNonUniformSeqposArrayAdd gt_bsGetNonUniformUInt32ArrayAdd
/** store array of Seqpos values in BitString */
#define gt_bsStoreUniformSeqposArray gt_bsStoreUniformUInt32Array
/** how many bits are required to store given Seqpos value */
#define requiredSeqposBits gt_requiredUInt32Bits
#else
/** retrieve Seqpos from BitString */
#define gt_bsGetSeqpos gt_bsGetUInt64
/** store Seqpos in BitString */
#define gt_bsStoreSeqpos gt_bsStoreUInt64
/** read back array of Seqpos values in BitString */
#define gt_bsGetUniformSeqposArray gt_bsGetUniformUInt64Array
/** read back array of Seqpos values in BitString */
#define gt_bsGetUniformSeqposArrayAdd gt_bsGetUniformUInt64ArrayAdd
/** read back array of Seqpos values in BitString */
#define gt_bsGetNonUniformSeqposArray gt_bsGetNonUniformUInt64Array
/** read back array of Seqpos values in BitString */
#define gt_bsGetNonUniformSeqposArrayAdd gt_bsGetNonUniformUInt64ArrayAdd
/** store array of Seqpos values in BitString  */
#define gt_bsStoreUniformSeqposArray gt_bsStoreUniformUInt64Array
/** how many bits are required to store given Seqpos value */
#define requiredSeqposBits gt_requiredUInt64Bits
#endif

#endif
