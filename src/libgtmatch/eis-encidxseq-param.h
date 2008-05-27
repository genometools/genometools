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

#ifndef EIS_ENCIDXSEQ_PARAM_H
#define EIS_ENCIDXSEQ_PARAM_H

#include "libgtcore/option.h"

/**
 * @file eis-encidxseq-param.h
 * @brief Parameter definitions for index creation routines.
 */

/**
 * Select features to use for index in-memory representation.
 */
enum EISFeatureBits
{
  EIS_FEATURE_NONE = 0,           /**< fallback value */
  EIS_FEATURE_REGION_SUMS = 1<<0, /**< if set construct sum tables for
                                   *   the special symbol ranges
                                   *   use this on index loading if
                                   *   you want to use this index in
                                   *   many queries, omit if memory
                                   *   is very tight (e.g. on construction) */
};

/**
 * Parameters to construct a block composition compressed index.
 */
struct blockEncParams
{
  unsigned blockSize;         /**< number of symbols to combine in
                               * one block a lookup-table
                               * containing
                               * $alphabetsize^{blockSize}$ entries is
                               * required so adjust with caution */
  unsigned bucketBlocks;      /**< number of blocks for which to
                               * store partial symbol sums (lower
                               * values increase index size and
                               * decrease computations for lookup) */
};

/**
 * Names the type of encoding used:
 */
enum seqBaseEncoding {
  BWT_BASETYPE_AUTOSELECT,      /**< automatic, load any index present
                                 *     (currently not implemented) */
  BWT_ON_RLE,                   /**< use original fmindex run-length
                                 * encoding  */
  BWT_ON_BLOCK_ENC,             /**< do block compression by dividing
                                 * sequence into strings of
                                 * composition and permutation
                                 * indices */
  BWT_ON_WAVELET_TREE_ENC,      /**< encode sequence with wavelet-trees */
};

/**
 * Stores information to construct the underlying sequence object.
 */
union seqBaseEncParam
{
  struct blockEncParams blockEnc; /**< parameters for block encoding  */
/*   struct  */
/*   { */
/*   } RLEParams; */
/*   struct  */
/*   { */
/*   } waveletTreeParams; */
};

struct seqBaseParam
{
  enum seqBaseEncoding encType;  /**< encType selects the encoding
                                   *   method of the sequence index *
                                   *   storing the BWT sequence (see
                                   *   enum seqBaseEncoding). */
  int EISFeatureSet;              /**< bitwise or of EIS_FEATURE_NONE
                                   * and other features selectable via enum
                                   * EISFeatureBits (see eis-encidxseq.h) */
  union seqBaseEncParam encParams;/**< a union holding extra parameter
                                   *   information specific to the
                                   *   type selected via parameter
                                   *   encType */
};

extern void
registerEncIdxSeqOptions(OptionParser *op, struct seqBaseParam *paramOutput);

/**
 * @brief Given the construction parameters for a sequence index,
 * estimate how many symbols will be encoded together.
 *
 * Extension information can be expected to be stored at a intervals
 * corresponding to the return value of this function.
 * @param params paramters for sequence index construction
 * @param encType select which type of index will be constructed
 * @param err
 * @return number of symbols stored consecutively
 */
extern unsigned
estimateSegmentSize(const struct seqBaseParam *params);

#endif
