/*
  Copyright (c) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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

#ifndef EIS_BWTSEQPARAM_H
#define EIS_BWTSEQPARAM_H

/**
 * @file eis-bwtseqparam.h
 * @brief Construction parameters for BWT index creation
 * routines.
 */

#include "libgtmatch/eis-encidxseqparam.h"

/**
 * all parameters for building the BWT sequence index
 */
struct bwtParam
{
  union seqBaseEncParam seqParams;/**< a union holding extra parameter
                                   *   information specific to the
                                   *   type selected via parameter
                                   *   baseType */
  enum seqBaseEncoding baseType;  /**< baseType selects the encoding
                                   *   method of the sequence index *
                                   *   storing the BWT sequence (see
                                   *   enum seqBaseEncoding). */
  int featureToggles;             /**< set of bittoggles composed
                                   *   from enum BWTFeatures */
  unsigned locateInterval;        /**< store locate information
                                   * (mapping from BWT sequence to
                                   * original sequence for every nth
                                   * position in original sequence) at
                                   * this frequency, unless
                                   * locateInterval = 0, which implies
                                   * storing no locate information at
                                   * all */
  const Str *projectName;         /**< base file name to derive name
                                   *   of suffixerator project from*/
};

/**
 * Select specifics of how
 * - the locate information is stored
 */
enum BWTFeatures
{
  BWTBaseFeatures      =      0,
  BWTLocateBitmap      = 1 << 0,
  BWTLocateCount       = 1 << 1,
  BWTProperlySorted    = 1 << 2, /**< unless special symbols are fully sorted
                                  * by suffixerator, there are two
                                  * restrictions for the use of the index:
                                  * - the bwt property is only true
                                  *   for subsequences without special
                                  *   symbols of the original sequence
                                  * - the regeneration of the original
                                  *   sequence is impossible and a
                                  *   reverse establishment of context
                                  *   impossible.
                                  */
};

#endif
