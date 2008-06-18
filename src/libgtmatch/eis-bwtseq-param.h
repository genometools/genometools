/*
  Copyright (c) 2007,2008 Thomas Jahns <Thomas.Jahns@gmx.net>

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
#ifndef EIS_BWTSEQ_PARAM_H
#define EIS_BWTSEQ_PARAM_H

/**
 * @file eis-bwtseq-param.h
 * @brief Routines to handle option processing for BWT index
 * construction and declaration of parameter types and corresponding
 * constant values.
 */

#include "libgtcore/option.h"

#include "libgtmatch/eis-encidxseq-param.h"
/* #include "libgtmatch/eis-bwtseq-context-param.h" not neccessary remove
 * Stefan */

/**
 * all parameters for building the BWT sequence index
 */
struct bwtParam
{
  struct seqBaseParam seqParams;  /**< holds extra parameter
                                   *   information specific to the
                                   *   base sequence storage */
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
  int sourceRankInterval;         /**< makes ranges which are sorted in
                                   * rank mode reversibly sorted:
                                   * -1:
                                   *     inactive
                                   * 0...sizeof(Seqpos)*CHAR_BIT:
                                   *     build accel table
                                   *     with 1<<SourceRankInterval
                                   *     sampling interval
                                   */
  int ctxMapILog;                 /**< according to the value of
                                   * ctxMapIlog:
                                   * 0..sizeof(Seqpos)*CHAR_BIT:
                                   *   a context map is produced with
                                   *   interval 1 << ctxMapILog
                                   * -1: use log(seqlen)
                                   * -2: inactive,
                                   * see enum ctxMapSize
                                   */
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
  BWTReversiblySorted  = 1 << 2, /**< unless special symbols are fully sorted
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

/**
 * Describes features to tune the constructed BWT sequence index
 * object for different applications.
 */
enum BWTOptionDefaultsOptimizationFlags
{
  BWTDEFOPT_LOW_RAM_OVERHEAD = 1 << 0, /**< try to minimize ram used */
  BWTDEFOPT_FAST_RANK        = 1 << 1, /**< build tables to speed up
                                        * rank queries for special symbols */
  /** optimize memory usage for construction process (i.e. don't
   * build lookup-tables only needed to speed-up queries) */
  BWTDEFOPT_CONSTRUCTION     = BWTDEFOPT_LOW_RAM_OVERHEAD,
  /** optimize for multiple queries */
  BWTDEFOPT_MULTI_QUERY      = BWTDEFOPT_FAST_RANK,
};

/**
 * Contains both the arguments finally passed on index construction
 * and the raw flags and values set via the option parser.
 */
struct bwtOptions
{
  struct bwtParam final;                /**< final flags, passed to
                                         * index construction */
  int defaultOptimizationFlags;         /**< set by user to one or
                                         * more of the values of enum
                                         * BWTOptionDefaultsOptimizationFlags,
                                         * determines ram usage of index */
  bool useLocateBitmap;                 /**< did the user request to
                                         * store locate flags for
                                         * every sequence position? */
  bool useSourceRank;                   /**< did the user request extra
                                         * information for sort reversing of
                                         * rank-sorted symbols? */
  Option *useLocateBitmapOption;        /**< used to query wether the
                                         * option was set or left
                                         * unspecified in which case a
                                         * reasonable default is
                                         * computed */
};

/**
 * @brief Add the options for BWT seqence index construction to an
 * option parser.
 * @param op
 * @param paramOutput used to hold options set by the user
 * @param defaultOptimizationFlags used to tune structures for
 * construction or matching
 * @param projectName reference to string which will hold the base
 * name (i.e. without extension) of the project once index creation begins
 */
extern void
registerPackedIndexOptions(OptionParser *op, struct bwtOptions *paramOutput,
                           int defaultOptimizationFlags,
                           const Str *projectName);

/**
 * @brief Compute the options for BWT seqence index construction from
 * the values specified by command-line.
 * @param paramOutput used to hold options set by the user
 * @param extraToggles add flags here which only become apparent at
 * after option processing is finished
 */
extern void
computePackedIndexDefaults(struct bwtOptions *paramOutput, int extraToggles);

/**
 * @brief Computes feature set of base index of type EISeq from
 * optimization flags set for BWT index.
 * @param BWTOptFlags
 * @return set of enum EISFeatureBits toggles
 */
extern int
convertBWTOptFlags2EISFeatures(int BWTOptFlags);

#endif
