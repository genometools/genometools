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
#ifndef EIS_BWTCONSTRUCT_PARAMS_H
#define EIS_BWTCONSTRUCT_PARAMS_H

/**
 * @file eis-bwtconstruct_params.h
 * @brief Routines to handle option processing for BWT index construction.
 */

#include "libgtcore/option.h"
#include "libgtmatch/eis-bwtseqparam.h"

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
 * @param err
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
 * @param err
 */
extern void
computePackedIndexDefaults(struct bwtOptions *paramOutput, int extraToggles,
                           Error *err);

#endif
