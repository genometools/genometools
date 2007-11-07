/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef SFX_OPTDEF_H
#define SFX_OPTDEF_H

#include <stdbool.h>
#include "libgtcore/str.h"
#include "libgtcore/strarray.h"
#include "readmode-def.h"

#define PREFIXLENGTH_AUTOMATIC 0

typedef struct
{
  unsigned int numofparts,
               prefixlength;
  Str *str_indexname,
      *str_smap,
      *str_sat;
  StrArray *filenametab;
  Readmode readmode;
  bool isdna,
       isprotein,
       isplain,
       beverbose,
       outtistab,
       outsuftab,
       outlcptab,
       outbwttab,
       outdestab;
  unsigned int blockSize;      /**< number of symbols to combine in
                                 * one block, a lookup-table
                                 * containing
                                 * $alphabetsize^{blockSize}$ entries is
                                 * required so adjust with caution */
  unsigned int bucketBlocks;   /**< number of blocks for which to
                                 * store partial symbol sums (lower
                                 * values increase index size and
                                 * decrease computations for lookup) */
  unsigned int locateInterval;  /**< sample suffix array values at
                                 * this interval, or not at all if 0 */
} Suffixeratoroptions;

#endif
