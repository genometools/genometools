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

#ifndef EIS_MRANGEALPHABET_PRIV_H
#define EIS_MRANGEALPHABET_PRIV_H

#include <limits.h>
#include <stddef.h>
#include <inttypes.h>

#include "match/eis-mrangealphabet.h"

struct multiRangeAlphabetEncoding
{
  enum sourceEncType encType;
  AlphabetRangeID numRanges;
  AlphabetRangeSize *rangeEndIndices, /*< maps to the last position + 1
                                       * for each range */
    *symbolsPerRange;   /*< gives number of symbols in each range */
};

struct multiRangeAlphabetEncodingUInt8
{
  struct multiRangeAlphabetEncoding baseClass;
  uint8_t mappings[UINT8_MAX+1], revMappings[UINT8_MAX+1];
};

typedef struct multiRangeAlphabetEncodingUInt8 MRAEncUInt8;

static inline MRAEncUInt8 *
MRAEnc2MRAEncUInt8(MRAEnc *mralpha)
{
  gt_assert(mralpha && mralpha->encType == sourceUInt8);
  return (MRAEncUInt8*)
    ((char *)mralpha - offsetof(MRAEncUInt8, baseClass));
}

static inline const MRAEncUInt8 *
constMRAEnc2MRAEncUInt8(const MRAEnc *mralpha)
{
  gt_assert(mralpha && mralpha->encType == sourceUInt8);
  return (const MRAEncUInt8*)
    ((char *)mralpha - offsetof(MRAEncUInt8, baseClass));
}

static inline MRAEnc*
MRAEncUInt82MRAEnc(MRAEncUInt8 *mralpha)
{
  gt_assert(mralpha);
  return &(mralpha->baseClass);
}

#endif
