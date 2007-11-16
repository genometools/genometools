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
#ifndef EIS_MRANGEALPHABETSIMPLEOP_H
#define EIS_MRANGEALPHABETSIMPLEOP_H

#include "libgtmatch/eis-mrangealphabet.h"
#include "libgtmatch/eis-mrangealphabetpriv.h"

static inline MRAEnc *
MRAEncUInt8New(int numRanges, int symbolsPerRange[],
               const uint8_t *mapping, Env *env)
{
  return newMultiRangeAlphabetEncodingUInt8(numRanges, symbolsPerRange,
                                            mapping, env);
}

static inline size_t
MRAEncGetRangeSize(const MRAEnc *mralpha, size_t range)
{
  assert(mralpha);
  assert(mralpha->numRanges > range);
  return mralpha->symbolsPerRange[range];
}

static inline size_t
MRAEncGetDomainSize(const MRAEnc *mralpha)
{
  assert(mralpha);
  switch (mralpha->encType)
  {
  case sourceUInt8:
    return UINT8_MAX + 1;
  default:
    abort();
  }
}

static inline Symbol
MRAEncMapSymbol(const MRAEnc *mralpha, Symbol sym)
{
  switch (mralpha->encType)
  {
  case sourceUInt8:
    return constMRAEnc2MRAEncUInt8(mralpha)->mappings[(uint8_t)sym];
  default:
    abort();
  }
}

static inline Symbol
MRAEncRevMapSymbol(const MRAEnc *mralpha, Symbol sym)
{
  switch (mralpha->encType)
  {
  case sourceUInt8:
    return constMRAEnc2MRAEncUInt8(mralpha)->revMappings[(uint8_t)sym];
  default:
    abort();
  }
}

static inline int
MRAEncSymbolHasValidMapping(const MRAEnc *mralpha, Symbol sym)
{
  switch (mralpha->encType)
  {
  case sourceUInt8:
    return mralpha->rangeEndIndices[mralpha->numRanges - 1] == UINT8_MAX
      || (constMRAEnc2MRAEncUInt8(mralpha)->mappings[(uint8_t)sym]
          != UNDEF_UCHAR);
  default:
    abort();
  }
}

#endif
