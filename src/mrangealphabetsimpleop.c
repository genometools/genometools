/*
** Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
**  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**  
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**  
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
**  
*/
#include "mrangealphabet.h"
#include "mrangealphabetpriv.h"

staticifinline inline MRAEnc *
MRAEncUInt8New(int numRanges, int symbolsPerRange[],
               const uint8_t *mapping, Env *env)
{
  return newMultiRangeAlphabetEncodingUInt8(numRanges, symbolsPerRange,
                                            mapping, env);
}

staticifinline inline size_t
MRAEncGetRangeSize(const MRAEnc *mralpha, size_t range)
{
  assert(mralpha);
  assert(mralpha->numRanges > range);
  return mralpha->symbolsPerRange[range];
}

staticifinline inline size_t
MRAEncGetDomainSize(const MRAEnc *mralpha)
{
  assert(mralpha);
  switch(mralpha->encType)
  {
  case sourceUInt8:
    return UINT8_MAX + 1;
  default:
    abort();
  }
}

staticifinline inline Symbol
MRAEncMapSymbol(const MRAEnc *mralpha, Symbol sym)
{
  switch(mralpha->encType)
  {
  case sourceUInt8:
    return constMRAEnc2MRAEncUInt8(mralpha)->mappings[(uint8_t)sym];
  default:
    abort();
  }  
}

staticifinline inline Symbol
MRAEncRevMapSymbol(const MRAEnc *mralpha, Symbol sym)
{
  switch(mralpha->encType)
  {
  case sourceUInt8:
    return constMRAEnc2MRAEncUInt8(mralpha)->revMappings[(uint8_t)sym];
  default:
    abort();
  }
}

staticifinline inline int
MRAEncSymbolHasValidMapping(const MRAEnc *mralpha, Symbol sym)
{
  switch(mralpha->encType)
  {
  case sourceUInt8:
    return mralpha->rangeEndIndices[mralpha->numRanges] == UINT8_MAX
      || (constMRAEnc2MRAEncUInt8(mralpha)->mappings[(uint8_t)sym]
          != UNDEF_UCHAR);
  default:
    abort();
  }  
}
