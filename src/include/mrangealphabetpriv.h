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
#ifndef MRANGEALPHABETPRIV_H_INCLUDED
#define MRANGEALPHABETPRIV_H_INCLUDED

#include <limits.h>
#include <stddef.h>
#include <stdint.h>

#include "mrangealphabet.h"

enum sourceEncType {
  sourceUInt8,
};

#define UNDEF_UCHAR ((unsigned char)~0)

struct multiRangeAlphabetEncoding
{
  enum sourceEncType encType;
  size_t numRanges;
  size_t *rangeEndIndices, /*< maps to the last position + 1 for each range */
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
  assert(mralpha && mralpha->encType == sourceUInt8);
  return (MRAEncUInt8*)
    ((char *)mralpha - offsetof(MRAEncUInt8, baseClass));
}

static inline const MRAEncUInt8 *
constMRAEnc2MRAEncUInt8(const MRAEnc *mralpha)
{
  assert(mralpha && mralpha->encType == sourceUInt8);
  return (const MRAEncUInt8*)
    ((char *)mralpha - offsetof(MRAEncUInt8, baseClass));
}

static inline MRAEnc*
MRAEncUInt82MRAEnc(MRAEncUInt8 *mralpha)
{
  assert(mralpha);
  return &(mralpha->baseClass);
}

#endif /* MRANGEALPHABETPRIV_H_INCLUDED */
