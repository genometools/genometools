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

#include <assert.h>
#include <limits.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>

#include "libgtcore/chardef.h"
#include "libgtcore/dataalign.h"
#include "libgtcore/ma.h"
#include "libgtcore/str.h"
#include "libgtcore/symboldef.h"
#include "libgtcore/strarray.h"
#include "libgtmatch/alphadef.h"

#include "libgtmatch/eis-mrangealphabet.h"
#include "libgtmatch/eis-mrangealphabet-priv.h"

MRAEnc *
newMultiRangeAlphabetEncodingUInt8(AlphabetRangeID numRanges,
                                   const AlphabetRangeSize symbolsPerRange[],
                                   const uint8_t *mappings)
{
  MRAEncUInt8 *newAlpha = NULL;
  size_t rEIOffset = offsetAlign(
    sizeof (MRAEncUInt8),
    sizeof (newAlpha->baseClass.rangeEndIndices[0]));
  size_t sPROffset = offsetAlign(
    rEIOffset  +  sizeof (newAlpha->baseClass.rangeEndIndices[0]) * numRanges,
    sizeof (newAlpha->baseClass.symbolsPerRange[0]));
  AlphabetRangeID i;
  if ((newAlpha = ma_calloc(
         sPROffset
         + sizeof (newAlpha->baseClass.symbolsPerRange[0]) * numRanges,
         1)))
  {
    newAlpha->baseClass.rangeEndIndices
      = (AlphabetRangeSize *)((char *)newAlpha + rEIOffset);
    newAlpha->baseClass.symbolsPerRange
      = (AlphabetRangeSize *)((char *)newAlpha + sPROffset);
    newAlpha->baseClass.encType = sourceUInt8;
    newAlpha->baseClass.numRanges = numRanges;
    memset(newAlpha->mappings, UNDEF_UCHAR, UINT8_MAX+1);
    memset(newAlpha->revMappings, UNDEF_UCHAR, UINT8_MAX+1);
    newAlpha->baseClass.rangeEndIndices[0] =
      newAlpha->baseClass.symbolsPerRange[0] = symbolsPerRange[0];
    for (i = 1; i < numRanges; ++i)
    {
      newAlpha->baseClass.rangeEndIndices[i] =
        newAlpha->baseClass.rangeEndIndices[i-1]
        + (newAlpha->baseClass.symbolsPerRange[i] = symbolsPerRange[i]);
    }
    for (i = 0; i <= UINT8_MAX; ++i)
    {
      newAlpha->mappings[i] = mappings[i];
      newAlpha->revMappings[mappings[i]] = i;
    }
  }
  else
  {
    if (newAlpha)
    {
      if (newAlpha->baseClass.symbolsPerRange)
        ma_free(newAlpha->baseClass.symbolsPerRange);
      if (newAlpha->baseClass.rangeEndIndices)
        ma_free(newAlpha->baseClass.rangeEndIndices);
      ma_free(newAlpha);
    }
    return NULL;
  }
  return &(newAlpha->baseClass);
}

MRAEnc *
MRAEncGTAlphaNew(const Alphabet *alpha)
{
  AlphabetRangeSize symsPerRange[2];
  uint8_t *mappings;
  MRAEnc *result;
  uint32_t numSyms = getmapsizeAlphabet(alpha);
  mappings = ma_malloc(sizeof (uint8_t) * (UINT8_MAX + 1));
  memset(mappings, UNDEF_UCHAR, UINT8_MAX+1);
  {
    /* handle regular symbols */
    int i;
    for (i = 0; i < numSyms - 1; ++i)
      mappings[i] = i;
    symsPerRange[0] = numSyms - 1;
  }
  /* handle special symbols */
  mappings[WILDCARD] = numSyms - 1;
  symsPerRange[1] = 1;
  result = newMultiRangeAlphabetEncodingUInt8(2, symsPerRange, mappings);
  ma_free(mappings);
  return result;
}

extern MRAEnc *
MRAEncCopy(const MRAEnc *alpha)
{
  assert(alpha);
  switch (alpha->encType)
  {
  case sourceUInt8:
    {
      MRAEncUInt8 *newAlpha = NULL;
      const MRAEncUInt8 *srcAlpha = constMRAEnc2MRAEncUInt8(alpha);
      int numRanges = alpha->numRanges;
      assert(numRanges > 0);
      if ((newAlpha = ma_calloc(sizeof (MRAEncUInt8), 1))
          && (newAlpha->baseClass.rangeEndIndices =
              ma_malloc(sizeof (newAlpha->baseClass.rangeEndIndices[0])
                        * numRanges))
          && (newAlpha->baseClass.symbolsPerRange =
              ma_malloc(sizeof (newAlpha->baseClass.rangeEndIndices[0])
                        * numRanges)))
      {
        newAlpha->baseClass.encType = sourceUInt8;
        newAlpha->baseClass.numRanges = srcAlpha->baseClass.numRanges;
        memcpy(newAlpha->mappings, srcAlpha->mappings, UINT8_MAX+1);
        memcpy(newAlpha->revMappings, srcAlpha->revMappings, UINT8_MAX+1);
        memcpy(newAlpha->baseClass.rangeEndIndices,
               srcAlpha->baseClass.rangeEndIndices,
               sizeof (newAlpha->baseClass.rangeEndIndices[0]) * numRanges);
        memcpy(newAlpha->baseClass.symbolsPerRange,
               srcAlpha->baseClass.symbolsPerRange,
               sizeof (newAlpha->baseClass.symbolsPerRange[0]) * numRanges);
        return &(newAlpha->baseClass);
      }
      else
      {
        if (newAlpha)
        {
          if (newAlpha->baseClass.symbolsPerRange)
            ma_free(newAlpha->baseClass.symbolsPerRange);
          if (newAlpha->baseClass.rangeEndIndices)
            ma_free(newAlpha->baseClass.rangeEndIndices);
          ma_free(newAlpha);
        }
        return NULL;
      }
    }
    break;
  default:
    return NULL;
    break;
  }
}

extern AlphabetRangeSize
MRAEncGetSize(const MRAEnc *mralpha)
{
  AlphabetRangeID range, numRanges = mralpha->numRanges, sumRanges = 0;
  for (range = 0; range < numRanges; ++range)
  {
    sumRanges += mralpha->symbolsPerRange[range];
  }
  return sumRanges;
}

extern MRAEnc *
MRAEncSecondaryMapping(const MRAEnc *srcAlpha, int selection,
                       const int *rangeSel, Symbol fallback)
{
  MRAEnc *newAlpha;
  switch (srcAlpha->encType)
  {
  case sourceUInt8:
    {
      const MRAEncUInt8 *ui8alpha;
      uint8_t *mappings, destSym;
      AlphabetRangeSize *newRanges, sym;
      AlphabetRangeID range, numRanges = MRAEncGetNumRanges(srcAlpha);
      ui8alpha = constMRAEnc2MRAEncUInt8(srcAlpha);
      mappings = ma_malloc(sizeof (uint8_t) * (UINT8_MAX + 1));
      memset(mappings, UNDEF_UCHAR, UINT8_MAX+1);
      newRanges = ma_malloc(sizeof (newRanges[0]) * numRanges);
      sym = 0;
      destSym = 0;
      for (range = 0; range < numRanges; ++range)
      {
        if (rangeSel[range] == selection)
        {
          for (; sym < srcAlpha->rangeEndIndices[range]; ++sym)
            mappings[sym] = destSym++;
          newRanges[range] = srcAlpha->symbolsPerRange[range];
        }
        else
        {
          for (; sym < srcAlpha->rangeEndIndices[range]; ++sym)
            mappings[sym] = fallback;
          newRanges[range] = 0;
        }
      }
      newAlpha = newMultiRangeAlphabetEncodingUInt8(numRanges, newRanges,
                                                    mappings);
      ma_free(mappings);
      ma_free(newRanges);
    }
    break;
  default:
    abort();
    break;
  }
  return newAlpha;
}

MRAEnc *
MRAEncAddSymbolToRange(MRAEnc *mralpha, Symbol sym, AlphabetRangeID range)
{
  Symbol insertPos, numSyms;
  assert(mralpha && range < mralpha->numRanges);
  insertPos = mralpha->rangeEndIndices[range];
  numSyms = mralpha->rangeEndIndices[mralpha->numRanges - 1];
  switch (mralpha->encType)
  {
  case sourceUInt8:
    {
      MRAEncUInt8 *ui8alpha;
      ui8alpha = MRAEnc2MRAEncUInt8(mralpha);
      assert(ui8alpha->mappings[sym] == UNDEF_UCHAR);
      /* first move all old mappings accordingly */
      {
        Symbol i;
        for (i = numSyms; i > insertPos; --i)
        {
          Symbol origSym = ui8alpha->revMappings[i - 1];
          ui8alpha->revMappings[i] = origSym;
          ui8alpha->mappings[origSym] += 1;
        }
      }
      /* do actual insertion */
      ui8alpha->mappings[sym] = insertPos;
      ui8alpha->revMappings[insertPos] = sym;
      /* adjust ranges */
      mralpha->symbolsPerRange[range] += 1;
      {
        AlphabetRangeID i;
        for (i = range; i < mralpha->numRanges; ++i)
        {
          mralpha->rangeEndIndices[i] += 1;
        }
      }
    }
    break;
  default:
    abort();
    break;
  }
  return mralpha;
}

/**
 * @return number of symbols actually read
 */
size_t
MRAEncReadAndTransform(const MRAEnc *mralpha, FILE *fp,
                       size_t numSyms, Symbol *dest)
{
  int retval = 0;
  switch (mralpha->encType)
  {
  case sourceUInt8:
    {
      const MRAEncUInt8 *ui8alpha;
      size_t i;
      ui8alpha = constMRAEnc2MRAEncUInt8(mralpha);
      for (i = 0; i < numSyms; ++i)
      {
        int c = getc(fp);
        if (c != EOF)
          dest[i] = ui8alpha->mappings[c];
        else
          break;
      }
      retval = i;
    }
    break;
  default:
    abort();
    break;
  }
  return retval;
}

void
MRAEncSymbolsTransform(const MRAEnc *mralpha, Symbol *symbols, size_t numSyms)
{
  switch (mralpha->encType)
  {
  case sourceUInt8:
    {
      const MRAEncUInt8 *ui8alpha;
      size_t i;
      ui8alpha = constMRAEnc2MRAEncUInt8(mralpha);
      for (i = 0; i < numSyms; ++i)
      {
        symbols[i] = ui8alpha->mappings[symbols[i]];
      }
    }
    break;
  default:
    abort();
    break;
  }
}

void
MRAEncSymbolsRevTransform(const MRAEnc *mralpha, Symbol *symbols,
                          size_t numSyms)
{
  switch (mralpha->encType)
  {
  case sourceUInt8:
    {
      const MRAEncUInt8 *ui8alpha;
      size_t i;
      ui8alpha = constMRAEnc2MRAEncUInt8(mralpha);
      for (i = 0; i < numSyms; ++i)
      {
        symbols[i] = ui8alpha->revMappings[symbols[i]];
      }
    }
    break;
  default:
    abort();
    break;
  }
}

extern int
MRAEncSymbolIsInSelectedRanges(const MRAEnc *mralpha, Symbol sym,
                               int selection, const int *rangeSel)
{
  AlphabetRangeID range = 0;
  assert(mralpha && rangeSel);
  while (range < mralpha->numRanges
         && sym >= mralpha->rangeEndIndices[range])
    ++range;
  if (range < mralpha->numRanges)
  {
    if (rangeSel[range] == selection
       && sym >= (mralpha->rangeEndIndices[range]
                  - mralpha->symbolsPerRange[range])
       /* implicitely: && sym < mralpha->rangeEndIndices[range] */)
      return 1;
    else
      return 0;
  }
  else
    return -1;
}

void
MRAEncDelete(struct multiRangeAlphabetEncoding *mralpha)
{
  assert(mralpha);
  switch (mralpha->encType)
  {
    MRAEncUInt8 *ui8alpha;
  case sourceUInt8:
    ui8alpha = MRAEnc2MRAEncUInt8(mralpha);
    ma_free(ui8alpha);
    break;
  default:
    abort();
    break;
  }
}
