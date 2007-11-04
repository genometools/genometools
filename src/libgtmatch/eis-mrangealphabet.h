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
#ifndef EIS_MRANGEALPHABET_H
#define EIS_MRANGEALPHABET_H

/**
 * \file mrangealphabet.h
 * \brief Methods for an alphabet mapping where the alphabet is mapped
 * to multiple contiguous ranges.
 *
 * The mapping is constructed so that an alphabet with n symbols
 * is mapped to the range 0..n-1 and any symbol can be queried for
 * what range it belongs to. 
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */

#include <inttypes.h>

#include "libgtcore/env.h"
#include "libgtmatch/alphadef.h"

typedef unsigned char Symbol;
typedef struct multiRangeAlphabetEncoding MRAEnc;

#define UNDEF_UCHAR ((unsigned char)~0)

enum sourceEncType {
  sourceUnknown = 0,
  sourceUInt8 = 1,
};

/**
 * \brief Create an alphabet to map a selection of unsigned characters
 * onto the resulting ranges.
 * @param numRanges number of distinct ranges the alphabet is to be
 * divided into.
 * @param symbolsPerRange gives the number of symbols for every range
 * @param mapping maps all uint8_t input symbols to values in the
 * range 0..sum(symbolsPerRange) or the input alphabet specific
 * value UNDEF_UCHAR
 * @param env 
 */
staticifinline inline MRAEnc *
MRAEncUInt8New(int numRanges, int symbolsPerRange[],
               const uint8_t *mapping, Env *env);

/**
 * \brief Creates a mapping to two ranges (regular and special
 * symbols) from a given genometools alphabet.
 * @param alpha original alphabet
 * @param env
 */
extern MRAEnc *
MRAEncGTAlphaNew(const Alphabet *alpha, Env *env);

/**
 * \brief alias of MRAEncUInt8New
 */
extern MRAEnc *
newMultiRangeAlphabetEncodingUInt8(int numRanges, const int symbolsPerRange[],
                                   const uint8_t *mappings, Env *env);

/**
 * \brief remap multirange alphabet by excluding some ranges
 *
 * Maps symbols from all included ranges to new values 0 to n where all
 * symbols from non-included ranges are mapped to fallback
 * @param srcAlpha alphabet to remap
 * @param rangeIncludeFlag array of integer flags, 0 => map all symbols in
 * range to fallback, otherwise append to already mapped symbols
 * @param fallback symbol to map not-included ranges to
 */
extern MRAEnc *
MRAEncSecondaryMapping(const MRAEnc *srcAlpha, int selection,
                       const int *rangeSel, Symbol fallback, Env *env);

/**
 * \brief Inserts a previously unmapped symbol into a range.
 *
 * All successor ranges will therefore be shifted by one in the
 * resulting mapping.
 * @param mralpha alphabet to add mapping to
 * @param sym input code to map
 * @param range number of range to insert new symbol into
 */
extern void
MRAEncAddSymbolToRange(MRAEnc *mralpha, Symbol sym, int range);

/**
 * \brief Query number of ranges in alphabet.
 * @param mralpha alphabet to query for number of ranges
 * @return number of ranges
 */
extern size_t
MRAEncGetNumRanges(const MRAEnc *mralpha);

/**
 * \brief Query number of symbols in given range of alphabet.
 * @param mralpha alphabet to get range from
 * @param range 
 */
staticifinline inline size_t
MRAEncGetRangeSize(const MRAEnc *mralpha, size_t range);

/**
 * @return number of symbols in alphabet
 */
extern size_t
MRAEncGetSize(const MRAEnc *mralpha);

/**
 * @return size of original value range of symbols in alphabet
 * (i.e. 256 for 8 bit mapping)
 */
staticifinline size_t
MRAEncGetDomainSize(const MRAEnc *mralpha);

/**
 * \brief Look up code of symbol from input domain in output range.
 * @return output code or input specific code for illegal symbol (in
 * which case MRAEncSymbolHasValidMapping would have returned false).
 */
staticifinline inline Symbol
MRAEncMapSymbol(const MRAEnc *mralpha, Symbol sym);

staticifinline inline int
MRAEncSymbolHasValidMapping(const MRAEnc *mralpha, Symbol sym);

staticifinline inline Symbol
MRAEncRevMapSymbol(const MRAEnc *mralpha, Symbol sym);

extern void
MRAEncSymbolsTransform(const MRAEnc *mralpha, Symbol *symbols, size_t numSyms);

extern void
MRAEncSymbolsRevTransform(const MRAEnc *mralpha, Symbol *symbols,
                          size_t numSyms);
/**
 * @param mralpha alphabet to query
 * @param sym symbol to look-up range for
 * @param selection value to test rangeSel for
 * @param rangeSel array of codes for every range in alphabet, to be
 * compared with selection
 * @return > 0 if sym is included in a range r which has rangeSel[r]
 * == selection, 0 if sym is not in any selected range, and <0 if sym
 * is out of the alphabets range of symbols
 */
extern int
MRAEncSymbolIsInSelectedRanges(const MRAEnc *mralpha, Symbol sym,
                               int selection, int *rangeSel);

extern int
MRAEncReadAndTransform(const MRAEnc *mralpha, FILE *fp,
                       size_t numSyms, Symbol *dest);

extern void
MRAEncDelete(struct multiRangeAlphabetEncoding *mralpha, Env *env);

#include "libgtmatch/eis-mrangealphabetsimpleop.h"

#endif /* EIS_MRANGEALPHABET_H */
