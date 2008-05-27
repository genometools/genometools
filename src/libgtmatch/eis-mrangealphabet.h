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
 * \file eis-mrangealphabet.h
 * \brief Methods for an alphabet mapping where the alphabet is mapped
 * to multiple contiguous ranges.
 *
 * The mapping is constructed so that an alphabet with n symbols
 * is mapped to the range 0..n-1 and any symbol can be queried for
 * what range it belongs to.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */

#include <inttypes.h>

#include "libgtmatch/alphadef.h"

/** symbols are stored in this scalar type */
typedef unsigned char Symbol;
#define PRIuSymbol "u"
/** let's assume 65536 ranges are enough */
typedef unsigned AlphabetRangeID;
/** Caution: sizeof (AlphabetRangeID) >= sizeof (AlphabetRangeSize)
 *           must be true, AlphabetRangeSize must be able to
 *           represent the total number of symbols */
typedef unsigned short AlphabetRangeSize;
/** retrieve symbol from BitString */
#define bsGetSymbol bsGetUInt8
/** store symbol in BitString */
#define bsStoreSymbol bsStoreUInt8
/** retrieve array of symbols from BitString  */
#define bsGetUniformSymbolArray bsGetUniformUInt8Array
/** store array of symbols in BitString  */
#define bsStoreUniformSymbolArray bsStoreUniformUInt8Array
/** how many bits are required to store given symbol value */
#define requiredSymbolBits requiredUInt8Bits

/**
 * Describes an alphabet i.e. a mapping function from a given input
 * type to a contiguous range of integers which are divided into
 * multiple ranges, which are again continuous.
 */
typedef struct multiRangeAlphabetEncoding MRAEnc;

/** used to describe a symbol not occurring in the alphabet */
#define UNDEF_UCHAR ((unsigned char)~0)

/** select width of symbol input */
enum sourceEncType {
  sourceUnknown = 0,            /**< invalid/undefined  */
  sourceUInt8 = 1,              /**< input in the range 0..255 */
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
 */
static inline MRAEnc *
MRAEncUInt8New(AlphabetRangeID numRanges, AlphabetRangeSize symbolsPerRange[],
               const uint8_t *mapping);

/**
 * \brief Creates a mapping to two ranges (regular and special
 * symbols) from a given genometools alphabet.
 * @param alpha original alphabet
 */
extern MRAEnc *
MRAEncGTAlphaNew(const Alphabet *alpha);

/**
 * \brief alias of MRAEncUInt8New
 */
extern MRAEnc *
newMultiRangeAlphabetEncodingUInt8(AlphabetRangeID numRanges,
                                   const AlphabetRangeSize symbolsPerRange[],
                                   const uint8_t *mappings);

/**
 * @brief Copy constructor for multi-range alphabets
 * @param alpha alphabet to copy
 * @return new alphabet object
 */
extern MRAEnc *
MRAEncCopy(const MRAEnc *alpha);

/**
 * \brief remap multirange alphabet by excluding some ranges
 *
 * Maps symbols from all included ranges to new values 0 to n where all
 * symbols from non-included ranges are mapped to fallback
 * @param srcAlpha alphabet to remap
 * @param selection ranges with this value in rangeSel are carried
 * over to new alphabet
 * @param rangeSel array of integer flags, if != selection for given
 * range => maps all symbols in range to fallback, otherwise append to
 * already mapped symbols
 * @param fallback symbol to map not-included ranges to
 */
extern MRAEnc *
MRAEncSecondaryMapping(const MRAEnc *srcAlpha, int selection,
                       const int *rangeSel, Symbol fallback);

/**
 * \brief Inserts a previously unmapped symbol into a range.
 *
 * All successor ranges will therefore be shifted by one in the
 * resulting mapping.
 * @param mralpha alphabet object reference to add mapping to
 * @param sym input code to map
 * @param range number of range to insert new symbol into
 * @return same reference as mralpha but after destructive change
 */
extern MRAEnc *
MRAEncAddSymbolToRange(MRAEnc *mralpha, Symbol sym, AlphabetRangeID range);

/**
 * \brief Query number of ranges in alphabet.
 * @param mralpha alphabet to query for number of ranges
 * @return number of ranges
 */
static inline AlphabetRangeID
MRAEncGetNumRanges(const MRAEnc *mralpha);

/**
 * \brief Query number of symbols in given range of alphabet.
 * @param mralpha alphabet to get range from
 * @param range
 */
static inline AlphabetRangeSize
MRAEncGetRangeSize(const MRAEnc *mralpha, AlphabetRangeID range);

/**
 * \brief Query symbol which starts range.
 * @param mralpha
 * @param range start symbol of range, symbols of range range from
 * MRAEncGetRangeBase(alph, range) to MRAEncGetRangeBase(alph, range + 1) - 1
 */
static inline Symbol
MRAEncGetRangeBase(const MRAEnc *mralpha, AlphabetRangeID range);

/**
 * @brief Get number of different symbols in alphabet.
 * @param mralpha
 * @return number of symbols in alphabet
 */
extern AlphabetRangeSize
MRAEncGetSize(const MRAEnc *mralpha);

/**
 * @brief Get range of input symbols.
 * @param mralpha
 * @return size of original value range of symbols in alphabet
 * (i.e. 256 for 8 bit mapping)
 */
static size_t
MRAEncGetDomainSize(const MRAEnc *mralpha);

/**
 * \brief Look up code of symbol from input domain in output range.
 * @return output code or input specific code for illegal symbol (in
 * which case MRAEncSymbolHasValidMapping would have returned false).
 */
static inline Symbol
MRAEncMapSymbol(const MRAEnc *mralpha, Symbol sym);

/**
 * @brief Find wether a symbol from input is accurately represented in
 * the alphabet or illegal input.
 * @param mralpha
 * @param sym symbol to map
 * @return 0 if alphabet has no valid mapping, !0 otherwise
 */
static inline int
MRAEncSymbolHasValidMapping(const MRAEnc *mralpha, Symbol sym);

/**
 * @brief Apply reverse mapping of transformed symbol to input alphabet.
 * @param mralpha
 * @param sym symbol to un-map
 * @return inverse mapping of sym
 */
static inline Symbol
MRAEncRevMapSymbol(const MRAEnc *mralpha, Symbol sym);

/**
 * @brief Apply mapping of input string.
 * @param mralpha
 * @param symbols symbols to convert
 * @param numSyms length of symbols string
 */
extern void
MRAEncSymbolsTransform(const MRAEnc *mralpha, Symbol *symbols, size_t numSyms);

/**
 * @brief Apply reverse mapping of string.
 * @param mralpha
 * @param symbols symbols to convert
 * @param numSyms length of symbols string
 */
extern void
MRAEncSymbolsRevTransform(const MRAEnc *mralpha, Symbol *symbols,
                          size_t numSyms);
/**
 * @brief Query wether a symbol belongs to a selected range
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
                               int selection, const int *rangeSel);

/**
 * @brief Query range for given symbol
 * @param mralpha alphabet to query
 * @param sym symbol to look-up range for (already transformed by
 * mapping corresponding to alphabet)
 * @return range id
 */
static inline AlphabetRangeID
MRAEncGetRangeOfSymbol(const MRAEnc *mralpha, Symbol sym);

/**
 * @brief Read symbols from file and transform according to
 * mapping.
 * @param mralpha
 * @param fp file pointer
 * @param numSyms read this many symbols
 * @param dest write converted symbols here
 * @return number of symbols actually read
 */
extern size_t
MRAEncReadAndTransform(const MRAEnc *mralpha, FILE *fp,
                       size_t numSyms, Symbol *dest);

/**
 * @brief Delete alphabet object.
 * @param mralpha
 */
extern void
MRAEncDelete(struct multiRangeAlphabetEncoding *mralpha);

#include "libgtmatch/eis-mrangealphabet-siop.h"

#endif
