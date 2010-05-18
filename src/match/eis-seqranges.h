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

#ifndef EIS_SEQRANGES_H
#define EIS_SEQRANGES_H

/**
 * @file eis-seqranges.h
 * @brief Store lists of ranges of symbols.
 *
 * This is intended to reduce the storage overhead introduced by the
 * presence of a special class of symbols in a sequence. These symbols
 * are generally expected to:
 * - be underrepresented in the sequence
 * - often form contiguous ranges
 * - require special treatment by algorithms
 */
#include <stdlib.h>

#include "core/error.h"
#include "core/bitpackstring.h"

#include "match/eis-mrangealphabet.h"

/**
 * Select wether certain internal data structures are constructed with
 * these flags.
 */
enum SRLFeatures {
  SRL_NO_FEATURES = 0,              /**< build no special structures */
  SRL_PARTIAL_SYMBOL_SUMS = 1 << 0, /**< for every range also keep
                                     * sums of other special symbols
                                     * encountered up to the
                                     * beginning of the range */
};

/**
 * Used to remember previous queries.
 */
typedef size_t seqRangeListSearchHint;

/**
 * @brief Constructor for list of sequence ranges.
 * @param rangesStartNum allocate space for this many ranges
 * @param alphabet map stored symbols according to this alphabet
 * @param features see enum SRLFeatures for a description
 * @return newly constructed list
 */
struct seqRangeList *
gt_newSeqRangeList(size_t rangesStartNum, const MRAEnc *alphabet,
                enum SRLFeatures features);
/**
 * @brief Shrink the list according to the exact number of ranges stored.
 * @param rangeList
 */
void
gt_SRLCompact(struct seqRangeList *rangeList);

/**
 * @brief Destructor for sequence range lists.
 * @param rangeList
 */
void
gt_deleteSeqRangeList(struct seqRangeList *rangeList);

/**
 * @brief Add a new range at the end of the current list.
 * @param rangeList
 * @param pos start of new range
 * @param len length of new range
 * @param sym the range is a contiguous sequence of this symbol
 */
void
gt_SRLAppendNewRange(struct seqRangeList *rangeList,
                  unsigned long pos,
                  unsigned long len,
                  Symbol sym);

/**
 * @brief Add a new range of length one
 * @param rangeList
 * @param pos at this position
 * @param sym this symbol occurs
 */
void
gt_SRLAddPosition(struct seqRangeList *rangeList, unsigned long pos,
               Symbol sym);

/**
 * @brief Initialize a search hint by this function.
 * @param rangeList
 * @param hint points to storage
 */
void
gt_SRLInitListSearchHint(struct seqRangeList *rangeList,
                      seqRangeListSearchHint *hint);

/**
 * @brief Find the range overlapping or if no such range exists
 * follows pos.
 * @param rangeList
 * @param pos
 * @param hint
 * @return NULL if no range overlaps or succeeds pos
 */
struct seqRange *
gt_SRLFindPositionNext(struct seqRangeList *rangeList, unsigned long pos,
                    seqRangeListSearchHint *hint);

/**
 * @brief This predicate is true if position is within one of the
 * ranges of the range list.
 * @param rangeList
 * @param pos
 * @param hint
 * @param symAtPos if not NULL, the symbol of the range overlapping
 * pos is written to this address, if such a range exists
 * @return true if an overlap exits, false if not
 */
int
gt_SRLOverlapsPosition(struct seqRangeList *rangeList, unsigned long pos,
                    seqRangeListSearchHint *hint, Symbol *symAtPos);

/**
 * @brief For the sequence region [start..end-1] count the number of
 * occurrences for each symbol represented in the overlapping ranges.
 * @param rangeList
 * @param start
 * @param end
 * @param occStore for each symbol occStore[MRAEncMapSymbol(alphabet,
 * sym)] is set to the number of occurrences of sym, where alphabet is
 * the alphabet originally used in the constructor
 * @param hint
 */
void
gt_SRLSymbolsInSeqRegion(struct seqRangeList *rangeList, unsigned long start,
                      unsigned long end, unsigned long *occStore,
                      seqRangeListSearchHint *hint);

/**
 * @brief Compute the occurrence count for one symbol in a given region.
 * @param rangeList
 * @param start
 * @param end
 * @param sym only account for ranges matching this symbol
 * @param hint
 */
unsigned long
gt_SRLSymbolCountInSeqRegion(struct seqRangeList *rangeList,
                          unsigned long start,
                          unsigned long end,
                          Symbol sym,
                          seqRangeListSearchHint *hint);

/**
 * @brief Sum over the occurrence counts for all symbols in a given region.
 * @param rangeList
 * @param start
 * @param end
 * @param hint
 */
unsigned long
gt_SRLAllSymbolsCountInSeqRegion(struct seqRangeList *rangeList,
                              unsigned long start,
                              unsigned long end,
                              seqRangeListSearchHint *hint);

/**
 * @brief Overwrite all positions in a string that coincide with
 * ranges in the list with the symbol for that range.
 *
 * @param rangeList
 * @param subString write symbols in ranges at subStringOffset+i to subString[i]
 * @param start start only use ranges overlapping [start..end-1]
 * @param len
 * @param subStringOffset offset of subString relative to the
 * underlying sequence
 * @param hint
 */
void
gt_SRLApplyRangesToSubString(struct seqRangeList *rangeList,
                          Symbol *subString,
                          unsigned long start,
                          unsigned long len,
                          unsigned long subStringOffset,
                          seqRangeListSearchHint *hint);

/**
 * @brief Print text description of ranges overlapping the given
 * region.
 *
 * @param rangeList
 * @param fp print description to this file
 * @param start
 * @param len
 * @param hint
 * @return <0 in case of I/O error, >=0 otherwise.
 */
int
gt_SRLPrintRangesInfo(struct seqRangeList *rangeList,
                   FILE *fp,
                   unsigned long start,
                   unsigned long len,
                   seqRangeListSearchHint *hint);
/**
 * @brief Save a range list structure to file.
 * @param rangeList
 * @param fp file to save rangeList to
 * @return <0 if an error occurred
 */
int
gt_SRLSaveToStream(struct seqRangeList *rangeList, FILE *fp);

/**
 * @brief Restore a sequence range list from file.
 * @param fp read from this file
 * @param alphabet symbols are interpreted under this alphabet
 * @param features see enum SRLFeatures
 * @param err
 */
struct seqRangeList *
gt_SRLReadFromStream(FILE *fp, const MRAEnc *alphabet,
                  enum SRLFeatures features, GtError *err);

#endif
