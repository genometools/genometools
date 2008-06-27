/*
  Copyright (c) 2007,2008 Thomas Jahns <Thomas.Jahns@gmx.net>

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
#ifndef EIS_BWTSEQ_H
#define EIS_BWTSEQ_H

/**
 * \file eis-bwtseq.h
 * Interface definitions for querying an indexed representation of the
 * BWT of a sequence as presented by Manzini and Ferragina (Compressed
 * Representations of Sequences and Full-Text Indexes, 2006)
 */

#include "libgtcore/error.h"
#include "libgtcore/str.h"
#include "libgtmatch/verbose-def.h"
#include "libgtmatch/seqpos-def.h"

#include "libgtmatch/eis-encidxseq.h"
#include "libgtmatch/eis-mrangealphabet.h"
#include "libgtmatch/eis-bwtseq-param.h"

/* TODO:
 * - implement other index types
 */

/**
 * Defines mode by which to derive sort of symbols in particular ranges
 */
enum rangeSortMode {
  SORTMODE_VALUE     = 0,
  SORTMODE_UNDEFINED = 1,
  SORTMODE_RANK      = 2,
};

/**
 * Stores column indices of the (virtual) matrix of rotations of the
 * input string used to construct the BWT, note that upper will
 * typically contain the lower value since rows are numbered from 0 at
 * the top to n-1 at the bottom.
 */
struct matchBound
{
  Seqpos start,                 /**< index of first boundary row */
    end;                        /**< index of second boundary row */
};

/** Object holding a BWT sequence index */
typedef struct BWTSeq BWTSeq;
/** Iterator of Matches produced by a search query from a BWT
 * sequence index */
typedef struct BWTSeqExactMatchesIterator BWTSeqExactMatchesIterator;

/**
 * \brief Creates or loads an encoded indexed sequence object of the
 * BWT transform.
 * @param params a struct holding parameter information for index construction
 * @param err genometools error object reference
 * @return reference to new BWT sequence object
 */
extern BWTSeq *
availBWTSeq(const struct bwtParam *params, Verboseinfo *verbosity, Error *err);

/**
 * \brief Creates an encoded indexed sequence object of the BWT
 * transform by translating the information in a pre-existing suffix
 * array.
 * @param params a struct holding parameter information for index construction
 * @param err genometools error object reference
 * @return reference to new BWT sequence object
 */
extern BWTSeq *
trSuftab2BWTSeq(const struct bwtParam *params, Verboseinfo *verbosity,
                Error *err);

/**
 * \brief Loads an encoded indexed sequence object of the
 * BWT transform.
 * @param projectName
 * @param BWTOptFlags Selects in-memory features of sequence index to
 * use, see enum BWTOptionDefaultsOptimizationFlags for possible settings.
 * @param err genometools error object reference
 * @return reference to new BWT sequence object
 */
extern BWTSeq *
loadBWTSeq(const Str *projectName, int BWTOptFlags, Verboseinfo *verbosity,
           Error *err);

/**
 * \brief Deallocate a previously loaded/created BWT sequence object.
 * @param bwtseq reference of object to delete
 */
extern void
deleteBWTSeq(BWTSeq *bwtseq);

/**
 * \brief Query BWT sequence object for availability of added
 * information to locate matches.
 * @param bwtSeq reference of object to query
 * @return 0 if no locate information is present, non-zero otherwise
 */
static inline bool
BWTSeqHasLocateInformation(const BWTSeq *bwtSeq);

/**
 * \brief Retrieve alphabet transformation from BWT sequence object
 * @param bwtSeq reference of object to query for alphabet
 * @return read-only reference of alphabet associated with sequence
 */
static inline const MRAEnc *
BWTSeqGetAlphabet(const BWTSeq *bwtSeq);

/**
 * \brief Retrieve sequence index in which the BWT is stored
 * @param bwtSeq reference of object to query for index
 * @return read-only reference of index containing the sequence
 */
static inline const EISeq *
BWTSeqGetEncIdxSeq(const BWTSeq *bwtSeq);

/**
 * \brief Retrieve Symbol at given position in BWT
 * @param bwtSeq reference of object to query for index
 * @return read-only reference of index containing the sequence
 */
static inline Symbol
BWTSeqGetSym(const BWTSeq *bwtSeq, Seqpos pos);

/**
 * \brief Query length of stored sequence.
 * @param bwtSeq reference of object to query
 * @return length of sequence
 */
static inline Seqpos
BWTSeqLength(const BWTSeq *bwtSeq);

/**
 * @brief Query position of rotation 0 in suffix array (aka
 * longest) which also is the position of the terminator symbol in the
 * BWT.
 * @param bwtSeq reference of object to query
 * @return position of terminator
 */
static inline Seqpos
BWTSeqTerminatorPos(const BWTSeq *bwtSeq);

/**
 * \brief Query BWT sequence for the number of occurences of a symbol in a
 * given prefix.
 * @param bwtSeq reference of object to query
 * @param tSym transformed symbol (as obtained by
 * MRAEncMapSymbol(BWTSeqGetAlphabet(bwtSeq), origSym)
 * @param pos right bound of BWT prefix queried
 * @return number of occurrences of symbol up to but not including pos
 */
static inline Seqpos
BWTSeqTransformedOcc(const BWTSeq *bwtSeq, Symbol tSym, Seqpos pos);

/**
 * \brief Query BWT sequence for the number of occurences of a symbol in a
 * given prefix.
 * @param bwtSeq reference of object to query
 * @param sym symbol (from original alphabet)
 * @param pos right bound of BWT prefix queried
 * @return number of occurrences of symbol up to but not including pos
 */
static inline Seqpos
BWTSeqOcc(const BWTSeq *bwtSeq, Symbol sym, Seqpos pos);

/**
 * \brief Query BWT sequence for the number of occurences of a symbol
 * in two given prefixes.
 * @param bwtSeq reference of object to query
 * @param tSym transformed symbol (as obtained by
 * MRAEncMapSymbol(BWTSeqGetAlphabet(bwtSeq), origSym)
 * @param posA right bound of first BWT prefix queried
 * @param posA right bound of second BWT prefix queried
 * @return number of occurrences of symbol up to but not including
 * posA and posB respectively in fields a and b of returned struct
 */
static inline struct SeqposPair
BWTSeqTransformedPosPairOcc(const BWTSeq *bwtSeq, Symbol tSym,
                            Seqpos posA, Seqpos posB);

/**
 * \brief Query BWT sequence for the number of occurences of a symbol
 * in two given prefixes.
 * @param bwtSeq reference of object to query
 * @param sym symbol (from original alphabet)
 * @param posA right bound of first BWT prefix queried
 * @param posB right bound of second BWT prefix queried
 * @return number of occurrences of symbol up to but not including
 * posA and posB respectively in fields a and b of returned struct
 */
static inline struct SeqposPair
BWTSeqPosPairOcc(const BWTSeq *bwtSeq, Symbol sym, Seqpos posA, Seqpos posB);

/**
 * \brief Query BWT sequence for the number of occurences of all symbols in a
 * given alphabet range and BWT sequence prefix.
 * @param bwtSeq reference of object to query
 * @param range range of symbols in alphabet to query
 * @param pos right bound of BWT prefix queried
 * @param rangeOccs occurrence counts for all symbols in range are written to
 * this array. The referenced memory must be sized appropriately to
 * accomodate as many symbols as are in range (MRAEncGetRangeSize if
 * in doubt) and rangeOccs[i] will hold the occurrence count of symbol
 * MRAEncRevMapSymbol(alphabet, i + MRAEncGetRangeBase(alphabet, range))
 */
static inline void
BWTSeqRangeOcc(const BWTSeq *bwtSeq, AlphabetRangeID range, Seqpos pos,
               Seqpos *rangeOccs);

/* XXx: range 0 for range = 0 for regular symbols;
   Seqpos rangeOcc[4]; */

/**
 * \brief Query BWT sequence for the number of occurences of all symbols in a
 * given alphabet range and two BWT sequence prefixes.
 * @param bwtSeq reference of object to query
 * @param range range of symbols in alphabet to query
 * @param posA right bound of first BWT prefix queried
 * @param posB right bound of second BWT prefix queried
 * @param rangeOccs occurrence counts for all symbols in range are written to
 * this array. The referenced memory must be sized appropriately to
 * accomodate two-times as many positions symbols as are in range
 * (use MRAEncGetRangeSize if in doubt) and rangeOccs[i] will hold the
 * occurrence count of symbol
 * MRAEncRevMapSymbol(alphabet, i + MRAEncGetRangeBase(alphabet, range))
 * up to position posA while the corresponding is true for
 * rangeOccs[rangeSize + i] concerning posB
 */
static inline void
BWTSeqPosPairRangeOcc(const BWTSeq *bwtSeq, AlphabetRangeID range, Seqpos posA,
                      Seqpos posB, Seqpos *rangeOccs);

/**
 * \brief Given a position in the L-column of the matrix of rotations,
 * find the corresponding row in the F-column.
 * @param bwtSeq reference of object to query
 * @param pos row index for L-column
 * @param extBits modified in intermediate queries
 * @return index of corresponding row F-column
 */
static inline Seqpos
BWTSeqLFMap(const BWTSeq *bwtSeq, Seqpos pos, struct extBitsRetrieval *extBits);

/**
 * \brief Given a symbol, query the aggregate count of symbols with
 * lower index, this corresponds to the first row in the C-column of
 * standard literature on the BWT on which the given symbol is found.
 * @param bwtSeq reference of object to query
 * @param sym symbol to query counts sum for
 * @return aggregate count
 */
static inline Seqpos
BWTSeqAggCount(const BWTSeq *bwtSeq, Symbol sym);

/**
 * \brief Given a symbol, query the aggregate count of symbols with
 * lower index.
 *
 * This corresponds to the first row in the C-column of
 * standard literature on the BWT on which the given symbol is found,
 * this function takes a symbol already transformed into the stored
 * alphabet as argument.
 *
 * @param bwtSeq reference of object to query
 * @param tSym symbol to query counts sum for
 * reference for core functions @return aggregate count
 */
static inline Seqpos
BWTSeqAggTransformedCount(const BWTSeq *bwtSeq, Symbol tSym);

/**
 * \brief Given a query string find number of matches in original
 * sequence (of which the sequence object is a BWT).
 * @param bwtSeq reference of object to query
 * @param query symbol string to search matches for
 * @param queryLen length of query string
 * @param forward direction of processing the query
 * @return number of matches
 */
extern Seqpos
BWTSeqMatchCount(const BWTSeq *bwtSeq, const Symbol *query, size_t queryLen,
                 bool forward);

/**
 * \brief Given a pair of limiting positions in the suffix array and a
 * symbol, compute the interval reached by matching one symbol further.
 * @param bwtSeq reference of sequence index to query
 * @param nextSym symbol by which to further restrict match
 * @param limits current restriction of match interval
 * @return limits, with bounds adjusted
 */
static inline struct matchBound *
BWTSeqIncrMatch(const BWTSeq *bwtSeq, struct matchBound *limits,
                Symbol nextSym);

/**
 * Error conditions encountered upon integrity check.
 */
enum verifyBWTSeqErrCode
{
  VERIFY_BWTSEQ_NO_ERROR = 0,       /**< every check completed okay */
  VERIFY_BWTSEQ_REFLOAD_ERROR = -1, /**< failed to load suffix array for
                                     *   reference comparisons */
  VERIFY_BWTSEQ_LENCOMPARE_ERROR = -2, /* lengths of bwt sequence
                                        * index and loaded suffix arry
                                        * don't match */
  VERIFY_BWTSEQ_SUFVAL_ERROR = -3, /**< a marked suffix array value
                                    * stored in the bwt sequence index
                                    * does not match the value read directly
                                    * from the suffix array table */
  VERIFY_BWTSEQ_LFMAPWALK_ERROR = -4, /**< while traversing the bwt
                                       * sequence in reverse original sequence
                                       * order, the symbol retrieved
                                       * does not match the
                                       * corresponding symbol in the
                                       * encoded sequence */
  VERIFY_BWTSEQ_LFMAPWALK_IMP_ERROR = -5, /**< original sequence
                                           * regeneration was
                                           * requested,
                                           * but is impossible */
  VERIFY_BWTSEQ_TERMPOS_ERROR = -6, /**< the position of the
                                     * 0-rotation does not match */
  VERIFY_BWTSEQ_CONTEXT_SYMFAIL = -7, /**< context regeneration
                                       * delivered an incorrect symbol */
  VERIFY_BWTSEQ_CONTEXT_LOADFAIL = -8, /**< context regeneration
                                        * is impossible because the
                                        * context failed to load */
};

enum verifyBWTSeqFlags
{
  VERIFY_BWTSEQ_SUFVAL    = 1 << 0, /**< check stored suffix arrays */
  VERIFY_BWTSEQ_LFMAPWALK = 1 << 1, /**< performs full backwards
                                     * regeneration of original
                                     * sequence (if possible)
                                     */
  VERIFY_BWTSEQ_CONTEXT   = 1 << 2, /**< try some random context
                                     * regenerations
                                     */
};

/**
 * \brief Perform various checks on the burrows wheeler transform
 *
 * - inspect all sampled suffix array values for equality with
 *   corresponding value of mapped reference suffix array
 * - check wether the last-to-first traversal of the BWT sequence
 *   index delivers the reversed encoded sequence
 *
 * @param bwtSeq index to check
 * @param projectName suffix array to load as reference
 * @param tickPrint print a dot every time tickPrint many symbols have
 *                  been processed
 * @param fp dots printed to this file
 */
extern enum verifyBWTSeqErrCode
BWTSeqVerifyIntegrity(BWTSeq *bwtSeq, const Str *projectName,
                      int checkFlags,
                      unsigned long tickPrint, FILE *fp,
                      Verboseinfo *verbosity, Error *err);

/**
 * \brief Given a query string produce iterator for all matches in
 * original sequence (of which the sequence object is a BWT).
 *
 * Warning: the iterator object will become invalid once the
 * corresponding bwt sequence object has been deleted.
 *
 * Warning: user must manage storage of iter manually
 *
 * @param iter points to storage for iterator
 * @param bwtSeq reference of bwt sequence object to use for matching
 * @param query symbol string to search matches for
 * @param queryLen length of query string
 * @param forward direction of processing the query
 * @return true if successfully initialized, false on error
 */
extern bool
initEMIterator(BWTSeqExactMatchesIterator *iter, const BWTSeq *bwtSeq,
               const Symbol *query, size_t queryLen, bool forward);

/**
 * \brief Only initializes empty iterator for given
 * sequence object.
 *
 * Warning: the iterator object will become invalid once the
 * corresponding bwt sequence object has been deleted.
 *
 * Warning: user must manage storage of iter manually
 * @param iter points to storage for iterator
 * @param bwtSeq reference of bwt sequence object to use for matching
 * @return true if successfully initialized, false on error
 */
extern bool
initEmptyEMIterator(BWTSeqExactMatchesIterator *iter, const BWTSeq *bwt);

/**
 * \brief Set up iterator for new query, iter must have been
 * initialized previously. Everything else is identical to
 * initEMIterator.
 *
 * Warning: the iterator object will become invalid once the
 * corresponding bwt sequence object has been deleted.
 * @param iter points to storage for iterator
 * @param bwtSeq reference of bwt sequence object to use for matching
 * @param query symbol string to search matches for
 * @param queryLen length of query string
 * @param forward direction of processing the query
 * @return true if successfully initialized, false on error
 */
extern bool
reinitEMIterator(BWTSeqExactMatchesIterator *iter, const BWTSeq *bwtSeq,
                 const Symbol *query, size_t queryLen, bool forward);
/**
 * \brief Destruct resources of matches iterator. Does not free the
 * storage of iterator itself.
 * Warning: user must manage storage of iter manually
 * @param iter reference to matches iterator
 */
extern void
destructEMIterator(struct BWTSeqExactMatchesIterator *iter);

/**
 * \brief Given a query string produce iterator for all matches in
 * original sequence (of which the sequence object is a BWT).
 *
 * Warning: the iterator object will become invalid once the
 * corresponding bwt sequence object has been deleted.
 * @param bwtSeq reference of bwt sequence object to use for matching
 * @param query symbol string to search matches for
 * @param queryLen length of query string
 * @param forward direction of processing the query
 * @return reference of iterator object, NULL on error
 */
extern BWTSeqExactMatchesIterator *
newEMIterator(const BWTSeq *bwtSeq, const Symbol *query, size_t queryLen,
              bool forward);

/**
 * \brief Deallocate an iterator object.
 * @param iter reference of iterator object
 */
extern void
deleteEMIterator(BWTSeqExactMatchesIterator *iter);

/**
 * \brief Get position of next match from an iterator.
 * @param iter reference of iterator object
 * @param bwtSeq reference of bwt sequence object to use for matching
 * @param pos location of match will be stored here if a further match
 * is available, not modified otherwise
 * @return true if another match could be found, false otherwise
 */
static inline bool
EMIGetNextMatch(BWTSeqExactMatchesIterator *iter, Seqpos *pos,
                const BWTSeq *bwtSeq);

/**
 * \brief Query an iterator for the total number of matches.
 * @param iter reference of iterator object
 * @return total number of matches
 */
extern Seqpos
EMINumMatchesTotal(const BWTSeqExactMatchesIterator *iter);

/**
 * \brief Query an iterator for the number of matches not yet
 * inspected via EMIGetNextMatch.
 * @param iter reference of iterator object
 * @return number of matches left
 */
extern Seqpos
EMINumMatchesLeft(const BWTSeqExactMatchesIterator *iter);

/**
 * @brief for packed index (given as void pointer), compute the longest
 * prefix of string in range between qstart and qend that occurs exactly
 * once in the index.
 * @param packed index (given as void pointer)
 * @param qstart points to memory area where query is found
 * @param qend points to memory area immediately after the query
 * @param err
 * @return 0 if not unique, otherwise length of minmum unique prefix.
 */
unsigned long packedindexuniqueforward(const BWTSeq *bwtseq,
                                       const Uchar *qstart,
                                       const Uchar *qend);

unsigned long packedindexmstatsforward(const BWTSeq *bwtseq,
                                       Seqpos *witnessposition,
                                       const Uchar *qstart,
                                       const Uchar *qend);

#include "libgtmatch/eis-bwtseq-siop.h"

#endif
