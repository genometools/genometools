/*
  Copyright (c) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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
#include "libgtmatch/seqpos-def.h"

#include "libgtmatch/eis-encidxseq.h"
#include "libgtmatch/eis-mrangealphabet.h"
#include "libgtmatch/eis-bwtseqparam.h"

/* TODO:
 * - implement other index types
 */

/**
 * Stores column indices of the (virtual) matrix of rotations of the
 * input string used to construct the BWT, note that upper will
 * typically contain the lower value since rows are numbered from 0 at
 * the top to n-1 at the bottom.
 */
struct matchBound
{
  Seqpos upper,                 /**< index of first boundary row */
    lower;                      /**< index of second boundary row */
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
availBWTSeq(const struct bwtParam *params, Error *err);

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
loadBWTSeq(const Str *projectName, int BWTOptFlags, Error *err);

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
 * \brief Query length of stored sequence.
 * @param bwtSeq reference of object to query
 * @return length of sequence
 */
static inline Seqpos
BWTSeqLength(const BWTSeq *bwtSeq);

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
 * @param tSym transformed symbol (as obtained by
 * MRAEncMapSymbol(BWTSeqGetAlphabet(bwtSeq), origSym)
 * @param pos right bound of BWT prefix queried
 * @return number of occurrences of symbol up to but not including pos
 */
static inline Seqpos
BWTSeqOcc(const BWTSeq *bwtSeq, Symbol tSym, Seqpos pos);

/**
 * \brief Given a position in the L-column of the matrix of rotations,
 * find the corresponding row in the F-column.
 * @param bwtSeq reference of object to query
 * @param pos row index for L-column
 * @return index of corresponding row F-column
 */
static inline Seqpos
BWTSeqLFMap(const BWTSeq *bwtSeq, Seqpos pos);

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
 * @return number of matches
 */
extern Seqpos
BWTSeqMatchCount(const BWTSeq *bwtSeq, const Symbol *query, size_t queryLen);

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
                      unsigned long tickPrint, FILE *fp, Error *err);

/**
 * \brief Given a query string produce iterator for all matches in
 * original sequence (of which the sequence object is a BWT).
 *
 * Warning: the iterator object will become invalid once the
 * corresponding bwt sequence object has been deleted.
 * @param bwtSeq reference of bwt sequence object to use for matching
 * @param query symbol string to search matches for
 * @param queryLen length of query string
 * @return reference of iterator object, NULL on error
 */
extern BWTSeqExactMatchesIterator *
newEMIterator(const BWTSeq *bwtSeq, const Symbol *query, size_t queryLen);

/**
 * \brief Deallocate an iterator object.
 * @param iter reference of iterator object
 */
extern void
deleteEMIterator(BWTSeqExactMatchesIterator *iter);

/**
 * location data corresponding to a match
 */
struct MatchData
{
  const char *dbFile;           /**< name of original sequence file
                                 *  the match is from  */
  Seqpos sfxArrayValue,         /**< position of match in concatenated
                                 *  encoded sequence representation
                                 *  of multiple database files */
    dbFilePos;                  /**< position of match in original
                                 *  sequence file dbFile */
};

/**
 * \brief Get position of next match from an iterator.
 * @param iter reference of iterator object
 * @param bwtSeq reference of bwt sequence object to use for matching
 * @return reference to a structure that specifies the location of a
 * match or NULL if no further match is available, the reference  will
 * become invalid  once the iterator has been queried again
 */
static inline struct MatchData *
EMIGetNextMatch(BWTSeqExactMatchesIterator *iter, const BWTSeq *bwtSeq);

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
 * \brief for packed index (given as void pointer), compute the longest
 * prefix of string in range between qstart and qend that occurs exactly
 * once in the index.
 * @param packed index (given as void pointer)
 * @param qstart points to memory area where query is found
 * @param qend points to memory area immediately after the query
 * @param err
 * @return 0 if not unique, otherwise length of minmum unique prefix.
 */
unsigned long packedindexuniqueforward(const void *genericindex,
                                       /*@unused@*/ unsigned long offset,
                                       /*@unused@*/ Seqpos left,
                                       /*@unused@*/ Seqpos right,
                                       /*@unused@*/ Seqpos *witnessposition,
                                       const Uchar *qstart,
                                       const Uchar *qend);

unsigned long packedindexmstatsforward(const void *genericindex,
                                       /*@unused@*/ unsigned long offset,
                                       /*@unused@*/ Seqpos left,
                                       /*@unused@*/ Seqpos right,
                                       Seqpos *witnessposition,
                                       const Uchar *qstart,
                                       const Uchar *qend);

#include "libgtmatch/eis-bwtseqsimpleop.h"

#endif
