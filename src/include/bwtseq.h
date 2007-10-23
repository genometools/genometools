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
#ifndef BWTSEQ_H_INCLUDED
#define BWTSEQ_H_INCLUDED

/**
 * \file bwtseq.h
 * Interface definitions for querying an indexed representation of the
 * BWT of a sequence as presented by Manzini and Ferragina (Compressed
 * Representations of Sequences and Full-Text Indexes, 2006)
 */

#include <libgtcore/env.h>
#include <libgtcore/str.h>
#include <libgtmatch/seqpos-def.h>

#include "encidxseq.h"
#include "mrangealphabet.h"

/* FIXME:
 * - implement other index types
 */

/**
 * Names the type of encoding used:
 */
enum seqBaseEncoding {
  BWT_BASETYPE_AUTOSELECT,      /**< automatic, load any index present
                                 *     (currently not implemented) */
  BWT_ON_RLE,                   /**< use original fmindex run-length
                                 * encoding  */
  BWT_ON_BLOCK_ENC,             /**< do block compression by dividing
                                 * sequence into strings of
                                 * composition and permutation
                                 * indices */
  BWT_ON_WAVELET_TREE_ENC,      /**< encode sequence with wavelet-trees */
};

/**
 * Stores information to construct the underlying sequence object of a
 * BWT sequence object.
 */
union bwtSeqParam
{
  struct 
  {
    unsigned blockSize;         /**< number of symbols to combine in
                                 * one block a lookup-table
                                 * containing
                                 * $alphabetsize^{blockSize}$ entries is
                                 * required so adjust with caution */
    unsigned bucketBlocks;      /**< number of blocks for which to
                                 * store partial symbol sums (lower
                                 * values increase index size and
                                 * decrease computations for lookup) */
  } blockEncParams;
/*   struct  */
/*   { */
/*   } RLEParams; */
/*   struct  */
/*   { */
/*   } waveletTreeParams; */
};

/**
 * Stores column indices of the (virtual) matrix of rotations of the
 * input string used to construct the BWT, note that upper will
 * typically contain the lower value since rows are numbered from 0 at
 * the top to n-1 at the bottom.
 */
struct matchBound
{
  Seqpos upper, lower;
};

typedef struct BWTSeq BWTSeq;

/**
 * \brief Creates or loads an encoded indexed sequence object of the
 * BWT transform.
 * @param baseType selects the encoding method of the sequence index
 * storing the BWT sequence (see enum seqBaseEncoding).
 * @param extraParams a union holding extra parameter information
 *   specific to the type selected via parameter baseType
 * @param projectName base file name to derive name of suffixerator
 *   project from
 * @param env genometools reference for core functions
 * @return reference to new BWT sequence object
 */
extern BWTSeq *
newBWTSeq(enum seqBaseEncoding baseType, unsigned locateInterval,
          union bwtSeqParam *extraParams, const Str *projectName, Env *env);

/**
 * \brief Creates or loads an encoded indexed sequence object of the
 * BWT transform.
 * @param baseType selects the encoding method of the sequence index
 * storing the BWT sequence (see enum seqBaseEncoding).
 * @param extraParams a union holding extra parameter information
 *   specific to the type selected via parameter baseType
 * @param sa Suffixarray data structure to build BWT index from
 * @param projectName base file name for index written (should be the
 * same as the one sa was read from
 * @param env genometools reference for core functions
 * @return reference to new BWT sequence object
 */
extern BWTSeq *
newBWTSeqFromSA(enum seqBaseEncoding baseType, unsigned locateInterval,
                union bwtSeqParam *extraParams, Suffixarray *sa,
                Seqpos totalLen, const Str *projectName, Env *env);

/**
 * \brief Deallocate a previously loaded/created BWT sequence object.
 * @param bwtseq reference of object to delete
 * @param env genometools reference for core functions
 */
extern void
deleteBWTSeq(BWTSeq *bwtseq, Env *env);

/**
 * \brief Query BWT sequence object for availability of added
 * information to locate matches.
 * @param bwtSeq reference of object to query
 * @return 0 if no locate information is present, non-zero otherwise
 */
extern int
BWTSeqHasLocateInformation(const BWTSeq *bwtSeq);

/**
 * \brief Retrieve alphabet transformation from BWT sequence object
 * @param bwtSeq reference of object to query for alphabet
 * @return read-only reference of alphabet associated with sequence
 */
staticifinline inline const MRAEnc *
BWTSeqGetAlphabet(const BWTSeq *bwtSeq);

/**
 * \brief Retrieve sequence index in which the BWT is stored
 * @param bwtSeq reference of object to query for index
 * @return read-only reference of index containing the sequence
 */
staticifinline inline const EISeq *
BWTSeqGetEncIdxSeq(const BWTSeq *bwtSeq);

/**
 * \brief Query length of stored sequence.
 * @param bwtseq reference of object to query
 * @return length of sequence
 */
staticifinline inline Seqpos
BWTSeqLength(const BWTSeq *seq);

/**
 * \brief Query BWT sequence for the number of occurences of a symbol in a
 * given prefix.
 * @param bwtSeq reference of object to query
 * @param tSym transformed symbol (as obtained by
 * MRAEncMapSymbol(BWTSeqGetAlphabet(bwtSeq), origSym)
 * @param pos right bound of BWT prefix queried
 * @param env genometools reference for core functions
 * @return number of occurrences of symbol up to but not including pos
 */
staticifinline inline Seqpos
BWTSeqTransformedOcc(const BWTSeq *bwtSeq, Symbol tsym, Seqpos pos, Env *env);

/**
 * \brief Query BWT sequence for the number of occurences of a symbol in a
 * given prefix.
 * @param bwtSeq reference of object to query
 * @param Sym symbol
 * @param pos right bound of BWT prefix queried
 * @param env genometools reference for core functions
 * @return number of occurrences of symbol up to but not including pos
 */
staticifinline inline Seqpos
BWTSeqOcc(const BWTSeq *bwtSeq, Symbol tsym, Seqpos pos, Env *env);

/**
 * \brief Given a position in the L-column of the matrix of rotations,
 * find the corresponding row in the F-column.
 * @param bwtSeq reference of object to query
 * @param pos row index for L-column
 * @param env genometools reference for core functions
 * @return index of corresponding row F-column
 */
staticifinline inline Seqpos
BWTSeqLFMap(const BWTSeq *bwtSeq, Seqpos pos, Env *env);

/**
 * \brief Given a symbol, query the aggregate count of symbols with
 * lower index, this corresponds to the first row in the C-column of
 * standard literature on the BWT on which the given symbol is found.
 * @param bwtSeq reference of object to query
 * @param sym symbol to query counts sum for
 * @param env genometools reference for core functions
 * @return aggregate count
 */
staticifinline inline Seqpos
BWTSeqAggCount(const BWTSeq *bwtSeq, Symbol sym, Env *env);

/**
 * \brief Given a symbol, query the aggregate count of symbols with
 * lower index, this corresponds to the first row in the C-column of
 * standard literature on the BWT on which the given symbol is found,
 * this function takes a symbol already transformed into the stored
 * alphabet as argument.  @param bwtSeq reference of object to query
 * @param sym symbol to query counts sum for @param env genometools
 * reference for core functions @return aggregate count
 */
staticifinline inline Seqpos
BWTSeqAggTransformedCount(const BWTSeq *bwtSeq, Symbol tSym, Env *env);

/**
 * \brief Given a query string find number of matches in original
 * sequence (of which the sequence object is a BWT).
 * @param bwtseq reference of object to query
 * @param query symbol string to search matches for
 * @param queryLen length of query string
 * @param env genometools reference for core functions
 * @return number of matches
 */
extern Seqpos
BWTSeqMatchCount(const BWTSeq *bwtseq, Symbol *query, size_t queryLen,
                 Env *env);

/**
 * \brief Given a pair of limiting .
 * @param bwtseq reference of object to query
 * @param query symbol string to search matches for
 * @param queryLen length of query string
 * @param env genometools reference for core functions
 * @return number of matches
 */
staticifinline inline struct matchBound *
BWTSeqIncrMatch(const BWTSeq *bwtSeq, struct matchBound *limits,
                Symbol nextSym, Env *env);

/**
 * \brief Given a query string produce iterator for all matches in
 * original sequence (of which the sequence object is a BWT).
 *
 * Warning: the iterator object will become invalid once the
 * corresponding bwt sequence object has been deleted.
 * @param bwtseq reference of bwt sequence object to use for matching
 * @param query symbol string to search matches for
 * @param queryLen length of query string
 * @param env genometools reference for core functions
 * @return reference of iterator object, NULL on error
 */
extern struct BWTSeqExactMatchesIterator *
newEMIterator(const BWTSeq *bwtSeq, Symbol *query, size_t queryLen, Env *env);

/**
 * \brief Deallocate an iterator object.
 * @param iter reference of iterator object
 * @param env genometools reference for core functions
 */
extern void
deleteEMIterator(struct BWTSeqExactMatchesIterator *iter, Env *env);

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
 * @param env genometools reference for core functions
 * @return reference to a structure that specifies the location of a
 * match or NULL if no further match is available, the reference  will
 * become invalid  once the iterator has been queried again
 */
extern struct MatchData *
EMIGetNextMatch(struct BWTSeqExactMatchesIterator *iter, const BWTSeq *bwtSeq,
                Env *env);

/**
 * \brief Query an iterator for the total number of matches.
 * @param iter reference of iterator object
 * @return total number of matches
 */
extern Seqpos
EMINumMatchesTotal(const struct BWTSeqExactMatchesIterator *iter);

/**
 * \brief Query an iterator for the number of matches not yet
 * inspected via EMIGetNextMatch.
 * @param iter reference of iterator object
 * @return number of matches left
 */
extern Seqpos
EMINumMatchesLeft(const struct BWTSeqExactMatchesIterator *iter);

#ifdef HAVE_WORKING_INLINE
#include "../bwtseqsimpleop.c"
#endif

#endif /* BWTSEQ_H_INCLUDED */
