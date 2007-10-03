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

#include "mrangealphabet.h"

/* FIXME:
 * - implement other index types
 */

enum seqBaseEncoding {
  BWT_BASETYPE_AUTOSELECT,
  BWT_ON_RLE,
  BWT_ON_BLOCK_ENC,
  BWT_ON_WAVELET_TREE_ENC,
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


typedef struct BWTSeq BWTSeq;

/**
 * \brief Creates or loads an encoded indexed sequence object of the
 * BWT transform.
 * @param baseType selects one of
 *   - BWT_BASETYPE_AUTOSELECT automatic, load any index present
 *     (currently not implemented)
 *   - BWT_ON_RLE use original fmindex run-length encoding
 *   - BWT_ON_BLOCK_ENC do block compression by dividing sequence into
 *     strings of composition and permutation indices
 *   - BWT_ON_WAVELET_TREE_ENC encode sequence with wavelet-trees
 * @param extraParams a union holding extra parameter information
 *   specific to the type selected via parameter baseType
 * @param projectName base file name to derive name of suffixerator
 *   project from
 * @param env genometools reference for core functions
 * @return reference to new BWT sequence object
 */
extern BWTSeq *
newBWTSeq(enum seqBaseEncoding baseType, union bwtSeqParam *extraParams,
          const Str *projectName, Env *env);

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
 * @param bwtseq reference of object to query
 * @return 0 if no locate information is present, non-zero otherwise
 */
extern int
BWTSeqHasLocateInformation(const BWTSeq *bwtseq);

/**
 * \brief Retrieve alphabet transformation from BWT sequence object
 * @param bwtseq reference of object to query for alphabet
 * @return read-only reference of alphabet associated with sequence
 */
staticifinline inline const MRAEnc *
BWTSeqGetAlphabet(const BWTSeq *bwtseq);

/**
 * \brief Query length of stored sequence.
 * @param bwtseq reference of object to query
 * @return length of sequence
 */
staticifinline inline Seqpos
BWTSeqLength(const BWTSeq *seq);

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
