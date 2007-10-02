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
#ifndef ENCIDXSEQ_H_INCLUDED
#define ENCIDXSEQ_H_INCLUDED

/**
 * \file encidxseq.h
 * Interface definitions for encoded indexed sequences.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */

#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif /* HAVE_STDINT_H */
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif /* HAVE_INTTYPES_H */

#include <libgtcore/bitpackstring.h>
#include <libgtcore/env.h>
#include <libgtcore/str.h>
#include <libgtmatch/seqpos-def.h>

#include "mrangealphabet.h"

/**
 * callback function to insert variable width data into encidx
 * construction
 * @param dest bitstring to append data to (can hold at least as many
 * bits as specified in the constructor)
 * @param sOffset position in dest at which the callee may start
 * writing
 * @param start the data will be written together with the encoding of
 * symbols at position start up to position start + len
 * @param len see parameter start for description
 * @param cbState is passed on every call back of bitInsertFunc to
 * pass information that is kept across individual calls
 * @param env genometools state, passes information about allocator etc
 * @return number of bits actually written, or (BitOffset)-1 if an
 * error occured
 */
typedef BitOffset (*bitInsertFunc)(BitString cwDest, BitOffset cwOffset,
                                   BitString varDest, BitOffset varOffset,
                                   Seqpos start, Seqpos len, void *cbState,
                                   Env *env);

typedef int (*headerWriteFunc)(FILE *fp, void *headerCBData);

/**
 * later inserted bits can be retrieved and will be presented in this
 * struct.
 * CAUTION: the bitstrings exposed in this manner become invalid if
 * delete is called for the corresponding sequence index.
 */
struct extBitsRetrieval
{
  BitOffset cwOffset,           /**< offset in constant-width string
                                 * at which data was earlier stored by
                                 * a bitInsertFunc */
    varOffset;                  /**< offset in variable-width string
                                 * at which data was earlier stored by
                                 * a bitInsertFunc  */
  Seqpos start,                 /**< start position of region for
                                 * which the retried information was
                                 * stored */
    len;                        /**< length of region */
  BitString cwPart,             /**< string containing
                                 * constant-width part of data  */
    varPart;                    /**< string containg variable-width data */
  int flags;                    /**< flags for internal use (Don't touch) */
};

enum extBitsRetrievalFlags
{
  EBRF_RETRIEVE_CWBITS = 0,
  EBRF_PERSISTENT_CWBITS = 1<<0,
  EBRF_RETRIEVE_VARBITS = 1<<1,
  EBRF_PERSISTENT_VARBITS = 1<<2,
};


enum rangeStoreMode {
  DIRECT_SYM_ENCODE,
  BLOCK_COMPOSITION_INCLUDE,
  REGIONS_LIST,
};

typedef union EISHint *EISHint;

extern struct encIdxSeq *
newBlockEncIdxSeq(const Str *projectName, unsigned blockSize,
                  unsigned superBlockBlocks,
                  size_t numExtHeaders, uint16_t *headerIDs,
                  uint32_t *extHeaderSizes, headerWriteFunc *extHeaderCallback,
                  void **headerCBData,
                  bitInsertFunc biFunc, BitOffset cwBitsPerPos,
                  BitOffset maxBitsPerPos, void *cbState, Env *env);
extern struct encIdxSeq *
loadBlockEncIdxSeq(const Str *projectName, Env *env);

extern void
deleteEncIdxSeq(struct encIdxSeq *seq, Env *env);

staticifinline inline const MRAEnc *
EISGetAlphabet(const struct encIdxSeq *seq);

/**
 * \brief Return number of occurrences of symbol sym in index up to
 * but not including given position.
 * @param seq sequence index object to query
 * @param sym original alphabet symbol to query occurrences of
 * @param pos occurences are counted up to this position
 * @param hint provides cache and direction information for queries
 * based on previous queries
 * @param env genometools state, passes information about allocator etc
 */
staticifinline inline Seqpos
EISRank(struct encIdxSeq *seq, Symbol sym, Seqpos pos, union EISHint *hint,
        Env *env);

/**
 * \brief Return number of occurrences of symbol sym in index up to
 * but not including given position.
 * @param seq sequence index object to query
 * @param sym symbol to query occurrences of, but already transformed
 * by input alphabet
 * @param pos occurences are counted up to this position
 * @param hint provides cache and direction information for queries
 * based on previous queries
 * @param env genometools state, passes information about allocator etc
 */
staticifinline inline Seqpos
EISSymTransformedRank(struct encIdxSeq *seq, Symbol msym, Seqpos pos,
                      union EISHint *hint, Env *env);

/**
 * Presents the bits previously stored by a bitInsertFunc callback.
 * @param seq
 * @param pos sequence position for which to retrieve corresponding
 * area
 * @param persistent if false, the retrieved BitString elements will
 * become invalid once the hint union is used in another call.
 * @param hint provides cache and direction information for queries
 * @param env genometools state, passes information about allocator etc
 */
staticifinline inline void
EISRetrieveExtraBits(struct encIdxSeq *seq, Seqpos pos, int flags,
                     struct extBitsRetrieval *retval, union EISHint *hint,
                     Env *env);

staticifinline inline void
initExtBitsRetrieval(struct extBitsRetrieval *r, Env *env);

staticifinline inline struct extBitsRetrieval *
newExtBitsRetrieval(Env *env);

staticifinline inline void
destructExtBitsRetrieval(struct extBitsRetrieval *r, Env *env);

staticifinline inline void
deleteExtBitsRetrieval(struct extBitsRetrieval *r, Env *env);

extern Seqpos
EISSelect(struct encIdxSeq *seq, Symbol sym, Seqpos count);

staticifinline inline Seqpos
EISLength(struct encIdxSeq *seq);

/**
 * Return symbol at specified position. Comparable to c[pos] if the
 * sequence was stored straight in an array c.
 * @param seq indexed sequence object to be queried
 * @param pos position to retrieve symbol for
 * @param hint optional caching/hinting structure (improves average
 * retrieval time)
 * @param env genometools state, passes information about allocator etc
 */
staticifinline inline Symbol
EISGetSym(struct encIdxSeq *seq, Seqpos pos, EISHint hint, Env *env);

staticifinline inline Symbol
EISGetTransformedSym(struct encIdxSeq *seq, Seqpos pos, EISHint hint, Env *env);

staticifinline inline EISHint
newEISHint(struct encIdxSeq *seq, Env *env);

staticifinline inline void
deleteEISHint(struct encIdxSeq *seq, EISHint hint, Env *env);

extern int
verifyIntegrity(struct encIdxSeq *seqIdx,
                Str *projectName, int tickPrint, FILE *fp, Env *env);

staticifinline inline FILE *
EISSeekToHeader(const struct encIdxSeq *seqIdx, uint16_t headerID,
                uint32_t *lenRet);

/**
 * Meant for testing purposes only.
 */
extern int
searchBlock2IndexPair(const struct encIdxSeq *seqIdx,
                      const Symbol *block,
                      size_t idxOutput[2], Env *env);

#ifdef HAVE_WORKING_INLINE
#include "../encidxseqsimpleop.c"
#endif

#endif /* ENCIDXSEQ_H_INCLUDED */
