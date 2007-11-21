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

#ifndef EIS_ENCIDXSEQ_H
#define EIS_ENCIDXSEQ_H

/**
 * \file encidxseq.h
 * Interface definitions for encoded indexed sequences.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */
#include <stdio.h>
#include <inttypes.h>

#include "libgtcore/bitpackstring.h"
#include "libgtcore/env.h"
#include "libgtcore/str.h"
#include "libgtmatch/seqpos-def.h"
#include "libgtmatch/eis-mrangealphabet.h"

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

struct seqStats
{
  Seqpos *symbolDistributionTable;
  enum sourceEncType sourceAlphaType;
};

typedef struct encIdxSeq EISeq;

typedef union EISHint *EISHint;

enum EISFeatureBits
{
  EIS_FEATURE_NONE = 0,
  EIS_FEATURE_REGION_SUMS = 1<<0, /**< if set construct sum tables for
                                   *   the special symbol ranges
                                   *   use this on index loading if
                                   *   you want to use this index in
                                   *   many queries, omit if memory
                                   *   is very tight (e.g. on construction) */
};

struct blockEncParams
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
  int EISFeatureSet;          /**< bitwise or of EIS_FEATURE_NONE
                               * and other features selectable via
                               * enum EISFeatureBits (see eis-encidxseq.h) */
};

/**
 * \brief Construct block-encoded indexed sequence object and write
 * corresponding representation to disk.
 * @param projectName base name of corresponding suffixerator project
 * @param blockSize number of symbol to combine in one block, a
 * lookup-table containing $alphabetsize^{blockSize}$ entries is
 * required so adjust with caution
 * @param bucketBlocks frequency at which to store partial symbol
 * sums, low values increase storage and decrease number of
 * computations required for rank queries on average
 * @param numExtHeaders number of extension headers to write via callbacks
 * @param headerIDs array of numExtHeaders ids to be used
 * for each extension header in turn
 * @param extHeaderSizes array of numExtHeaders sizes
 * representing the length of each extension header
 * @param extHeaderCallback array of numExtHeaders function pointers
 * each of which will be called once upon writing the header
 * @param headerCBData array of pointers passed as argument when the
 * corresponding header writing function is called
 * @param biFunc function to be called when a chunk of data has been
 * accumulated for a given region of sequence data
 * @param cwBitsPerPos exactly this many bits will be appended by
 * biFunc for each symbol of the input sequence
 * @param maxBitsPerPos at most this many bits will be appended to the
 * variable width part of the data
 * @param cbState will be passed on each call of biFunc
 * @param env genometools reference for core functions
 */
extern EISeq *
newBlockEncIdxSeq(const Str *projectName, const struct blockEncParams *params,
                  size_t numExtHeaders, uint16_t *headerIDs,
                  uint32_t *extHeaderSizes, headerWriteFunc *extHeaderCallbacks,
                  void **headerCBData,
                  bitInsertFunc biFunc, BitOffset cwBitsPerPos,
                  BitOffset maxBitsPerPos, void *cbState, Env *env);

/**
 * \brief Load previously written block encoded sequence
 * representation.
 * @param projectName base name of corresponding suffixerator project
 * @param env genometools reference for core functions
 */
extern EISeq *
loadBlockEncIdxSeq(const Str *projectName, int features, Env *env);

/**
 * \brief Deallocate a previously loaded/created sequence object.
 * @param seq reference of object to delete
 * @param env genometools reference for core functions
 */
extern void
deleteEncIdxSeq(EISeq *seq, Env *env);

/**
 * \brief Retrieve alphabet transformation from sequence object
 * @param seq reference of object to query for alphabet
 * @return read-only reference of alphabet associated with sequence
 */
static inline const MRAEnc *
EISGetAlphabet(const EISeq *seq);

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
static inline Seqpos
EISRank(EISeq *seq, Symbol sym, Seqpos pos, union EISHint *hint,
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
static inline Seqpos
EISSymTransformedRank(EISeq *seq, Symbol msym, Seqpos pos,
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
static inline void
EISRetrieveExtraBits(EISeq *seq, Seqpos pos, int flags,
                     struct extBitsRetrieval *retval, union EISHint *hint,
                     Env *env);

/**
 * \brief Initialize empty structure to hold retrieval results later.
 * @param r struct to initialize
 * @param env genometools state, passes information about allocator etc
 */
static inline void
initExtBitsRetrieval(struct extBitsRetrieval *r, Env *env);

/**
 * \brief Allocate and initialize empty structure to hold retrieval
 * results later.
 * @param env genometools state, passes information about allocator
 * etc.
 * @return reference of new retrieval object
 */
static inline struct extBitsRetrieval *
newExtBitsRetrieval(Env *env);

/**
 * \brief Destruct structure holding retrieval data, deallocates
 * referenced data as necessary.
 * @param r struct to destruct
 * @param env genometools state, passes information about allocator etc.
 */
static inline void
destructExtBitsRetrieval(struct extBitsRetrieval *r, Env *env);

/**
 * \brief Destruct structure holding retrieval data, deallocates
 * referenced data as necessary and the structure itself.
 * @param r struct to delete
 * @param env genometools state, passes information about allocator etc.
 */
static inline void
deleteExtBitsRetrieval(struct extBitsRetrieval *r, Env *env);

/**
 * \brief Find positions of nth symbol occurrence. TODO: NOT IMPLEMENTED
 */
extern Seqpos
EISSelect(EISeq *seq, Symbol sym, Seqpos count);

/**
 * \brief Query length of stored sequence.
 * @param seq indexed sequence object to query
 * @return length of sequence
 */
static inline Seqpos
EISLength(const EISeq *seq);

/**
 * \brief Return symbol at specified position. Comparable to c[pos] if the
 * sequence was stored straight in an array c.
 * @param seq indexed sequence object to be queried
 * @param pos position to retrieve symbol for
 * @param hint optional caching/hinting structure (improves average
 * retrieval time)
 * @param env genometools state, passes information about allocator etc
 * @return value of symbol as encoded in original alphabet
 */
static inline Symbol
EISGetSym(EISeq *seq, Seqpos pos, EISHint hint, Env *env);

/**
 * \brief Return symbol at specified position. Comparable to c[pos] if the
 * sequence was stored straight in an array c. Returns value of
 * encoding, i.e. does not retranslate to original alphabet.
 * @param seq indexed sequence object to be queried
 * @param pos position to retrieve symbol for
 * @param hint optional caching/hinting structure (improves average
 * retrieval time)
 * @param env genometools state, passes information about allocator etc
 * @return value of symbol as processed by encoding alphabet
 */
static inline Symbol
EISGetTransformedSym(EISeq *seq, Seqpos pos, EISHint hint, Env *env);

/**
 * \brief Construct new hinting structure to accelerate operations on
 * related positions.
 * @param seq reference of sequence object to use
 * @param env genometools state, passes information about allocator etc

 */
static inline EISHint
newEISHint(EISeq *seq, Env *env);

static inline void
deleteEISHint(EISeq *seq, EISHint hint, Env *env);

enum EISIntegrityCheckResults
{
  EIS_INTEGRITY_CHECK_NO_ERROR = 0,
  EIS_INTEGRITY_CHECK_INVALID_SYMBOL, /**< reading from the index
                                       *   produced an incorrect symbol  */
  EIS_INTEGRITY_CHECK_BWT_READ_ERROR, /**< reading from the BWT
                                       *   reference file failed or
                                       *   delivered a symbol not in
                                       *   the alphabet (bwt file
                                       *   corrupt?) */
  EIS_INTEGRITY_CHECK_SA_LOAD_ERROR,  /**< loading/mapping of the
                                       *   suffix array project failed
                                       *   (did you generate the BWT) */
  EIS_INTEGRITY_CHECK_RANK_FAILED     /**< the rank operation
                                       *   delivered a wrong count */
};

extern const char *EISIntegrityCheckResultStrings[];

enum EISIntegrityCheckFlags
{
  EIS_VERIFY_BASIC = 0,
  EIS_VERIFY_EXT_RANK = 1 << 0,
};

extern enum EISIntegrityCheckResults
EISVerifyIntegrity(EISeq *seqIdx, const Str *projectName, Seqpos skip,
                   unsigned long tickPrint, FILE *fp, int chkFlags, Env *env);

static inline FILE *
EISSeekToHeader(const EISeq *seqIdx, uint16_t headerID,
                uint32_t *lenRet);

static inline int
EISPrintDiagsForPos(const EISeq *seqIdx, Seqpos pos, FILE *fp, EISHint hint,
                    Env *env);

#include "libgtmatch/eis-encidxseqsimpleop.h"

#endif
