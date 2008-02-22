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
 * \file eis-encidxseq.h
 * Interface definitions for encoded indexed sequences.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */
#include <stdio.h>
#include <inttypes.h>

#include "libgtcore/bitpackstring.h"
#include "libgtcore/error.h"
#include "libgtcore/str.h"
#include "libgtmatch/seqpos-def.h"
#include "libgtmatch/eis-mrangealphabet.h"
#include "libgtmatch/eis-encidxseqparam.h"

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
 * @param err genometools error object reference
 * @return number of bits actually written, or (BitOffset)-1 if an
 * error occured
 */
typedef BitOffset (*bitInsertFunc)(BitString cwDest, BitOffset cwOffset,
                                   BitString varDest, BitOffset varOffset,
                                   Seqpos start, Seqpos len, void *cbState,
                                   Error *err);

/**
 * Callback to insert one header field. The data written can later be
 * retrieved via EISSeekToHeader.
 */
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

/**
 * Select which parts of the stored data to retrieve: only constant width
 * part or also from the variable part. Join with bitwise or to form a
 * request for extBitsRetrieval.
 */
enum extBitsRetrievalFlags
{
  EBRF_RETRIEVE_CWBITS = 0,        /**< retrieve the constant-width part */
  EBRF_PERSISTENT_CWBITS = 1<<0,   /**< make constant width string persistent,
                                    * i.e. survive further calls to
                                    * query functions */
  EBRF_RETRIEVE_VARBITS = 1<<1,    /**< retrieve variable-width string */
  EBRF_PERSISTENT_VARBITS = 1<<2,  /**< make variable-width string persistent,
                                    * i.e. survive further calls to
                                    * query functions */
};

/**
 * Stores the number of occurrences for every symbol of the input alphabet.
 */
struct seqStats
{
  Seqpos *symbolDistributionTable;          /**< indexed by symbol value */
  enum sourceEncType sourceAlphaType;       /**< defines alphabet width */
};

/** holds encoded indexed sequence object */
typedef struct encIdxSeq EISeq;
/** hints to speed up retrievals when accessing positions in sequence
 * (or at least close to one another) */
typedef union EISHint *EISHint;

/**
 * \brief Construct block-encoded indexed sequence object and write
 * corresponding representation to disk.
 * @param projectName base name of corresponding suffixerator project
 * @param params parameters for index construction
 * @param numExtHeaders number of extension headers to write via callbacks
 * @param headerIDs array of numExtHeaders ids to be used
 * for each extension header in turn
 * @param extHeaderSizes array of numExtHeaders sizes
 * representing the length of each extension header
 * @param extHeaderCallbacks array of numExtHeaders function pointers
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
 * @param err genometools error object reference
 */
extern EISeq *
newBlockEncIdxSeq(const Str *projectName, const struct blockEncParams *params,
                  size_t numExtHeaders, uint16_t *headerIDs,
                  uint32_t *extHeaderSizes, headerWriteFunc *extHeaderCallbacks,
                  void **headerCBData,
                  bitInsertFunc biFunc, BitOffset cwBitsPerPos,
                  BitOffset maxBitsPerPos, void *cbState, Error *err);

/**
 * \brief Load previously written block encoded sequence
 * representation.
 * @param projectName base name of corresponding suffixerator project
 * @param features select optional in-memory data structures for speed-up
 * @param err genometools error object reference
 */
extern EISeq *
loadBlockEncIdxSeq(const Str *projectName, int features, Error *err);

/**
 * \brief Deallocate a previously loaded/created sequence object.
 * @param seq reference of object to delete
 */
extern void
deleteEncIdxSeq(EISeq *seq);

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
 */
static inline Seqpos
EISRank(EISeq *seq, Symbol sym, Seqpos pos, union EISHint *hint);

/**
 * \brief Return number of occurrences of symbol tSym in index up to
 * but not including given position.
 * @param seq sequence index object to query
 * @param tSym symbol to query occurrences of, but already transformed
 * by input alphabet
 * @param pos occurences are counted up to this position
 * @param hint provides cache and direction information for queries
 * based on previous queries
 */
static inline Seqpos
EISSymTransformedRank(EISeq *seq, Symbol tSym, Seqpos pos,
                      union EISHint *hint);

struct SeqposPair
{
  Seqpos a,b;
};

/**
 * \brief Return number of occurrences of symbol sym in index up to
 * but not including given positions.
 *
 * precondition: posA <= posB
 * @param seq sequence index object to query
 * @param sym original alphabet symbol to query occurrences of
 * @param posA occurences are counted up to this position
 * @param posB as for posA occurences are counted up to this position
 * @param hint provides cache and direction information for queries
 * based on previous queries
 * @return members a and b of returned struct contain Occ results for
 * posA and posB respectively
 */
static inline struct SeqposPair
EISPosPairRank(EISeq *seq, Symbol sym, Seqpos posA, Seqpos posB,
               union EISHint *hint);

/**
 * \brief Return number of occurrences of symbol sym in index up to
 * but not including given positions.
 *
 * precondition: posA <= posB
 * @param seq sequence index object to query
 * @param tSym symbol to query occurrences of, but already transformed
 * by input alphabet
 * @param posA occurences are counted up to this position
 * @param posB as for posA occurences are counted up to this position
 * @param hint provides cache and direction information for queries
 * based on previous queries
 * @return members a and b of returned struct contain Occ results for
 * posA and posB respectively
 */
static inline struct SeqposPair
EISSymTransformedPosPairRank(EISeq *seq, Symbol tSym, Seqpos posA, Seqpos posB,
                             union EISHint *hint);

/**
 * Presents the bits previously stored by a bitInsertFunc callback.
 * @param seq
 * @param pos sequence position for which to retrieve corresponding
 * area
 * @param flags select which part of the extension bits to query
 * @param retval store information of retrieved bits here
 * @param hint provides cache and direction information for queries
 */
static inline void
EISRetrieveExtraBits(EISeq *seq, Seqpos pos, int flags,
                     struct extBitsRetrieval *retval, union EISHint *hint);

/**
 * \brief Initialize empty structure to hold retrieval results later.
 * @param r struct to initialize
 */
static inline void
initExtBitsRetrieval(struct extBitsRetrieval *r);

/**
 * \brief Allocate and initialize empty structure to hold retrieval
 * results later.
 * @return reference of new retrieval object
 */
static inline struct extBitsRetrieval *
newExtBitsRetrieval();

/**
 * \brief Destruct structure holding retrieval data, deallocates
 * referenced data as necessary.
 * @param r struct to destruct
 */
static inline void
destructExtBitsRetrieval(struct extBitsRetrieval *r);

/**
 * \brief Destruct structure holding retrieval data, deallocates
 * referenced data as necessary and the structure itself.
 * @param r struct to delete
 */
static inline void
deleteExtBitsRetrieval(struct extBitsRetrieval *r);

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
 * @return value of symbol as encoded in original alphabet
 */
static inline Symbol
EISGetSym(EISeq *seq, Seqpos pos, EISHint hint);

/**
 * \brief Return symbol at specified position. Comparable to c[pos] if the
 * sequence was stored straight in an array c. Returns value of
 * encoding, i.e. does not retranslate to original alphabet.
 * @param seq indexed sequence object to be queried
 * @param pos position to retrieve symbol for
 * @param hint optional caching/hinting structure (improves average
 * retrieval time)
 * @return value of symbol as processed by encoding alphabet
 */
static inline Symbol
EISGetTransformedSym(EISeq *seq, Seqpos pos, EISHint hint);

/**
 * \brief Construct new hinting structure to accelerate operations on
 * related positions.
 * @param seq reference of sequence object to use
 * @return new EISHint handle
 */
static inline EISHint
newEISHint(const EISeq *seq);

/**
 * Deallocate hinting data.
 * @param seq sequence associated
 * @param hint hint to free
 */
static inline void
deleteEISHint(EISeq *seq, EISHint hint);

/**
 * Possible outcome of index integrity check.
 */
enum EISIntegrityCheckResults
{
  EIS_INTEGRITY_CHECK_NO_ERROR = 0,   /**< all tests completed
                                       * without error */
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

/** table of error descriptions indexed by return values of
 * EISVerifyIntegrity  */
extern const char *EISIntegrityCheckResultStrings[];

/**
 * Select potentially expensive checks to perform on index.
 */
enum EISIntegrityCheckFlags
{
  EIS_VERIFY_BASIC = 0,         /**< just basic checks  */
  EIS_VERIFY_EXT_RANK = 1 << 0, /**< do rank queries for all
                                 * positions for all symbols in the
                                 * alphabet, not only for ranks changed */
};

/**
 * @brief Check wether index contains same symbols as reference
 * sequence and delivers correct rank counts.
 * @param seqIdx index to check
 * @param projectName name of corresponding suffix array project to
 * use as reference
 * @param skip omit this many symbols at the beginning
 * @param tickPrint print a dot every tickPrint symbols processed
 * @param fp print dots to this file pointer
 * @param chkFlags select additional tests (see enum EISIntegrityCheckFlags)
 */
extern enum EISIntegrityCheckResults
EISVerifyIntegrity(EISeq *seqIdx, const Str *projectName, Seqpos skip,
                   unsigned long tickPrint, FILE *fp, int chkFlags, Error *err);

/**
 * @brief Position file pointer at header written by upper layer.
 * @param seqIdx index to search header in
 * @param headerID number of header to search for
 * @param lenRet write length of header here
 * @return appropriate file pointer or NULL if header was not found
 */
static inline FILE *
EISSeekToHeader(const EISeq *seqIdx, uint16_t headerID,
                uint32_t *lenRet);

/**
 * Given a position write debugging output for surrounding sequence.
 * @param seqIdx sequence index to query
 * @param pos position for which to print context
 * @param fp print diagnostics to this file pointer
 * @param hint use this structure for hinting
 * @return 0 if an I/O error occured wrt fp
 */
static inline int
EISPrintDiagsForPos(const EISeq *seqIdx, Seqpos pos, FILE *fp, EISHint hint);

#include "libgtmatch/eis-encidxseqsimpleop.h"

#endif
