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

#ifndef EIS_ENCIDXSEQCONSTRUCT_H
#define EIS_ENCIDXSEQCONSTRUCT_H

/**
 * \file eis-encidxseqconstruct.h
 * \brief Methods to construct encoded indexed sequence objects from a
 * variety of sources.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */

#include "libgtcore/error.h"
#include "libgtcore/str.h"
#include "libgtmatch/eis-encidxseq.h"
#include "libgtmatch/eis-suffixerator-interface.h"
#include "libgtmatch/eis-suffixarray-interface.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/seqpos-def.h"

/**
 * \brief Construct block-encoded indexed sequence object and write
 * corresponding representation to disk.
 * @param si reference of suffixerator interface used to read
 * construction from
 * @param totalLen length of indexed sequence (including terminator
 * and separators)
 * @param projectName base name of corresponding suffixerator project
 * @param params parameters describing the composition of blocks,
 * buckets and construction of in memory structures
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
 * @param cwExtBitsPerPos exactly this many bits will be appended by
 * biFunc for each symbol of the input sequence
 * @param maxVarExtBitsPerPos at most this many bits will be appended to the
 * variable width part of the data
 * @param cbState will be passed on each call of biFunc
 * @param err genometools reference for core functions
 */
extern EISeq *
newBlockEncIdxSeqFromSfxI(sfxInterface *si, Seqpos totalLen,
                          const Str *projectName,
                          const struct blockEncParams *params,
                          size_t numExtHeaders, uint16_t *headerIDs,
                          uint32_t *extHeaderSizes,
                          headerWriteFunc *extHeaderCallbacks,
                          void **headerCBData, bitInsertFunc biFunc,
                          BitOffset cwExtBitsPerPos,
                          varExtBitsEstimator biVarBits, void *cbState,
                          Error *err);

/**
 * \brief Construct block-encoded indexed sequence object and write
 * corresponding representation to disk.
 * @param sa reference of suffix-array data structure to read
 * construction from
 * @param totalLen length of indexed sequence (including terminator
 * and separators)
 * @param projectName base name of corresponding suffixerator project
 * @param params parameters describing the composition of blocks,
 * buckets and construction of in memory structures
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
 * @param cwExtBitsPerPos exactly this many bits will be appended by
 * biFunc for each symbol of the input sequence
 * @param maxVarExtBitsPerPos at most this many bits will be appended to the
 * variable width part of the data
 * @param cbState will be passed on each call of biFunc
 * @param err genometools reference for core functions
 */
extern EISeq *
newBlockEncIdxSeqFromSA(Suffixarray *sa, Seqpos totalLen,
                        const Str *projectName,
                        const struct blockEncParams *params,
                        size_t numExtHeaders, uint16_t *headerIDs,
                        uint32_t *extHeaderSizes,
                        headerWriteFunc *extHeaderCallbacks,
                        void **headerCBData,
                        bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                        varExtBitsEstimator biVarBits, void *cbState,
                        Error *err);

/**
 * \brief Construct block-encoded indexed sequence object and write
 * corresponding representation to disk.
 * @param sai reference of state struct holding reference to
 * suffix-array data structure to read construction from and progress
 * information
 * @param alphabet encoding to use for the built sequence, "ownership"
 * of alphabet will pass to the returned object if no error occured.
 * @param totalLen length of indexed sequence (including terminator
 * and separators)
 * @param projectName base name of corresponding suffixerator project
 * @param params parameters describing the composition of blocks,
 * buckets and construction of in memory structures
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
 * @param cwExtBitsPerPos exactly this many bits will be appended by
 * biFunc for each symbol of the input sequence
 * @param maxVarExtBitsPerPos at most this many bits will be appended to the
 * variable width part of the data
 * @param cbState will be passed on each call of biFunc
 * @param err genometools reference for core functions
 */
extern EISeq *
newBlockEncIdxSeqFromSAI(SuffixarrayFileInterface *sai,
                         Seqpos totalLen, const Str *projectName,
                         const struct blockEncParams *params,
                         size_t numExtHeaders, uint16_t *headerIDs,
                         uint32_t *extHeaderSizes,
                         headerWriteFunc *extHeaderCallbacks,
                         void **headerCBData,
                         bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                         varExtBitsEstimator biVarBits, void *cbState,
                         Error *err);

/**
 * \brief Load previously written block encoded sequence
 * representation for already opened suffix-array structure.
 * @param sa reference of suffix-array data structure to read
 * construction from
 * @param totalLen length of indexed sequence (including terminator
 * and separators)
 * @param features select optional in-memory structures
 * @param projectName base name of corresponding suffixerator project
 * @param err genometools reference for core functions
 */
struct encIdxSeq *
loadBlockEncIdxSeqForSA(const Suffixarray *sa, Seqpos totalLen,
                        const Str *projectName, int features, Error *err);

#endif
