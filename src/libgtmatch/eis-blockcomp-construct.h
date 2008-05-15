/*
  Copyright (C) 2008 Thomas Jahns <Thomas.Jahns@gmx.net>

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

#ifndef EIS_BLOCKCOMP_CONSTRUCT_H
#define EIS_BLOCKCOMP_CONSTRUCT_H

#include "libgtmatch/eis-encidxseq.h"
#include "libgtmatch/sarr-def.h"

/**
 * @param alphabet ownership of alphabet is transferred to the sequence
 * index produced unless NULL is returned
 */
extern EISeq *
newGenBlockEncIdxSeq(Seqpos totalLen, const Str *projectName,
                     MRAEnc *alphabet, const struct seqStats *stats,
                     SeqDataReader BWTGenerator,
                     const struct seqBaseParam *params,
                     size_t numExtHeaders, const uint16_t *headerIDs,
                     const uint32_t *extHeaderSizes,
                     headerWriteFunc *extHeaderCallbacks,
                     void **headerCBData,
                     bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                     varExtBitsEstimator biVarBits, void *cbState, Error *err);

/**
 * @brief Load previously written block encoded sequence
 * representation.
 * @param alphabet
 * @param totalLen
 * Caution: ownership of alphabet changes to returned index.
 * @param projectName base name of corresponding suffixerator project
 * @param features select optional in-memory data structures for speed-up
 * @param err genometools error object reference
 * @return new encoded indexed sequence object reference
 */
extern EISeq *
loadBlockEncIdxSeqGen(MRAEnc *alphabet, Seqpos totalLen,
                      const Str *projectName, int features, Error *err);

#endif
