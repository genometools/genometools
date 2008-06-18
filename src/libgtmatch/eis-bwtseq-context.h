/*
  Copyright (c) 2008 Thomas Jahns <Thomas.Jahns@gmx.net>

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
#ifndef EIS_BWTSEQ_CONTEXT_H
#define EIS_BWTSEQ_CONTEXT_H

/**
 * @file eis-bwtseq-context.h
 * @brief interface for generating arbitrary context from a
 * packedindex object
 */

#include "libgtmatch/eis-seqdatasrc.h"
#include "libgtmatch/seqpos-def.h"
#include "libgtmatch/eis-bwtseq-context-param.h"

typedef struct BWTSeqContextRetriever BWTSeqContextRetriever;

typedef struct BWTSeqContextRetrieverFactory BWTSeqContextRetrieverFactory;

/**
 * @param seqLen length of sequence to build context map for
 * @param mapIntervalLog2 unless equal to CTX_MAP_ILOG_NOMAP or
 * CTX_MAP_ILOG_AUTOSIZE compute map with interval 1<<mapIntervalLog2.
 * If equal CTX_MAP_ILOG_AUTOSIZE, uses @f$log_2\mathrm{seqLen}@f$.
 */
extern BWTSeqContextRetrieverFactory *
newBWTSeqContextRetrieverFactory(Seqpos seqLen,
                                 short mapIntervalLog2);

extern void
deleteBWTSeqContextRetrieverFactory(BWTSeqContextRetrieverFactory *factory);

extern Seqpos
BWTSCRFReadAdvance(BWTSeqContextRetrieverFactory *factory, Seqpos chunkSize,
                   SeqDataReader readSfxIdx, Error *err);

extern size_t
BWTSCRFMapAdvance(BWTSeqContextRetrieverFactory *factory, const Seqpos *src,
                  size_t chunkSize);

extern bool
BWTSCRFFinished(const BWTSeqContextRetrieverFactory *factory);

/**
 * @param bwtSeq reference to companion BWT sequence object, can be
 * NULL if the retriever object is not used for queries
 * (i.e. immediately destructed).
 */
extern BWTSeqContextRetriever *
BWTSCRFGet(BWTSeqContextRetrieverFactory *factory, const BWTSeq *bwtSeq,
           const Str *projectName);

/**
 * @brief Load context retriever that can be used to generate arbitrary
 * subsequences of the original sequence.
 * @param projectName base filename for which to load the retriever
 * table
 * @param bwtSeq reference to previously loaded BWT index object
 * @param mapIntervalLog2 unless equal to CTX_MAP_ILOG_AUTOSIZE
 * load map with interval 1<<mapIntervalLog2.  If equal
 * CTX_MAP_ILOG_AUTOSIZE, uses largest table that can be successfully
 * mapped.
 * @return NULL in case of error (e.g. no corresponding table file or
 * not enought memory for mmap).
 */
extern BWTSeqContextRetriever *
BWTSeqCRLoad(const BWTSeq *bwtSeq, const Str *projectName,
             short mapIntervalLog2);

extern void
deleteBWTSeqCR(BWTSeqContextRetriever *bwtSeqCR);

struct SeqMark
{
  Seqpos textPos, bwtPos;
};

/**
 * Writes the retrieved symbols to subseq which must accomodate for
 * len or more symbols.
 */
extern void
BWTSeqCRAccessSubseq(const BWTSeqContextRetriever *bwtSeqCR,
                     Seqpos start, size_t len, Symbol subseq[]);

/**
 * @brief Compute next position in original sequence following pos
 * that is marked for efficient retrieval.
 */
static inline struct SeqMark
BWTSeqCRNextMark(const BWTSeqContextRetriever *bwtSeqCR, Seqpos pos);

#include "libgtmatch/eis-bwtseq-context-siop.h"

#endif
