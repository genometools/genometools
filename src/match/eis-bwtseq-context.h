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

#include "match/eis-seqdatasrc.h"

#include "match/eis-bwtseq-context-param.h"

typedef struct BWTSeqContextRetriever BWTSeqContextRetriever;

typedef struct BWTSeqContextRetrieverFactory BWTSeqContextRetrieverFactory;

/**
 * @param seqLen length of sequence to build context map for
 * @param mapIntervalLog2 unless equal to CTX_MAP_ILOG_NOMAP or
 * CTX_MAP_ILOG_AUTOSIZE compute map with interval 1<<mapIntervalLog2.
 * If equal CTX_MAP_ILOG_AUTOSIZE, uses @f$gt_log_2\mathrm{seqLen}@f$.
 */
BWTSeqContextRetrieverFactory *
gt_newBWTSeqContextRetrieverFactory(unsigned long seqLen,
                                 short mapIntervalLog2);

void
gt_deleteBWTSeqContextRetrieverFactory(BWTSeqContextRetrieverFactory *factory);

unsigned long
gt_BWTSCRFReadAdvance(BWTSeqContextRetrieverFactory *factory,
                   unsigned long chunkSize,
                   SeqDataReader readSfxIdx);

size_t
gt_BWTSCRFMapAdvance(BWTSeqContextRetrieverFactory *factory,
                  const unsigned long *src,
                  size_t chunkSize);

bool
gt_BWTSCRFFinished(const BWTSeqContextRetrieverFactory *factory);

/**
 * @param bwtSeq reference to companion BWT sequence object, can be
 * NULL if the retriever object is not used for queries
 * (i.e. immediately destructed).
 */
BWTSeqContextRetriever *
gt_BWTSCRFGet(BWTSeqContextRetrieverFactory *factory, const BWTSeq *bwtSeq,
           const char *projectName);

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
BWTSeqContextRetriever *
gt_BWTSeqCRLoad(const BWTSeq *bwtSeq, const char *projectName,
             short mapIntervalLog2);

void
gt_deleteBWTSeqCR(BWTSeqContextRetriever *bwtSeqCR);

struct SeqMark
{
  unsigned long textPos, bwtPos;
};

/**
 * Writes the retrieved symbols to subseq which must accomodate for
 * len or more symbols.
 */
void
gt_BWTSeqCRAccessSubseq(const BWTSeqContextRetriever *bwtSeqCR,
                     unsigned long start, size_t len, Symbol subseq[]);

/**
 * @brief Compute next position in original sequence following pos
 * that is marked for efficient retrieval.
 */
static inline struct SeqMark
BWTSeqCRNextMark(const BWTSeqContextRetriever *bwtSeqCR, unsigned long pos);

#include "match/eis-bwtseq-context-siop.h"

#endif
