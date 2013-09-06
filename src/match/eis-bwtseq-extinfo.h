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

#ifndef EIS_BWTSEQ_EXTINFO_H
#define EIS_BWTSEQ_EXTINFO_H

/**
 * @file eis-bwtseq-extinfo.h
 * @brief Routines to add and query extension information embedded in
 * the base sequence index.
 * - generic BWT index creation routine
 */

#include "match/eis-specialsrank.h"
#include "match/eis-bwtseq.h"
#include "match/eis-encidxseq.h"
#include "match/eis-mrangealphabet.h"
#include "match/eis-random-seqaccess.h"
#include "match/eis-sa-common.h"
#include "match/eis-seqdatasrc.h"

/** The constructor for the base index must conform to this
 * signature. */
typedef EISeq *(*indexCreateFunc)(
  GtUword totalLen, const char *projectName, MRAEnc *alphabet,
  const struct seqStats *stats, SeqDataReader BWTGenerator,
  const struct seqBaseParam *params, size_t numExtHeaders,
  const uint16_t *headerIDs, const uint32_t *extHeaderSizes,
  headerWriteFunc *extHeaderCallbacks,
  void **headerCBData, bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
  varExtBitsEstimator biVarBits, void *cbState, GtError *err);

/**
 * To enrich a base index with the information required for a BWT
 * index, wrap the base index constructor in a call to gt_createBWTSeqGeneric.
 * @param params holds all parameters for both, the BWT sequence
 * object and the base index
 * @param createIndex wrapped constructor
 * @param src passed to createIndex and used as source of suffix array
 * @param totalLen length of the sorted sequence plus terminator symbol
 * @param alphabet encoding to use for symbols of the input sequence
 * @param specialRanges one value describing for each range of
 * alphabet how symbols in this range are sorted
 * @param readOrigSeq makes the original sequence available
 * @param origSeqState opaque sequence object to pass to readOrigSeq
 * @param readNextUlong stream reader for the values of the suffix
 * array
 * @param spReadState opaque suffix array source object to pass to
 * readNextUlong
 * @param lrepFunc reports the position of the null-rotation
 * @param lrepState passed to lrepFunc
 * @param err
 */
EISeq *
gt_createBWTSeqGeneric(const struct bwtParam *params,
                       indexCreateFunc createIndex,
                       SASeqSrc *src,
                       const enum rangeSortMode rangeSort[],
                       const SpecialsRankLookup *sprTable,
                       GtError *err);

int
gt_BWTSeqPosHasLocateInfo(const BWTSeq *bwtSeq, GtUword pos,
                       struct extBitsRetrieval *extBits);

GtUword
gt_BWTSeqGetRankSort(const BWTSeq *bwtSeq, GtUword pos,
                  AlphabetRangeID range, struct extBitsRetrieval *extBits);

void
gt_BWTSeqInitLocateHandling(BWTSeq *bwtSeq,
                         const enum rangeSortMode *defaultRangeSort);

#endif
