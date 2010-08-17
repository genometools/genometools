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

/**
 * \file eis-suffixerator-interface.h
 * \brief Methods to call suffixerator functions through one object,
 * but have the same data available to multiple listeners.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */
#ifndef EIS_SA_COMMON_H
#define EIS_SA_COMMON_H

#include <stdlib.h>
#include "core/defined-types.h"
#include "core/encseq.h"
#include "core/types_api.h"
#include "match/eis-mrangealphabet.h"
#include "match/eis-random-seqaccess.h"
#include "match/eis-seqdatasrc.h"
#include "match/sfx-suffixgetset.h"

/**
 * Describes what kind of information will be read by a requestor:
 */
enum sfxDataRequest {
  SFX_REQUEST_NONE = 0,         /**< empty request, used for special purposes */
  SFX_REQUEST_SUFTAB = 1<<0,    /**< request for suffix array entries */
  SFX_REQUEST_BWTTAB = 1<<1,    /**< request for bwt table */
  SFX_REQUEST_ALL = SFX_REQUEST_SUFTAB |
                    SFX_REQUEST_BWTTAB, /**< used as bitmask  */
  SFX_REQUEST_ANY = SFX_REQUEST_ALL,             /**< used as bitmask  */
};

/**
 * @return Symbol at position sufIdx or UNDEFBWTCHAR i.e. the terminator
 */
static inline GtUchar
sfxIdx2BWTSym(unsigned long sufIdx, const GtEncseq *encseq,
              GtReadmode readmode);

static inline size_t
EncSeqGetSubSeq(const GtEncseq *encseq, GtReadmode readmode,
                unsigned long pos, size_t len, GtUchar *subStr);

struct encSeqTrState
{
  const GtEncseq *encseq;
  GtReadmode readmode;
};

/**
 * @brief Produce given length of symbols from the BWT by translating
 * an array of suffix indices, this version uses the values of
 * suffix array and encoded sequence instead of reading the BWT file.
 * @param state reference of a SuffixarrayFileInterface
 * @param dest write symbols here
 * @param src read suffix indices from here
 * @param len length of string to read
 */

size_t gt_translateSuftab2BWT(void *data,
                              void *voiddest,
                              const unsigned long *src,
                              size_t len);

size_t gt_translateSuftab2BWTSuffixsortspace(
                                       void *translator,
                                       void *voiddest,
                                       const GtSuffixsortspace *suffixsortspace,
                                       unsigned long offset,
                                       size_t len);

struct saTaggedXltorState
{
  enum sfxDataRequest typeTag;
  struct encSeqTrState state;
};

struct saTaggedXltorStateList
{
  size_t numXltors;
  struct saTaggedXltorStateLE *stateList;
};

void
gt_initSATaggedXltorStateList(struct saTaggedXltorStateList *saXltorStateList);

void
gt_destructSATaggedXltorStateList(
  struct saTaggedXltorStateList *saXltorStateList);

struct saTaggedXltorState *
gt_addSuffixarrayXltor(struct saTaggedXltorStateList *saXltorStateList,
                    enum sfxDataRequest request,
                    struct encSeqTrState state);

typedef struct SASeqSrc SASeqSrc;

static inline SeqDataReader
SASSCreateReader(SASeqSrc *src, enum sfxDataRequest request);

static inline Definedunsignedlong
SASSGetRot0Pos(const SASeqSrc *src);

static inline unsigned long
SASSGetLength(const SASeqSrc *src);

static inline MRAEnc *
SASSNewMRAEnc(const SASeqSrc *src);

static inline const MRAEnc *
SASSGetMRAEnc(SASeqSrc *src);

static inline const struct seqStats *
SASSGetSeqStats(const SASeqSrc *src);

static inline size_t
SASSAccessSequence(const SASeqSrc *src,
                   Symbol *dest,
                   unsigned long pos,
                   size_t len);

static inline RandomSeqAccessor
SASSGetOrigSeqAccessor(const SASeqSrc *src);

static inline void
SASSDelete(SASeqSrc *src);

#include "match/eis-sa-common-siop.h"

#endif
