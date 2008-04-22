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

#include "libgtcore/symboldef.h"
#include "libgtmatch/encseq-def.h"
#include "libgtmatch/seqpos-def.h"

/**
 * Describes what kind of information will be read by a requestor:
 */
enum sfxDataRequest {
  SFX_REQUEST_NONE = 0,         /**< empty request, used for special purposes */
  SFX_REQUEST_SUFTAB = 1<<0,    /**< request for suffix array entries */
  SFX_REQUEST_LCPTAB = 1<<1,    /**< request for lcp table entries */
  SFX_REQUEST_BWTTAB = 1<<2,    /**< request for bwt table */
  SFX_REQUEST_ALL = SFX_REQUEST_SUFTAB | SFX_REQUEST_LCPTAB
                  | SFX_REQUEST_BWTTAB,          /**< used as bitmask  */
  SFX_REQUEST_ANY = SFX_REQUEST_ALL,             /**< used as bitmask  */
};

static inline Uchar
sfxIdx2BWTSym(Seqpos sufIdx, const Encodedsequence *encseq, Readmode readmode)
{
  return sufIdx != 0
    ? getencodedchar(encseq, sufIdx - 1, readmode)
    : (Uchar) UNDEFBWTCHAR;
}

struct encSeqTrState
{
  const Encodedsequence *encseq;
  Readmode readmode;
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
extern size_t
translateSuftab2BWT(struct encSeqTrState *trState, Uchar *dest, Seqpos *src,
                    size_t len);

struct encSeqLCPState
{
  Seqpos lastSufIdx;
  const Encodedsequence *encseq;
  Readmode readmode;
};

extern size_t
translateSuftab2LCP(struct encSeqLCPState *lcpState, Seqpos *dest, Seqpos *src,
                    size_t len);

union saXltorState
{
  struct encSeqTrState encSeqTr;
  struct encSeqLCPState lcpState;
};

struct saTaggedXltorState
{
  enum sfxDataRequest typeTag;
  union saXltorState state;
};

struct saTaggedXltorStateList
{
  size_t numXltors;
  struct saTaggedXltorStateLE *stateList;
};

extern void
initSATaggedXltorStateList(struct saTaggedXltorStateList *saXltorStateList);

extern void
destructSATaggedXltorStateList(
  struct saTaggedXltorStateList *saXltorStateList);

extern struct saTaggedXltorState *
addSuffixarrayXltor(struct saTaggedXltorStateList *saXltorStateList,
                    enum sfxDataRequest request,
                    union saXltorState saXltorState);

#endif
