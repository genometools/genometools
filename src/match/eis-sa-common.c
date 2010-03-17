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

#include <stdlib.h>

#include "core/ma_api.h"
#include "match/eis-sa-common.h"
#include "match/encodedsequence.h"

extern size_t
translateSuftab2BWT(struct encSeqTrState *trState, GtUchar *dest, Seqpos *src,
                    size_t len)
{
  size_t i;
  gt_assert(trState);
  for (i = 0; i < len; ++i)
  {
    dest[i] = sfxIdx2BWTSym(src[i], trState->encseq, trState->readmode);
  }
  return len * sizeof (GtUchar);
}

static inline void
writeLCPVal(const GtEncodedsequence *encseq, GtReadmode readmode,
            Seqpos *dest, Seqpos a, Seqpos b)
{
#ifndef NDEBUG
  int cmp =
#endif /* NDEBUG */
    comparetwosuffixes(encseq,
                       readmode,
                       dest,
                       false,
                       false,
                       0,
                       a,
                       b,
                       NULL,  /* XXX FIXME */
                       NULL);  /* XXX FIXME */
#ifndef NDEBUG
  if (cmp > 0)
  {
    fprintf(stderr, ": cmp " FormatSeqpos " " FormatSeqpos " = %d",
            PRINTSeqposcast(a), PRINTSeqposcast(b), cmp);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
#endif /* NDEBUG */
}

extern size_t
translateSuftab2LCP(struct encSeqLCPState *lcpState, Seqpos *dest, Seqpos *src,
                    size_t len)
{
  size_t elemsLeft = len;
  Seqpos lastSufIdx;
  gt_assert(lcpState && dest && src);
  lastSufIdx = lcpState->lastSufIdx;
  if (elemsLeft)
  {
    Seqpos nextSufIdx = *src;
    if (lastSufIdx != -1)
    {
      writeLCPVal(lcpState->encseq, lcpState->readmode, dest,
                  lastSufIdx, nextSufIdx);
    }
    else
    {
      *dest = 0;
    }
    while (--elemsLeft)
    {
      lastSufIdx = nextSufIdx;
      nextSufIdx = *(++src);
      writeLCPVal(lcpState->encseq, lcpState->readmode, ++dest,
                  lastSufIdx, nextSufIdx);
    }
    lcpState->lastSufIdx = nextSufIdx;
  }
  return len * sizeof (dest[0]);
}

struct saTaggedXltorStateLE
{
  struct saTaggedXltorStateLE *next;
  struct saTaggedXltorState state;
};

extern void
initSATaggedXltorStateList(struct saTaggedXltorStateList *saXltorStateList)
{
  saXltorStateList->numXltors = 0;
  saXltorStateList->stateList = NULL;
}

extern void
destructSATaggedXltorStateList(
  struct saTaggedXltorStateList *saXltorStateList)
{
  struct saTaggedXltorStateLE *next;
  gt_assert(saXltorStateList);
  next = saXltorStateList->stateList;
  while (next)
  {
    struct saTaggedXltorStateLE *prev = next;
    next = next->next;
    gt_free(prev);
  }
}

extern struct saTaggedXltorState *
addSuffixarrayXltor(struct saTaggedXltorStateList *saXltorStateList,
                    enum sfxDataRequest request,
                    union saXltorState saXltorState)
{
  struct saTaggedXltorStateLE *newSAXltorState;
  newSAXltorState = gt_malloc(sizeof (*newSAXltorState));
  newSAXltorState->state.typeTag = request;
  newSAXltorState->state.state = saXltorState;
  newSAXltorState->next = saXltorStateList->stateList;
  saXltorStateList->stateList = newSAXltorState;
  ++(saXltorStateList->numXltors);
  return &newSAXltorState->state;
}

extern SeqDataReader
SASSGenericCreateReader(SASeqSrc *src, enum sfxDataRequest request)
{
  return seqReaderSetRegisterConsumer(&src->readerSet, request,
                                      src->createTranslator(src, request));
}
