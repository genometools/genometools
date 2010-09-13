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
#include "core/encseq.h"
#include "match/eis-sa-common.h"

size_t gt_translateSuftab2BWT(void *translator,
                              void *voiddest,
                              const unsigned long *src,
                              size_t len)
{
  struct encSeqTrState *trState = (struct encSeqTrState *) translator;
  GtUchar *dest = (GtUchar *) voiddest;
  size_t idx;

  gt_assert(trState);
  for (idx = 0; idx < len; ++idx)
  {
    dest[idx] = sfxIdx2BWTSym(src[idx], trState->encseq,
                              trState->readmode);
  }
  return len * sizeof (GtUchar);
}

size_t gt_translateSuftab2BWTSuffixsortspace(
                                       void *translator,
                                       void *voiddest,
                                       const GtSuffixsortspace *suffixsortspace,
                                       unsigned long offset,
                                       size_t len)
{
  struct encSeqTrState *trState = (struct encSeqTrState *) translator;
  GtUchar *dest = (GtUchar *) voiddest;
  size_t idx;

  gt_assert(trState);
  for (idx = 0; idx < len; ++idx)
  {
    dest[idx]
      = sfxIdx2BWTSym(gt_suffixsortspace_getdirect(suffixsortspace,offset+idx),
                      trState->encseq, trState->readmode);
  }
  return len * sizeof (GtUchar);
}

struct saTaggedXltorStateLE
{
  struct saTaggedXltorStateLE *next;
  struct saTaggedXltorState state;
};

void
gt_initSATaggedXltorStateList(struct saTaggedXltorStateList *saXltorStateList)
{
  saXltorStateList->numXltors = 0;
  saXltorStateList->stateList = NULL;
}

void
gt_destructSATaggedXltorStateList(
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

struct saTaggedXltorState *
gt_addSuffixarrayXltor(struct saTaggedXltorStateList *saXltorStateList,
                    enum sfxDataRequest request,
                    struct encSeqTrState state)
{
  struct saTaggedXltorStateLE *newSAXltorState;
  newSAXltorState = gt_malloc(sizeof (*newSAXltorState));
  newSAXltorState->state.typeTag = request;
  newSAXltorState->state.state = state;
  newSAXltorState->next = saXltorStateList->stateList;
  saXltorStateList->stateList = newSAXltorState;
  ++(saXltorStateList->numXltors);
  return &newSAXltorState->state;
}

SeqDataReader
gt_SASSGenericCreateReader(SASeqSrc *src, enum sfxDataRequest request)
{
  return gt_seqReaderSetRegisterConsumer(&src->readerSet, request,
                                      src->createTranslator(src, request));
}
