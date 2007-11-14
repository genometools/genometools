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
#include <string.h>

#include "libgtcore/minmax.h"

#include "libgtmatch/eis-seqranges.h"

struct seqRangeList *
newSeqRangeList(size_t rangesStartNum, const MRAEnc *alphabet,
                enum SRLFeatures features, Env *env)
{
  assert(env);
  struct seqRangeList *newList;
  newList = env_ma_malloc(env, sizeof (struct seqRangeList));
  newList->numRanges = 0;
  newList->numRangesStorable = rangesStartNum;
  newList->ranges = env_ma_malloc(env,
                                  sizeof (newList->ranges[0]) * rangesStartNum);
  if (features & SRL_PARTIAL_SYMBOL_SUMS)
    newList->partialSymSums = env_ma_malloc(env, sizeof (Seqpos)
                                            * MRAEncGetSize(alphabet)
                                            * rangesStartNum);
  else
    newList->partialSymSums = NULL;
  newList->alphabet = alphabet;
  return newList;
}

void
SRLCompact(struct seqRangeList *rangeList, Env *env)
{
  assert(rangeList && env);
  rangeList->ranges = env_ma_realloc(env, rangeList->ranges,
                                     sizeof (rangeList->ranges[0])
                                     * rangeList->numRanges);
  if (rangeList->partialSymSums)
    rangeList->partialSymSums =
      env_ma_realloc(env, rangeList->partialSymSums, sizeof (Seqpos)
                     * MRAEncGetSize(rangeList->alphabet)
                     * rangeList->numRanges);
  rangeList->numRangesStorable = rangeList->numRanges;
}

void
deleteSeqRangeList(struct seqRangeList *rangeList, Env *env)
{
  assert(rangeList && env);
  if (rangeList->ranges)
    env_ma_free(rangeList->ranges, env);
  if (rangeList->partialSymSums)
    env_ma_free(rangeList->partialSymSums, env);
  env_ma_free(rangeList, env);
}

void
SRLAppendNewRange(struct seqRangeList *rangeList, Seqpos pos, Seqpos len,
                  Symbol esym, Env *env)
{
  assert(rangeList && env);
  if (len)
  {
    Symbol sym = MRAEncMapSymbol(rangeList->alphabet, esym);
    struct seqRange *p;
    Seqpos *pSums;
    size_t numRanges = rangeList->numRanges,
      numSyms = MRAEncGetSize(rangeList->alphabet),
      numNewRanges = len/MAX_SEQRANGE_LEN + ((len%MAX_SEQRANGE_LEN)?1:0);
    if (numRanges + numNewRanges > rangeList->numRangesStorable)
    {
      size_t newSize = numRanges + 2 * numNewRanges;
      rangeList->ranges = env_ma_realloc(env, rangeList->ranges,
                                         sizeof (struct seqRange) * newSize);
      if (rangeList->partialSymSums)
        rangeList->partialSymSums =
          env_ma_realloc(env, rangeList->partialSymSums, sizeof (Seqpos)
                         * numSyms * newSize);
      rangeList->numRangesStorable = newSize;
    }
    p = rangeList->ranges + numRanges;
    if (rangeList->partialSymSums)
      pSums = rangeList->partialSymSums + numSyms * numRanges;
    else
      pSums = NULL;
    while (len > MAX_SEQRANGE_LEN)
    {
      p->startPos = pos;
      p->sym = sym;
      p->len = MAX_SEQRANGE_LEN;
      pos += MAX_SEQRANGE_LEN;
      len -= MAX_SEQRANGE_LEN;
      if (numRanges && pSums)
      {
        Seqpos *spDest = pSums, *spSrc = pSums - numSyms;
        size_t i;
        for (i = 0; i < numSyms; ++i)
        {
          *spDest++ = *spSrc++;
        }
        pSums[(p - 1)->sym] += (p-1)->len;
        pSums += numSyms;
      }
      else if (pSums)
      {
        memset(pSums, 0, sizeof (Seqpos) * numSyms);
        pSums += numSyms;
      }
      ++p;
      ++numRanges;
    }
    if (len)
    {
      p->startPos = pos;
      p->len = len;
      p->sym = sym;
      if (numRanges && pSums)
      {
        Seqpos *spDest = pSums, *spSrc = pSums - numSyms;
        size_t i;
        for (i = 0; i < numSyms; ++i)
        {
          *spDest++ = *spSrc++;
        }
        pSums[(p - 1)->sym] += (p-1)->len;
        pSums += numSyms;
      }
      else if (pSums)
      {
        memset(pSums, 0, sizeof (Seqpos) * numSyms);
        pSums += numSyms;
      }
      ++numRanges;
    }
    assert(numRanges == rangeList->numRanges + numNewRanges);
    rangeList->numRanges = numRanges;
  }
}

void
SRLinsertNewRange(struct seqRangeList *rangeList, Seqpos pos, Seqpos len,
                  Symbol esym, Env *env)
{
  assert(rangeList && env);
  abort();
/*   { */
/*     Symbol sym = MRAEncMapSymbol(rangeList->alphabet, esym); */
/*   } */
  /* currently not implemented because only append is currently
   * needed */
}

void
SRLAddPosition(struct seqRangeList *rangeList, Seqpos pos,
               Symbol esym, Env *env)
{
  size_t numRanges;
  struct seqRange *lastRange;
  Symbol sym;
  assert(rangeList && env);
  sym = MRAEncMapSymbol(rangeList->alphabet, esym);
  numRanges = rangeList->numRanges;
  lastRange = rangeList->ranges + numRanges - 1;
  if (numRanges && lastRange->startPos + lastRange->len > pos)
  {
    /* TODO: search for range */
    SRLinsertNewRange(rangeList, pos, 1, esym, env);
  }
  else if (numRanges
          && (lastRange->sym == sym)
          && (lastRange->startPos + lastRange->len == pos)
          && (lastRange->len < MAX_SEQRANGE_LEN))
    ++(lastRange->len);
  else
    SRLAppendNewRange(rangeList, pos, 1, esym, env);
}

void
SRLInitListSearchHint(struct seqRangeList *rangeList,
                      seqRangeListSearchHint *hint)
{
  assert(rangeList && hint);
  *hint = rangeList->numRanges / 2;
}

static int
posSeqRangeOverlapCompare(const void *key, const void *elem)
{
  Seqpos pos;
  struct seqRange *range;
  pos = *(Seqpos *)key;
  range = (struct seqRange *)elem;
  if (pos < range->startPos)
    return -1;
  else if (pos >= range->startPos + range->len)
    return 1;
  else
    return 0;
}

int
SRLOverlapsPosition(struct seqRangeList *rangeList, Seqpos pos,
                    seqRangeListSearchHint *hint, Symbol *symAtPos)
{
  size_t rangeIdx;
  struct seqRange *p;
  assert(rangeList);
  if (hint)
    rangeIdx = *hint;
  else
    rangeIdx = rangeList->numRanges / 2;
  p = rangeList->ranges + rangeIdx;
  if (posIsInSeqRange(p, pos))
  {
    if (symAtPos)
      *symAtPos = MRAEncRevMapSymbol(rangeList->alphabet, p->sym);
    return 1;
  }
  /* TODO: implement "else if" case to also look at next range */
  else
  {
    Seqpos searchPos = pos;
    struct seqRange *searchRes =
      bsearch(&searchPos, rangeList->ranges, rangeList->numRanges,
              sizeof (struct seqRange), posSeqRangeOverlapCompare);
    if (searchRes)
    {
      if (symAtPos)
        *symAtPos = searchRes->sym;
      if (hint)
        *hint = searchRes - rangeList->ranges;
      return 1;
    }
    else
      return 0;
  }
}

static int
posSeqRangeNextCompare(const void *key, const void *elem)
{
  Seqpos pos;
  struct seqRange *range;
  pos = *(Seqpos *)key;
  range = (struct seqRange *)elem;
  if (pos < range->startPos)
    if (pos < range[-1].startPos + range[-1].len)
      return -1;
    else
      return 0;
  else if (pos >= range->startPos + range->len)
    return 1;
  else /* pos < range->startPos + range->len  && pos >= range->startPos */
    return 0;
}

#define rangeNextMatch(idx)                             \
  ((rangeList->ranges[idx].startPos >= pos              \
    || rangeList->ranges[idx].startPos                  \
    + rangeList->ranges[idx].len > pos)                 \
   && rangeList->ranges[idx - 1].startPos               \
   + rangeList->ranges[idx - 1].len <= pos)

/**
 * \brief Given a position and an optional search hint, find the next
 * range either overlapping or ahead of the position.
 * @param rangeList list of ranges to be searched
 * @param pos position to search range for
 * @param hint pointer to position of possibly near range
 * @return pointer to sequence range or NULL if no corresponding range
 * was found
 */
struct seqRange *
SRLFindPositionNext(struct seqRangeList *rangeList, Seqpos pos,
                    seqRangeListSearchHint *hint)
{
  /* no ranges no matches in ranges */
  seqRangeListSearchHint hintCopy;
  size_t numRanges;
  assert(rangeList);
  if (hint)
    hintCopy = *hint;
  else
    SRLInitListSearchHint(rangeList, &hintCopy);
  if ((numRanges = rangeList->numRanges) == 0)
  {
    return NULL;
  }
  if (rangeList->ranges[0].startPos >= pos
     || pos < rangeList->ranges[0].startPos + rangeList->ranges[0].len)
  {
    return rangeList->ranges + 0;
  }
  else if (hintCopy && rangeNextMatch(hintCopy))
  {
    return rangeList->ranges + hintCopy;
  }
  else if ((numRanges > hintCopy + 1) && rangeNextMatch(hintCopy + 1))
  {
    ++hintCopy;
    if (hint)
      *hint = hintCopy;
    return rangeList->ranges + hintCopy;
  }
  else if (numRanges > 2)
  {
    Seqpos searchPos = pos;
    struct seqRange *searchRes =
      bsearch(&searchPos, rangeList->ranges + 1, rangeList->numRanges - 1,
              sizeof (struct seqRange), posSeqRangeNextCompare);
    if (searchRes && hint)
      *hint = searchRes - rangeList->ranges;
    return searchRes;
  }
  else /* numRanges <= 2 */
    return NULL;
}

static int
posSeqRangeLastCompare(const void *key, const void *elem)
{
  Seqpos pos;
  struct seqRange *range;
  pos = *(Seqpos *)key;
  range = (struct seqRange *)elem;
  if (pos >= range->startPos)
    if (pos < range[1].startPos)
      return 0;
    else
      return 1;
  else
    return -1;
}

#define rangeLastMatch(idx)                             \
  (rangeList->ranges[idx].startPos <= pos               \
   && rangeList->ranges[idx + 1].startPos > pos)

/**
 * \brief Given a position and an optional search hint, find the next
 * range either overlapping or ahead of the position.
 * @param rangeList list of ranges to be searched
 * @param pos position to search range for
 * @param hint pointer to position of possibly near range
 * @return pointer to sequence range or NULL if no corresponding range
 * was found
 */
struct seqRange *
SRLFindPositionLast(struct seqRangeList *rangeList, Seqpos pos,
                    seqRangeListSearchHint *hint)
{
  /* no ranges no matches in ranges */
  seqRangeListSearchHint hintCopy;
  size_t numRanges;
  assert(rangeList);
  if (hint)
    hintCopy = *hint;
  else
    SRLInitListSearchHint(rangeList, &hintCopy);
  if ((numRanges = rangeList->numRanges) == 0)
  {
    return NULL;
  }
  if (rangeList->ranges[0].startPos > pos)
    return NULL;
  else if (rangeList->ranges[numRanges - 1].startPos <= pos)
  {
    return rangeList->ranges + numRanges - 1;
  }
  else if ((hintCopy < numRanges - 1) && rangeLastMatch(hintCopy))
  {
    return rangeList->ranges + hintCopy;
  }
  else if ((hintCopy < numRanges - 2) && rangeLastMatch(hintCopy + 1))
  {
    ++hintCopy;
    if (hint)
      *hint = hintCopy;
    return rangeList->ranges + hintCopy;
  }
  else if (hintCopy && rangeLastMatch(hintCopy - 1))
  {
    --hintCopy;
    if (hint)
      *hint = hintCopy;
    return rangeList->ranges + hintCopy;
  }
  else if (numRanges > 2)
  {
    Seqpos searchPos = pos;
    struct seqRange *searchRes =
      bsearch(&searchPos, rangeList->ranges, rangeList->numRanges - 1,
              sizeof (struct seqRange), posSeqRangeLastCompare);
    if (searchRes && hint)
      *hint = searchRes - rangeList->ranges;
    return searchRes;
  }
  else /* numRanges <= 2 */
    return NULL;
}

/**
 * @return 0 on error, >0 on success
 */
/* FIXME: convert to platform-independent variant */
int
SRLSaveToStream(struct seqRangeList *rangeList, FILE *fp)
{
  size_t numRanges;
  assert(rangeList && fp);
  numRanges = rangeList->numRanges;
  if (!fwrite(&(rangeList->numRanges), sizeof (rangeList->numRanges), 1, fp))
    return 0;
  if (fwrite(rangeList->ranges, sizeof (struct seqRange),
            numRanges, fp) != numRanges)
    return 0;
  return 1;
}

/* FIXME: convert to platform-independent variant */
struct seqRangeList *
SRLReadFromStream(FILE *fp, const MRAEnc *alphabet,
                  enum SRLFeatures features, Env *env)
{
  struct seqRangeList *newRangeList;
  size_t numRanges;
  assert(fp && env);
  newRangeList = env_ma_malloc(env, sizeof (struct seqRangeList));
  if (!fread(&(newRangeList->numRanges),
            sizeof (newRangeList->numRanges), 1, fp))
  {
    env_ma_free(newRangeList, env);
    return NULL;
  }
  numRanges = newRangeList->numRanges;
  newRangeList->partialSymSums = NULL;
  newRangeList->ranges = env_ma_malloc(env, sizeof (struct seqRange) *
                                       (newRangeList->numRangesStorable
                                        = numRanges));
  if (fread(newRangeList->ranges, sizeof (struct seqRange),
           numRanges, fp) != numRanges)
  {
    deleteSeqRangeList(newRangeList, env);
    return NULL;
  }
  if (features & SRL_PARTIAL_SYMBOL_SUMS)
  {
    Seqpos *partialSymSums;
    size_t numSyms = MRAEncGetSize(alphabet), i;
    newRangeList->partialSymSums = partialSymSums =
      env_ma_malloc(env, sizeof (Seqpos) * MRAEncGetSize(alphabet) * numRanges);
    memset(partialSymSums, 0, sizeof (Seqpos) * numSyms);
    for (i = 1; i < numRanges; ++i)
    {
      struct seqRange *lastRange = newRangeList->ranges + i - 1;
      Symbol lastSym = lastRange->sym;
      memcpy(partialSymSums + i * numSyms, partialSymSums + (i - 1) * numSyms,
             sizeof (Seqpos) * numSyms);
      partialSymSums[i * numSyms + lastSym] += lastRange->len;
    }
  }
  newRangeList->alphabet = alphabet;
  return newRangeList;
}

/**
 * @param occStore must be pre-initialized (either to some previous
 * count or zero) table of count for every symbol in alphabet,
 * i.e. indexing like occStore[sym] for every value of sym stored in
 * rangeList must be a valid memory reference.
 */
void
SRLSymbolsInSeqRegion(struct seqRangeList *rangeList, Seqpos start,
                      Seqpos end, Seqpos *occStore,
                      seqRangeListSearchHint *hint)
{
  struct seqRange *p;
  /* find range next to start position */
  /* first solve special case where pos precedes very first list
   * entry, because we want to be able to look at the preceding range
   * in binary search */
  /* no ranges => no symbols in ranges */
  if (rangeList->numRanges == 0)
    return;
  p = SRLFindPositionNext(rangeList, start, hint);
  /* no range overlapping or between start and end of sequence? */
  if (!p)
    return;
  /* iterate over ranges left */
  {
    Seqpos s = MAX(start, p->startPos);
    struct seqRange *maxRange = rangeList->ranges + rangeList->numRanges - 1;
    while (s <= end)
    {
      Seqpos overlap = MIN(p->startPos + p->len, end + 1) - s;
      occStore[MRAEncRevMapSymbol(rangeList->alphabet, p->sym)] += overlap;
      if (p == maxRange)
        break;
      s = (++p)->startPos;
    }
  }
}

Seqpos
SRLSymbolCountInSeqRegion(struct seqRangeList *rangeList, Seqpos start,
                          Seqpos end, Symbol esym, seqRangeListSearchHint *hint)
{
  const struct seqRange *p;
  if (rangeList->numRanges == 0)
    return 0;
  p = SRLFindPositionNext(rangeList, start, NULL);
  if (p)
  {
    if (rangeList->partialSymSums)
    {
      const struct seqRange *q = SRLFindPositionLast(rangeList, end, hint);
      if (!q)
      {
        /* in case
         * - end precedes first region
         * or
         * - start is after last region
         * no symbols can be encoded in range
         */
        return 0;
      }
      {
        size_t pOff = p - rangeList->ranges, qOff = q - rangeList->ranges;
        size_t numSyms = MRAEncGetSize(rangeList->alphabet);
        Seqpos symCount;
        Symbol sym = MRAEncMapSymbol(rangeList->alphabet, esym);
        symCount = rangeList->partialSymSums[qOff * numSyms + sym]
          - rangeList->partialSymSums[pOff * numSyms + sym];
        /* two special cases: start might be inside p and end inside q */
        if (sym == p->sym && start >= p->startPos)
          symCount -= start - p->startPos;
        if (sym == q->sym)
        {
          regionLength lastOverlap = MIN(end - q->startPos, q->len);
          symCount += lastOverlap;
        }
        return symCount;
      }
    }
    else
    {
      Seqpos symCount = 0;
      Seqpos s = MAX(start, p->startPos);
      Symbol sym = MRAEncMapSymbol(rangeList->alphabet, esym);
      struct seqRange *maxRange = rangeList->ranges + rangeList->numRanges - 1;
      while (s < end)
      {
        if (p->sym == sym)
          symCount += MIN(p->startPos + p->len, end) - s;
        if (p == maxRange)
          break;
        s = (++p)->startPos;
      }
      return symCount;
    }
  }
  else
    return 0;
}

Seqpos
SRLAllSymbolsCountInSeqRegion(struct seqRangeList *rangeList, Seqpos start,
                              Seqpos end, seqRangeListSearchHint *hint)
{
  const struct seqRange *p;
  if (rangeList->numRanges == 0)
    return 0;
  p = SRLFindPositionNext(rangeList, start, hint);
  if (p)
  {
    if (rangeList->partialSymSums)
    {
      const struct seqRange *q = SRLFindPositionLast(rangeList, end, hint);
      if (!q)
      {
        /* in case
         * - end precedes first region
         * or
         * - start is after last region
         * no symbols can be encoded in range
         */
        return 0;
      }
      {
        size_t pOff = p - rangeList->ranges, qOff = q - rangeList->ranges;
        size_t numSyms = MRAEncGetSize(rangeList->alphabet);
        Seqpos symCount = 0;
        Symbol sym;
        for (sym = 0; sym < numSyms; ++sym)
        {
          symCount += rangeList->partialSymSums[qOff * numSyms + sym]
            - rangeList->partialSymSums[pOff * numSyms + sym];
          /* two special cases: start might be inside p and end inside q */
        }
        if (start >= p->startPos)
          symCount -= start - p->startPos;

        {
          regionLength lastOverlap = MIN(end - q->startPos, q->len);
          symCount += lastOverlap;
        }
        return symCount;
      }
    }
    else
    {
      Seqpos symCount = 0;
      Seqpos s = MAX(start, p->startPos);
      struct seqRange *maxRange = rangeList->ranges + rangeList->numRanges - 1;
      while (s < end)
      {
        symCount += MIN(p->startPos + p->len, end) - s;
        if (p == maxRange)
          break;
        s = (++p)->startPos;
      }
      return symCount;
    }
  }
  else
    return 0;
}

void
SRLapplyRangesToSubString(struct seqRangeList *rangeList, MRAEnc *alphabet,
                          Symbol *subString, Seqpos start, Seqpos len,
                          Seqpos subStringOffset, seqRangeListSearchHint *hint)
{
  struct seqRange *nextRange;
  Seqpos inSeqPos = start;
  nextRange = SRLFindPositionNext(rangeList, inSeqPos, hint);
  do {
    if (inSeqPos < nextRange->startPos)
    {
      inSeqPos = nextRange->startPos;
    }
    else
    {
      size_t i;
      unsigned maxSubstPos =
        MIN(nextRange->startPos + nextRange->len,
            start + len) - subStringOffset;
      Symbol sym = MRAEncRevMapSymbol(alphabet,
                                      nextRange->sym);
      for (i = inSeqPos - subStringOffset; i < maxSubstPos; ++i)
        subString[i] = sym;
      inSeqPos = subStringOffset + maxSubstPos;
      ++nextRange;
    }
  } while (inSeqPos < start + len);
}

