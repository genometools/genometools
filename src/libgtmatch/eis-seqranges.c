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

#include "libgtcore/bsearch.h"
#include "libgtcore/error.h"
#include "libgtcore/ma.h"
#include "libgtcore/minmax.h"
#include "libgtcore/unused.h"
#include "libgtcore/xansi.h"

#include "libgtmatch/eis-seqranges.h"
#include "libgtmatch/eis-seqranges-priv.h"

struct seqRangeList *
newSeqRangeList(size_t rangesStartNum, const MRAEnc *alphabet,
                enum SRLFeatures features)
{
  struct seqRangeList *newList;
  newList = ma_malloc(sizeof (struct seqRangeList));
  newList->numRanges = 0;
  newList->numRangesStorable = rangesStartNum;
  newList->ranges = ma_malloc(sizeof (newList->ranges[0]) * rangesStartNum);
  if (features & SRL_PARTIAL_SYMBOL_SUMS)
    newList->partialSymSums = ma_malloc(sizeof (Seqpos)
                                        * MRAEncGetSize(alphabet)
                                        * rangesStartNum);
  else
    newList->partialSymSums = NULL;
  newList->alphabet = alphabet;
  newList->symBits = requiredSymbolBits(MRAEncGetSize(alphabet) - 1);
  if (newList->symBits)
    newList->maxRangeLen =
      (((Seqpos)1) << (symLenStrBits - newList->symBits)) - 1;
  else
    newList->maxRangeLen = ~(Seqpos)0;
  return newList;
}

void
SRLCompact(struct seqRangeList *rangeList)
{
  assert(rangeList);
  rangeList->ranges = ma_realloc(rangeList->ranges,
                                 sizeof (rangeList->ranges[0])
                                 * rangeList->numRanges);
  if (rangeList->partialSymSums)
    rangeList->partialSymSums =
      ma_realloc(rangeList->partialSymSums, sizeof (Seqpos)
                 * MRAEncGetSize(rangeList->alphabet)
                 * rangeList->numRanges);
  rangeList->numRangesStorable = rangeList->numRanges;
}

void
deleteSeqRangeList(struct seqRangeList *rangeList)
{
  assert(rangeList);
  if (rangeList->ranges)
    ma_free(rangeList->ranges);
  if (rangeList->partialSymSums)
    ma_free(rangeList->partialSymSums);
  ma_free(rangeList);
}

void
SRLAppendNewRange(struct seqRangeList *rangeList, Seqpos pos, Seqpos len,
                  Symbol esym)
{
  assert(rangeList);
  if (len)
  {
    Symbol sym = MRAEncMapSymbol(rangeList->alphabet, esym);
    struct seqRange *p;
    Seqpos *pSums;
    size_t numRanges = rangeList->numRanges,
      numSyms = MRAEncGetSize(rangeList->alphabet),
      numNewRanges = len/rangeList->maxRangeLen
      + ((len%rangeList->maxRangeLen)?1:0);
    unsigned symBits = rangeList->symBits;
    if (numRanges + numNewRanges > rangeList->numRangesStorable)
    {
      size_t newSize = numRanges + 2 * numNewRanges;
      rangeList->ranges = ma_realloc(rangeList->ranges,
                                     sizeof (struct seqRange) * newSize);
      if (rangeList->partialSymSums)
        rangeList->partialSymSums =
          ma_realloc(rangeList->partialSymSums, sizeof (Seqpos)
                     * numSyms * newSize);
      rangeList->numRangesStorable = newSize;
    }
    p = rangeList->ranges + numRanges;
    if (rangeList->partialSymSums)
      pSums = rangeList->partialSymSums + numSyms * numRanges;
    else
      pSums = NULL;
    while (len > rangeList->maxRangeLen)
    {
      p->startPos = pos;
      seqRangeSetSym(p, sym, symBits);
      seqRangeSetLen(p, rangeList->maxRangeLen, symBits);
      pos += rangeList->maxRangeLen;
      len -= rangeList->maxRangeLen;
      if (numRanges && pSums)
      {
        Seqpos *spDest = pSums, *spSrc = pSums - numSyms;
        size_t i;
        for (i = 0; i < numSyms; ++i)
        {
          *spDest++ = *spSrc++;
        }
        pSums[seqRangeSym(p - 1, symBits)] += seqRangeLen(p - 1, symBits);
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
      seqRangeSetSym(p, sym, symBits);
      seqRangeSetLen(p, len, symBits);
      if (numRanges && pSums)
      {
        Seqpos *spDest = pSums, *spSrc = pSums - numSyms;
        size_t i;
        for (i = 0; i < numSyms; ++i)
        {
          *spDest++ = *spSrc++;
        }
        pSums[seqRangeSym(p - 1, symBits)] += seqRangeLen(p - 1, symBits);
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
SRLinsertNewRange(UNUSED struct seqRangeList *rangeList, UNUSED Seqpos pos,
                  UNUSED Seqpos len, UNUSED Symbol esym)
{
  assert(rangeList);
  abort();
/*   { */
/*     Symbol sym = MRAEncMapSymbol(rangeList->alphabet, esym); */
/*   } */
  /* currently not implemented because only append is currently
   * needed */
}

void
SRLAddPosition(struct seqRangeList *rangeList, Seqpos pos, Symbol esym)
{
  size_t numRanges;
  struct seqRange *lastRange;
  Seqpos lastRangeLen;
  unsigned symBits;
  Symbol sym;
  assert(rangeList);
  sym = MRAEncMapSymbol(rangeList->alphabet, esym);
  numRanges = rangeList->numRanges;
  lastRange = rangeList->ranges + numRanges - 1;
  symBits = rangeList->symBits;
  if (numRanges)
    lastRangeLen = seqRangeLen(lastRange, symBits);

  if (numRanges && lastRange->startPos + lastRangeLen > pos)
  {
    /* TODO: search for range */
    SRLinsertNewRange(rangeList, pos, 1, esym);
  }
  else if (numRanges
           && (seqRangeSym(lastRange, symBits) == sym)
           && (lastRange->startPos + lastRangeLen == pos)
           && (lastRangeLen < rangeList->maxRangeLen))
    seqRangeSetLen(lastRange, ++lastRangeLen, symBits);
  else
    SRLAppendNewRange(rangeList, pos, 1, esym);
}

void
SRLInitListSearchHint(struct seqRangeList *rangeList,
                      seqRangeListSearchHint *hint)
{
  assert(rangeList && hint);
  *hint = rangeList->numRanges / 2;
}

static int
posSeqRangeOverlapCompare(const void *key, const void *elem, const void *data)
{
  Seqpos pos;
  const struct seqRange *range;
  const struct seqRangeList *rangeList;
  pos = *(const Seqpos *)key;
  range = elem;
  rangeList = data;
  if (pos < range->startPos)
    return -1;
  else if (pos >= range->startPos + seqRangeLen(range, rangeList->symBits))
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
  if (posIsInSeqRange(p, pos, rangeList->symBits))
  {
    if (symAtPos)
      *symAtPos = MRAEncRevMapSymbol(rangeList->alphabet,
                                     seqRangeSym(p, rangeList->symBits));
    return 1;
  }
  /* TODO: implement "else if" case to also look at next range */
  else
  {
    Seqpos searchPos = pos;
    struct seqRange *searchRes =
      bsearch_data(&searchPos, rangeList->ranges, rangeList->numRanges,
                   sizeof (struct seqRange), posSeqRangeOverlapCompare,
                   rangeList);
    if (searchRes)
    {
      if (symAtPos)
        *symAtPos = seqRangeSym(searchRes, rangeList->symBits);
      if (hint)
        *hint = searchRes - rangeList->ranges;
      return 1;
    }
    else
      return 0;
  }
}

static int
posSeqRangeNextCompare(const void *key, const void *elem, const void *data)
{
  Seqpos pos;
  const struct seqRange *range;
  const struct seqRangeList *rangeList;
  pos = *(const Seqpos *)key;
  range = elem;
  rangeList = data;
  if (pos < range->startPos)
    if (pos < range[-1].startPos + seqRangeLen(range - 1, rangeList->symBits))
      return -1;
    else
      return 0;
  else if (pos >= range->startPos + seqRangeLen(range, rangeList->symBits))
    return 1;
  else /* pos < range->startPos + range->len  && pos >= range->startPos */
    return 0;
}

#define rangeNextMatch(idx)                                  \
  ((ranges[idx].startPos >= pos                              \
    || ranges[idx].startPos + seqRangeLen(ranges + idx, symBits) > pos) \
   && ranges[idx - 1].startPos + seqRangeLen(ranges + idx - 1, symBits) <= pos)

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
  struct seqRange *ranges;
  unsigned symBits;
  assert(rangeList);
  ranges = rangeList->ranges;
  symBits = rangeList->symBits;
  if (hint)
    hintCopy = *hint;
  else
    SRLInitListSearchHint(rangeList, &hintCopy);
  if ((numRanges = rangeList->numRanges) == 0)
  {
    return NULL;
  }
  if (ranges[0].startPos >= pos
      || pos < ranges[0].startPos + seqRangeLen(ranges + 0, symBits))
  {
    return ranges + 0;
  }
  else if (hintCopy && rangeNextMatch(hintCopy))
  {
    return ranges + hintCopy;
  }
  else if ((numRanges > hintCopy + 1) && rangeNextMatch(hintCopy + 1))
  {
    ++hintCopy;
    if (hint)
      *hint = hintCopy;
    return rangeList->ranges + hintCopy;
  }
  else if (hintCopy && rangeNextMatch(hintCopy - 1))
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
      bsearch_data(&searchPos, rangeList->ranges + 1, rangeList->numRanges - 1,
                   sizeof (struct seqRange), posSeqRangeNextCompare, rangeList);
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

#define rangeLastMatch(idx)                                             \
  (ranges[idx].startPos <= pos && ranges[idx + 1].startPos > pos)

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
  struct seqRange *ranges;
  unsigned symBits;
  assert(rangeList);
  ranges = rangeList->ranges;
  symBits = rangeList->symBits;
  if (hint)
    hintCopy = *hint;
  else
    SRLInitListSearchHint(rangeList, &hintCopy);
  if ((numRanges = rangeList->numRanges) == 0)
  {
    return NULL;
  }
  if (ranges[0].startPos > pos)
    return NULL;
  else if (ranges[numRanges - 1].startPos <= pos)
  {
    return ranges + numRanges - 1;
  }
  else if ((hintCopy < numRanges - 1) && rangeLastMatch(hintCopy))
  {
    return ranges + hintCopy;
  }
  else if ((hintCopy < numRanges - 2) && rangeLastMatch(hintCopy + 1))
  {
    ++hintCopy;
    if (hint)
      *hint = hintCopy;
    return ranges + hintCopy;
  }
  else if (hintCopy && rangeLastMatch(hintCopy - 1))
  {
    --hintCopy;
    if (hint)
      *hint = hintCopy;
    return ranges + hintCopy;
  }
  else if (numRanges > 2)
  {
    Seqpos searchPos = pos;
    struct seqRange *searchRes =
      bsearch(&searchPos, ranges, rangeList->numRanges - 1,
              sizeof (struct seqRange), posSeqRangeLastCompare);
    if (searchRes && hint)
      *hint = searchRes - ranges;
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
  xfwrite(&(rangeList->numRanges), sizeof (rangeList->numRanges), 1, fp);
  xfwrite(rangeList->ranges, sizeof (struct seqRange), numRanges, fp);
  return 1;
}

/* FIXME: convert to platform-independent variant */
struct seqRangeList *
SRLReadFromStream(FILE *fp, const MRAEnc *alphabet,
                  enum SRLFeatures features, UNUSED Error *err)
{
  struct seqRangeList *newRangeList;
  size_t numRanges;
  assert(fp && err);
  newRangeList = ma_malloc(sizeof (struct seqRangeList));
  newRangeList->alphabet = alphabet;
  newRangeList->symBits = requiredSymbolBits(MRAEncGetSize(alphabet) - 1);
  if (newRangeList->symBits)
    newRangeList->maxRangeLen =
      (((Seqpos)1) << (symLenStrBits - newRangeList->symBits)) - 1;
  else
    newRangeList->maxRangeLen = ~(Seqpos)0;
  xfread(&(newRangeList->numRanges), sizeof (newRangeList->numRanges), 1, fp);
  numRanges = newRangeList->numRanges;
  newRangeList->partialSymSums = NULL;
  newRangeList->ranges = ma_malloc(sizeof (struct seqRange) *
                                   (newRangeList->numRangesStorable
                                   = numRanges));
  xfread(newRangeList->ranges, sizeof (struct seqRange), numRanges, fp);
  if (features & SRL_PARTIAL_SYMBOL_SUMS)
  {
    Seqpos *partialSymSums;
    size_t numSyms = MRAEncGetSize(alphabet), i;
    newRangeList->partialSymSums = partialSymSums =
      ma_malloc(sizeof (Seqpos) * MRAEncGetSize(alphabet) * numRanges);
    memset(partialSymSums, 0, sizeof (Seqpos) * numSyms);
    for (i = 1; i < numRanges; ++i)
    {
      struct seqRange *lastRange = newRangeList->ranges + i - 1;
      Symbol lastSym = seqRangeSym(lastRange, newRangeList->symBits);
      memcpy(partialSymSums + i * numSyms, partialSymSums + (i - 1) * numSyms,
             sizeof (Seqpos) * numSyms);
      partialSymSums[i * numSyms + lastSym] +=
        seqRangeLen(lastRange, newRangeList->symBits);
    }
  }
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
      Seqpos overlap = MIN(p->startPos + seqRangeLen(p, rangeList->symBits),
                           end + 1) - s;
      occStore[MRAEncRevMapSymbol(rangeList->alphabet,
                                  seqRangeSym(p, rangeList->symBits))]
               += overlap;
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
        unsigned symBits = rangeList->symBits;
        symCount = rangeList->partialSymSums[qOff * numSyms + sym]
          - rangeList->partialSymSums[pOff * numSyms + sym];
        /* two special cases: start might be inside p and end inside q */
        if (sym == seqRangeSym(p, symBits) && start >= p->startPos)
          symCount -= start - p->startPos;
        if (sym == seqRangeSym(q, symBits))
        {
          regionLength lastOverlap = MIN(end - q->startPos,
                                         seqRangeLen(q, symBits));
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
      unsigned symBits = rangeList->symBits;
      struct seqRange *maxRange = rangeList->ranges + rangeList->numRanges - 1;
      while (s < end)
      {
        if (seqRangeSym(p, symBits) == sym)
          symCount += MIN(p->startPos + seqRangeLen(p, symBits), end) - s;
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
        unsigned symBits = rangeList->symBits;
        for (sym = 0; sym < numSyms; ++sym)
        {
          symCount += rangeList->partialSymSums[qOff * numSyms + sym]
            - rangeList->partialSymSums[pOff * numSyms + sym];
          /* two special cases: start might be inside p and end inside q */
        }
        if (start >= p->startPos)
          symCount -= start - p->startPos;

        {
          regionLength lastOverlap = MIN(end - q->startPos,
                                         seqRangeLen(q, symBits));
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
      unsigned symBits = rangeList->symBits;
      while (s < end)
      {
        symCount += MIN(p->startPos + seqRangeLen(p, symBits), end) - s;
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
SRLApplyRangesToSubString(struct seqRangeList *rangeList,
                          Symbol *subString, Seqpos start, Seqpos len,
                          Seqpos subStringOffset, seqRangeListSearchHint *hint)
{
  struct seqRange *nextRange;
  Seqpos inSeqPos = start;
  unsigned symBits = rangeList->symBits;
  assert(rangeList);
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
        MIN(nextRange->startPos + seqRangeLen(nextRange, symBits),
            start + len) - subStringOffset;
      Symbol sym = MRAEncRevMapSymbol(rangeList->alphabet,
                                      seqRangeSym(nextRange, symBits));
      for (i = inSeqPos - subStringOffset; i < maxSubstPos; ++i)
        subString[i] = sym;
      inSeqPos = subStringOffset + maxSubstPos;
      ++nextRange;
    }
  } while (inSeqPos < start + len);
}

extern int
SRLPrintRangesInfo(struct seqRangeList *rangeList, FILE *fp, Seqpos start,
                   Seqpos len, seqRangeListSearchHint *hint)
{
  struct seqRange *nextRange;
  Seqpos end = start + len;
  unsigned symBits = rangeList->symBits;
  int result = 0;
  assert(rangeList);
  nextRange = SRLFindPositionNext(rangeList, start, hint);
  while (nextRange->startPos < end)
  {
    if (rangeList->partialSymSums)
    {
      size_t numSyms = MRAEncGetSize(rangeList->alphabet);
      size_t pOff = nextRange - rangeList->ranges;
      Symbol sym;
      fputs("# range partial sums:", fp);
      for (sym = 0; sym < numSyms; ++sym)
        fprintf(fp, " sum[%u]="FormatSeqpos,
                MRAEncRevMapSymbol(rangeList->alphabet, sym),
                rangeList->partialSymSums[pOff * numSyms + sym]);
      fputs("\n", fp);
    }
    if (result +=
        fprintf(fp, "# range overlap: symbol %u, start="FormatSeqpos", length="
                FormatSeqpos"\n",
                MRAEncRevMapSymbol(rangeList->alphabet,
                                   seqRangeSym(nextRange, symBits)),
                nextRange->startPos, seqRangeLen(nextRange, symBits)) < 0)
    {
      result = -1;
      break;
    }
    ++nextRange;
  }
  return result;
}
