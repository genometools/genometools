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

#include "core/bsearch.h"
#include "core/error.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"

#include "match/eis-seqranges.h"
#include "match/eis-seqranges-priv.h"

struct seqRangeList *
gt_newSeqRangeList(size_t rangesStartNum, const MRAEnc *alphabet,
                enum SRLFeatures features)
{
  struct seqRangeList *newList;
  newList = gt_malloc(sizeof (struct seqRangeList));
  newList->numRanges = 0;
  newList->numRangesStorable = rangesStartNum;
  newList->ranges = gt_malloc(sizeof (newList->ranges[0]) * rangesStartNum);
  if (features & SRL_PARTIAL_SYMBOL_SUMS)
    newList->partialSymSums = gt_malloc(sizeof (GtUword)
                                        * gt_MRAEncGetSize(alphabet)
                                        * rangesStartNum);
  else
    newList->partialSymSums = NULL;
  newList->alphabet = alphabet;
  newList->symBits = requiredSymbolBits(gt_MRAEncGetSize(alphabet) - 1);
  if (newList->symBits)
    newList->maxRangeLen =
      (((GtUword)1) << (symLenStrBits - newList->symBits)) - 1;
  else
    newList->maxRangeLen = ~(GtUword)0;
  return newList;
}

void
gt_SRLCompact(struct seqRangeList *rangeList)
{
  gt_assert(rangeList);
  rangeList->ranges = gt_realloc(rangeList->ranges,
                                 sizeof (rangeList->ranges[0])
                                 * rangeList->numRanges);
  if (rangeList->partialSymSums)
    rangeList->partialSymSums =
      gt_realloc(rangeList->partialSymSums, sizeof (GtUword)
                 * gt_MRAEncGetSize(rangeList->alphabet)
                 * rangeList->numRanges);
  rangeList->numRangesStorable = rangeList->numRanges;
}

void
gt_deleteSeqRangeList(struct seqRangeList *rangeList)
{
  gt_assert(rangeList);
  if (rangeList->ranges)
    gt_free(rangeList->ranges);
  if (rangeList->partialSymSums)
    gt_free(rangeList->partialSymSums);
  gt_free(rangeList);
}

void
gt_SRLAppendNewRange(struct seqRangeList *rangeList,
                  GtUword pos,
                  GtUword len,
                  Symbol esym)
{
  gt_assert(rangeList);
  if (len)
  {
    Symbol sym = MRAEncMapSymbol(rangeList->alphabet, esym);
    struct seqRange *p;
    GtUword *pSums;
    size_t numRanges = rangeList->numRanges,
      numSyms = gt_MRAEncGetSize(rangeList->alphabet),
      numNewRanges = len/rangeList->maxRangeLen
      + ((len%rangeList->maxRangeLen)?1:0);
    unsigned symBits = rangeList->symBits;
    if (numRanges + numNewRanges > rangeList->numRangesStorable)
    {
      size_t newSize = numRanges + 2 * numNewRanges;
      rangeList->ranges = gt_realloc(rangeList->ranges,
                                     sizeof (struct seqRange) * newSize);
      if (rangeList->partialSymSums)
        rangeList->partialSymSums =
          gt_realloc(rangeList->partialSymSums, sizeof (GtUword)
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
        GtUword *spDest = pSums, *spSrc = pSums - numSyms;
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
        memset(pSums, 0, sizeof (GtUword) * numSyms);
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
        GtUword *spDest = pSums, *spSrc = pSums - numSyms;
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
        memset(pSums, 0, sizeof (GtUword) * numSyms);
        pSums += numSyms;
      }
      ++numRanges;
    }
    gt_assert(numRanges == rangeList->numRanges + numNewRanges);
    rangeList->numRanges = numRanges;
  }
}

void
gt_SRLinsertNewRange(GT_UNUSED struct seqRangeList *rangeList,
                  GT_UNUSED GtUword pos, GT_UNUSED GtUword len,
                  GT_UNUSED Symbol esym)
{
  gt_assert(rangeList);
  abort();
/*   { */
/*     Symbol sym = MRAEncMapSymbol(rangeList->alphabet, esym); */
/*   } */
  /* currently not implemented because only append is currently
   * needed */
}

void
gt_SRLAddPosition(struct seqRangeList *rangeList, GtUword pos,
                  Symbol esym)
{
  size_t numRanges;
  struct seqRange *lastRange;
  GtUword lastRangeLen;
  unsigned symBits;
  Symbol sym;
  gt_assert(rangeList);
  sym = MRAEncMapSymbol(rangeList->alphabet, esym);
  numRanges = rangeList->numRanges;
  lastRange = rangeList->ranges + numRanges - 1;
  symBits = rangeList->symBits;
  if (numRanges)
    lastRangeLen = seqRangeLen(lastRange, symBits);

  if (numRanges && lastRange->startPos + lastRangeLen > pos)
  {
    /* TODO: search for range */
    gt_SRLinsertNewRange(rangeList, pos, 1, esym);
  }
  else if (numRanges
           && (seqRangeSym(lastRange, symBits) == sym)
           && (lastRange->startPos + lastRangeLen == pos)
           && (lastRangeLen < rangeList->maxRangeLen))
    seqRangeSetLen(lastRange, ++lastRangeLen, symBits);
  else
    gt_SRLAppendNewRange(rangeList, pos, 1, esym);
}

void
gt_SRLInitListSearchHint(struct seqRangeList *rangeList,
                      seqRangeListSearchHint *hint)
{
  gt_assert(rangeList && hint);
  *hint = rangeList->numRanges / 2;
}

static int
posSeqRangeOverlapCompare(const void *key, const void *elem, void *data)
{
  GtUword pos;
  const struct seqRange *range;
  const struct seqRangeList *rangeList;
  pos = *(const GtUword *)key;
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
gt_SRLOverlapsPosition(struct seqRangeList *rangeList, GtUword pos,
                    seqRangeListSearchHint *hint, Symbol *symAtPos)
{
  size_t rangeIdx;
  struct seqRange *p;
  gt_assert(rangeList);
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
    GtUword searchPos = pos;
    struct seqRange *searchRes =
      gt_bsearch_data(&searchPos, rangeList->ranges, rangeList->numRanges,
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
posSeqRangeNextCompare(const void *key, const void *elem, void *data)
{
  GtUword pos;
  const struct seqRange *range;
  const struct seqRangeList *rangeList;
  pos = *(const GtUword *)key;
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
gt_SRLFindPositionNext(struct seqRangeList *rangeList, GtUword pos,
                    seqRangeListSearchHint *hint)
{
  /* no ranges no matches in ranges */
  seqRangeListSearchHint hintCopy;
  size_t numRanges;
  struct seqRange *ranges;
  unsigned symBits;
  gt_assert(rangeList);
  ranges = rangeList->ranges;
  symBits = rangeList->symBits;
  if (hint)
    hintCopy = *hint;
  else
    gt_SRLInitListSearchHint(rangeList, &hintCopy);
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
    GtUword searchPos = pos;
    struct seqRange *searchRes =
      gt_bsearch_data(&searchPos, rangeList->ranges + 1,
                      rangeList->numRanges - 1, sizeof (struct seqRange),
                      posSeqRangeNextCompare, rangeList);
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
  GtUword pos;
  struct seqRange *range;
  pos = *(GtUword *)key;
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
gt_SRLFindPositionLast(struct seqRangeList *rangeList, GtUword pos,
                    seqRangeListSearchHint *hint)
{
  /* no ranges no matches in ranges */
  seqRangeListSearchHint hintCopy;
  size_t numRanges;
  struct seqRange *ranges;
  GT_UNUSED unsigned symBits;
  gt_assert(rangeList);
  ranges = rangeList->ranges;
  symBits = rangeList->symBits;
  if (hint)
    hintCopy = *hint;
  else
    gt_SRLInitListSearchHint(rangeList, &hintCopy);
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
    GtUword searchPos = pos;
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
gt_SRLSaveToStream(struct seqRangeList *rangeList, FILE *fp)
{
  size_t numRanges;
  gt_assert(rangeList && fp);
  numRanges = rangeList->numRanges;
  gt_xfwrite(&(rangeList->numRanges), sizeof (rangeList->numRanges), 1, fp);
  gt_xfwrite(rangeList->ranges, sizeof (struct seqRange), numRanges, fp);
  return 1;
}

/* FIXME: convert to platform-independent variant */
struct seqRangeList *
gt_SRLReadFromStream(FILE *fp, const MRAEnc *alphabet,
                  enum SRLFeatures features, GT_UNUSED GtError *err)
{
  struct seqRangeList *newRangeList;
  size_t numRanges;
  gt_assert(fp && err);
  newRangeList = gt_malloc(sizeof (struct seqRangeList));
  newRangeList->alphabet = alphabet;
  newRangeList->symBits = requiredSymbolBits(gt_MRAEncGetSize(alphabet) - 1);
  if (newRangeList->symBits)
    newRangeList->maxRangeLen =
      (((GtUword)1) << (symLenStrBits - newRangeList->symBits)) - 1;
  else
    newRangeList->maxRangeLen = ~(GtUword)0;
  gt_xfread(&(newRangeList->numRanges), sizeof (newRangeList->numRanges), 1,
            fp);
  numRanges = newRangeList->numRanges;
  newRangeList->partialSymSums = NULL;
  newRangeList->ranges = gt_malloc(sizeof (struct seqRange) *
                                   (newRangeList->numRangesStorable
                                   = numRanges));
  gt_xfread(newRangeList->ranges, sizeof (struct seqRange), numRanges, fp);
  if (features & SRL_PARTIAL_SYMBOL_SUMS)
  {
    GtUword *partialSymSums;
    size_t numSyms = gt_MRAEncGetSize(alphabet), i;
    newRangeList->partialSymSums = partialSymSums =
     gt_malloc(sizeof (GtUword) * gt_MRAEncGetSize(alphabet) * numRanges);
    memset(partialSymSums, 0, sizeof (GtUword) * numSyms);
    for (i = 1; i < numRanges; ++i)
    {
      struct seqRange *lastRange = newRangeList->ranges + i - 1;
      Symbol lastSym = seqRangeSym(lastRange, newRangeList->symBits);
      memcpy(partialSymSums + i * numSyms, partialSymSums + (i - 1) * numSyms,
             sizeof (GtUword) * numSyms);
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
gt_SRLSymbolsInSeqRegion(struct seqRangeList *rangeList, GtUword start,
                      GtUword end, GtUword *occStore,
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
  p = gt_SRLFindPositionNext(rangeList, start, hint);
  /* no range overlapping or between start and end of sequence? */
  if (!p)
    return;
  /* iterate over ranges left */
  {
    GtUword s = MAX(start, p->startPos);
    struct seqRange *maxRange = rangeList->ranges + rangeList->numRanges - 1;
    while (s <= end)
    {
      GtUword overlap = MIN(p->startPos + seqRangeLen(p,
                                                            rangeList->symBits),
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

GtUword
gt_SRLSymbolCountInSeqRegion(struct seqRangeList *rangeList,
                          GtUword start,
                          GtUword end,
                          Symbol esym,
                          seqRangeListSearchHint *hint)
{
  const struct seqRange *p;
  if (rangeList->numRanges == 0)
    return 0;
  p = gt_SRLFindPositionNext(rangeList, start, NULL);
  if (p)
  {
    if (rangeList->partialSymSums)
    {
      const struct seqRange *q = gt_SRLFindPositionLast(rangeList, end, hint);
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
        size_t numSyms = gt_MRAEncGetSize(rangeList->alphabet);
        GtUword symCount;
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
      GtUword symCount = 0;
      GtUword s = MAX(start, p->startPos);
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

GtUword
gt_SRLAllSymbolsCountInSeqRegion(struct seqRangeList *rangeList,
                              GtUword start,
                              GtUword end,
                              seqRangeListSearchHint *hint)
{
  const struct seqRange *p;
  if (rangeList->numRanges == 0)
    return 0;
  p = gt_SRLFindPositionNext(rangeList, start, hint);
  if (p)
  {
    if (rangeList->partialSymSums)
    {
      const struct seqRange *q = gt_SRLFindPositionLast(rangeList, end, hint);
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
        size_t numSyms = gt_MRAEncGetSize(rangeList->alphabet);
        GtUword symCount = 0;
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
      GtUword symCount = 0;
      GtUword s = MAX(start, p->startPos);
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
gt_SRLApplyRangesToSubString(struct seqRangeList *rangeList,
                          Symbol *subString,
                          GtUword start,
                          GtUword len,
                          GtUword subStringOffset,
                          seqRangeListSearchHint *hint)
{
  struct seqRange *nextRange;
  GtUword inSeqPos = start;
  unsigned symBits = rangeList->symBits;
  gt_assert(rangeList);
  nextRange = gt_SRLFindPositionNext(rangeList, inSeqPos, hint);
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

int
gt_SRLPrintRangesInfo(struct seqRangeList *rangeList,
                   FILE *fp,
                   GtUword start,
                   GtUword len,
                   seqRangeListSearchHint *hint)
{
  struct seqRange *nextRange;
  GtUword end = start + len;
  unsigned symBits = rangeList->symBits;
  int result = 0;
  gt_assert(rangeList);
  nextRange = gt_SRLFindPositionNext(rangeList, start, hint);
  while (nextRange->startPos < end)
  {
    if (rangeList->partialSymSums)
    {
      size_t numSyms = gt_MRAEncGetSize(rangeList->alphabet);
      size_t pOff = nextRange - rangeList->ranges;
      Symbol sym;
      fputs("# range partial sums:", fp);
      for (sym = 0; sym < numSyms; ++sym)
        fprintf(fp, " sum[%u]="GT_WU"",
                MRAEncRevMapSymbol(rangeList->alphabet, sym),
                rangeList->partialSymSums[pOff * numSyms + sym]);
      fputs("\n", fp);
    }
    if (result +=
        fprintf(fp, "# range overlap: symbol %u, start="GT_WU", length="
                GT_WU"\n", MRAEncRevMapSymbol(rangeList->alphabet,
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
