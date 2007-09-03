/*
** Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
**  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**  
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**  
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
**  
*/

#include <libgtcore/minmax.h>

#include "seqranges.h"

struct seqRangeList *
newSeqRangeList(size_t rangesStartNum, Env *env)
{
  assert(env);
  struct seqRangeList *newList;
  newList = env_ma_malloc(env, sizeof(struct seqRangeList));
  newList->numRanges = 0;
  newList->numRangesStorable = rangesStartNum;
  newList->ranges = env_ma_malloc(env,
                                  sizeof(newList->ranges[0]) * rangesStartNum);
  return newList;
}

void
SRLCompact(struct seqRangeList *rangeList, Env *env)
{
  assert(rangeList && env);
  rangeList->ranges = env_ma_realloc(env, rangeList->ranges,
                                     sizeof(rangeList->ranges[0])
                                     * rangeList->numRanges);
  rangeList->numRangesStorable = rangeList->numRanges;
}

void
deleteSeqRangeList(struct seqRangeList *rangeList, Env *env)
{
  assert(rangeList && env);
  if(rangeList->ranges)
    env_ma_free(rangeList->ranges, env);
  env_ma_free(rangeList, env);
}

void
SRLAppendNewRange(struct seqRangeList *rangeList, Seqpos pos, Seqpos len,
                  Symbol sym, Env *env)
{
  assert(rangeList && env);
  if(len)
  {
    struct seqRange *p;
    size_t numRanges = rangeList->numRanges,
      numNewRanges = len/MAX_SEQRANGE_LEN + ((len%MAX_SEQRANGE_LEN)?1:0);
    if(numRanges + numNewRanges > rangeList->numRangesStorable)
      rangeList->ranges =
        env_ma_realloc(env, rangeList->ranges,
                       sizeof(struct seqRange)
                       * (rangeList->numRangesStorable
                          = (numRanges + 2 * numNewRanges)));
    p = rangeList->ranges + numRanges;
    while(len > MAX_SEQRANGE_LEN)
    {
      p->startPos = pos;
      p->sym = sym;
      p->len = MAX_SEQRANGE_LEN;
      pos += MAX_SEQRANGE_LEN;
      len -= MAX_SEQRANGE_LEN;
      ++p;
      ++numRanges;
    }
    if(len)
    {
      p->startPos = pos;
      p->len = len;
      p->sym = sym;
      ++numRanges;
    }
    assert(numRanges == rangeList->numRanges + numNewRanges);
    rangeList->numRanges = numRanges;
  }
}

void
SRLinsertNewRange(struct seqRangeList *rangeList, Seqpos pos, Seqpos len,
                  Symbol sym, Env *env)
{
  assert(rangeList && env);
  /* currently not implemented because only append is currently
   * needed */
  abort();
}

void
SRLAddPosition(struct seqRangeList *rangeList, Seqpos pos,
               Symbol sym, Env *env)
{
  size_t numRanges;
  struct seqRange *lastRange;
  assert(rangeList && env);
  numRanges = rangeList->numRanges;
  lastRange = rangeList->ranges + numRanges - 1;
  if(numRanges && lastRange->startPos > pos)
  {
    /* TODO: search for range */
    SRLinsertNewRange(rangeList, pos, 1, sym, env);
  }
  else if(numRanges
          && (lastRange->sym == sym)
          && (lastRange->startPos + lastRange->len == pos)
          && (lastRange->len < MAX_SEQRANGE_LEN))
    ++(lastRange->len);
  else
    SRLAppendNewRange(rangeList, pos, 1, sym, env);
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
  if(pos < range->startPos)
    return -1;
  else if(pos >= range->startPos + range->len)
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
  if(hint)
    rangeIdx = *hint;
  else
    rangeIdx = rangeList->numRanges / 2;
  p = rangeList->ranges + rangeIdx;
  if(posIsInSeqRange(p, pos))
  {
    if(symAtPos)
      *symAtPos = p->sym;
    return 1;
  }
  /* TODO: implement "else if" case to also look at next range */
  else
  {
    Seqpos searchPos = pos;
    struct seqRange *searchRes =
      bsearch(&searchPos, rangeList->ranges, rangeList->numRanges,
              sizeof(struct seqRange), posSeqRangeOverlapCompare);
    if(searchRes)
    {
      if(symAtPos)
        *symAtPos = searchRes->sym;
      if(hint)
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
  if(pos < range->startPos)
    if(pos < range[-1].startPos + range[-1].len)
      return -1;
    else
      return 0;
  else if(pos >= range->startPos + range->len)
    return 1;
  else /* pos < range->startPos + range->len  && pos >= range->startPos */
    return 0;
}


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
  if(hint)
    hintCopy = *hint;
  else
    SRLInitListSearchHint(rangeList, &hintCopy);
  if((numRanges = rangeList->numRanges) == 0)
  {
    return NULL;
  }
  if(rangeList->ranges[0].startPos >= pos
     || pos < rangeList->ranges[0].startPos + rangeList->ranges[0].len)
  {
    return rangeList->ranges + 0;
  }
  else if(hintCopy
          && (rangeList->ranges[hintCopy].startPos >= pos
              || rangeList->ranges[hintCopy].startPos
              + rangeList->ranges[hintCopy].len > pos)
          && rangeList->ranges[hintCopy - 1].startPos
          + rangeList->ranges[hintCopy - 1].len <= pos)
  {
    return rangeList->ranges + hintCopy;
  }
  else if((numRanges > hintCopy + 1)
          && (rangeList->ranges[hintCopy + 1].startPos >= pos
              || pos < rangeList->ranges[hintCopy + 1].startPos
              + rangeList->ranges[hintCopy + 1].len)
          && (rangeList->ranges[hintCopy].startPos
              + rangeList->ranges[hintCopy].len <= pos))
  {
    ++hintCopy;
    if(hint)
      *hint = hintCopy;
    return rangeList->ranges + hintCopy;
  }
  else if(numRanges > 2)
  {
    Seqpos searchPos = pos;
    struct seqRange *searchRes =
      bsearch(&searchPos, rangeList->ranges + 1, rangeList->numRanges - 1,
              sizeof(struct seqRange), posSeqRangeNextCompare);
    if(searchRes && hint)
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
  if(!fwrite(&(rangeList->numRanges), sizeof(rangeList->numRanges), 1, fp))
    return 0;
  if(fwrite(rangeList->ranges, sizeof(struct seqRange),
            numRanges, fp) != numRanges)
    return 0;
  return 1;
}

/* FIXME: convert to platform-independent variant */
struct seqRangeList *
SRLReadFromStream(FILE *fp, Env *env)
{
  struct seqRangeList *newRangeList;
  size_t numRanges;
  assert(fp && env);
  newRangeList = env_ma_malloc(env, sizeof(struct seqRangeList));
  if(!fread(&(newRangeList->numRanges),
            sizeof(newRangeList->numRanges), 1, fp))
  {
    env_ma_free(newRangeList, env);
    return NULL;
  }
  numRanges = newRangeList->numRanges;
  newRangeList->ranges = env_ma_malloc(env, sizeof(struct seqRange) *
                                       (newRangeList->numRangesStorable
                                        = numRanges));
  if(fread(newRangeList->ranges, sizeof(struct seqRange),
           numRanges, fp) != numRanges)
  {
    deleteSeqRangeList(newRangeList, env);
    return NULL;
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
  if(rangeList->numRanges == 0)
    return;
  p = SRLFindPositionNext(rangeList, start, hint);
  /* no range overlapping or between start and end of sequence? */
  if(!p)
    return;
  /* iterate over ranges left */
  {
    Seqpos s = MAX(start, p->startPos);
    struct seqRange *maxRange = rangeList->ranges + rangeList->numRanges - 1;
    while(s <= end)
    {
      Seqpos overlap = MIN(p->startPos + p->len, end + 1) - s;
      occStore[p->sym] += overlap;
      if(p == maxRange)
        break;
      s = (++p)->startPos;
    }
  }
}

extern Seqpos
SRLSymbolCountInSeqRegion(struct seqRangeList *rangeList, Seqpos start,
                      Seqpos end, Symbol sym, seqRangeListSearchHint *hint)
{
  struct seqRange *p;
  if(rangeList->numRanges == 0)
    return 0;
  p = SRLFindPositionNext(rangeList, start, hint);
  if(p)
  {
    Seqpos symCount = 0;
    Seqpos s = MAX(start, p->startPos);
    struct seqRange *maxRange = rangeList->ranges + rangeList->numRanges - 1;
    while(s <= end)
    {
      if(p->sym == sym)
        symCount += MIN(p->startPos + p->len, end + 1) - s;
      if(p == maxRange)
        break;
      s = (++p)->startPos;
    }
    return symCount;
  }
  else
    return 0;
}

extern Seqpos
SRLAllSymbolsCountInSeqRegion(struct seqRangeList *rangeList, Seqpos start,
                              Seqpos end, seqRangeListSearchHint *hint)
{
  struct seqRange *p;
  if(rangeList->numRanges == 0)
    return 0;
  p = SRLFindPositionNext(rangeList, start, hint);
  if(p)
  {
    Seqpos symCount = 0;
    Seqpos s = MAX(start, p->startPos);
    struct seqRange *maxRange = rangeList->ranges + rangeList->numRanges - 1;
    while(s <= end)
    {
      symCount += MIN(p->startPos + p->len, end + 1) - s;
      if(p == maxRange)
        break;
      s = (++p)->startPos;
    }
    return symCount;
  }
  else
    return 0;
}

