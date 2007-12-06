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

#include <errno.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "libgtcore/dataalign.h"
#include "libgtcore/fa.h"
#include "libgtcore/filelengthvalues.h"
#include "libgtcore/minmax.h"
#include "libgtcore/seqiterator.h"
#include "libgtcore/str.h"
#include "libgtcore/strarray.h"
#include "libgtmatch/seqpos-def.h"
#include "libgtcore/symboldef.h"
#include "libgtmatch/verbose-def.h"
#include "libgtmatch/esafileend.h"
#include "libgtmatch/sfx-optdef.h"
#include "libgtmatch/intcode-def.h"
#include "libgtmatch/encseq-def.h"
#include "libgtmatch/sfx-cmpsuf.pr"

#include "libgtmatch/eis-encidxseq.h"
#include "libgtmatch/eis-suffixerator-interface.h"

struct sfxIReaderState
{
  Seqpos nextReadPos;
  int readFlag;
};

struct sfxInterface
{
  Suffixeratoroptions so;
  Measuretime *mtime;
  Seqpos length;
  const Alphabet *alpha;
  const Encodedsequence *encseq;
  struct seqStats *stats;
  Sfxiterator *sfi;
  DefinedSeqpos longest;
  bool specialsuffixes;
  /* data relevant to holding portions of the suffix array */
  Seqpos lastGeneratedLen, lastGeneratedStart;
  const Seqpos *lastGeneratedSufTabSegment;
  /**< pointer to part of suffix array
   * returned by last iterator call*/
  Seqpos *prevGeneratedSufTabSegments; /**< holds cache of suffix array
                                  * parts returned by previous to last
                                  * iterator calls as far back as
                                  * required by the readers */
  Seqpos prevGeneratedStart;
  size_t prevGeneratedSize,     /**< maximum number of Seqpos values
                                 * the cache may hold  */
    prevGeneratedLen;           /**< number of Seqpos values the
                                 * cache holds currently */
  /* data relevant to listeners */
  int allRequests;
  size_t numReaders;
  struct sfxIReaderState *readers;
};

static Seqpos
sfxIReadAdvance(sfxInterface *iface,
                Seqpos requestMaxPos,
                Error *err);

extern sfxInterface *
newSfxInterface(Suffixeratoroptions *so,
                const Encodedsequence *encseq,
                const Specialcharinfo *specialcharinfo,
                unsigned long numofsequences,
                Measuretime *mtime,
                Seqpos length,
                const Alphabet *alpha,
                const unsigned long *characterdistribution,
                Verboseinfo *verbosity,
                Error *err)
{
  return newSfxInterfaceWithReaders(so, 0, NULL, NULL, encseq,
                                    specialcharinfo, numofsequences, mtime,
                                    length, alpha, characterdistribution,
                                    verbosity, err);
}

static struct seqStats *
newSeqStatsFromCharDist(const Alphabet *alpha, Seqpos len, unsigned numOfSeqs,
                        const unsigned long *characterdistribution, Error *err)
{
  struct seqStats *stats = NULL;
  unsigned i, mapSize;
  Seqpos regularSymsSum = 0;
  stats = ma_malloc(offsetAlign(sizeof (*stats), sizeof (Seqpos))
                    + (UINT8_MAX + 1) * sizeof (Seqpos));
  stats->sourceAlphaType = sourceUInt8;
  stats->symbolDistributionTable =
    (Seqpos *)((char *)stats + offsetAlign(sizeof (*stats), sizeof (Seqpos)));
  memset(stats->symbolDistributionTable, 0, sizeof (Seqpos) * (UINT8_MAX + 1));
  mapSize = getmapsizeAlphabet(alpha);
  for (i = 0; i < mapSize - 1; ++i)
    regularSymsSum +=
      (stats->symbolDistributionTable[i] = characterdistribution[i]);
  stats->symbolDistributionTable[WILDCARD] = len - regularSymsSum - numOfSeqs;
  stats->symbolDistributionTable[SEPARATOR] += numOfSeqs;
  stats->symbolDistributionTable[UNDEFBWTCHAR] += 1;
  return stats;
}

static void
deleteSeqStats(struct seqStats *stats, Error *err)
{
  ma_free(stats);
}

#define INITOUTFILEPTR(PTR,FLAG,SUFFIX)                                 \
  ((FLAG)?(PTR = opensfxfile(so->str_indexname,SUFFIX,"wb",err)):((void*)1))

#define newSfxInterfaceWithReadersErrRet()        \
  do {                                            \
    if (iface->stats)                             \
      deleteSeqStats(iface->stats, err);          \
    if (iface->readers)                           \
      ma_free(iface->readers);                    \
    if (iface) ma_free(iface);                    \
  } while (0)

extern sfxInterface *
newSfxInterfaceWithReaders(Suffixeratoroptions *so,
                           size_t numReaders,
                           enum sfxDataRequest *requests,
                           listenerID *ids,
                           const Encodedsequence *encseq,
                           const Specialcharinfo *specialcharinfo,
                           unsigned long numofsequences,
                           Measuretime *mtime,
                           Seqpos length,
                           const Alphabet *alpha,
                           const unsigned long *characterdistribution,
                           Verboseinfo *verbosity, Error *err)
{
  sfxInterface *iface = NULL;

  error_check(err);

  iface = ma_calloc(1, sizeof (*iface));
  memcpy(&iface->so, so, sizeof (*so));
  iface->mtime = mtime;
  iface->length = length;
  iface->alpha = alpha;
  iface->encseq = encseq;
  iface->stats = newSeqStatsFromCharDist(iface->alpha, iface->length,
                                         numofsequences,
                                         characterdistribution, err);
  if (!(iface->sfi = newSfxiterator(specialcharinfo->specialcharacters,
                                    specialcharinfo->specialranges,
                                    encseq, so->readmode,
                                    getnumofcharsAlphabet(alpha),
                                    so->prefixlength,
                                    so->numofparts, 
                                    NULL,
                                    iface->mtime,
                                    verbosity, err)))
    newSfxInterfaceWithReadersErrRet();
  iface->longest.defined = false;
  iface->specialsuffixes = false;

  iface->prevGeneratedSize = iface->lastGeneratedStart
    = iface->prevGeneratedStart = 0;
  iface->lastGeneratedSufTabSegment = iface->prevGeneratedSufTabSegments = NULL;

  iface->allRequests = SFX_REQUEST_NONE;
  {
    size_t i;
    iface->numReaders = 0;
    iface->readers = NULL;
    for (i = 0; i < numReaders; ++i)
      if (!SfxIRegisterReader(iface, ids + i, requests[i], err))
        newSfxInterfaceWithReadersErrRet();
  }
  return iface;
}

extern void
deleteSfxInterface(sfxInterface *iface, Error *err)
{
  ma_free(iface->prevGeneratedSufTabSegments);
  freeSfxiterator(&iface->sfi);
  deleteSeqStats(iface->stats, err);
  ma_free(iface->readers);
  ma_free(iface);
}

const Uchar *
SfxIReadESQRange(sfxInterface *iface, Seqpos start, Seqpos len,
                 Uchar *dest)
{
  size_t i;
  assert(dest);
  for (i = 0; i < (size_t)len; ++i)
  {
    dest[i] = getencodedchar(iface->encseq, start + i, iface->so.readmode);
  }
  return dest;
}

extern const Alphabet *
getSfxIAlphabet(const sfxInterface *si)
{
  return si->alpha;
}

extern MRAEnc *
newMRAEncFromSfxI(const sfxInterface *si, Error *err)
{
  MRAEnc *alphabet;
  assert(si && err);
  alphabet = MRAEncGTAlphaNew(getSfxIAlphabet(si), err);
  MRAEncAddSymbolToRange(alphabet, SEPARATOR, 1);
  return alphabet;
}

Seqpos
getSfxILength(const sfxInterface *si)
{
  assert(si);
  return si->length;
}

extern const struct seqStats *
getSfxISeqStats(const sfxInterface *si)
{
  return si->stats;
}

extern DefinedSeqpos
getSfxILongestPos(const struct sfxInterface *si)
{
  return si->longest;
}

extern const Encodedsequence *
getSfxIEncSeq(const sfxInterface *si)
{
  return si->encseq;
}

int
SfxIRegisterReader(sfxInterface *iface, listenerID *id,
                   enum sfxDataRequest request, Error *err)
{
  size_t availId = iface->numReaders++;
  iface->readers = ma_realloc(iface->readers,
                              sizeof (iface->readers[0]) * (availId + 1));
  iface->readers[availId].nextReadPos = 0;
  iface->readers[availId].readFlag = request;
  iface->allRequests |= request;
  *id = availId;
  return 1;
}

static Seqpos
getSufTabVal(sfxInterface *iface, Seqpos pos, Error *err)
{
  assert(iface && pos < iface->length);
  while (1)
  {
    if (pos >= iface->lastGeneratedStart
       && pos < iface->lastGeneratedStart + iface->lastGeneratedLen)
    {
      /* pos is in last block */
      return iface->lastGeneratedSufTabSegment[pos
                                               - iface->lastGeneratedStart];
    }
    else if (pos >= iface->prevGeneratedStart
            && pos < iface->prevGeneratedStart + iface->prevGeneratedLen)
    {
      /* pos is in prev block */
      return iface->prevGeneratedSufTabSegments[pos
                                                - iface->prevGeneratedStart];
    }
    else
    {
      while (pos >= iface->lastGeneratedStart + iface->lastGeneratedLen)
        if (!sfxIReadAdvance(iface, pos, err))
          break;
    }
  }
}

extern int
SfxIGetOrigSeq(void *state, Symbol *dest, Seqpos pos, size_t len)
{
  struct sfxInterface *iface;
  size_t i;
  assert(state);
  iface = state;
  for (i = 0; i < len; ++i)
    dest[i] = getencodedchar(iface->encseq, pos + i, iface->so.readmode);
  return len;
}

extern size_t
readSfxIBWTRange(sfxInterface *iface, listenerID id, size_t len,
                 Uchar *dest, Error *err)
{
  size_t i, effLen;
  Seqpos start;
  assert(iface && id < iface->numReaders && dest);
  assert(iface->readers[id].readFlag == SFX_REQUEST_BWTTAB);
  start = iface->readers[id].nextReadPos;
  effLen = MIN((size_t)len, (size_t)(iface->length - start));
  for (i = 0; i < (size_t)effLen; ++i)
  {
    Seqpos sufTabVal = getSufTabVal(iface, start + i, err);
    if (sufTabVal != 0)
      dest[i] = getencodedchar(iface->encseq, sufTabVal - 1,
                               iface->so.readmode);
    else
      dest[i] = (Uchar) UNDEFBWTCHAR;
  }
  iface->readers[id].nextReadPos = start + effLen;
  return effLen;
}

extern size_t
readSfxILCPRange(sfxInterface *iface, listenerID id, size_t len,
                 Seqpos *dest, Error *err)
{
  size_t i, effLen;
  Seqpos start;
  assert(iface && id < iface->numReaders && dest);
  assert(iface->readers[id].readFlag == SFX_REQUEST_LCPTAB);
  start = iface->readers[id].nextReadPos;
  effLen = MIN(len, iface->length - start);
  if (start == 0)
    dest[0] = 0, i = 1;
  else
    i = 0;
  for (; i < (size_t)len; ++i)
  {
#ifndef NDEBUG
    int cmp =
#endif /* NDEBUG */
      comparetwosuffixes(iface->encseq,
                         iface->so.readmode,
                         dest + i,
                         false,
                         false,
                         0,
                         getSufTabVal(iface, start + i - 1, err),
                         getSufTabVal(iface, start + i, err),
                         NULL,  /* XXX FIXME */
                         NULL);  /* XXX FIXME */
#ifndef NDEBUG
    if (cmp > 0)
    {
      fprintf(stderr,"pos = " FormatSeqpos
              ": cmp " FormatSeqpos
              " " FormatSeqpos " = %d",
              PRINTSeqposcast(start + (Seqpos)i - 1),
              PRINTSeqposcast(start + (Seqpos)i),
              PRINTSeqposcast(getSufTabVal(iface, start + (Seqpos)i, err)),
              cmp);
      exit(EXIT_FAILURE);
    }
#endif /* NDEBUG */
  }
  iface->readers[id].nextReadPos = start + len - 1;
  return effLen;
}

size_t
readSfxISufTabRange(sfxInterface *iface, listenerID id, size_t len,
                    Seqpos *dest, Error *err)
{
  Seqpos start;
  assert(iface && id < iface->numReaders && dest);
  start = iface->readers[id].nextReadPos;
  while (1)
  {
    if (start >= iface->lastGeneratedStart
       && len <= iface->lastGeneratedLen)
    {
      /* chunk is completely contained in last block */
      iface->readers[id].nextReadPos = start + len;
      memcpy(dest,
             iface->lastGeneratedSufTabSegment
             + start - iface->lastGeneratedStart,
             sizeof (dest[0]) * len);
      return len;
    }
    else if (start >= iface->prevGeneratedStart
            && len <= iface->prevGeneratedLen)
    {
      /* chunk is completely contained in prev block */
      iface->readers[id].nextReadPos = start + len;
      memcpy(dest,
             iface->prevGeneratedSufTabSegments
             + start - iface->prevGeneratedStart,
             sizeof (dest[0]) * len);
      return len;
    }
    else if (start + len < iface->lastGeneratedStart + iface->lastGeneratedLen)
    {
      /* chunk overlaps both segments, but does not extend beyond
       * last segment */
      iface->readers[id].nextReadPos = start + len;
      memcpy(dest,
             iface->prevGeneratedSufTabSegments
             + start - iface->prevGeneratedStart,
             sizeof (dest[0])
             * (iface->prevGeneratedLen + iface->prevGeneratedStart - start));
      memcpy(dest,
             iface->lastGeneratedSufTabSegment,
             sizeof (dest[0])
             * (start + len - iface->lastGeneratedStart));
      return len;
    }
    else
    {
      while (start + len > iface->prevGeneratedStart + iface->prevGeneratedLen)
        if (!sfxIReadAdvance(iface, start + len, err))
        {
          /* Caution: sneakily updates the length parameter to obtain a
           * request we can fulfill */
          len = iface->lastGeneratedStart + iface->lastGeneratedLen - start;
          break;
        }
    }
  }
}

static Seqpos
findMinOpenRequest(sfxInterface *iface, int reqType)
{
  Seqpos min;
  size_t i;
  assert(iface);
  min = iface->length;
  for (i = 0; i < iface->numReaders; ++i)
    if (iface->readers[i].readFlag & reqType)
      min = MIN(min, iface->readers[i].nextReadPos);
  return min;
}

static Seqpos
sfxIReadAdvance(sfxInterface *iface,
                Seqpos requestMaxPos,
                Error *err)
{
  Seqpos requestMinPos, lengthOfExtension = 0;
  assert(iface && err);
  /* 1. find first position that still needs to be read */
  requestMinPos = findMinOpenRequest(iface, SFX_REQUEST_ANY);
  /* 2. move still unread old values as far as possible to head of copy */
  assert(requestMinPos >= iface->prevGeneratedStart);
  assert(requestMinPos <= iface->lastGeneratedStart + iface->lastGeneratedLen);
  if (requestMinPos < iface->prevGeneratedStart + iface->prevGeneratedLen)
  {
    size_t prevSegmentUnread
      = iface->prevGeneratedLen - requestMinPos + iface->prevGeneratedStart;
    memmove(iface->prevGeneratedSufTabSegments,
            iface->prevGeneratedSufTabSegments + requestMinPos
            - iface->prevGeneratedStart,
            sizeof (iface->prevGeneratedSufTabSegments[0]) *
            prevSegmentUnread);
    iface->prevGeneratedLen = prevSegmentUnread;
    iface->prevGeneratedStart = requestMinPos;
  }
  else
  {
    iface->prevGeneratedLen = 0;
    iface->prevGeneratedStart = requestMinPos;
  }
  do
  {
    Seqpos copyStartPos = MAX(requestMinPos, iface->lastGeneratedStart);
    size_t copyLen = iface->lastGeneratedLen
      - copyStartPos + iface->lastGeneratedStart;
    /* 3. extend cache to also accept all values in last region read */
    if (copyLen)
    {
      size_t prevGeneratedSizeLeft = iface->prevGeneratedSize
        - iface->prevGeneratedLen;
      if (copyLen > prevGeneratedSizeLeft)
      {
        size_t newSize = iface->prevGeneratedLen + copyLen;
        iface->prevGeneratedSufTabSegments =
          ma_realloc(iface->prevGeneratedSufTabSegments,
                     sizeof (iface->lastGeneratedSufTabSegment[0]) * newSize);
        iface->prevGeneratedSize = newSize;
      }
      memcpy(iface->prevGeneratedSufTabSegments + iface->prevGeneratedLen,
             iface->lastGeneratedSufTabSegment
             + copyStartPos - iface->lastGeneratedStart,
             sizeof (iface->lastGeneratedSufTabSegment[0]) * copyLen);
      iface->prevGeneratedLen += copyLen;
    }
    /* 4. read next region of sequence by calling nextSfxIterator */
    iface->lastGeneratedStart += iface->lastGeneratedLen;
    if ((iface->lastGeneratedSufTabSegment =
        nextSfxiterator(&iface->lastGeneratedLen, &iface->specialsuffixes,
                        iface->mtime, iface->sfi, err)))
    {
      size_t pos;
      /* size_t because the current approach cannot generate more
       * than memory will hold anyway */
      size_t lastGeneratedLen = iface->lastGeneratedLen;
      const Seqpos *suftab = iface->lastGeneratedSufTabSegment;
      lengthOfExtension += iface->lastGeneratedLen;
      if (!iface->longest.defined)
        for (pos=0; pos < lastGeneratedLen; pos++)
        {
          if (suftab[pos] == 0)
          {
            iface->longest.defined = true;
            iface->longest.valueseqpos = iface->lastGeneratedStart + pos;
          }
        }
    /* uncomment this to reenable synchronous writing of tables */
/*if (iface->lastGeneratedSufTabSegment == NULL */
/*    || suftab2file(&iface->outfileinfo, iface->lastGeneratedSufTabSegment, */
/*                   iface->so.readmode, iface->lastGeneratedLen, err) != 0) */
/*       break; */
    }
    /* 5. if positions in region don't suffice go back to step 3. */
  } while (requestMaxPos
           >= iface->lastGeneratedStart + iface->lastGeneratedLen);
  return lengthOfExtension;
}
