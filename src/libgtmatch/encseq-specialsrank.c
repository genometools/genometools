/*
  Copyright (c) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <assert.h>
#include <limits.h>

#include "libgtcore/dataalign.h"
#include "libgtcore/ma.h"
#include "libgtcore/minmax.h"

#include "libgtmatch/encseq-specialsrank.h"
#include "libgtmatch/encseq-specialsrank-priv.h"

typedef struct specialsRankTable SpecialsRankTable;

static Seqpos
specialsRankFromSampleTable(const SpecialsRankLookup *rankTable, Seqpos pos);

static Seqpos
specialsRankFromTermPos(const SpecialsRankLookup *rankTable, Seqpos pos);

static inline struct specialsRankLookup *
allocSpecialsRankTable(const Encodedsequence *encseq, Seqpos lastSeqPos,
                       unsigned sampleIntervalLog2, Readmode readmode)
{
  struct specialsRankTable *rankTable;
  struct specialsRankLookup *ranker;
  Seqpos numSamples = (lastSeqPos >> sampleIntervalLog2) + 1;

  ranker = ma_malloc(offsetAlign(sizeof (*ranker), sizeof (Seqpos))
                     + numSamples * sizeof (Seqpos));
  rankTable = &ranker->implementationData.sampleTable;
  rankTable->rankSumSamples
    = (Seqpos *)((char *)ranker + offsetAlign(sizeof (*ranker),
                                              sizeof (Seqpos)));
  rankTable->sampleIntervalLog2 = sampleIntervalLog2;
  rankTable->sampleInterval = ((Seqpos)1) << sampleIntervalLog2;
  rankTable->readmode = readmode;
  rankTable->numSamples = numSamples;
  rankTable->scanState = newEncodedsequencescanstate();
  ranker->encseq = encseq;
  ranker->rankFunc = specialsRankFromSampleTable;
  return ranker;
}

static inline struct specialsRankLookup *
allocEmptySpecialsRankLookup(const Encodedsequence *encseq, Seqpos lastSeqPos)
{
  struct specialsRankLookup *ranker;
  ranker = ma_malloc(sizeof (*ranker));
  ranker->implementationData.lastSeqPos = lastSeqPos;
  ranker->rankFunc = specialsRankFromTermPos;
  ranker->encseq = encseq;
  return ranker;
}

static inline bool
nextRange(Sequencerange *range, Specialrangeiterator *sri,
          Readmode readmode, Seqpos seqLastPos)
{
  bool hasNextRange = nextspecialrangeiterator(range, sri);
  if (hasNextRange)
  {
    if (ISDIRREVERSE(readmode))
    {
      Seqpos temp = range->rightpos;
      range->rightpos = seqLastPos - range->leftpos;
      range->leftpos = seqLastPos - temp;
    }
  }
  else
    range->rightpos = seqLastPos + 1, range->leftpos = seqLastPos;
  return hasNextRange;
}

extern SpecialsRankLookup *
newSpecialsRankLookup(const Encodedsequence *encseq, Readmode readmode,
                     unsigned sampleIntervalLog2)
{
  struct specialsRankLookup *ranker;
  Seqpos seqLastPos, seqLen;
  Seqpos sampleInterval = ((Seqpos)1) << sampleIntervalLog2;
  assert(encseq);
  assert(sampleIntervalLog2 < sizeof (Seqpos) * CHAR_BIT);
  seqLastPos = getencseqtotallength(encseq);
  seqLen = seqLastPos + 1;
  if (hasspecialranges(encseq))
  {
    /* this sequence has some special characters */
    struct specialsRankTable *rankTable;
    Specialrangeiterator *sri;
    Seqpos *sample, *maxSample, sum = 0, pos = 0, nextSamplePos;
    Sequencerange range = { 0, 0 };
    ranker = allocSpecialsRankTable(encseq, seqLen, sampleIntervalLog2,
                                    readmode);
    rankTable = &ranker->implementationData.sampleTable;
    sri = newspecialrangeiterator(encseq, !ISDIRREVERSE(readmode));
    sample = rankTable->rankSumSamples;
    maxSample = sample + rankTable->numSamples;
    *sample++ = sum;
    nextSamplePos = sampleInterval;
    nextRange(&range, sri, readmode, seqLastPos);
    while (sample < maxSample)
    {
      while (pos < nextSamplePos)
      {
        pos = MIN(MAX(pos, range.leftpos), nextSamplePos);
        sum += MIN(range.rightpos - pos, nextSamplePos - pos);
        pos = MIN(range.rightpos, nextSamplePos);
        if (pos < nextSamplePos)
        {
          nextRange(&range, sri, readmode, seqLastPos);
        }
      }
      *sample++ = sum;
      nextSamplePos += sampleInterval;
    }
    freespecialrangeiterator(&sri);
  }
  else
  {
    /* While there is no special characters in this sequence, there
     * is of course the terminator */
    ranker = allocEmptySpecialsRankLookup(encseq, seqLastPos);
  }
  return ranker;
}

extern void
deleteSpecialsRankLookup(SpecialsRankLookup *ranker)
{
  if (ranker->rankFunc == specialsRankFromSampleTable)
    freeEncodedsequencescanstate(
      &ranker->implementationData.sampleTable.scanState);
  ma_free(ranker);
}

static Seqpos
specialsRankFromSampleTable(const SpecialsRankLookup *ranker, Seqpos pos)
{
  const SpecialsRankTable *rankTable = &ranker->implementationData.sampleTable;
  Seqpos rankCount, samplePos, encSeqLen;
  assert(ranker);
  encSeqLen = getencseqtotallength(ranker->encseq);
  assert(pos <= encSeqLen + 1);
  samplePos = pos & ~(rankTable->sampleInterval - 1);
  {
    size_t sampleIdx = pos >> rankTable->sampleIntervalLog2;
    rankCount = rankTable->rankSumSamples[sampleIdx];
  }
  {
    const Encodedsequence *encseq = ranker->encseq;
    Encodedsequencescanstate *esr = rankTable->scanState;
    Readmode readmode = rankTable->readmode;
    Seqpos encseqQueryMax = MIN(pos, encSeqLen);
    if (samplePos < encseqQueryMax)
    {
      initEncodedsequencescanstate(esr, encseq, readmode, samplePos);
      do {
        if (ISSPECIAL(sequentialgetencodedchar(encseq, esr, samplePos++,
                                               readmode)))
          ++rankCount;
      } while (samplePos < encseqQueryMax);
    }
    if (pos == encSeqLen + 1)
      ++rankCount;
  }
  return rankCount;
}

static Seqpos
specialsRankFromTermPos(const SpecialsRankLookup *ranker, Seqpos pos)
{
  assert(pos <= ranker->implementationData.lastSeqPos + 1);
  return ((pos == ranker->implementationData.lastSeqPos + 1)?1:0);
}

extern const Encodedsequence *
SPRTGetOrigEncseq(const SpecialsRankLookup *ranker)
{
  return ranker->encseq;
}
