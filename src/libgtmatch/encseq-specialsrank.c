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

struct specialsRankTable
{
  Encodedsequence *encseq;
  Encodedsequencescanstate *scanState;
  Seqpos *rankSumSamples, numSamples, sampleInterval;
  Readmode readmode;
  unsigned sampleIntervalLog2;
};

static inline struct specialsRankTable *
allocSpecialsRankTable(Encodedsequence *encseq, Seqpos seqLen,
                       unsigned sampleIntervalLog2, Readmode readmode)
{
  struct specialsRankTable *rankTable;
  Seqpos numSamples = (seqLen >> sampleIntervalLog2) + 1;

  rankTable = ma_malloc(offsetAlign(sizeof (SpecialsRankTable),
                                    sizeof (Seqpos))
                        + numSamples * sizeof (Seqpos));
  rankTable->rankSumSamples
    = (Seqpos *)((char *)rankTable + offsetAlign(sizeof (SpecialsRankTable),
                                                 sizeof (Seqpos)));
  rankTable->sampleIntervalLog2 = sampleIntervalLog2;
  rankTable->sampleInterval = ((Seqpos)1) << sampleIntervalLog2;
  rankTable->readmode = readmode;
  rankTable->encseq = encseq;
  rankTable->numSamples = numSamples;
  rankTable->scanState = newEncodedsequencescanstate();
  return rankTable;
}

extern SpecialsRankTable *
newSpecialsRankTable(Encodedsequence *encseq, Readmode readmode,
                     unsigned sampleIntervalLog2)
{
  struct specialsRankTable *rankTable;
  Seqpos seqLen;
  Seqpos sampleInterval = ((Seqpos)1) << sampleIntervalLog2;
  assert(encseq);
  assert(sampleIntervalLog2 < sizeof (Seqpos) * CHAR_BIT);
  seqLen = getencseqtotallength(encseq);
  if (hasspecialranges(encseq))
  {
    Specialrangeiterator *sri;
    Seqpos *sample, *maxSample, sum = 0, pos = 0, nextSamplePos;
    Sequencerange range = { 0, 0 };
    rankTable
      = allocSpecialsRankTable(encseq, seqLen, sampleIntervalLog2, readmode);
    sri = newspecialrangeiterator(encseq, !ISDIRREVERSE(readmode));
    sample = rankTable->rankSumSamples;
    maxSample = sample + rankTable->numSamples;
    *sample++ = sum;
    nextSamplePos = sampleInterval;
    nextspecialrangeiterator(&range, sri);
    while (sample < maxSample)
    {
      while (pos < nextSamplePos)
      {
        pos = MIN(MAX(pos, range.leftpos), nextSamplePos);
        sum += MIN(range.rightpos - pos, nextSamplePos - pos);
        pos = MIN(range.rightpos, nextSamplePos);
        if (pos < nextSamplePos)
        {
          if (!nextspecialrangeiterator(&range, sri))
            range.rightpos = range.leftpos = seqLen;
        }
      }
      *sample++ = sum;
      nextSamplePos += sampleInterval;
    }
  }
  else
  {
    /* implementation pending */
    abort();
  }
  return rankTable;
}

extern void
deleteSpecialsRankTable(SpecialsRankTable *table)
{
  freeEncodedsequencescanstate(&table->scanState);
  ma_free(table);
}

extern Seqpos
specialsRank(const SpecialsRankTable *rankTable, Seqpos pos)
{
  Seqpos rankCount, samplePos;
  assert(rankTable);
  assert(pos <= getencseqtotallength(rankTable->encseq));
  samplePos = pos & ~(rankTable->sampleInterval - 1);
  {
    size_t sampleIdx = pos >> rankTable->sampleIntervalLog2;
    rankCount = rankTable->rankSumSamples[sampleIdx];
  }
  {
    Encodedsequence *encseq = rankTable->encseq;
    Encodedsequencescanstate *esr = rankTable->scanState;
    initEncodedsequencescanstate(esr, encseq,
                                 rankTable->readmode, samplePos);
    while (samplePos < pos)
      if (ISSPECIAL(sequentialgetencodedchar(encseq, esr, samplePos++)))
        ++rankCount;
  }
  return rankCount;
}

extern const Encodedsequence *
SPRTGetOrigEncseq(const SpecialsRankTable *rankTable)
{
  return rankTable->encseq;
}
