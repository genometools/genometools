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

#include <limits.h>

#include "core/assert_api.h"
#include "match/dataalign.h"
#include "core/ma_api.h"
#include "core/minmax.h"

#include "match/eis-specialsrank.h"
#include "match/eis-specialsrank-priv.h"

typedef struct specialsRankTable SpecialsRankTable;

static GtUword
specialsRankFromSampleTable(const SpecialsRankLookup *rankTable,
                            GtUword pos);

static GtUword
specialsRankFromTermPos(const SpecialsRankLookup *rankTable,
                        GtUword pos);

static inline struct specialsRankLookup *
allocSpecialsRankTable(const GtEncseq *encseq,
                       GtUword lastSeqPos,
                       unsigned sampleIntervalLog2,
                       GtReadmode readmode)
{
  struct specialsRankTable *rankTable;
  struct specialsRankLookup *ranker;
  GtUword numSamples = (lastSeqPos >> sampleIntervalLog2) + 1;

  ranker = gt_malloc(offsetAlign(sizeof (*ranker), sizeof (GtUword))
                     + numSamples * sizeof (GtUword));
  rankTable = &ranker->implementationData.sampleTable;
  rankTable->rankSumSamples
    = (GtUword *)((char *)ranker + offsetAlign(sizeof (*ranker),
                                              sizeof (GtUword)));
  rankTable->sampleIntervalLog2 = sampleIntervalLog2;
  rankTable->sampleInterval = ((GtUword)1) << sampleIntervalLog2;
  rankTable->readmode = readmode;
  rankTable->numSamples = numSamples;
  rankTable->scanState = gt_encseq_create_reader_with_readmode(encseq,
                                                               readmode, 0);
  ranker->encseq = encseq;
  ranker->rankFunc = specialsRankFromSampleTable;
  return ranker;
}

static inline struct specialsRankLookup *
allocEmptySpecialsRankLookup(const GtEncseq *encseq,
                             GtUword lastSeqPos)
{
  struct specialsRankLookup *ranker;
  ranker = gt_malloc(sizeof (*ranker));
  ranker->implementationData.lastSeqPos = lastSeqPos;
  ranker->rankFunc = specialsRankFromTermPos;
  ranker->encseq = encseq;
  return ranker;
}

static inline bool
nextRange(GtRange *range, GtSpecialrangeiterator *sri,
          GtReadmode readmode, GtUword seqLastPos)
{
  bool hasNextRange = gt_specialrangeiterator_next(sri, range);
  if (hasNextRange)
  {
    if (GT_ISDIRREVERSE(readmode))
    {
      GtUword temp = range->end;
      range->end = seqLastPos - range->start;
      range->start = seqLastPos - temp;
    }
  }
  else
    range->end = seqLastPos + 1, range->start = seqLastPos;
  return hasNextRange;
}

SpecialsRankLookup *
gt_newSpecialsRankLookup(const GtEncseq *encseq, GtReadmode readmode,
                         unsigned sampleIntervalLog2)
{
  struct specialsRankLookup *ranker;
  GtUword seqLastPos, seqLen;
  GtUword sampleInterval = ((GtUword)1) << sampleIntervalLog2;
  gt_assert(encseq);
  gt_assert(sampleIntervalLog2 < sizeof (GtUword) * CHAR_BIT);
  seqLastPos = gt_encseq_total_length(encseq);
  seqLen = seqLastPos + 1;
  if (gt_encseq_has_specialranges(encseq))
  {
    /* this sequence has some special characters */
    struct specialsRankTable *rankTable;
    GtSpecialrangeiterator *sri;
    GtUword *sample, *maxSample, sum = 0, pos = 0, nextSamplePos;
    GtRange range = { 0, 0 };
    ranker = allocSpecialsRankTable(encseq, seqLen, sampleIntervalLog2,
                                    readmode);
    rankTable = &ranker->implementationData.sampleTable;
    sri = gt_specialrangeiterator_new(encseq, !GT_ISDIRREVERSE(readmode));
    sample = rankTable->rankSumSamples;
    maxSample = sample + rankTable->numSamples;
    *sample++ = sum;
    nextSamplePos = sampleInterval;
    nextRange(&range, sri, readmode, seqLastPos);
    while (sample < maxSample)
    {
      while (pos < nextSamplePos)
      {
        pos = MIN(MAX(pos, range.start), nextSamplePos);
        sum += MIN(range.end - pos, nextSamplePos - pos);
        pos = MIN(range.end, nextSamplePos);
        if (pos < nextSamplePos)
        {
          nextRange(&range, sri, readmode, seqLastPos);
        }
      }
      *sample++ = sum;
      nextSamplePos += sampleInterval;
    }
    gt_specialrangeiterator_delete(sri);
    sri = NULL;
  }
  else
  {
    /* While there is no special characters in this sequence, there
     * is of course the terminator */
    ranker = allocEmptySpecialsRankLookup(encseq, seqLastPos);
  }
  return ranker;
}

void
gt_deleteSpecialsRankLookup(SpecialsRankLookup *ranker)
{
  if (ranker->rankFunc == specialsRankFromSampleTable)
    gt_encseq_reader_delete(ranker->implementationData.sampleTable.scanState);
  gt_free(ranker);
}

static GtUword
specialsRankFromSampleTable(const SpecialsRankLookup *ranker, GtUword pos)
{
  const SpecialsRankTable *rankTable = &ranker->implementationData.sampleTable;
  GtUword rankCount, samplePos, encSeqLen;
  gt_assert(ranker);
  encSeqLen = gt_encseq_total_length(ranker->encseq);
  gt_assert(pos <= encSeqLen + 1);
  samplePos = pos & ~(rankTable->sampleInterval - 1);
  {
    size_t sampleIdx = pos >> rankTable->sampleIntervalLog2;
    rankCount = rankTable->rankSumSamples[sampleIdx];
  }
  {
    const GtEncseq *encseq = ranker->encseq;
    GtEncseqReader *esr = rankTable->scanState;
    GtReadmode readmode = rankTable->readmode;
    GtUword encseqQueryMax = MIN(pos, encSeqLen);
    if (samplePos < encseqQueryMax)
    {
      gt_encseq_reader_reinit_with_readmode(esr, encseq, readmode, samplePos);
      do {
        samplePos++;
        if (ISSPECIAL(gt_encseq_reader_next_encoded_char(esr)))
          ++rankCount;
      } while (samplePos < encseqQueryMax);
    }
    if (pos == encSeqLen + 1)
      ++rankCount;
  }
  return rankCount;
}

static GtUword
specialsRankFromTermPos(const SpecialsRankLookup *ranker, GtUword pos)
{
  gt_assert(pos <= ranker->implementationData.lastSeqPos + 1);
  return ((pos == ranker->implementationData.lastSeqPos + 1)?1:0);
}

const GtEncseq *
gt_SPRTGetOrigEncseq(const SpecialsRankLookup *ranker)
{
  return ranker->encseq;
}
