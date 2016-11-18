/*
  Copyright (c) 2016-2016 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2016-2016 Center for Bioinformatics, University of Hamburg

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

#include "core/range_api.h"
#include "match/seed_extend_parts.h"

typedef struct
{
  GtUword start,      /* index number of first sequence in part */
          end,        /* index number of last sequence in part */
          max_length; /* length of longest sequence in range */
} GtSequenceRangeWithMaxLength;

struct GtSequencePartsInfo
{
  GtSequenceRangeWithMaxLength *ranges;
  GtUword *ssptab,
          same_length,
          totallength,
          numofsequences,
          parts_number;
};

GtUword gt_sequence_parts_info_seqstartpos(const GtSequencePartsInfo *spi,
                                           GtUword seqnum)
{
  if (spi->ssptab != NULL)
  {
    return seqnum > 0 ? (spi->ssptab[seqnum - 1] + 1) : 0;
  } else
  {
    return seqnum * (spi->same_length + 1);
  }
}

GtUword gt_sequence_parts_info_seqendpos(const GtSequencePartsInfo *spi,
                                         GtUword seqnum)
{
  if (spi->ssptab != NULL)
  {
    return seqnum < spi->numofsequences - 1 ? (spi->ssptab[seqnum] - 1)
                                            : (spi->totallength - 1);
  } else
  {
    return (seqnum + 1) * (spi->same_length + 1) - 2;
  }
}

GtUword gt_sequence_parts_info_partlength(const GtSequencePartsInfo *spi,
                                          GtUword fromseq,
                                          GtUword toseq)
{
  return gt_sequence_parts_info_seqendpos(spi,toseq) -
         gt_sequence_parts_info_seqstartpos(spi,fromseq) + 1;
}

static GtUword gt_encseq_next_larger_width(const GtSequencePartsInfo *spi,
                                           GtUword startseqnum,
                                           GtUword width,
                                           GtUword numofsequences)
{
  GtUword left, right, found = GT_UWORD_MAX,
          start_segment = gt_sequence_parts_info_seqstartpos(spi,startseqnum);

  left = startseqnum;
  gt_assert(numofsequences > 0);
  right = numofsequences - 1;
  while (left <= right)
  {
    GtUword mid = left + (right - left + 1)/2, mid_end, this_width;

    gt_assert(mid < numofsequences);
    mid_end = gt_sequence_parts_info_seqendpos(spi,mid);
    gt_assert(mid_end > start_segment);
    this_width = mid_end - start_segment;
    if (this_width > width)
    {
      found = mid;
      if (right == 0)
      {
        break;
      }
      right = mid - 1;
    } else
    {
      if (left == numofsequences - 1)
      {
        break;
      }
      left = mid + 1;
    }
  }
  return found;
}

GtSequencePartsInfo *gt_sequence_parts_info_new(const GtEncseq *encseq,
                                                GtUword numofsequences,
                                                GtUword numparts)
{
  GtSequencePartsInfo *spi = gt_malloc(sizeof *spi);
  spi->ranges = gt_malloc(sizeof *spi->ranges * numparts);
  spi->ssptab = gt_all_sequence_separators_get(encseq);
  spi->totallength = gt_encseq_total_length(encseq);
  spi->numofsequences = gt_encseq_num_of_sequences(encseq);

  if (spi->ssptab == NULL)
  {
    spi->same_length = gt_encseq_seqlength(encseq,0);
  } else
  {
    spi->same_length = 0;
  }
  if (numparts >= numofsequences)
  { /* assign one seq for each part */
    GtUword idx;

    for (idx = 0; idx < numofsequences; ++idx) {
      spi->ranges[idx].start = spi->ranges[idx].end = idx;
      spi->ranges[idx].max_length
        = gt_sequence_parts_info_partlength(spi,idx,idx);
    }
    spi->parts_number = numofsequences;
  } else
  {
    GtUword seqnum, idx, effective_num_parts;
    const GtUword partwidth = spi->totallength/numparts;

    for (idx = 0, seqnum = 0; idx < numparts && seqnum < numofsequences; idx++)
    {
      const GtUword seqnum_next_width
        = gt_encseq_next_larger_width(spi,seqnum,partwidth,
                                      numofsequences);
      spi->ranges[idx].start = seqnum;
      if (seqnum_next_width == GT_UWORD_MAX)
      {
        spi->ranges[idx].end = numofsequences - 1;
        idx++;
        break;
      }
      spi->ranges[idx].end = seqnum_next_width;
      seqnum = seqnum_next_width + 1;
    }
    gt_assert(idx > 0 &&
              spi->ranges[idx-1].end == numofsequences - 1);
    effective_num_parts = idx;
    if (effective_num_parts == 1)
    {
      spi->ranges[0].max_length = gt_encseq_max_seq_length(encseq);
    } else
    {
      if (spi->ssptab == NULL)
      {
        for (idx = 0; idx < effective_num_parts; idx++)
        {
          spi->ranges[idx].max_length = spi->same_length;
        }
      } else
      {
        GtUword currentstart = 0, idx = 0, currentlength, maxlength = 0;

        gt_assert(numofsequences > 1);
        for (seqnum = 0; seqnum < numofsequences - 1; seqnum++)
        {
          gt_assert(currentstart < spi->ssptab[seqnum]);
          currentlength = spi->ssptab[seqnum] - currentstart;
          if (maxlength < currentlength)
          {
            maxlength = currentlength;
          }
          if (seqnum == spi->ranges[idx].end)
          {
            spi->ranges[idx].max_length = maxlength;
            idx++;
            maxlength = 0;
          }
          currentstart = spi->ssptab[seqnum] + 1;
        }
        currentlength = spi->totallength - currentstart;
        if (maxlength < currentlength)
        {
          maxlength = currentlength;
        }
        gt_assert(idx + 1 == effective_num_parts &&
                  seqnum == spi->ranges[idx].end);
        spi->ranges[idx].max_length = maxlength;
      }
    }
    spi->parts_number = effective_num_parts;
  }
  return spi;
}

GtUword gt_sequence_parts_info_number(const GtSequencePartsInfo *spi)
{
  gt_assert(spi != NULL);
  return spi->parts_number;
}

void gt_sequence_parts_info_delete(GtSequencePartsInfo *spi)
{
  if (spi != NULL)
  {
    gt_free(spi->ssptab);
    gt_free(spi->ranges);
    gt_free(spi);
  }
}

GtUword gt_sequence_parts_info_start_get(const GtSequencePartsInfo *spi,
                                         GtUword idx)
{
  gt_assert(spi != NULL && idx < spi->parts_number);
  return spi->ranges[idx].start;
}

GtUword gt_sequence_parts_info_end_get(const GtSequencePartsInfo *spi,
                                       GtUword idx)
{
  gt_assert(spi != NULL && idx < spi->parts_number);
  return spi->ranges[idx].end;
}

GtUword gt_sequence_parts_info_numofsequences_get(
                        const GtSequencePartsInfo *spi,GtUword idx)
{
  gt_assert(spi != NULL && idx < spi->parts_number);
  return gt_sequence_parts_info_end_get(spi,idx) -
         gt_sequence_parts_info_start_get(spi,idx) + 1;
}

GtUword gt_sequence_parts_info_max_length_get(const GtSequencePartsInfo *spi,
                                              GtUword idx)
{
  gt_assert(spi != NULL && idx < spi->parts_number);
  return spi->ranges[idx].max_length;
}

bool gt_sequence_parts_info_overlap(const GtSequencePartsInfo *spia,
                                    GtUword aidx,
                                    const GtSequencePartsInfo *spib,
                                    GtUword bidx)
{
  GtRange range_a, range_b;

  range_a.start = gt_sequence_parts_info_start_get(spia,aidx);
  range_a.end = gt_sequence_parts_info_end_get(spia,aidx);
  range_b.start = gt_sequence_parts_info_start_get(spib,bidx);
  range_b.end = gt_sequence_parts_info_end_get(spib,bidx);
  return gt_range_overlap(&range_a,&range_b);
}

bool gt_sequence_parts_info_equal(const GtSequencePartsInfo *spia,
                                  GtUword aidx,
                                  const GtSequencePartsInfo *spib,
                                  GtUword bidx)
{
  return (gt_sequence_parts_info_start_get(spia,aidx) ==
          gt_sequence_parts_info_start_get(spib,bidx)) &&
         (gt_sequence_parts_info_end_get(spia,aidx) ==
          gt_sequence_parts_info_end_get(spib,bidx)) ? true : false;

}

void gt_sequence_parts_info_variance_show(const GtSequencePartsInfo *spi)
{
  GtUword variance_sum = 0, idx, avgpartlength;

  avgpartlength = spi->totallength/spi->parts_number;
  for (idx = 0; idx < spi->parts_number; idx++)
  {
    const GtUword partlength
      = gt_sequence_parts_info_partlength(
             spi,
             gt_sequence_parts_info_start_get(spi,idx),
             gt_sequence_parts_info_end_get(spi,idx));
    if (partlength > avgpartlength)
    {
      GtUword diff = partlength - avgpartlength;
      variance_sum += diff * diff;
    } else
    {
      GtUword diff = avgpartlength - partlength;
      variance_sum += diff * diff;
    }
    printf("# Part " GT_WU ": sequence " GT_WU "..." GT_WU ", total length="
           GT_WU ", max_length=" GT_WU "\n",idx+1,
           gt_sequence_parts_info_start_get(spi,idx),
           gt_sequence_parts_info_end_get(spi,idx),
           partlength,
           gt_sequence_parts_info_max_length_get(spi,idx));
  }
  printf("# Variance of parts is %.2e\n",
           (double) variance_sum/spi->parts_number);
}

GtUchar *gt_sequence_parts_info_seq_extract(const GtEncseq *encseq,
                                            const GtSequencePartsInfo *spi,
                                            GtUword idx)
{
  GtUchar *byte_sequence;
  const GtUword
    firstseqnum = gt_sequence_parts_info_start_get(spi,idx),
    lastseqnum = gt_sequence_parts_info_end_get(spi,idx),
    firstpos = gt_sequence_parts_info_seqstartpos(spi,firstseqnum),
    lastpos = gt_sequence_parts_info_seqendpos(spi,lastseqnum);

  gt_assert(firstpos <= lastpos);
  byte_sequence = gt_malloc(sizeof *byte_sequence * (lastpos - firstpos + 1));
  gt_encseq_extract_encoded(encseq,byte_sequence,firstpos,lastpos);
  return byte_sequence;
}
