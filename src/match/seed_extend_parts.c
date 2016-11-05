#include "match/seed_extend_parts.h"

struct GtSequencePartsInfo
{
  GtSequenceRangeWithMaxLength *ranges;
  GtUword numofranges;
};

static GtUword gt_encseq_next_larger_width(const GtEncseq *encseq,
                                           GtUword startseqnum,
                                           GtUword width,
                                           GtUword numofsequences)
{
  GtUword left, right, found = GT_UWORD_MAX,
          start_segment = gt_encseq_seqstartpos(encseq,startseqnum);

  left = startseqnum;
  gt_assert(numofsequences > 0);
  right = numofsequences - 1;
  while (left <= right)
  {
    GtUword mid = left + (right - left + 1)/2, mid_end, this_width;

    gt_assert(mid < numofsequences);
    mid_end = gt_encseq_seqstartpos(encseq,mid) +
              gt_encseq_seqlength(encseq,mid) - 1;
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

GtSequencePartsInfo *gt_seed_extend_parts_new(const GtEncseq *encseq,
                                              GtUword numofsequences,
                                              GtUword numparts)
{
  GtSequencePartsInfo *seqrangeinfo = gt_malloc(sizeof *seqrangeinfo);
  seqrangeinfo->ranges = gt_malloc(sizeof *seqrangeinfo->ranges * numparts);

  if (numparts >= numofsequences) { /* assign one seq for each part */
    GtUword idx;

    for (idx = 0; idx < numofsequences; ++idx) {
      seqrangeinfo->ranges[idx].start = seqrangeinfo->ranges[idx].end = idx;
      seqrangeinfo->ranges[idx].max_length = gt_encseq_seqlength(encseq,idx);
      seqrangeinfo->ranges[idx].partnum = idx;
    }
    seqrangeinfo->numofranges = numofsequences;
  } else
  {
    GtUword seqnum, idx, effective_num_parts;
    const GtUword totallength = gt_encseq_total_length(encseq),
                  partwidth = totallength/numparts;

    for (idx = 0, seqnum = 0; idx < numparts && seqnum < numofsequences; idx++)
    {
      const GtUword seqnum_next_width
        = gt_encseq_next_larger_width(encseq,seqnum,partwidth,numofsequences);
      seqrangeinfo->ranges[idx].start = seqnum;
      seqrangeinfo->ranges[idx].partnum = idx;
      if (seqnum_next_width == GT_UWORD_MAX)
      {
        seqrangeinfo->ranges[idx].end = numofsequences - 1;
        idx++;
        break;
      }
      seqrangeinfo->ranges[idx].end = seqnum_next_width;
      seqnum = seqnum_next_width + 1;
    }
    gt_assert(idx > 0 &&
              seqrangeinfo->ranges[idx-1].end == numofsequences - 1);
    effective_num_parts = idx;
    if (effective_num_parts == 1)
    {
      seqrangeinfo->ranges[0].max_length = gt_encseq_max_seq_length(encseq);
    } else
    {
      GtUword *ssptab = gt_all_sequence_separators_get(encseq);
      if (ssptab == NULL)
      {
        GtUword samelength = gt_encseq_seqlength(encseq,0);
        for (idx = 0; idx < effective_num_parts; idx++)
        {
          seqrangeinfo->ranges[idx].max_length = samelength;
        }
      } else
      {
        GtUword currentstart = 0, idx = 0, currentlength, maxlength = 0;

        gt_assert(numofsequences > 1);
        for (seqnum = 0; seqnum < numofsequences - 1; seqnum++)
        {
          gt_assert(currentstart < ssptab[seqnum]);
          currentlength = ssptab[seqnum] - currentstart;
          if (maxlength < currentlength)
          {
            maxlength = currentlength;
          }
          if (seqnum == seqrangeinfo->ranges[idx].end)
          {
            seqrangeinfo->ranges[idx].max_length = maxlength;
            idx++;
            maxlength = 0;
          }
          currentstart = ssptab[seqnum] + 1;
        }
        currentlength = totallength - currentstart;
        if (maxlength < currentlength)
        {
          maxlength = currentlength;
        }
        gt_assert(idx + 1 == effective_num_parts &&
                  seqnum == seqrangeinfo->ranges[idx].end);
        seqrangeinfo->ranges[idx].max_length = maxlength;
        gt_free(ssptab);
      }
    }
    seqrangeinfo->numofranges = effective_num_parts;
  }
  return seqrangeinfo;
}

void gt_seed_extend_parts_delete(GtSequencePartsInfo *seqrangeinfo)
{
  if (seqrangeinfo != NULL)
  {
    gt_free(seqrangeinfo->ranges);
    gt_free(seqrangeinfo);
  }
}

GtUword gt_seed_extend_parts_number(const GtSequencePartsInfo *seqrangeinfo)
{
  gt_assert(seqrangeinfo != NULL);
  return seqrangeinfo->numofranges;
}

const GtSequenceRangeWithMaxLength *gt_seed_extend_parts_get(
                        const GtSequencePartsInfo *seqrangeinfo,
                        GtUword idx)
{
  gt_assert(seqrangeinfo != NULL && idx < seqrangeinfo->numofranges);
  return seqrangeinfo->ranges + idx;
}

void gt_seed_extend_parts_variance_show(const GtSequencePartsInfo *seqrangeinfo,
                                        const GtEncseq *encseq)
{
  GtUword variance_sum = 0, idx, avgpartsize;

  gt_assert(encseq != NULL);
  avgpartsize = gt_encseq_total_length(encseq)/seqrangeinfo->numofranges;
  for (idx = 0; idx < seqrangeinfo->numofranges; idx++)
  {
    GtUword segmentstart, segmentend, width;

    segmentstart = gt_encseq_seqstartpos(encseq,
                                         seqrangeinfo->ranges[idx].start);
    segmentend = gt_encseq_seqstartpos(encseq,seqrangeinfo->ranges[idx].end) +
                 gt_encseq_seqlength(encseq,seqrangeinfo->ranges[idx].end) - 1;
    gt_assert(segmentstart < segmentend);
    width = segmentend - segmentstart + 1;
    if (width > avgpartsize)
    {
      GtUword diff = width - avgpartsize;
      variance_sum += diff * diff;
    } else
    {
      GtUword diff = avgpartsize - width;
      variance_sum += diff * diff;
    }
    printf("# Part " GT_WU ": sequence " GT_WU "..." GT_WU ", total length="
           GT_WU ", max_length=" GT_WU "\n",idx+1,
                                            seqrangeinfo->ranges[idx].start,
           seqrangeinfo->ranges[idx].end,segmentend - segmentstart + 1,
           seqrangeinfo->ranges[idx].max_length);
  }
  printf("# Variance of parts is %.2e\n",
           (double) variance_sum/seqrangeinfo->numofranges);
}
