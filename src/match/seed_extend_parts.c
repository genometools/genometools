#include "match/seed_extend_parts.h"

struct GtSequencePartsInfo
{
  GtSequenceRangeWithMaxLength *ranges;
  GtUword *ssptab,
          same_length,
          totallength,
          numofsequences,
          parts_number;
};

static GtUword gt_sspbased_seqstartpos(const GtSequencePartsInfo *seqrangeinfo,
                                       GtUword seqnum)
{
  if (seqrangeinfo->ssptab != NULL)
  {
    return seqnum > 0 ? (seqrangeinfo->ssptab[seqnum - 1] + 1) : 0;
  } else
  {
    return seqnum * (seqrangeinfo->same_length + 1);
  }
}

static GtUword gt_sspbased_seqendpos(const GtSequencePartsInfo *seqrangeinfo,
                                     GtUword seqnum)
{
  if (seqrangeinfo->ssptab != NULL)
  {
    return seqnum < seqrangeinfo->numofsequences - 1
             ? (seqrangeinfo->ssptab[seqnum] - 1)
             : (seqrangeinfo->totallength - 1);
  } else
  {
    return (seqnum + 1) * (seqrangeinfo->same_length + 1) - 2;
  }
}

static GtUword gt_sspbased_seqlength(const GtSequencePartsInfo *seqrangeinfo,
                                     GtUword seqnum)
{
  return gt_sspbased_seqendpos(seqrangeinfo,seqnum) -
         gt_sspbased_seqstartpos(seqrangeinfo,seqnum) + 1;
}

static GtUword gt_encseq_next_larger_width(const GtSequencePartsInfo
                                             *seqrangeinfo,
                                           GtUword startseqnum,
                                           GtUword width,
                                           GtUword numofsequences)
{
  GtUword left, right, found = GT_UWORD_MAX,
          start_segment = gt_sspbased_seqstartpos(seqrangeinfo,startseqnum);

  left = startseqnum;
  gt_assert(numofsequences > 0);
  right = numofsequences - 1;
  while (left <= right)
  {
    GtUword mid = left + (right - left + 1)/2, mid_end, this_width;

    gt_assert(mid < numofsequences);
    mid_end = gt_sspbased_seqendpos(seqrangeinfo,mid);
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
  seqrangeinfo->ssptab = gt_all_sequence_separators_get(encseq);
  seqrangeinfo->totallength = gt_encseq_total_length(encseq);
  seqrangeinfo->numofsequences = gt_encseq_num_of_sequences(encseq);

  if (seqrangeinfo->ssptab == NULL)
  {
    seqrangeinfo->same_length = gt_encseq_seqlength(encseq,0);
  } else
  {
    seqrangeinfo->same_length = 0;
  }
  if (numparts >= numofsequences) { /* assign one seq for each part */
    GtUword idx;

    for (idx = 0; idx < numofsequences; ++idx) {
      seqrangeinfo->ranges[idx].start = seqrangeinfo->ranges[idx].end = idx;
      seqrangeinfo->ranges[idx].max_length
        = gt_sspbased_seqlength(seqrangeinfo,idx);
      seqrangeinfo->ranges[idx].partnum = idx;
    }
    seqrangeinfo->parts_number = numofsequences;
  } else
  {
    GtUword seqnum, idx, effective_num_parts;
    const GtUword partwidth = seqrangeinfo->totallength/numparts;

    for (idx = 0, seqnum = 0; idx < numparts && seqnum < numofsequences; idx++)
    {
      const GtUword seqnum_next_width
        = gt_encseq_next_larger_width(seqrangeinfo,seqnum,partwidth,
                                      numofsequences);
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
      if (seqrangeinfo->ssptab == NULL)
      {
        for (idx = 0; idx < effective_num_parts; idx++)
        {
          seqrangeinfo->ranges[idx].max_length = seqrangeinfo->same_length;
        }
      } else
      {
        GtUword currentstart = 0, idx = 0, currentlength, maxlength = 0;

        gt_assert(numofsequences > 1);
        for (seqnum = 0; seqnum < numofsequences - 1; seqnum++)
        {
          gt_assert(currentstart < seqrangeinfo->ssptab[seqnum]);
          currentlength = seqrangeinfo->ssptab[seqnum] - currentstart;
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
          currentstart = seqrangeinfo->ssptab[seqnum] + 1;
        }
        currentlength = seqrangeinfo->totallength - currentstart;
        if (maxlength < currentlength)
        {
          maxlength = currentlength;
        }
        gt_assert(idx + 1 == effective_num_parts &&
                  seqnum == seqrangeinfo->ranges[idx].end);
        seqrangeinfo->ranges[idx].max_length = maxlength;
      }
    }
    seqrangeinfo->parts_number = effective_num_parts;
  }
  return seqrangeinfo;
}

void gt_seed_extend_parts_delete(GtSequencePartsInfo *seqrangeinfo)
{
  if (seqrangeinfo != NULL)
  {
    gt_free(seqrangeinfo->ranges);
    gt_free(seqrangeinfo->ssptab);
    gt_free(seqrangeinfo);
  }
}

GtUword gt_seed_extend_parts_number(const GtSequencePartsInfo *seqrangeinfo)
{
  gt_assert(seqrangeinfo != NULL);
  return seqrangeinfo->parts_number;
}

const GtSequenceRangeWithMaxLength *gt_seed_extend_parts_get(
                        const GtSequencePartsInfo *seqrangeinfo,
                        GtUword idx)
{
  gt_assert(seqrangeinfo != NULL && idx < seqrangeinfo->parts_number);
  return seqrangeinfo->ranges + idx;
}

void gt_seed_extend_parts_variance_show(const GtSequencePartsInfo *seqrangeinfo,
                                        const GtEncseq *encseq)
{
  GtUword variance_sum = 0, idx, avgpartsize;

  gt_assert(encseq != NULL);
  avgpartsize = seqrangeinfo->totallength/seqrangeinfo->parts_number;
  for (idx = 0; idx < seqrangeinfo->parts_number; idx++)
  {
    GtUword segmentstart, segmentend, width;

    segmentstart = gt_sspbased_seqstartpos(seqrangeinfo,
                                           seqrangeinfo->ranges[idx].start);
    segmentend = gt_sspbased_seqendpos(seqrangeinfo,
                                       seqrangeinfo->ranges[idx].end);
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
           (double) variance_sum/seqrangeinfo->parts_number);
}
