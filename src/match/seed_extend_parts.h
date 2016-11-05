#ifndef SEED_EXTEND_PARTS_H
#define SEED_EXTEND_PARTS_H

#include "core/ma_api.h"
#include "core/encseq.h"

typedef struct
{
  GtUword partnum,    /* order number of part */
          start,      /* index number of first sequence in part */
          end,        /* index number of last sequence in part */
          max_length; /* length of longest sequence in range */
} GtSequenceRangeWithMaxLength;

typedef struct GtSequencePartsInfo GtSequencePartsInfo;

GtSequencePartsInfo *gt_seed_extend_parts_new(const GtEncseq *encseq,
                                              GtUword numofsequences,
                                              GtUword numparts);

void gt_seed_extend_parts_delete(GtSequencePartsInfo *seqrangeinfo);

GtUword gt_seed_extend_parts_number(const GtSequencePartsInfo *seqrangeinfo);

const GtSequenceRangeWithMaxLength *gt_seed_extend_parts_get(
                        const GtSequencePartsInfo *seqrangeinfo,
                        GtUword idx);

void gt_seed_extend_parts_variance_show(const GtSequencePartsInfo *seqrangeinfo,
                                        const GtEncseq *encseq);

#endif
