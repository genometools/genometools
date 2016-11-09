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

#ifndef SEED_EXTEND_PARTS_H
#define SEED_EXTEND_PARTS_H

#include "core/ma_api.h"
#include "core/encseq.h"

typedef struct GtSequencePartsInfo GtSequencePartsInfo;

GtSequencePartsInfo *gt_sequence_parts_info_new(const GtEncseq *encseq,
                                                GtUword numofsequences,
                                                GtUword numparts);

void gt_sequence_parts_info_delete(GtSequencePartsInfo *spi);

GtUword gt_sequence_parts_info_number(const GtSequencePartsInfo *spi);

GtUword gt_sequence_parts_info_start_get(const GtSequencePartsInfo *spi,
                                         GtUword idx);

GtUword gt_sequence_parts_info_end_get(const GtSequencePartsInfo *spi,
                                       GtUword idx);

GtUword gt_sequence_parts_info_numofsequences_get(
                                       const GtSequencePartsInfo *spi,
                                       GtUword idx);

GtUword gt_sequence_parts_info_max_length_get(const GtSequencePartsInfo *spi,
                                              GtUword idx);

bool gt_sequence_parts_info_overlap(const GtSequencePartsInfo *spia,
                                    GtUword aidx,
                                    const GtSequencePartsInfo *spib,
                                    GtUword bidx);

bool gt_sequence_parts_info_equal(const GtSequencePartsInfo *spia,
                                  GtUword aidx,
                                  const GtSequencePartsInfo *spib,
                                  GtUword bidx);

GtUword gt_sequence_parts_info_seqstartpos(const GtSequencePartsInfo *spi,
                                           GtUword seqnum);

GtUword gt_sequence_parts_info_seqendpos(const GtSequencePartsInfo *spi,
                                         GtUword seqnum);

GtUword gt_sequence_parts_info_partlength(const GtSequencePartsInfo *spi,
                                          GtUword fromseq,
                                          GtUword toseq);

void gt_sequence_parts_info_variance_show(const GtSequencePartsInfo *spi);

GtUchar *gt_sequence_parts_info_seq_extract(const GtSequencePartsInfo *spi,
                                            GtUword idx);

const GtEncseq *gt_seuence_part_info_encseq_get(const GtSequencePartsInfo *spi);

#endif
