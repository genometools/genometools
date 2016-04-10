/*
  Copyright (c) 2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#ifndef QUERYMATCH_ALIGN_H
#define QUERYMATCH_ALIGN_H

#include "core/types_api.h"
#include "match/ft-front-prune.h"
#include "match/seq_or_encseq.h"
#include "extended/alignment.h"

typedef struct GtQuerymatchoutoptions GtQuerymatchoutoptions;

GtQuerymatchoutoptions *gt_querymatchoutoptions_new(bool generatealignment,
                                                    bool showeoplist,
                                                    GtUword alignmentwidth);

void gt_querymatchoutoptions_extend(
                  GtQuerymatchoutoptions *querymatchoutoptions,
                  GtUword errorpercentage,
                  GtUword maxalignedlendifference,
                  GtUword history,
                  GtUword perc_mat_history,
                  GtExtendCharAccess extend_char_access,
                  bool weakends,
                  GtUword sensitivity,
                  double matchscore_bias,
                  bool always_polished_ends,
                  unsigned int display_flag);

void gt_querymatchoutoptions_for_align_only(
                  GtQuerymatchoutoptions *querymatchoutoptions,
                  GtUword errorpercentage,
                  double matchscore_bias,
                  GtUword history_size,
                  bool always_polished_ends,
                  unsigned int display_flag);

void gt_querymatchoutoptions_delete(
        GtQuerymatchoutoptions *querymatchoutoptions);

bool gt_querymatchoutoptions_alignment_prepare(
                                     GtQuerymatchoutoptions
                                       *querymatchoutoptions,
                                     const GtEncseq *encseq,
                                     const GtSeqorEncseq *query,
                                     GtReadmode query_readmode,
                                     GtUword query_seqstartpos,
                                     GtUword query_totallength,
                                     GtUword dbstart,
                                     GtUword dblen,
                                     GtUword querystart,
                                     GtUword querystart_fwdstrand,
                                     GtUword querylen,
                                     GtUword edist,
                                     GtUword seedpos1,
                                     GtUword seedpos2,
                                     GtUword seedlen,
                                     GT_UNUSED bool greedyextension);

void gt_querymatchoutoptions_alignment_show(const GtQuerymatchoutoptions
                                              *querymatchoutoptions,
                                            GtUword distance,
                                            bool verify_alignment,
                                            FILE *fp);

typedef struct
{
  GtUword uoffset, voffset, ulen, vlen, sumdist;
} GtSeqpaircoordinates;

const GtSeqpaircoordinates *gt_querymatchoutoptions_correction_get(
              const GtQuerymatchoutoptions *querymatchoutoptions);

const GtAlignment *gt_querymatchoutoptions_alignment_get(
              const GtQuerymatchoutoptions *querymatchoutoptions);

#endif
