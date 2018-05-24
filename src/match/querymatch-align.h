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
#include "core/error_api.h"
#include "match/ft-front-prune.h"
#include "match/seq_or_encseq.h"
#include "match/querymatch-display.h"
#include "match/affine_dband.h"

typedef struct GtQuerymatchoutoptions GtQuerymatchoutoptions;

GtQuerymatchoutoptions *gt_querymatchoutoptions_new(const
                                                     GtSeedExtendDisplayFlag
                                                      *out_display_flag,
                                                    const char *indexname,
                                                    GtError *err);

void gt_querymatchoutoptions_reset(GtQuerymatchoutoptions
                                     *querymatchoutoptions);

void gt_querymatchoutoptions_extend(
                  GtQuerymatchoutoptions *querymatchoutoptions,
                  GtUword errorpercentage,
                  double evalue_threshold,
                  GtUword maxalignedlendifference,
                  GtUword history,
                  GtUword perc_mat_history,
                  GtExtendCharAccess a_extend_char_access,
                  GtExtendCharAccess b_extend_char_access,
                  bool cam_generic,
                  bool weakends,
                  GtUword sensitivity,
                  double matchscore_bias,
                  bool always_polished_ends,
                  const GtSeedExtendDisplayFlag *out_display_flag);

void gt_querymatchoutoptions_for_align_only(
                  GtQuerymatchoutoptions *querymatchoutoptions,
                  GtUword errorpercentage,
                  double matchscore_bias,
                  GtUword history_size,
                  bool always_polished_ends,
                  GtExtendCharAccess a_extend_char_access,
                  GtExtendCharAccess b_extend_char_access,
                  const GtSeedExtendDisplayFlag *out_display_flag);

void gt_querymatchoutoptions_delete(
        GtQuerymatchoutoptions *querymatchoutoptions);

void gt_frontprune2eoplist(GtQuerymatchoutoptions *querymatchoutoptions,
                           const GtSeqorEncseq *dbes,
                           GtUword dbstart,
                           GtUword dblen,
                           const GtSeqorEncseq *queryes,
                           GtReadmode query_readmode,
                           GtUword query_seqstart,
                           GtUword query_seqlen,
                           GtUword querystart,
                           GtUword querylen,
                           bool verify_alignment);

void gt_querymatchoutoptions_cigar_show(const GtQuerymatchoutoptions
                                              *querymatchoutoptions,
                                        bool distinguish_mismatch_match,
                                        FILE *fp);

void gt_querymatchoutoptions_exact_cigar_show(bool distinguish_mismatch_match,
                                              GtUword matchlength,
                                              FILE *fp);

void gt_querymatchoutoptions_trace_show(const GtQuerymatchoutoptions
                                              *querymatchoutoptions,
                                        bool dtrace,
                                        FILE *fp);

void gt_querymatchoutoptions_alignment_show(const GtQuerymatchoutoptions
                                              *querymatchoutoptions,
                                            GtUword subject_seqlength,
                                            GtUword query_reference,
                                            GtUword one_off,
                                            GtUword distance,
                                            bool distinguish_mismatch_match,
                                            bool verify_alignment,
                                            bool subject_first,
                                            bool alignment_show_forward,
                                            bool show_complement_characters,
                                            FILE *fp);

typedef struct
{
  GtUword uoffset, voffset, ulen, vlen, sumdist, sum_max_mismatches;
} GtSeqpaircoordinates;

const GtSeqpaircoordinates *gt_querymatchoutoptions_correction_get(
              const GtQuerymatchoutoptions *querymatchoutoptions);

void gt_querymatch_column_header_output(const GtSeedExtendDisplayFlag
                                         *out_display_flag,FILE *stream);

void gt_querymatchoutoptions_extract_seq(GtQuerymatchoutoptions
                                           *querymatchoutoptions,
                                         const GtSeqorEncseq *dbes,
                                         GtUword dbstart,
                                         GtUword dblen,
                                         GtReadmode query_readmode,
                                         const GtSeqorEncseq *queryes,
                                         GtUword abs_querystart_fwdstrand,
                                         GtUword querylen);

GtUword gt_querymatchoutoptions_affine_alignment(GtQuerymatchoutoptions
                                                   *querymatchoutoptions,
                                                 GtUword dblen,
                                                 GtUword querylen,
                                                 GtWord score,
                                                 GtAffineDPreservoir *adpr);

void gt_querymatchoutoptions_set_sequences(GtQuerymatchoutoptions
                                             *querymatchoutoptions,
                                           GtUword dbstart_relative,
                                           GtUword dblen,
                                           GtUword querystart,
                                           GtUword querylen,
                                           bool withcorrection);

void gt_querymatchoutoptions_seededmatch2eoplist(
                                GtQuerymatchoutoptions *querymatchoutoptions,
                                const GtSeqorEncseq *dbes,
                                GtUword dbstart_relative,
                                GtUword db_seqstart,
                                GtUword dblen,
                                GtReadmode query_readmode,
                                const GtSeqorEncseq *queryes,
                                GtUword query_seqstart,
                                GtUword query_seqlen,
                                GtUword querystart,
                                GtUword querylen,
                                GtUword db_seedpos,
                                GtUword query_seedpos,
                                GtUword seedlen,
                                bool verify_alignment,
                                bool greedyextension);

GtEoplist *gt_querymatchoutoptions_eoplist(const GtQuerymatchoutoptions
                                             *querymatchoutoptions);

#endif
