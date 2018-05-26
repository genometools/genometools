/*
  Copyright (c) 2007-2018 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2018 Center for Bioinformatics, University of Hamburg

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

#ifndef QUERYMATCH_H
#define QUERYMATCH_H

#include <inttypes.h>
#include <stdio.h>
#include "core/readmode.h"
#include "core/encseq.h"
#include "core/arraydef.h"
#include "querymatch-align.h"
#include "karlin_altschul_stat.h"
#include "querymatch-display.h"
#include "seq_or_encseq.h"

typedef struct GtQuerymatch GtQuerymatch;

GtQuerymatch *gt_querymatch_new(void);

void gt_querymatch_set_order(GtQuerymatch *querymatch,bool force);

void gt_querymatch_file_set(GtQuerymatch *querymatch, FILE *fp);

void gt_querymatch_db_keyvalues_set(GtQuerymatch *querymatch,
                                    GtUword db_totallength,
                                    GtUword db_numofsequences);

void gt_querymatch_init(GtQuerymatch *querymatch,
                        GtUword dblen,
                        GtUword dbseqnum,
                        GtUword dbstart_relative,
                        GtUword db_seqstart,
                        GtUword dbseqlen,
                        GtUword distance,
                        GtUword mismatches,
                        bool selfmatch,
                        GtUword queryseqnum,
                        GtUword querylen,
                        GtUword querystart,
                        GtUword query_seqstart,
                        GtUword query_seqlen,
                        const char *db_desc,
                        const char *query_desc);

void gt_querymatch_read_line(GtQuerymatch *querymatch,
                             double *evalue_ptr,
                             double *bitscore_ptr,
                             const char *line_ptr,
                             const GtSeedExtendDisplayFlag *in_display_flag,
                             bool selfmatch,
                             const GtEncseq *dbencseq,
                             const GtEncseq *queryencseq);

void gt_querymatch_delete(GtQuerymatch *querymatch);

bool gt_querymatch_complete(GtQuerymatch *querymatch,
                            const GtSeedExtendDisplayFlag *out_display_flag,
                            GtUword dblen,
                            GtUword dbseqnum,
                            GtUword dbstart_relative,
                            GtUword db_seqstart,
                            GtUword dbseqlen,
                            GtUword distance,
                            GtUword mismatches,
                            bool selfmatch,
                            GtUword queryseqnum,
                            GtUword querylen,
                            GtUword querystart,
                            const GtSeqorEncseq *dbes,
                            const GtSeqorEncseq *queryes,
                            GtUword query_seqstart,
                            GtUword query_seqlen,
                            GtUword db_seedpos_rel,
                            GtUword query_seedpos,
                            GtUword seedlen,
                            GtAffineDPreservoir *adpr,
                            bool greedyextension);

void gt_querymatch_outoptions_set(GtQuerymatch *querymatch,
                GtQuerymatchoutoptions *querymatchoutoptions);

GtUword gt_querymatch_dblen(const GtQuerymatch *querymatch);

GtUword gt_querymatch_querylen(const GtQuerymatch *querymatch);

GtUword gt_querymatch_dbstart(const GtQuerymatch *querymatch);

GtUword gt_querymatch_dbstart_relative(const GtQuerymatch *querymatch);

GtUword gt_querymatch_querystart(const GtQuerymatch *querymatch);

void gt_querymatch_db_coordinates(GtUword *db_seqnum,GtUword *db_seqstart,
                                  GtUword *db_seqlen,
                                  const GtQuerymatch *querymatch);

void gt_querymatch_query_coordinates(GtUword *query_seqnum,
                                     GtUword *query_seqstart,
                                     GtUword *query_seqlen,
                                     const GtQuerymatch *querymatch);

double gt_querymatch_error_percentage(GtUword distance,GtUword alignedlen);

double gt_querymatch_error_rate(GtUword distance,GtUword alignedlen);

GtUword gt_querymatch_distance(const GtQuerymatch *querymatch);

void gt_querymatch_gfa2_edge(const GtQuerymatch *querymatch,GtUword edgenum);

void gt_querymatch_prettyprint(double evalue,double bit_score,
                               const GtSeedExtendDisplayFlag *out_display_flag,
                               const GtQuerymatch *querymatch);

void gt_querymatch_show_failed_seed(const GtSeedExtendDisplayFlag
                                       *out_display_flag,
                                    const GtQuerymatch *querymatch);

void gt_querymatch_stat_failed_seed(GtUword userdefinedlieastlength,
                                    const GtQuerymatch *querymatch);

bool gt_querymatch_ordered(const GtQuerymatch *querymatch);

bool gt_querymatch_check_final(double *evalue_ptr,
                               double *bit_score_ptr,
                               const GtKarlinAltschulStat *karlin_altschul_stat,
                               const GtQuerymatch *querymatch,
                               GtUword userdefinedleastlength,
                               GtUword errorpercentage,
                               double evalue_threshold);

bool gt_querymatch_check_final_generic(
                               double *evalue_ptr,
                               double *bit_score_ptr,
                               const GtKarlinAltschulStat *karlin_altschul_stat,
                               GtUword query_seqlen,
                               GtUword aligned_len,
                               GtUword distance,
                               GtUword mismatches,
                               GtUword userdefinedleastlength,
                               GtUword errorpercentage,
                               double evalue_threshold);

void gt_querymatch_query_readmode_set(GtQuerymatch *querymatch,
                                      GtReadmode query_readmode);

void gt_querymatch_verify_alignment_set(GtQuerymatch *querymatch);

GtReadmode gt_querymatch_query_readmode(const GtQuerymatch *querymatch);

GtWord gt_querymatch_distance2unit_score(GtUword distance,GtUword alignedlen);

GT_DECLAREARRAYSTRUCT(GtQuerymatch);

void gt_querymatch_table_add(GtArrayGtQuerymatch *querymatch_table,
                             const GtQuerymatch *querymatch);

void gt_querymatch_table_sort(GtArrayGtQuerymatch *querymatch_table,
                              bool ascending);

GtQuerymatch *gt_querymatch_table_get(const GtArrayGtQuerymatch
                                        *querymatch_table,GtUword idx);

void gt_querymatch_recompute_alignment(GtQuerymatch *querymatch,
                                       const GtSeedExtendDisplayFlag
                                         *out_display_flag,
                                       bool match_has_cigar,
                                       bool dtrace,
                                       GtUword trace_delta,
                                       bool match_has_seed,
                                       GtAffineDPreservoir *adpr,
                                       const GtEncseq *db_encseq,
                                       const GtEncseq *query_encseq,
                                       const GtKarlinAltschulStat
                                         *karlin_altschul_stat,
                                       double evalue,
                                       double bitscore);

static inline GtUword gt_querymatch_score2distance(GtWord score,
                                                   GtUword alignedlen)
{
  gt_assert(score < 0 || alignedlen >= score);
  return (score >= 0) ? (((GtWord) alignedlen - score)/3)
                      : -(((GtWord) alignedlen + score)/3);
}

#endif
