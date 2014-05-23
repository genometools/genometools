/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef SEQ_COL_H
#define SEQ_COL_H

#include "core/error_api.h"
#include "core/str_api.h"

typedef struct GtSeqCol GtSeqCol;

void        gt_seq_col_delete(GtSeqCol*);
void        gt_seq_col_enable_match_desc_start(GtSeqCol*);
int         gt_seq_col_grep_desc(GtSeqCol*, char **seq,
                                 GtUword start, GtUword end,
                                 GtStr *seqid, GtError*);
int         gt_seq_col_grep_desc_md5(GtSeqCol*, const char **md5,
                                     GtStr *seqid, GtError*);
int         gt_seq_col_grep_desc_sequence_length(GtSeqCol *sc,
                                                 GtUword *length,
                                                 GtStr *seqid,
                                                 GtError *err);
int         gt_seq_col_md5_to_seq(GtSeqCol*, char **seq,
                                  GtUword start, GtUword end,
                                  GtStr *md5_seqid, GtError *err);
int         gt_seq_col_md5_to_description(GtSeqCol*, GtStr *desc,
                                          GtStr *md5_seqid, GtError *err);
int         gt_seq_col_md5_to_sequence_length(GtSeqCol*, GtUword *len,
                                              GtStr *md5_seqid, GtError *err);
GtUword     gt_seq_col_num_of_files(const GtSeqCol*);
GtUword     gt_seq_col_num_of_seqs(const GtSeqCol*, GtUword filenum);
const char* gt_seq_col_get_md5_fingerprint(const GtSeqCol*,
                                           GtUword filenum,
                                           GtUword seqnum);
char*       gt_seq_col_get_sequence(const GtSeqCol*,
                                    GtUword filenum,
                                    GtUword seqnum,
                                    GtUword start,
                                    GtUword end);
char*       gt_seq_col_get_description(const GtSeqCol*,
                                       GtUword filenum,
                                       GtUword seqnum);
GtUword     gt_seq_col_get_sequence_length(const GtSeqCol*,
                                           GtUword filenum,
                                           GtUword seqnum);

#endif
