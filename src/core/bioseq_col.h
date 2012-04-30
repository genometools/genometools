/*
  Copyright (c) 2010, 2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#ifndef BIOSEQ_COL_H
#define BIOSEQ_COL_H

#include "core/str_array_api.h"

typedef struct GtBioseqCol GtBioseqCol;

GtBioseqCol*  gt_bioseq_col_new(GtStrArray *sequence_files, GtError *err);
void          gt_bioseq_col_delete(GtBioseqCol*);
int           gt_bioseq_col_grep_desc(GtBioseqCol*, const char **rawseq,
                                      unsigned long *length, GtStr *seqid,
                                      GtError*);
int           gt_bioseq_col_grep_desc_md5(GtBioseqCol*, const char **md5,
                                          GtStr *seqid, GtError*);
int           gt_bioseq_col_md5_to_seq(GtBioseqCol*, const char **seq,
                                       unsigned long *length, GtStr *md5_seqid,
                                       GtError *err);
int           gt_bioseq_col_md5_to_description(GtBioseqCol*, GtStr *desc,
                                               GtStr *md5_seqid, GtError *err);
unsigned long gt_bioseq_col_num_of_files(const GtBioseqCol*);
unsigned long gt_bioseq_col_num_of_seqs(const GtBioseqCol*,
                                        unsigned long filenum);
const char*   gt_bioseq_col_get_md5_fingerprint(const GtBioseqCol*,
                                                unsigned long filenum,
                                                unsigned long seqnum);
const char*   gt_bioseq_col_get_sequence(const GtBioseqCol*,
                                         unsigned long filenum,
                                         unsigned long seqnum);
unsigned long gt_bioseq_col_get_sequence_length(const GtBioseqCol*,
                                                unsigned long filenum,
                                                unsigned long seqnum);

#endif
