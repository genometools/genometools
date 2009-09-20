/*
  Copyright (c) 2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#ifndef SEQ_COL_REP_H
#define SEQ_COL_REP_H

#include "gth/seq_col.h"

typedef GtUchar* (*GthSeqColGetOrigSeqFunc)(GthSeqCol*, unsigned long seq_num);
typedef GtUchar* (*GthSeqColGetTranSeqFunc)(GthSeqCol*, unsigned long seq_num);
typedef GtUchar* (*GthSeqColGetOrigSeqRCFunc)(GthSeqCol*,
                                              unsigned long seq_num);
typedef GtUchar* (*GthSeqColGetTranSeqRCFunc)(GthSeqCol*,
                                              unsigned long seq_num);
typedef void     (*GthSeqColGetDescriptionFunc)(GthSeqCol*,
                                                unsigned long seq_num,
                                                GtStr *desc);
typedef void     (*GthSeqColEchoDescriptionFunc)(GthSeqCol*,
                                                 unsigned long seq_num,
                                                 GtFile *outfp);
typedef unsigned long (*GthSeqColNumOfSeqsFunc)(GthSeqCol*);
typedef unsigned long (*GthSeqColTotalLengthFunc)(GthSeqCol*);
typedef GtRange       (*GthSeqColGetRangeFunc)(GthSeqCol*,
                                               unsigned long seq_num);
typedef GtAlphabet*   (*GthSeqColGetAlphabetFunc)(GthSeqCol*);
typedef void          (*GthSeqColFreeFunc)(GthSeqCol*);

typedef struct GthSeqColMembers GthSeqColMembers;

struct GthSeqCol {
  const GthSeqColClass *c_class;
  GthSeqColMembers *pvt;
};

const GthSeqColClass* gth_seq_col_class_new(size_t size,
                                            GthSeqColGetOrigSeqFunc
                                            get_orig_seq,
                                            GthSeqColGetTranSeqFunc
                                            get_tran_seq,
                                            GthSeqColGetOrigSeqRCFunc
                                            get_orig_seq_rc,
                                            GthSeqColGetTranSeqRCFunc
                                            get_tran_seq_rc,
                                            GthSeqColGetDescriptionFunc
                                            get_description,
                                            GthSeqColEchoDescriptionFunc
                                            echo_description,
                                            GthSeqColNumOfSeqsFunc num_of_seqs,
                                            GthSeqColTotalLengthFunc
                                            total_length,
                                            GthSeqColGetRangeFunc get_range,
                                            GthSeqColGetAlphabetFunc
                                            get_alphabet,
                                            GthSeqColFreeFunc free);
GthSeqCol*           gth_seq_col_create(const GthSeqColClass*);
void*                gth_seq_col_cast(const GthSeqColClass*, GthSeqCol*);

#endif
