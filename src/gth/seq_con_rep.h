/*
  Copyright (c) 2009-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#ifndef SEQ_CON_REP_H
#define SEQ_CON_REP_H

#include "gth/seq_con.h"

typedef void     (*GthSeqConDemandOrigSeqFunc)(GthSeqCon*);
typedef GtUchar* (*GthSeqConGetOrigSeqFunc)(GthSeqCon*, unsigned long seq_num);
typedef GtUchar* (*GthSeqConGetTranSeqFunc)(GthSeqCon*, unsigned long seq_num);
typedef GtUchar* (*GthSeqConGetOrigSeqRCFunc)(GthSeqCon*,
                                              unsigned long seq_num);
typedef GtUchar* (*GthSeqConGetTranSeqRCFunc)(GthSeqCon*,
                                              unsigned long seq_num);
typedef void     (*GthSeqConGetDescriptionFunc)(GthSeqCon*,
                                                unsigned long seq_num,
                                                GtStr *desc);
typedef void     (*GthSeqConEchoDescriptionFunc)(GthSeqCon*,
                                                 unsigned long seq_num,
                                                 GtFile *outfp);
typedef unsigned long (*GthSeqConNumOfSeqsFunc)(GthSeqCon*);
typedef unsigned long (*GthSeqConTotalLengthFunc)(GthSeqCon*);
typedef GtRange       (*GthSeqConGetRangeFunc)(GthSeqCon*,
                                               unsigned long seq_num);
typedef GtAlphabet*   (*GthSeqConGetAlphabetFunc)(GthSeqCon*);
typedef void          (*GthSeqConFreeFunc)(GthSeqCon*);

typedef struct GthSeqConMembers GthSeqConMembers;

struct GthSeqCon {
  const GthSeqConClass *c_class;
  GthSeqConMembers *pvt;
};

const GthSeqConClass* gth_seq_con_class_new(size_t size,
                                            GthSeqConDemandOrigSeqFunc
                                            demand_orig_seq,
                                            GthSeqConGetOrigSeqFunc
                                            get_orig_seq,
                                            GthSeqConGetTranSeqFunc
                                            get_tran_seq,
                                            GthSeqConGetOrigSeqRCFunc
                                            get_orig_seq_rc,
                                            GthSeqConGetTranSeqRCFunc
                                            get_tran_seq_rc,
                                            GthSeqConGetDescriptionFunc
                                            get_description,
                                            GthSeqConEchoDescriptionFunc
                                            echo_description,
                                            GthSeqConNumOfSeqsFunc num_of_seqs,
                                            GthSeqConTotalLengthFunc
                                            total_length,
                                            GthSeqConGetRangeFunc get_range,
                                            GthSeqConGetAlphabetFunc
                                            get_alphabet,
                                            GthSeqConFreeFunc free);
GthSeqCon*           gth_seq_con_create(const GthSeqConClass*);
void*                gth_seq_con_cast(const GthSeqConClass*, GthSeqCon*);

#endif
