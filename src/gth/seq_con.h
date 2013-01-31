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

#ifndef SEQ_CON_H
#define SEQ_CON_H

#include "core/alphabet.h"
#include "core/file.h"
#include "core/range.h"
#include "core/str_api.h"
#include "core/types_api.h"

/* The sequence container class */
typedef struct GthSeqConClass GthSeqConClass;
typedef struct GthSeqCon GthSeqCon;

typedef GthSeqCon* (*GthSeqConConstructor)(const char *indexname,
                                           bool assign_rc, bool orig_seq,
                                           bool tran_seq);

void          gth_seq_con_delete(GthSeqCon*);
void          gth_seq_con_demand_orig_seq(GthSeqCon *seq_con);
GtUchar*      gth_seq_con_get_orig_seq(GthSeqCon *seq_con,
                                       unsigned long seq_num);
GtUchar*      gth_seq_con_get_tran_seq(GthSeqCon *seq_con,
                                       unsigned long seq_num);
GtUchar*      gth_seq_con_get_orig_seq_rc(GthSeqCon *seq_con,
                                          unsigned long seq_num);
GtUchar*      gth_seq_con_get_tran_seq_rc(GthSeqCon *seq_con,
                                          unsigned long seq_num);
void          gth_seq_con_get_description(GthSeqCon *seq_con,
                                          unsigned long seq_num, GtStr *desc);
void          gth_seq_con_echo_description(GthSeqCon *seq_con,
                                           unsigned long seq_num,
                                           GtFile *outfp);
unsigned long gth_seq_con_num_of_seqs(GthSeqCon *seq_con);
unsigned long gth_seq_con_total_length(GthSeqCon *seq_con);
GtRange       gth_seq_con_get_range(GthSeqCon *seq_con, unsigned long seq_num);

GtRange       gth_seq_con_get_relative_range(GthSeqCon *seq_con,
                                             unsigned long seq_num);
unsigned long gth_seq_con_get_length(GthSeqCon *seq_con, unsigned long seq_num);
GtAlphabet*   gth_seq_con_get_alphabet(GthSeqCon *seq_con);

#endif
