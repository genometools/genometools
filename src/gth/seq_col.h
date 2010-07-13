/*
  Copyright (c) 2009-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include "core/alphabet.h"
#include "core/file.h"
#include "core/range.h"
#include "core/str_api.h"
#include "core/types_api.h"

typedef struct GthSeqColClass GthSeqColClass;
typedef struct GthSeqCol GthSeqCol;

typedef GthSeqCol* (*GthSeqColConstructor)(const char *indexname,
                                           bool assign_rc, bool translate);

void          gth_seq_col_delete(GthSeqCol*);
GtUchar*      gth_seq_col_get_orig_seq(GthSeqCol *seq_col,
                                       unsigned long seq_num);
GtUchar*      gth_seq_col_get_tran_seq(GthSeqCol *seq_col,
                                       unsigned long seq_num);
GtUchar*      gth_seq_col_get_orig_seq_rc(GthSeqCol *seq_col,
                                          unsigned long seq_num);
GtUchar*      gth_seq_col_get_tran_seq_rc(GthSeqCol *seq_col,
                                          unsigned long seq_num);
void          gth_seq_col_get_description(GthSeqCol *seq_col,
                                          unsigned long seq_num, GtStr *desc);
void          gth_seq_col_echo_description(GthSeqCol *seq_col,
                                           unsigned long seq_num,
                                           GtFile *outfp);
unsigned long gth_seq_col_num_of_seqs(GthSeqCol *seq_col);
unsigned long gth_seq_col_total_length(GthSeqCol *seq_col);
GtRange       gth_seq_col_get_range(GthSeqCol *seq_col, unsigned long seq_num);

GtRange       gth_seq_col_get_relative_range(GthSeqCol *seq_col,
                                             unsigned long seq_num);
unsigned long gth_seq_col_get_length(GthSeqCol *seq_col, unsigned long seq_num);
GtAlphabet*   gth_seq_col_get_alphabet(GthSeqCol *seq_col);

#endif
