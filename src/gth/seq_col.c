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

#include "core/class_alloc.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "gth/seq_col_rep.h"

struct GthSeqColClass {
  size_t size;
  GthSeqColGetOrigSeqFunc get_orig_seq;
  GthSeqColGetTranSeqFunc get_tran_seq;
  GthSeqColGetOrigSeqRCFunc get_orig_seq_rc;
  GthSeqColGetTranSeqRCFunc get_tran_seq_rc;
  GthSeqColGetDescriptionFunc get_description;
  GthSeqColEchoDescriptionFunc echo_description;
  GthSeqColNumOfSeqsFunc num_of_seqs;
  GthSeqColTotalLengthFunc total_length;
  GthSeqColGetRangeFunc get_range;
  GthSeqColGetAlphabetFunc get_alphabet;
  GthSeqColFreeFunc free;
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
                                            GthSeqColFreeFunc free)
{
  GthSeqColClass *c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->get_orig_seq = get_orig_seq;
  c_class->get_tran_seq = get_tran_seq;
  c_class->get_orig_seq_rc = get_orig_seq_rc;
  c_class->get_tran_seq_rc = get_tran_seq_rc;
  c_class->get_description = get_description;
  c_class->echo_description = echo_description;
  c_class->num_of_seqs = num_of_seqs;
  c_class->total_length = total_length;
  c_class->get_range = get_range;
  c_class->get_alphabet = get_alphabet;
  c_class->free = free;
  return c_class;
}

GthSeqCol* gth_seq_col_create(const GthSeqColClass *scc)
{
  GthSeqCol *sc;
  gt_assert(scc && scc->size);
  sc = gt_calloc(1, scc->size);
  sc->c_class = scc;
  sc->pvt = NULL; /* XXX */
  return sc;
}

void* gth_seq_col_cast(GT_UNUSED const GthSeqColClass *scc, GthSeqCol *sc)
{
  gt_assert(scc && sc && sc->c_class == scc);
  return sc;
}

void gth_seq_col_delete(GthSeqCol *sc)
{
  if (!sc) return;
  if (sc->c_class->free)
    sc->c_class->free(sc);
  gt_free(sc->pvt);
  gt_free(sc);
}

GtUchar* gth_seq_col_get_orig_seq(GthSeqCol *seq_col, unsigned long seq_num)
{
  gt_assert(seq_col && seq_col->c_class && seq_col->c_class->get_orig_seq);
  return seq_col->c_class->get_orig_seq(seq_col, seq_num);
}

GtUchar* gth_seq_col_get_tran_seq(GthSeqCol *seq_col, unsigned long seq_num)
{
  gt_assert(seq_col && seq_col->c_class && seq_col->c_class->get_tran_seq);
  return seq_col->c_class->get_tran_seq(seq_col, seq_num);
}

GtUchar* gth_seq_col_get_orig_seq_rc(GthSeqCol *seq_col, unsigned long seq_num)
{
  gt_assert(seq_col && seq_col->c_class && seq_col->c_class->get_orig_seq_rc);
  return seq_col->c_class->get_orig_seq_rc(seq_col, seq_num);
}

GtUchar* gth_seq_col_get_tran_seq_rc(GthSeqCol *seq_col, unsigned long seq_num)
{
  gt_assert(seq_col && seq_col->c_class && seq_col->c_class->get_tran_seq_rc);
  return seq_col->c_class->get_tran_seq_rc(seq_col, seq_num);
}

void gth_seq_col_get_description(GthSeqCol *seq_col, unsigned long seq_num,
                                 GtStr *desc)
{
  gt_assert(seq_col && seq_col->c_class && seq_col->c_class->get_description);
  seq_col->c_class->get_description(seq_col, seq_num, desc);
}

void gth_seq_col_echo_description(GthSeqCol *seq_col, unsigned long seq_num,
                                  GtFile *outfp)
{
  gt_assert(seq_col && seq_col->c_class && seq_col->c_class->echo_description);
  seq_col->c_class->echo_description(seq_col, seq_num, outfp);
}

unsigned long gth_seq_col_num_of_seqs(GthSeqCol *seq_col)
{
  gt_assert(seq_col && seq_col->c_class && seq_col->c_class->num_of_seqs);
  return seq_col->c_class->num_of_seqs(seq_col);
}

unsigned long gth_seq_col_total_length(GthSeqCol *seq_col)
{
  gt_assert(seq_col && seq_col->c_class && seq_col->c_class->total_length);
  return seq_col->c_class->total_length(seq_col);
}

GtRange gth_seq_col_get_range(GthSeqCol *seq_col, unsigned long seq_num)
{
  gt_assert(seq_col && seq_col->c_class && seq_col->c_class->get_range);
  return seq_col->c_class->get_range(seq_col, seq_num);
}

GtRange gth_seq_col_get_relative_range(GthSeqCol *seq_col,
                                       unsigned long seq_num)
{
  GtRange relative_range, range;
  gt_assert(seq_col);
  gt_assert(seq_num < gth_seq_col_num_of_seqs(seq_col));
  range = gth_seq_col_get_range(seq_col, seq_num);
  relative_range.start = 0;
  relative_range.end = range.end - range.start;
  return relative_range;
}

unsigned long gth_seq_col_get_length(GthSeqCol *seq_col, unsigned long seq_num)
{
  GtRange range;
  gt_assert(seq_col);
  range = gth_seq_col_get_range(seq_col, seq_num);
  return gt_range_length(&range);
}

GtAlphabet* gth_seq_col_get_alphabet(GthSeqCol *seq_col)
{
  gt_assert(seq_col && seq_col->c_class && seq_col->c_class->get_alphabet);
  return seq_col->c_class->get_alphabet(seq_col);
}
