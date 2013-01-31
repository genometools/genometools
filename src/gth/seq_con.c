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

#include "core/class_alloc.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "gth/seq_con_rep.h"

struct GthSeqConClass {
  size_t size;
  GthSeqConDemandOrigSeqFunc demand_orig_seq;
  GthSeqConGetOrigSeqFunc get_orig_seq;
  GthSeqConGetTranSeqFunc get_tran_seq;
  GthSeqConGetOrigSeqRCFunc get_orig_seq_rc;
  GthSeqConGetTranSeqRCFunc get_tran_seq_rc;
  GthSeqConGetDescriptionFunc get_description;
  GthSeqConEchoDescriptionFunc echo_description;
  GthSeqConNumOfSeqsFunc num_of_seqs;
  GthSeqConTotalLengthFunc total_length;
  GthSeqConGetRangeFunc get_range;
  GthSeqConGetAlphabetFunc get_alphabet;
  GthSeqConFreeFunc free;
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
                                            GthSeqConFreeFunc free)
{
  GthSeqConClass *c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->demand_orig_seq = demand_orig_seq;
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

GthSeqCon* gth_seq_con_create(const GthSeqConClass *scc)
{
  GthSeqCon *sc;
  gt_assert(scc && scc->size);
  sc = gt_calloc(1, scc->size);
  sc->c_class = scc;
  sc->pvt = NULL; /* XXX */
  return sc;
}

void* gth_seq_con_cast(GT_UNUSED const GthSeqConClass *scc, GthSeqCon *sc)
{
  gt_assert(scc && sc && sc->c_class == scc);
  return sc;
}

void gth_seq_con_delete(GthSeqCon *sc)
{
  if (!sc) return;
  if (sc->c_class->free)
    sc->c_class->free(sc);
  gt_free(sc->pvt);
  gt_free(sc);
}

void gth_seq_con_demand_orig_seq(GthSeqCon *seq_con)
{
  gt_assert(seq_con && seq_con->c_class && seq_con->c_class->demand_orig_seq);
  seq_con->c_class->demand_orig_seq(seq_con);
}
GtUchar* gth_seq_con_get_orig_seq(GthSeqCon *seq_con, unsigned long seq_num)
{
  gt_assert(seq_con && seq_con->c_class && seq_con->c_class->get_orig_seq);
  return seq_con->c_class->get_orig_seq(seq_con, seq_num);
}

GtUchar* gth_seq_con_get_tran_seq(GthSeqCon *seq_con, unsigned long seq_num)
{
  gt_assert(seq_con && seq_con->c_class && seq_con->c_class->get_tran_seq);
  return seq_con->c_class->get_tran_seq(seq_con, seq_num);
}

GtUchar* gth_seq_con_get_orig_seq_rc(GthSeqCon *seq_con, unsigned long seq_num)
{
  gt_assert(seq_con && seq_con->c_class && seq_con->c_class->get_orig_seq_rc);
  return seq_con->c_class->get_orig_seq_rc(seq_con, seq_num);
}

GtUchar* gth_seq_con_get_tran_seq_rc(GthSeqCon *seq_con, unsigned long seq_num)
{
  gt_assert(seq_con && seq_con->c_class && seq_con->c_class->get_tran_seq_rc);
  return seq_con->c_class->get_tran_seq_rc(seq_con, seq_num);
}

void gth_seq_con_get_description(GthSeqCon *seq_con, unsigned long seq_num,
                                 GtStr *desc)
{
  gt_assert(seq_con && seq_con->c_class && seq_con->c_class->get_description);
  seq_con->c_class->get_description(seq_con, seq_num, desc);
}

void gth_seq_con_echo_description(GthSeqCon *seq_con, unsigned long seq_num,
                                  GtFile *outfp)
{
  gt_assert(seq_con && seq_con->c_class && seq_con->c_class->echo_description);
  seq_con->c_class->echo_description(seq_con, seq_num, outfp);
}

unsigned long gth_seq_con_num_of_seqs(GthSeqCon *seq_con)
{
  gt_assert(seq_con && seq_con->c_class && seq_con->c_class->num_of_seqs);
  return seq_con->c_class->num_of_seqs(seq_con);
}

unsigned long gth_seq_con_total_length(GthSeqCon *seq_con)
{
  gt_assert(seq_con && seq_con->c_class && seq_con->c_class->total_length);
  return seq_con->c_class->total_length(seq_con);
}

GtRange gth_seq_con_get_range(GthSeqCon *seq_con, unsigned long seq_num)
{
  gt_assert(seq_con && seq_con->c_class && seq_con->c_class->get_range);
  return seq_con->c_class->get_range(seq_con, seq_num);
}

GtRange gth_seq_con_get_relative_range(GthSeqCon *seq_con,
                                       unsigned long seq_num)
{
  GtRange relative_range, range;
  gt_assert(seq_con);
  gt_assert(seq_num < gth_seq_con_num_of_seqs(seq_con));
  range = gth_seq_con_get_range(seq_con, seq_num);
  relative_range.start = 0;
  relative_range.end = range.end - range.start;
  return relative_range;
}

unsigned long gth_seq_con_get_length(GthSeqCon *seq_con, unsigned long seq_num)
{
  GtRange range;
  gt_assert(seq_con);
  range = gth_seq_con_get_range(seq_con, seq_num);
  return gt_range_length(&range);
}

GtAlphabet* gth_seq_con_get_alphabet(GthSeqCon *seq_con)
{
  gt_assert(seq_con && seq_con->c_class && seq_con->c_class->get_alphabet);
  return seq_con->c_class->get_alphabet(seq_con);
}
