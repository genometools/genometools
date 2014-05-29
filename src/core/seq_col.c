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

#include "core/class_alloc.h"
#include "core/seq_col.h"
#include "core/seq_col_rep.h"
#include "core/ma.h"
#include "core/undef_api.h"
#include "core/unused_api.h"

const GtSeqColClass* gt_seq_col_class_new(size_t size,
                                          GtSeqColFreeFunc free,
                                          GtSeqColEnableMatchDescStartFunc
                                                        enable_match_desc_start,
                                          GtSeqColGrepDescFunc grep_desc,
                                          GtSeqColGrepDescMD5Func grep_desc_md5,
                                          GtSeqColGrepDescSeqlenFunc
                                                               grep_desc_seqlen,
                                          GtSeqColMD5ToSeqFunc md5_to_seq,
                                          GtSeqColMD5ToDescFunc md5_to_desc,
                                          GtSeqColMD5ToSeqlenFunc md5_to_seqlen,
                                          GtSeqColNumFilesFunc num_files,
                                          GtSeqColNumSeqsFunc num_seqs,
                                          GtSeqColGetMD5Func get_md5,
                                          GtSeqColGetSeqFunc get_seq,
                                          GtSeqColGetDescFunc get_desc,
                                          GtSeqColGetSeqlenFunc get_seqlen)
{
  GtSeqColClass *c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->free = free;
  c_class->enable_match_desc_start = enable_match_desc_start;
  c_class->grep_desc = grep_desc;
  c_class->grep_desc_md5 = grep_desc_md5;
  c_class->grep_desc_seqlen = grep_desc_seqlen;
  c_class->md5_to_desc = md5_to_desc;
  c_class->md5_to_seq = md5_to_seq;
  c_class->md5_to_seqlen = md5_to_seqlen;
  c_class->num_files = num_files;
  c_class->num_seqs = num_seqs;
  c_class->get_md5 = get_md5;
  c_class->get_seq = get_seq;
  c_class->get_desc = get_desc;
  c_class->get_seqlen = get_seqlen;
  return c_class;
}

GtSeqCol* gt_seq_col_create(const GtSeqColClass *cc)
{
  GtSeqCol *c;
  gt_assert(cc && cc->size);
  c = gt_calloc(1, cc->size);
  c->c_class = cc;
  return c;
}

void gt_seq_col_delete(GtSeqCol *sc)
{
  if (!sc) return;
  gt_assert(sc->c_class);
  if (sc->c_class->free)
    sc->c_class->free(sc);
  gt_free(sc);
}

void* gt_seq_col_cast(GT_UNUSED const GtSeqColClass *scc, const GtSeqCol *sc)
{
  gt_assert(scc && sc && sc->c_class == scc);
  return (void*) sc;
}

void gt_seq_col_enable_match_desc_start(GtSeqCol *sc)
{
  gt_assert(sc);
  if (sc->c_class->enable_match_desc_start)
    sc->c_class->enable_match_desc_start(sc);
}

int gt_seq_col_grep_desc(GtSeqCol *sc, char **seq, GtUword start,
                         GtUword end, GtStr *seqid, GtError *err)
{
  gt_assert(sc && seq && seqid);
  if (sc->c_class->grep_desc)
    return sc->c_class->grep_desc(sc, seq, start, end, seqid, err);
  return 0;
}

int gt_seq_col_grep_desc_md5(GtSeqCol *sc, const char **md5, GtStr *seqid,
                             GtError *err)
{
  gt_assert(sc && md5 && seqid);
  if (sc->c_class->grep_desc_md5)
    return sc->c_class->grep_desc_md5(sc, md5, seqid, err);
  return 0;
}

int gt_seq_col_grep_desc_sequence_length(GtSeqCol *sc, GtUword *length,
                                         GtStr *seqid, GtError *err)
{
  gt_assert(sc && length && seqid);
  if (sc->c_class->grep_desc_seqlen)
    return sc->c_class->grep_desc_seqlen(sc, length, seqid, err);
  return 0;
}

int gt_seq_col_md5_to_seq(GtSeqCol *sc, char **seq,
                          GtUword start, GtUword end,
                          GtStr *md5_seqid, GtError *err)
{
  gt_assert(sc && seq && md5_seqid);
  if (sc->c_class->md5_to_seq)
    return sc->c_class->md5_to_seq(sc, seq, start, end, md5_seqid, err);
  return 0;
}

int gt_seq_col_md5_to_description(GtSeqCol *sc, GtStr *desc, GtStr *md5_seqid,
                                  GtError *err)
{
  gt_assert(sc && desc && md5_seqid);
  if (sc->c_class->md5_to_desc)
    return sc->c_class->md5_to_desc(sc, desc, md5_seqid, err);
  return 0;
}

int gt_seq_col_md5_to_sequence_length(GtSeqCol *sc, GtUword *len,
                                      GtStr *md5_seqid, GtError *err)
{
  gt_assert(sc && len && md5_seqid);
  if (sc->c_class->md5_to_seqlen)
    return sc->c_class->md5_to_seqlen(sc, len, md5_seqid, err);
  return 0;
}

GtUword gt_seq_col_num_of_files(const GtSeqCol *sc)
{
  gt_assert(sc);
  if (sc->c_class->num_files)
    return sc->c_class->num_files(sc);
  return 0;
}

GtUword gt_seq_col_num_of_seqs(const GtSeqCol *sc, GtUword filenum)
{
 gt_assert(sc);
  if (sc->c_class->num_seqs)
    return sc->c_class->num_seqs(sc, filenum);
  return 0;
}

const char* gt_seq_col_get_md5_fingerprint(const GtSeqCol *sc,
                                           GtUword filenum,
                                           GtUword seqnum)
{
  gt_assert(sc);
  if (sc->c_class->get_md5)
    return sc->c_class->get_md5(sc, filenum, seqnum);
  return 0;
}

char* gt_seq_col_get_sequence(const GtSeqCol *sc, GtUword filenum,
                              GtUword seqnum,GtUword start,
                              GtUword end)
{
  gt_assert(sc);
  if (sc->c_class->get_seq)
    return sc->c_class->get_seq(sc, filenum, seqnum, start, end);
  return 0;
}

char* gt_seq_col_get_description(const GtSeqCol *sc, GtUword filenum,
                                 GtUword seqnum)
{
  gt_assert(sc);
  if (sc->c_class->get_seq)
    return sc->c_class->get_desc(sc, filenum, seqnum);
  return 0;
}

GtUword gt_seq_col_get_sequence_length(const GtSeqCol *sc,
                                             GtUword filenum,
                                             GtUword seqnum)
{
  gt_assert(sc);
  if (sc->c_class->get_seqlen)
    return sc->c_class->get_seqlen(sc, filenum, seqnum);
  return 0;
}
