/*
  Copyright (c) 2011-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include "core/ma_api.h"
#include "core/md5_tab.h"
#include "core/str_cache.h"
#include "gth/gthdef.h"
#include "gth/md5_cache.h"

struct GthMD5Cache {
  GtMD5Tab *md5_tab;
  GtStrCache *str_cache;
};

static const char* seq_con_get_seq(void *seqs, unsigned long index)
{
  GthSeqCon *seq_con = seqs;
  const char *seq;
  gt_assert(seq_con);
  seq = (const char*) gth_seq_con_get_orig_seq(seq_con, index);
  if (!seq) {
    gth_seq_con_demand_orig_seq(seq_con);
    seq = (const char*) gth_seq_con_get_orig_seq(seq_con, index);
    gt_assert(seq);
  }
  return seq;
}

static unsigned long seq_con_get_seq_len(void *seqs, unsigned long index)
{
  GthSeqCon *seq_con = seqs;
  gt_assert(seq_con);
  return gth_seq_con_get_length(seq_con, index);
}

static GtStr* get_md5_str(void *str_source, unsigned long index)
{
  const GtMD5Tab *md5_tab = str_source;
  gt_assert(md5_tab);
  return gt_str_new_cstr(gt_md5_tab_get(md5_tab, index));
}

GthMD5Cache* gth_md5_cache_new(const char *indexname, GthSeqCon *seq_con)
{
  GthMD5Cache *md5_cache;
  gt_assert(indexname && seq_con);
  md5_cache = gt_malloc(sizeof *md5_cache);
  md5_cache->md5_tab = gt_md5_tab_new(indexname, seq_con, seq_con_get_seq,
                                      seq_con_get_seq_len,
                                      gth_seq_con_num_of_seqs(seq_con), true,
                                      !getenv(GTHNOFLOCKENVNAME));
  md5_cache->str_cache = gt_str_cache_new(md5_cache->md5_tab, get_md5_str,
                                          gt_md5_tab_size(md5_cache->md5_tab));
  return md5_cache;
}

void gth_md5_cache_delete(GthMD5Cache *md5_cache)
{
  if (!md5_cache) return;
  gt_str_cache_delete(md5_cache->str_cache);
  gt_md5_tab_delete(md5_cache->md5_tab);
  gt_free(md5_cache);
}

GtStr* gth_md5_cache_get(GthMD5Cache *md5_cache, unsigned long seq_num)
{
  gt_assert(md5_cache);
  return gt_str_cache_get(md5_cache->str_cache, seq_num);
}
