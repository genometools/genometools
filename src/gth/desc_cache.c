/*
  Copyright (c) 2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "core/str_cache.h"
#include "gth/desc_cache.h"

struct GthDescCache {
  GtStrCache *str_cache;
};

static GtStr* get_desc_str(void *str_source, unsigned long index)
{
  GtStr *str;
  GthSeqCol *seq_col = str_source;
  gt_assert(seq_col);
  str = gt_str_new();
  gth_seq_col_get_description(seq_col, index, str);
  return str;
}

GthDescCache* gth_desc_cache_new(GthSeqCol *seq_col)
{
  GthDescCache *desc_cache;
  gt_assert(seq_col);
  desc_cache = gt_malloc(sizeof *desc_cache);
  desc_cache->str_cache = gt_str_cache_new(seq_col, get_desc_str,
                                           gth_seq_col_num_of_seqs(seq_col));
  return desc_cache;
}

void gth_desc_cache_delete(GthDescCache *desc_cache)
{
  if (!desc_cache) return;
  gt_str_cache_delete(desc_cache->str_cache);
  gt_free(desc_cache);
}

GtStr* gth_desc_cache_get(GthDescCache *desc_cache, unsigned long seq_num)
{
  gt_assert(desc_cache);
  return gt_str_cache_get(desc_cache->str_cache, seq_num);
}
