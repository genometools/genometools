/*
  Copyright (c) 2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include "core/cstr_api.h"
#include "core/hashmap_api.h"
#include "core/ma.h"
#include "core/seq_info_cache.h"

struct GtSeqInfoCache {
  GtHashmap *cache;
};

GtSeqInfoCache* gt_seq_info_cache_new(void)
{
  GtSeqInfoCache *sic = gt_malloc(sizeof *sic);
  sic->cache = gt_hashmap_new(GT_HASH_STRING, gt_free_func, gt_free_func);
  return sic;
}

void gt_seq_info_cache_delete(GtSeqInfoCache *sic)
{
  if (!sic) return;
  gt_hashmap_delete(sic->cache);
  gt_free(sic);
}

const GtSeqInfo* gt_seq_info_cache_get(GtSeqInfoCache *sic, const char *key)
{
  gt_assert(sic);
  return gt_hashmap_get(sic->cache, key);
}

void gt_seq_info_cache_add(GtSeqInfoCache *sic, const char *key,
                           const GtSeqInfo *si)
{
  GtSeqInfo *si_dup;
  gt_assert(sic && key && si);
  gt_assert(!gt_seq_info_cache_get(sic, key));
  si_dup = gt_malloc(sizeof *si_dup);
  si_dup->filenum = si->filenum;
  si_dup->seqnum = si->seqnum;
  gt_hashmap_add(sic->cache, gt_cstr_dup(key), si_dup);
}
