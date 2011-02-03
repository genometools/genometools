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

struct GtStrCache {
  GtStr **cache;
  void *source;
  GtStrConstructorFunc constructor;
  unsigned long num_of_strings;
};

GtStrCache* gt_str_cache_new(void *str_source,
                             GtStrConstructorFunc str_constructor,
                             unsigned long num_of_strings)
{
  GtStrCache *str_cache;
  gt_assert(str_source && str_constructor && num_of_strings);
  str_cache = gt_malloc(sizeof *str_cache);
  str_cache->cache = gt_calloc(num_of_strings, sizeof (GtStr*));
  str_cache->source = str_source;
  str_cache->constructor = str_constructor;
  str_cache->num_of_strings = num_of_strings;
  return str_cache;
}

void gt_str_cache_delete(GtStrCache *str_cache)
{
  unsigned long i;
  if (!str_cache) return;
  for (i = 0; i < str_cache->num_of_strings; i++)
    gt_str_delete(str_cache->cache[i]);
  gt_free(str_cache->cache);
  gt_free(str_cache);
}

GtStr* gt_str_cache_get(GtStrCache *str_cache, unsigned long index)
{
  gt_assert(str_cache && index < str_cache->num_of_strings);
  if (!str_cache->cache[index])
    str_cache->cache[index] = str_cache->constructor(str_cache->source, index);
  gt_assert(str_cache->cache[index]);
  return gt_str_ref(str_cache->cache[index]);
}
