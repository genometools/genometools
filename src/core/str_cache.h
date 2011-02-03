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

#ifndef STR_CACHE_H
#define STR_CACHE_H

#include "core/str_api.h"

typedef struct GtStrCache GtStrCache;

typedef GtStr* (*GtStrConstructorFunc)(void *str_source, unsigned long index);

/* Create a new string cache object for <num_of_strings> many strings creatable
   from <str_source> with the function <str_constructor>.
   That is, the first time a certain string is requested with the method
   <gt_str_cache_get()>, a new string object is created via <str_constructor>
   and then cached and returned. Subsequent calls to <gt_str_cache_get()> for
   the same string return a new reference made from the cached string. */
GtStrCache* gt_str_cache_new(void *str_source,
                             GtStrConstructorFunc str_constructor,
                             unsigned long num_of_strings);
void        gt_str_cache_delete(GtStrCache*);
/* Return a new (i.e., the caller is responsible to free it) <GtStr*> object
   for string with given <index>. The mechanics of the cache are described in
   detail in the documentation of <gt_str_cache_new()>. */
GtStr*      gt_str_cache_get(GtStrCache*, unsigned long index);

#endif
