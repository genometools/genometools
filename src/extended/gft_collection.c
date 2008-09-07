/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "core/cstr.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/unused.h"
#include "extended/genome_feature_type_imp.h"
#include "extended/gft_collection.h"

struct GFTCollection {
  Hashmap *genome_feature_types;
};

GFTCollection* gft_collection_new(void)
{
  GFTCollection *gftc = ma_malloc(sizeof (GFTCollection));
  gftc->genome_feature_types = hashmap_new(
    HASH_STRING, ma_free_func, (FreeFunc)genome_feature_type_delete);
  return gftc;
}

void gft_collection_delete(GFTCollection *gftc)
{
  if (!gftc) return;
  hashmap_delete(gftc->genome_feature_types);
  ma_free(gftc);
}

void gft_collection_add(GFTCollection *gftc, const char *type,
                         GT_GenomeFeatureType *gft)
{
  assert(gftc && type && gft);
  hashmap_add(gftc->genome_feature_types, cstr_dup(type), gft);
}

GT_GenomeFeatureType* gft_collection_get(GFTCollection *gftc, const char *type)
{
  assert(gftc && type);
  return hashmap_get(gftc->genome_feature_types, type);
}

static int store_type(void *key, UNUSED void *value, void *data,
                      UNUSED GT_Error *err)
{
  GT_StrArray *types = data;
  assert(key && types);
  gt_strarray_add_cstr(types, key);
  return 0;
}

GT_StrArray* gft_collection_get_types(const GFTCollection *gftc)
{
  GT_StrArray *types;
  int had_err;
  assert(gftc);
  types = gt_strarray_new();
  had_err = hashmap_foreach_in_key_order(gftc->genome_feature_types,
                                         store_type, types, NULL);
  assert(!had_err);
  return types;
}
