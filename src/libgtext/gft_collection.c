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

#include "libgtcore/cstr.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"
#include "libgtcore/unused.h"
#include "libgtext/genome_feature_type_imp.h"
#include "libgtext/gft_collection.h"

struct GFTCollection {
  Hashtable *genome_feature_types;
};

GFTCollection* gft_collection_new(void)
{
  GFTCollection *gftc = ma_malloc(sizeof (GFTCollection));
  gftc->genome_feature_types = hashtable_new(HASH_STRING, ma_free_func,
                                             (FreeFunc)
                                             genome_feature_type_delete);
  return gftc;
}

void gft_collection_delete(GFTCollection *gftc)
{
  if (!gftc) return;
  hashtable_delete(gftc->genome_feature_types);
  ma_free(gftc);
}

void gft_collection_add(GFTCollection *gftc, const char *type,
                         GenomeFeatureType *gft)
{
  assert(gftc && type && gft);
  hashtable_add(gftc->genome_feature_types, cstr_dup(type), gft);
}

GenomeFeatureType* gft_collection_get(GFTCollection *gftc, const char *type)
{
  assert(gftc && type);
  return hashtable_get(gftc->genome_feature_types, type);
}

static int store_type(void *key, UNUSED void *value, void *data,
                      UNUSED Error *err)
{
  StrArray *types = data;
  assert(key && types);
  strarray_add_cstr(types, key);
  return 0;
}

StrArray* gft_collection_get_types(const GFTCollection *gftc)
{
  StrArray *types;
  int had_err;
  assert(gftc);
  types = strarray_new();
  had_err = hashtable_foreach_ao(gftc->genome_feature_types, store_type, types,
                                 NULL);
  assert(!had_err);
  return types;
}
