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
#include "libgtext/genome_feature_type_imp.h"
#include "libgtext/feature_type_factory_any.h"
#include "libgtext/feature_type_factory_rep.h"

struct FeatureTypeFactoryAny {
  const FeatureTypeFactory parent_instance;
  Hashtable *genome_feature_types;
};

#define feature_type_factory_any_cast(FTF)\
        feature_type_factory_cast(feature_type_factory_any_class(), FTF)

static void feature_type_factory_any_free(FeatureTypeFactory *ftf)
{
  FeatureTypeFactoryAny *ftfa = feature_type_factory_any_cast(ftf);
  hashtable_delete(ftfa->genome_feature_types);
}

static GenomeFeatureType*
feature_type_factory_any_create_gft(FeatureTypeFactory *ftf,
                                        const char *type)
{
  FeatureTypeFactoryAny *ftfa;
  GenomeFeatureType *gft = NULL;
  assert(ftf && type);
  ftfa = feature_type_factory_any_cast(ftf);
  if (!(gft = hashtable_get(ftfa->genome_feature_types, type))) {
    char *type_dup = cstr_dup(type);
    gft = genome_feature_type_construct(ftf, type_dup);
    hashtable_add(ftfa->genome_feature_types, type_dup, gft);
  }
  return gft;
}

const FeatureTypeFactoryClass* feature_type_factory_any_class(void)
{
  static const FeatureTypeFactoryClass feature_type_factory_class =
    { sizeof (FeatureTypeFactoryAny),
      feature_type_factory_any_create_gft,
      feature_type_factory_any_free };
  return &feature_type_factory_class;
}

FeatureTypeFactory* feature_type_factory_any_new(void)
{
  FeatureTypeFactoryAny *ftfa;
  FeatureTypeFactory *ftf;
  ftf = feature_type_factory_create(feature_type_factory_any_class());
  ftfa = feature_type_factory_any_cast(ftf);
  ftfa->genome_feature_types = hashtable_new(HASH_STRING, ma_free_func,
                                             (FreeFunc)
                                             genome_feature_type_delete);
  return ftf;
}
