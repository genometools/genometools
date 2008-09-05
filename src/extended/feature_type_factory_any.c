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
#include "core/hashtable.h"
#include "core/ma.h"
#include "extended/genome_feature_type_imp.h"
#include "extended/feature_type_factory_any.h"
#include "extended/feature_type_factory_rep.h"

struct FeatureTypeFactoryAny {
  const FeatureTypeFactory parent_instance;
};

#define feature_type_factory_any_cast(FTF)\
        feature_type_factory_cast(feature_type_factory_any_class(), FTF)

static GenomeFeatureType*
feature_type_factory_any_create_gft(FeatureTypeFactory *ftf,
                                        const char *type)
{
  FeatureTypeFactoryAny *ftfa;
  GenomeFeatureType *gft = NULL;
  assert(ftf && type);
  ftfa = feature_type_factory_any_cast(ftf);
  if (!(gft = gft_collection_get(ftf->used_types, type))) {
    gft = genome_feature_type_construct(ftf, type);
    gft_collection_add(ftf->used_types, type, gft);
  }
  return gft;
}

const FeatureTypeFactoryClass* feature_type_factory_any_class(void)
{
  static const FeatureTypeFactoryClass feature_type_factory_class =
    { sizeof (FeatureTypeFactoryAny),
      feature_type_factory_any_create_gft,
      NULL };
  return &feature_type_factory_class;
}

FeatureTypeFactory* feature_type_factory_any_new(void)
{
  FeatureTypeFactoryAny *ftfa;
  FeatureTypeFactory *ftf;
  ftf = feature_type_factory_create(feature_type_factory_any_class());
  ftfa = feature_type_factory_any_cast(ftf);
  return ftf;
}
