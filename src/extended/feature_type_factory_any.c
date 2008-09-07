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

struct GT_FeatureTypeFactoryAny {
  const GT_FeatureTypeFactory parent_instance;
};

#define gt_feature_type_factory_any_cast(FTF)\
        gt_feature_type_factory_cast(gt_feature_type_factory_any_class(), FTF)

static GT_GenomeFeatureType*
gt_feature_type_factory_any_create_gft(GT_FeatureTypeFactory *ftf,
                                        const char *type)
{
  GT_FeatureTypeFactoryAny *ftfa;
  GT_GenomeFeatureType *gft = NULL;
  assert(ftf && type);
  ftfa = gt_feature_type_factory_any_cast(ftf);
  if (!(gft = gft_collection_get(ftf->used_types, type))) {
    gft = gt_genome_feature_type_construct(ftf, type);
    gft_collection_add(ftf->used_types, type, gft);
  }
  return gft;
}

const GT_FeatureTypeFactoryClass* gt_feature_type_factory_any_class(void)
{
  static const GT_FeatureTypeFactoryClass gt_feature_type_factory_class =
    { sizeof (GT_FeatureTypeFactoryAny),
      gt_feature_type_factory_any_create_gft,
      NULL };
  return &gt_feature_type_factory_class;
}

GT_FeatureTypeFactory* gt_feature_type_factory_any_new(void)
{
  GT_FeatureTypeFactoryAny *ftfa;
  GT_FeatureTypeFactory *ftf;
  ftf = gt_feature_type_factory_create(gt_feature_type_factory_any_class());
  ftfa = gt_feature_type_factory_any_cast(ftf);
  return ftf;
}
