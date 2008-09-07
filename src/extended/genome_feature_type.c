/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include <assert.h>
#include "core/cstr.h"
#include "core/ma.h"
#include "extended/feature_type_factory.h"
#include "extended/genome_feature_type.h"

struct GT_GenomeFeatureType {
  GT_FeatureTypeFactory *feature_type_factory;
  char *type;
};

GT_GenomeFeatureType* gt_genome_feature_type_construct(GT_FeatureTypeFactory
                                                 *feature_type_factory,
                                                 const char *type)
{
  GT_GenomeFeatureType *gft;
  assert(feature_type_factory && type);
  gft = ma_calloc(1, sizeof *gft);
  gft->feature_type_factory = feature_type_factory;
  gft->type = cstr_dup(type);
  return gft;
}

GT_GenomeFeatureType* gt_genome_feature_type_create_gft(GT_GenomeFeatureType *gft,
                                                  const char *type)
{
  assert(gft && type);
  return gt_feature_type_factory_create_gft(gft->feature_type_factory, type);
}

void gt_genome_feature_type_delete(GT_GenomeFeatureType *gft)
{
  if (!gft) return;
  ma_free(gft->type);
  ma_free(gft);
}

bool gt_genome_feature_type_is(GT_GenomeFeatureType *gft, const char *type)
{
  if (gft == gt_feature_type_factory_create_gft(gft->feature_type_factory, type))
    return true;
  return false;
}

const char* gt_genome_feature_type_get_cstr(const GT_GenomeFeatureType *gft)
{
  assert(gft);
  return gft->type;
}

GT_FeatureTypeFactory* gt_genome_feature_type_get_ftf(const GT_GenomeFeatureType *gft)
{
  assert(gft);
  return gft->feature_type_factory;
}
