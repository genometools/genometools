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
#include "extended/feature_type.h"
#include "extended/type_factory.h"

struct GT_FeatureType {
  GT_TypeFactory *feature_type_factory;
  char *type;
};

GT_FeatureType* gt_feature_type_construct(GT_TypeFactory
                                                 *feature_type_factory,
                                                 const char *type)
{
  GT_FeatureType *gft;
  assert(feature_type_factory && type);
  gft = gt_calloc(1, sizeof *gft);
  gft->feature_type_factory = feature_type_factory;
  gft->type = gt_cstr_dup(type);
  return gft;
}

GT_FeatureType* gt_feature_type_create_gft(GT_FeatureType *gft,
                                                  const char *type)
{
  assert(gft && type);
  return gt_type_factory_create_gft(gft->feature_type_factory, type);
}

void gt_feature_type_delete(GT_FeatureType *gft)
{
  if (!gft) return;
  gt_free(gft->type);
  gt_free(gft);
}

bool gt_feature_type_is(GT_FeatureType *gft, const char *type)
{
  if (gft == gt_type_factory_create_gft(gft->feature_type_factory, type))
    return true;
  return false;
}

const char* gt_feature_type_get_cstr(const GT_FeatureType *gft)
{
  assert(gft);
  return gft->type;
}

GT_TypeFactory* gt_feature_type_get_ftf(const GT_FeatureType *gft)
{
  assert(gft);
  return gft->feature_type_factory;
}
