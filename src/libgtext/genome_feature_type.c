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
#include "libgtcore/cstr.h"
#include "libgtcore/ma.h"
#include "libgtext/feature_type_factory.h"
#include "libgtext/genome_feature_type.h"

struct GenomeFeatureType {
  FeatureTypeFactory *feature_type_factory;
  char *type;
};

GenomeFeatureType* genome_feature_type_construct(FeatureTypeFactory
                                                 *feature_type_factory,
                                                 const char *type)
{
  GenomeFeatureType *gft;
  assert(feature_type_factory && type);
  gft = ma_calloc(1, sizeof *gft);
  gft->feature_type_factory = feature_type_factory;
  gft->type = cstr_dup(type);
  return gft;
}

GenomeFeatureType* genome_feature_type_create_gft(GenomeFeatureType *gft,
                                                  const char *type)
{
  assert(gft && type);
  return feature_type_factory_create_gft(gft->feature_type_factory, type);
}

void genome_feature_type_delete(GenomeFeatureType *gft)
{
  if (!gft) return;
  ma_free(gft->type);
  ma_free(gft);
}

bool genome_feature_type_is(GenomeFeatureType *gft, const char *type)
{
  if (gft == feature_type_factory_create_gft(gft->feature_type_factory, type))
    return true;
  return false;
}

const char* genome_feature_type_get_cstr(const GenomeFeatureType *gft)
{
  assert(gft);
  return gft->type;
}
