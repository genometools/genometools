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

#include <assert.h>
#include "libgtcore/ma.h"
#include "libgtcore/unused.h"
#include "libgtext/feature_type_factory_rep.h"

FeatureTypeFactory* feature_type_factory_create(const FeatureTypeFactoryClass
                                                *ftfc)
{
  FeatureTypeFactory *ftf;
  assert(ftfc && ftfc->size);
  ftf = ma_calloc(1, ftfc->size);
  ftf->c_class = ftfc;
  ftf->used_types = gft_collection_new();
  return ftf;
}

FeatureTypeFactory* feature_type_factory_ref(FeatureTypeFactory *ftf)
{
  assert(ftf);
  ftf->reference_count++;
  return ftf;
}

GenomeFeatureType* feature_type_factory_create_gft(FeatureTypeFactory *ftf,
                                                   const char *type)
{
  assert(ftf && ftf->c_class && ftf->c_class->create_gft);
  return ftf->c_class->create_gft(ftf, type);
}

StrArray* feature_type_factory_get_used_types(const FeatureTypeFactory *ftf)
{
  return gft_collection_get_types(ftf->used_types);
}

void feature_type_factory_delete(FeatureTypeFactory *ftf)
{
  if (!ftf) return;
  if (ftf->reference_count) {
    ftf->reference_count--;
    return;
  }
  assert(ftf->c_class);
  if (ftf->c_class->free)
    ftf->c_class->free(ftf);
  gft_collection_delete(ftf->used_types);
  ma_free(ftf);
}

void* feature_type_factory_cast(UNUSED const FeatureTypeFactoryClass *ftfc,
                                FeatureTypeFactory *ftf)
{
  assert(ftfc && ftf && ftf->c_class == ftfc);
  return ftf;
}
