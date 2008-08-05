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

#ifndef FEATURE_TYPE_FACTORY_REP_H
#define FEATURE_TYPE_FACTORY_REP_H

#include <stdio.h>
#include "libgtext/feature_type_factory.h"
#include "libgtext/genome_feature_type.h"
#include "libgtext/gft_collection.h"

struct FeatureTypeFactoryClass {
  size_t size;
  GenomeFeatureType* (*create_gft)(FeatureTypeFactory*, const char *type);
  void               (*free)(FeatureTypeFactory*);
};

struct FeatureTypeFactory {
  const FeatureTypeFactoryClass *c_class;
  GFTCollection *used_types;
  unsigned int reference_count;
};

FeatureTypeFactory* feature_type_factory_create(const FeatureTypeFactoryClass*);
void*               feature_type_factory_cast(const FeatureTypeFactoryClass*,
                                              FeatureTypeFactory*);

#endif
