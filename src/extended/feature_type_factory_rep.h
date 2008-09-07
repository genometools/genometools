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
#include "extended/feature_type_factory.h"
#include "extended/genome_feature_type.h"
#include "extended/gft_collection.h"

struct GT_FeatureTypeFactoryClass {
  size_t size;
  GT_GenomeFeatureType* (*create_gft)(GT_FeatureTypeFactory*, const char *type);
  void               (*free)(GT_FeatureTypeFactory*);
};

struct GT_FeatureTypeFactory {
  const GT_FeatureTypeFactoryClass *c_class;
  GFTCollection *used_types;
  unsigned int reference_count;
};

GT_FeatureTypeFactory* gt_feature_type_factory_create(const GT_FeatureTypeFactoryClass*);
void*               gt_feature_type_factory_cast(const GT_FeatureTypeFactoryClass*,
                                              GT_FeatureTypeFactory*);

#endif
