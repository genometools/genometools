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

#ifndef FEATURE_TYPE_FACTORY_OBO_H
#define FEATURE_TYPE_FACTORY_OBO_H

#include "extended/feature_type_factory.h"

/* Implements the GT_FeatureTypeFactory interface with types from an OBO file. */
typedef struct GT_FeatureTypeFactoryOBO GT_FeatureTypeFactoryOBO;

const GT_FeatureTypeFactoryClass* gt_feature_type_factory_obo_class(void);
GT_FeatureTypeFactory*            gt_feature_type_factory_obo_new(const char
                                                            *obo_file_path,
                                                            GT_Error*);

#endif
