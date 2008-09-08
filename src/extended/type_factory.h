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

#ifndef FEATURE_TYPE_FACTORY_H
#define FEATURE_TYPE_FACTORY_H

/* The GT_TypeFactory interface. Implementations of the GT_TypeFactory
   interface are used to create GT_GenomeFeatureTypes.
   Since a GT_TypeFactory owns all GT_GenomeFeatureTypes it creates, you have
   to make sure to keep it around until all references to the created
   GT_GenomeFeatureTypes have been removed.
*/

typedef struct GT_TypeFactoryClass GT_TypeFactoryClass;
typedef struct GT_TypeFactory GT_TypeFactory;

#include "core/strarray.h"
#include "extended/genome_feature_type.h"

/* Return a new reference to <feature_type_factory>. */
GT_TypeFactory* gt_type_factory_ref(GT_TypeFactory
                                             *feature_type_factory);
/* Uses the factory to create a new genome feature type object of the given
   <type>. Returns NULL, if <type> is not a valid type. */
GT_GenomeFeatureType*  gt_type_factory_create_gft(GT_TypeFactory*,
                                                    const char *type);
/* Returns a GT_StrArray which contains all type names in alphabetical order
   which have been created by this factory.
   The caller is responsible to free it! */
GT_StrArray*           gt_type_factory_get_used_types(const
                                                        GT_TypeFactory*);
void                gt_type_factory_delete(GT_TypeFactory*);

#endif
