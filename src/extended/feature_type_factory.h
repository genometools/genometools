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

/* The FeatureTypeFactory interface. Implementations of the FeatureTypeFactory
   interface are used to create GenomeFeatureTypes.
   Since a FeatureTypeFactory owns all GenomeFeatureTypes it creates, you have
   to make sure to keep it around until all references to the created
   GenomeFeatureTypes have been removed.
*/

typedef struct FeatureTypeFactoryClass FeatureTypeFactoryClass;
typedef struct FeatureTypeFactory FeatureTypeFactory;

#include "core/strarray.h"
#include "extended/genome_feature_type.h"

/* Return a new reference to <feature_type_factory>. */
FeatureTypeFactory* feature_type_factory_ref(FeatureTypeFactory
                                             *feature_type_factory);
/* Uses the factory to create a new genome feature type object of the given
   <type>. Returns NULL, if <type> is not a valid type. */
GenomeFeatureType*  feature_type_factory_create_gft(FeatureTypeFactory*,
                                                    const char *type);
/* Returns a StrArray which contains all type names in alphabetical order which
   have been created by this factory. The caller is responsible to free it! */
StrArray*           feature_type_factory_get_used_types(const
                                                        FeatureTypeFactory*);
void                feature_type_factory_delete(FeatureTypeFactory*);

#endif
