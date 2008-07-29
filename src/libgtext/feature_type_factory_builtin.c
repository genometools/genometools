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
#include <stdlib.h>
#include <string.h>
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"
#include "libgtext/compare.h"
#include "libgtext/genome_feature_type_imp.h"
#include "libgtext/feature_type_factory_builtin.h"
#include "libgtext/feature_type_factory_rep.h"

struct FeatureTypeFactoryBuiltin {
  const FeatureTypeFactory parent_instance;
};

#define feature_type_factory_builtin_cast(FTF)\
        feature_type_factory_cast(feature_type_factory_builtin_class(), FTF)

static const char *genome_feature_type_strings[] = { "CDS",
                                                     "EST_match",
                                                     "LTR_retrotransposon",
                                                     "SNP",
                                                     "TF_binding_site",
                                                     "cDNA_match",
                                                     "exon",
                                                     "five_prime_UTR",
                                                     "five_prime_splice_site",
                                                     "gene",
                                                     "intron",
                                                     "inverted_repeat",
                                                     "long_terminal_repeat",
                                                     "mRNA",
                                                     "protein_match",
                                                     "repeat_region",
                                                     "target_site_duplication",
                                                     "three_prime_UTR",
                                                     "three_prime_splice_site",
                                                     "transcript",
                                                     "undefined"
                                                   };

static const char* find_type(const char *gft_string)
{
  void *result;
  assert(gft_string);
  assert(strcmp(gft_string, "undefined")); /* do not convert undefined string */
  result = bsearch(&gft_string,
                   genome_feature_type_strings,
                   sizeof (genome_feature_type_strings) /
                   sizeof (genome_feature_type_strings[0]),
                   sizeof (char*),
                   compare);
  if (result)
    return *(char**) result;
  return NULL;
}

static GenomeFeatureType*
feature_type_factory_builtin_create_gft(FeatureTypeFactory *ftf,
                                        const char *type)
{
  FeatureTypeFactoryBuiltin *ftfb;
  GenomeFeatureType *gft = NULL;
  assert(ftf && type);
  ftfb = feature_type_factory_builtin_cast(ftf);
  if (!(gft = gft_collection_get(ftf->used_types, type))) {
    if ((find_type(type))) {
      gft = genome_feature_type_construct(ftf, type);
      gft_collection_add(ftf->used_types, type, gft);
    }
  }
  return gft;
}

const FeatureTypeFactoryClass* feature_type_factory_builtin_class(void)
{
  static const FeatureTypeFactoryClass feature_type_factory_class =
    { sizeof (FeatureTypeFactoryBuiltin),
      feature_type_factory_builtin_create_gft,
      NULL };
  return &feature_type_factory_class;
}

FeatureTypeFactory* feature_type_factory_builtin_new(void)
{
  FeatureTypeFactoryBuiltin *ftfb;
  FeatureTypeFactory *ftf;
  ftf = feature_type_factory_create(feature_type_factory_builtin_class());
  ftfb = feature_type_factory_builtin_cast(ftf);
  return ftf;
}
