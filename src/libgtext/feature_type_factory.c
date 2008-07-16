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
#include "libgtext/genome_feature_type.h"
#include "libgtext/feature_type_factory.h"

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
  result =  bsearch(&gft_string,
                    genome_feature_type_strings,
                    sizeof (genome_feature_type_strings) /
                    sizeof (genome_feature_type_strings[0]),
                    sizeof (char*),
                    compare);
  if (result)
    return *(char**) result;
  return NULL;
}

struct FeatureTypeFactory {
  Hashtable *genome_feature_types;
};

FeatureTypeFactory* feature_type_factory_new(void)
{
  FeatureTypeFactory *ftf = ma_calloc(1, sizeof *ftf);
  ftf->genome_feature_types = hashtable_new(HASH_STRING, NULL,
                                            (FreeFunc)
                                            genome_feature_type_delete);
  return ftf;
}

void feature_type_factory_delete(FeatureTypeFactory *ftf)
{
  if (!ftf) return;
  hashtable_delete(ftf->genome_feature_types);
  ma_free(ftf);
}

GenomeFeatureType* feature_type_factory_create_gft(FeatureTypeFactory *ftf,
                                                   const char *type)
{
  GenomeFeatureType *gft = NULL;
  assert(ftf && type);
  if (!(gft = hashtable_get(ftf->genome_feature_types, type))) {
    const char *static_type;
    if ((static_type = find_type(type))) {
      gft = genome_feature_type_construct(ftf, static_type);
      hashtable_add(ftf->genome_feature_types, (char*) static_type, gft);
    }
  }
  return gft;
}
