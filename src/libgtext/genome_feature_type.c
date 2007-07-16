/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <libgtext/compare.h>
#include <libgtext/genome_feature_type.h>

/*
  Keep in sync with ``genome_feature_type'' in ``genome_feature.h''!
  The feature types have to be sorted. We assume ASCII encoding.
*/

static const char *genome_feature_type_strings[] = { "CDS",
                                                     "EST_match",
                                                     "LTR_retrotransposon",
                                                     "TF_binding_site",
                                                     "cDNA_match",
                                                     "exon",
                                                     "five_prime_UTR",
                                                     "gene",
                                                     "intron",
                                                     "inverted_repeat",
                                                     "long_terminal_repeat",
                                                     "mRNA",
                                                     "protein_match",
                                                     "repeat_region",
                                                     "target_site_duplication",
                                                     "three_prime_UTR",
                                                     "transcript"
                                                   };

int genome_feature_type_get(GenomeFeatureType *type, char *gft_string)
{
  void *result;

  assert(type && gft_string);

  result = bsearch(&gft_string,
                   genome_feature_type_strings,
                   sizeof (genome_feature_type_strings) /
                   sizeof (genome_feature_type_strings[0]),
                   sizeof (char*),
                   compare);

  if (result) {
    *type = (GenomeFeatureType)
            ((char**) result - (char**) genome_feature_type_strings);
    return 0;
  }
  /* else type not found */
  return -1;
}

const char* genome_feature_type_get_cstr(GenomeFeatureType type)
{
  return genome_feature_type_strings[type];
}

unsigned long genome_feature_type_num_of_features(void)
{
  return sizeof (genome_feature_type_strings) /
         sizeof (genome_feature_type_strings[0]);
}
