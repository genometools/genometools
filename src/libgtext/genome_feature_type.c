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
#include <stdlib.h>
#include <string.h>
#include "libgtext/compare.h"
#include "libgtext/genome_feature_type.h"

/*
  Keep in sync with ``genome_feature_type'' in ``genome_feature.h''!
  The feature types have to be sorted. We assume ASCII encoding.
*/

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

int genome_feature_type_get(GenomeFeatureType *type, const char *gft_string)
{
  void *result;

  assert(type && gft_string);
  assert(strcmp(gft_string, "undefined")); /* do not convert undefined string */

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
  assert(type != undefined); /* do not convert undefined type */
  return genome_feature_type_strings[type];
}

unsigned long genome_feature_type_num_of_features(void)
{
  return sizeof (genome_feature_type_strings) /
         sizeof (genome_feature_type_strings[0]);
}
