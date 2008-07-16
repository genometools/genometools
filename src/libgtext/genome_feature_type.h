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

#ifndef GENOME_FEATURE_TYPE_H
#define GENOME_FEATURE_TYPE_H

#include <stdbool.h>

/* The GenomeFeatureType represents the feature type mainly used in
   GenomeFeatures and corresponds to the type column in GFF3 files.
   To create new GenomeFeatureTypes a FeatureTypeFactory or an already existing
   GenomeFeatureType (which internally uses the FeatureTypeFactory it was
   created from) should be used!
*/

/* Some predefined (genome feature) type strings. */
#define gft_CDS                      "CDS"
#define gft_EST_match                "EST_match"
#define gft_LTR_retrotransposon      "LTR_retrotransposon"
#define gft_SNP                      "SNP"
#define gft_TF_binding_site          "TF_binding_site"
#define gft_cDNA_match               "cDNA_match"
#define gft_exon                     "exon"
#define gft_five_prime_UTR           "five_prime_UTR"
#define gft_five_prime_splice_site   "five_prime_splice_site"
#define gft_gene                     "gene"
#define gft_intron                   "intron"
#define gft_inverted_repeat          "inverted_repeat"
#define gft_long_terminal_repeat     "long_terminal_repeat"
#define gft_mRNA                     "mRNA"
#define gft_protein_match            "protein_match"
#define gft_repeat_region            "repeat_region"
#define gft_target_site_duplication  "target_site_duplication"
#define gft_three_prime_UTR          "three_prime_UTR"
#define gft_three_prime_splice_site  "three_prime_splice_site"
#define gft_transcript               "transcript"

typedef struct GenomeFeatureType GenomeFeatureType;

#include "libgtext/feature_type_factory.h"

GenomeFeatureType* genome_feature_type_construct(FeatureTypeFactory*,
                                                 const char *type);
GenomeFeatureType* genome_feature_type_create_gft(GenomeFeatureType*,
                                                  const char *type);
void               genome_feature_type_delete(GenomeFeatureType*);
bool               genome_feature_type_is(GenomeFeatureType*, const char *type);
const char*        genome_feature_type_get_cstr(const GenomeFeatureType*);

#endif
