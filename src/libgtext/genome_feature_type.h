/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_FEATURE_TYPE_H
#define GENOME_FEATURE_TYPE_H

/*
  Keep in sync with `genome_feature_type_strings' in `genome_feature_type.c'!
  The feature types have to be sorted. We assume ASCII encoding.
*/

typedef enum {
  gft_CDS,
  gft_EST_match,
  gft_LTR_retrotransposon,
  gft_TF_binding_site,
  gft_cDNA_match,
  gft_exon,
  gft_five_prime_UTR,
  gft_gene,
  gft_intron,
  gft_inverted_repeat,
  gft_long_terminal_repeat,
  gft_mRNA,
  gft_protein_match,
  gft_repeat_region,
  gft_target_site_duplication,
  gft_three_prime_UTR,
  gft_transcript,
} GenomeFeatureType;

/*
  Determine a genome feature type ``gft'' from the string ``gft_string''.
  If such a feature does not exits, -1 is returned.
*/
int           genome_feature_type_get(GenomeFeatureType*, char *gft_string);
const char*   genome_feature_type_get_cstr(GenomeFeatureType);
unsigned long genome_feature_type_num_of_features(void);

#endif
