/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
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
  gft_TF_binding_site,
  gft_exon,
  gft_gene,
  gft_intron,
  gft_mRNA
} Genome_feature_type;

/*
  Determine a genome feature type ``gft'' from the string ``gft_string''.
  If such a feature does not exits, -1 is returned.
*/
int           genome_feature_type_get(Genome_feature_type*, char *gft_string);
const char*   genome_feature_type_get_cstr(Genome_feature_type);
unsigned long genome_feature_type_num_of_features(void);

#endif
