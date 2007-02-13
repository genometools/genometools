/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_FEATURE_H
#define GENOME_FEATURE_H

/* implements the ``genome node'' interface */
typedef struct Genome_feature Genome_feature;

#include "genome_node.h"
#include "genome_feature_type.h"
#include "phase.h"
#include "range.h"
#include "strand.h"

const GenomeNodeClass* genome_feature_class(void);
GenomeNode*        genome_feature_new(GenomeFeatureType, Range, Strand,
                                       const char *filename,
                                       unsigned long line_number);
const char*         genome_feature_get_source(Genome_feature*);
GenomeFeatureType genome_feature_get_type(Genome_feature*);
double              genome_feature_get_score(Genome_feature*);
Strand              genome_feature_get_strand(Genome_feature*);
Phase               genome_feature_get_phase(Genome_feature*);
void                genome_feature_get_exons(Genome_feature*,
                                             Array *exon_features);
void                genome_feature_set_end(Genome_feature*, unsigned long);
/* XXX: move to genome_node */
void                genome_feature_set_score(Genome_feature*, double);

#endif
