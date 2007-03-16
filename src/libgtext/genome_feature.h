/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_FEATURE_H
#define GENOME_FEATURE_H

/* implements the ``genome node'' interface */
typedef struct GenomeFeature GenomeFeature;

#include <gtcore.h>
#include <libgtext/genome_node.h>
#include <libgtext/genome_feature_type.h>
#include <libgtext/transcript_feature_type.h>

typedef int (*AttributeIterFunc)(const char *attr_name, const char *attr_value,
                                 void *data, Env*);

const GenomeNodeClass* genome_feature_class(void);
GenomeNode*            genome_feature_new(GenomeFeatureType, Range, Strand,
                                          const char *filename,
                                          unsigned long line_number, Env*);
const char*            genome_feature_get_source(GenomeFeature*);
GenomeFeatureType      genome_feature_get_type(GenomeFeature*);
double                 genome_feature_get_score(GenomeFeature*);
Strand                 genome_feature_get_strand(GenomeFeature*);
Phase                  genome_feature_get_phase(GenomeFeature*);
void                   genome_feature_get_exons(GenomeFeature*,
                                                Array *exon_features, Env*);
void                   genome_feature_determine_transcripttypes(GenomeFeature*,
                                                                Env*);
TranscriptFeatureType  genome_feature_get_transcriptfeaturetype(GenomeFeature*);
void                   genome_feature_set_end(GenomeFeature*, unsigned long);
void                   genome_feature_set_score(GenomeFeature*, double);
void                   genome_feature_add_attribute(GenomeFeature*,
                                                    const char *attr_name,
                                                    const char *attr_value,
                                                    Env*);
bool                   genome_feature_has_attribute(const GenomeFeature*);
int                    genome_feature_foreach_attribute(GenomeFeature*,
                                                        AttributeIterFunc,
                                                        void *data, Env*);

#endif
