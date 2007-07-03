/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Malte Mader <mmader@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FEATUREINDEX_H
#define FEATUREINDEX_H

#include <libgtcore/array.h>
#include <libgtcore/str.h>
#include <libgtext/sequence_region.h>
#include <libgtext/genome_feature.h>

typedef struct FeatureIndex FeatureIndex;

FeatureIndex* feature_index_new(Env*);
int           feature_index_add_sequence_region(FeatureIndex*, char*, Env*);
void          feature_index_add_genome_feature_for_seqid(FeatureIndex*,
                                                         GenomeFeature*, Env*);
Array*        feature_index_get_features_for_seqid(FeatureIndex*,
                                             char*);
int           feature_index_get_features_for_range(FeatureIndex*, Array*, char*,
                                                   Range, Env*);
char*         feature_index_get_first_seqid(FeatureIndex*);
Range         feature_index_get_range_for_seqid(FeatureIndex*, char*);
void          feature_index_print_contents(FeatureIndex*, Env*);
bool          feature_index_has_seqid(FeatureIndex*, char* seqid, Env* env);

/* XXX: move this function to proper place */
int           genome_node_print_feature_children(GenomeNode*, void*, Env*);

int           feature_index_unit_test(Env*);
void          feature_index_delete(FeatureIndex*, Env*);

#endif
