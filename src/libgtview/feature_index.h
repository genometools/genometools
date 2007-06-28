#ifndef FEATUREINDEX_H
#define FEATUREINDEX_H

#include <libgtcore/array.h>
#include <libgtcore/str.h>
#include <libgtext/sequence_region.h>
#include <libgtext/genome_feature.h>

typedef struct FeatureIndex FeatureIndex;

FeatureIndex* 	feature_index_new(Env*);
void	  feature_index_delete(FeatureIndex*,
                             Env*);
int     feature_index_add_sequence_region(FeatureIndex*,
                                             char*,
                                             Env*);
void	  feature_index_add_genome_feature_for_seqid(FeatureIndex*,
                                                   GenomeFeature*,
                                                   Env*);
Array*  feature_index_get_features_for_seqid(FeatureIndex*,
                                             char*);
int  	  feature_index_get_features_for_range(FeatureIndex*,
                                             Array*,
                                             char*,
                                             Range,
                                             Env*);
void 	  feature_index_print_contents(FeatureIndex*,
                                     Env*);
bool feature_index_has_seqid(FeatureIndex*,
                             char* seqid,
                             Env* env);

int genome_node_print_feature_children(GenomeNode*, void*, Env*);
int feature_index_unit_test(Env*);

#endif
