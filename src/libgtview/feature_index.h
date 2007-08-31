/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>,
  Copyright (c) 2007 Malte Mader <mmader@stud.zbh.uni-hamburg.de>,
  Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FEATUREINDEX_H
#define FEATUREINDEX_H

#include "libgtcore/array.h"
#include "libgtcore/str.h"
#include "libgtext/sequence_region.h"
#include "libgtext/genome_feature.h"

typedef struct FeatureIndex FeatureIndex;

FeatureIndex* feature_index_new(Env*);
FeatureIndex* feature_index_ref(FeatureIndex*);
void          feature_index_add_sequence_region(FeatureIndex*, SequenceRegion*,
                                                Env*);
/* Add a GenomeFeature to the index, associating it with a sequence region
   denoted by its identifier string. */
void          feature_index_add_genome_feature(FeatureIndex*, GenomeFeature*,
                                               Env*);
/* Returns an array of GenomeFeatures associated with a given sequence region
   identifier. */
Array*        feature_index_get_features_for_seqid(FeatureIndex*, const char*);
/* Look up genome features for sequence region <seqid> in <range> and store them
   in <results> */
int           feature_index_get_features_for_range(FeatureIndex*,
                                                   Array *results,
                                                   const char *seqid,
                                                   Range range, Env*);
/* Returns the first sequence region identifier added to the index. */
const char*   feature_index_get_first_seqid(FeatureIndex*);
Range         feature_index_get_range_for_seqid(FeatureIndex*, const char*);
bool          feature_index_has_seqid(const FeatureIndex*, const char*, Env*);
int           feature_index_unit_test(Env*);
void          feature_index_delete(FeatureIndex*, Env*);

#endif
