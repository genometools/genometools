/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>,
  Copyright (c) 2007 Malte Mader <mmader@stud.zbh.uni-hamburg.de>,
  Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef FEATURE_INDEX_H
#define FEATURE_INDEX_H

#include "libgtcore/array.h"
#include "libgtcore/str.h"
#include "libgtext/sequence_region.h"
#include "libgtext/genome_feature.h"

typedef struct FeatureIndex FeatureIndex;

FeatureIndex* feature_index_new(Env*);
FeatureIndex* feature_index_ref(FeatureIndex*);
void          feature_index_add_sequence_region(FeatureIndex*, SequenceRegion*,
                                                Env*);
/* add a GenomeFeature to the index, associating it with a sequence region
   denoted by its identifier string */
void          feature_index_add_genome_feature(FeatureIndex*, GenomeFeature*,
                                               Env*);
/* returns an array of GenomeFeatures associated with a given sequence region
   identifier */
Array*        feature_index_get_features_for_seqid(FeatureIndex*, const char*);
/* look up genome features for sequence region <seqid> in <range> and store them
   in <results> */
int           feature_index_get_features_for_range(FeatureIndex*,
                                                   Array *results,
                                                   const char *seqid, Range,
                                                   Env*);
/* returns the first sequence region identifier added to the index */
const char*   feature_index_get_first_seqid(FeatureIndex*);
Range         feature_index_get_range_for_seqid(FeatureIndex*, const char*);
bool          feature_index_has_seqid(const FeatureIndex*, const char*, Env*);
int           feature_index_unit_test(Env*);
void          feature_index_delete(FeatureIndex*, Env*);

#endif
