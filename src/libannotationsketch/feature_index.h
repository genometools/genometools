/*
  Copyright (c) 2007-2008 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>,
  Copyright (c) 2007      Malte Mader <mmader@stud.zbh.uni-hamburg.de>,
  Copyright (c) 2007      Chr. Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/strarray.h"
#include "libgtext/sequence_region.h"
#include "libgtext/genome_feature.h"

typedef struct FeatureIndex FeatureIndex;

FeatureIndex* feature_index_new(void);
FeatureIndex* feature_index_ref(FeatureIndex*);
void          feature_index_add_sequence_region(FeatureIndex*, SequenceRegion*);
/* Add <genome_feature> to <feature_index>, associating it with a sequence
   region denoted by its identifier string. */
void          feature_index_add_genome_feature(FeatureIndex *feature_index,
                                               GenomeFeature *genome_feature);
/* Returns an array of GenomeFeatures associated with a given sequence region
   identifier <seqid>. */
Array*        feature_index_get_features_for_seqid(FeatureIndex*,
                                                   const char *seqid);
/* Look up genome features in <feature_index> for sequence region <seqid> in
   <range> and store them in <results>. */
int           feature_index_get_features_for_range(FeatureIndex *feature_index,
                                                   Array *results,
                                                   const char *seqid, Range,
                                                   Error*);
/* Returns the first sequence region identifier added to <feature_index>. */
const char*   feature_index_get_first_seqid(const FeatureIndex *feature_index);
/* Returns a StrArray of all sequence region identifiers contained in
   <feature_index> (in alphabetical order). */
StrArray*     feature_index_get_seqids(const FeatureIndex *feature_index);
Range         feature_index_get_range_for_seqid(FeatureIndex*, const char*);
/* Similar to previous function. Necessary for Ruby bindings, because
   apparently 'dl/import' cannot handle returned structs. */
void          feature_index_get_rangeptr_for_seqid(FeatureIndex*, Range*,
                                                   const char *);
bool          feature_index_has_seqid(const FeatureIndex*, const char*);
int           feature_index_unit_test(Error*);
void          feature_index_delete(FeatureIndex*);

#endif
