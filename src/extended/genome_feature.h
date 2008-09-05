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

#ifndef GENOME_FEATURE_H
#define GENOME_FEATURE_H

/* implements the ``genome node'' interface */
typedef struct GenomeFeature GenomeFeature;

#include "core/range.h"
#include "core/phase.h"
#include "core/strand.h"
#include "extended/feature_type_factory.h"
#include "extended/genome_node.h"
#include "extended/genome_feature_type.h"
#include "extended/transcript_feature_type.h"

typedef void (*AttributeIterFunc)(const char *attr_name, const char *attr_value,
                                  void *data);

const GenomeNodeClass* genome_feature_class(void);
GenomeNode*            genome_feature_new(Str *seqid, GenomeFeatureType*, Range,
                                          Strand);
GenomeNode*            genome_feature_new_pseudo(GenomeFeature*);
/* Return the ``standard gene'' (mainly for testing purposes). */
GenomeNode*            genome_feature_new_standard_gene(FeatureTypeFactory*);
const char*            genome_feature_get_source(GenomeFeature*);
const char*            genome_feature_get_attribute(GenomeNode *gn,
                                                    const char *attr_name);
/* Return a GT_StrArray containing the used attribute names. */
GT_StrArray*              genome_feature_get_attribute_list(GenomeFeature*);
GenomeFeatureType*     genome_feature_get_type(GenomeFeature*);
GenomeFeatureType*     genome_feature_create_gft(GenomeFeature*, const char*);
bool                   genome_feature_has_type(GenomeFeature*, const char*);
bool                   genome_feature_score_is_defined(const GenomeFeature*);
bool                   genome_feature_is_multi(const GenomeFeature*);
bool                   genome_feature_is_pseudo(const GenomeFeature*);
void                   genome_feature_make_multi_representative(const
                                                                GenomeFeature*);
void                   genome_feature_set_multi_representative(GenomeFeature*,
                                                               GenomeFeature*);
GenomeFeature*         genome_feature_get_multi_representative(GenomeFeature*);
float                  genome_feature_get_score(GenomeFeature*);
Strand                 genome_feature_get_strand(GenomeFeature*);
Phase                  genome_feature_get_phase(GenomeFeature*);
void                   genome_feature_get_exons(GenomeFeature*,
                                                Array *exon_features);
void                   genome_feature_determine_transcripttypes(GenomeFeature*);
TranscriptFeatureType  genome_feature_get_transcriptfeaturetype(GenomeFeature*);
void                   genome_feature_set_source(GenomeNode*, Str *source);
void                   genome_feature_set_phase(GenomeNode*, Phase);
void                   genome_feature_set_end(GenomeFeature*, unsigned long);
void                   genome_feature_set_score(GenomeFeature*, float);
void                   genome_feature_unset_score(GenomeFeature*);
void                   genome_feature_add_attribute(GenomeFeature*,
                                                    const char *attr_name,
                                                    const char *attr_value);
void                   genome_feature_foreach_attribute(GenomeFeature*,
                                                        AttributeIterFunc,
                                                        void *data);
bool                   genome_feature_has_CDS(const GenomeFeature*);
bool                   genome_feature_has_splice_site(const GenomeFeature*);
double                 genome_feature_average_splice_site_prob(const
                                                               GenomeFeature*);
/* Returns true, if the given features have the same seqid, feature type, range,
   strand, and phase. */
bool                   genome_features_are_similar(GenomeFeature*,
                                                   GenomeFeature*);
int                    genome_feature_unit_test(Error*);

#endif
