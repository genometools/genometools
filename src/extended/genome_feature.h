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
typedef struct GtGenomeFeature GtGenomeFeature;

#include "core/range.h"
#include "core/phase.h"
#include "core/strand.h"
#include "core/strarray.h"
#include "extended/feature_type.h"
#include "extended/genome_node.h"
#include "extended/transcript_feature_type.h"

typedef void (*AttributeIterFunc)(const char *attr_name, const char *attr_value,
                                  void *data);

const GtGenomeNodeClass* gt_genome_feature_class(void);
/* Create an new <GtGenomeFeature*> on sequence with ID <seqid> and type <type>
   which lies from <start> to <end> on strand <strand>.
   <start> and <end> always refer to the forward strand, therefore <start> has
   to be smaller or equal than <end>. */
GtGenomeNode*        gt_genome_feature_new(GtStr *seqid, const char *type,
                                            unsigned long start,
                                            unsigned long end,
                                            GtStrand strand);
GtGenomeNode*        gt_genome_feature_new_pseudo(GtGenomeFeature*);
/* Return the ``standard gene'' (mainly for testing purposes). */
GtGenomeNode*        gt_genome_feature_new_standard_gene(void);
const char*           gt_genome_feature_get_source(GtGenomeFeature*);
const char*           gt_genome_feature_get_attribute(GtGenomeNode *gn,
                                                      const char *attr_name);
/* Return a GtStrArray containing the used attribute names. */
GtStrArray*          gt_genome_feature_get_attribute_list(GtGenomeFeature*);
const char*           gt_genome_feature_get_type(GtGenomeFeature*);
bool                  gt_genome_feature_has_type(GtGenomeFeature*,
                                                 const char*);
bool                  gt_genome_feature_score_is_defined(const
                                                         GtGenomeFeature*);
bool                  gt_genome_feature_is_multi(const GtGenomeFeature*);
bool                  gt_genome_feature_is_pseudo(const GtGenomeFeature*);
void                  gt_genome_feature_make_multi_representative(const
                                                             GtGenomeFeature*);
void                  gt_genome_feature_set_multi_representative(
                                                             GtGenomeFeature*,
                                                             GtGenomeFeature*);
GtGenomeFeature*     gt_genome_feature_get_multi_representative(
                                                             GtGenomeFeature*);
float                 gt_genome_feature_get_score(GtGenomeFeature*);
GtStrand             gt_genome_feature_get_strand(GtGenomeFeature*);
Phase                 gt_genome_feature_get_phase(GtGenomeFeature*);
void                  gt_genome_feature_get_exons(GtGenomeFeature*,
                                                  GtArray *exon_features);
void                  gt_genome_feature_determine_transcripttypes(
                                                             GtGenomeFeature*);
TranscriptFeatureType gt_genome_feature_get_transcriptfeaturetype(
                                                             GtGenomeFeature*);
void                  gt_genome_feature_set_source(GtGenomeNode*,
                                                   GtStr *source);
void                  gt_genome_feature_set_phase(GtGenomeNode*, Phase);
void                  gt_genome_feature_set_end(GtGenomeFeature*,
                                                unsigned long);
void                  gt_genome_feature_set_score(GtGenomeFeature*, float);
void                  gt_genome_feature_unset_score(GtGenomeFeature*);
void                  gt_genome_feature_add_attribute(GtGenomeFeature*,
                                                      const char *attr_name,
                                                      const char *attr_value);
void                  gt_genome_feature_foreach_attribute(GtGenomeFeature*,
                                                          AttributeIterFunc,
                                                          void *data);
bool                  gt_genome_feature_has_CDS(const GtGenomeFeature*);
bool                  gt_genome_feature_has_splice_site(const
                                                        GtGenomeFeature*);
double                gt_genome_feature_average_splice_site_prob(const
                                                             GtGenomeFeature*);
/* Returns true, if the given features have the same seqid, feature type, range,
   strand, and phase. */
bool                  gt_genome_features_are_similar(GtGenomeFeature*,
                                                     GtGenomeFeature*);
int                   gt_genome_feature_unit_test(GtError*);

#define gt_genome_feature_cast(genome_node) \
        gt_genome_node_cast(gt_genome_feature_class(), genome_node)

#define gt_genome_feature_try_cast(genome_node) \
        gt_genome_node_try_cast(gt_genome_feature_class(), genome_node)

#endif
