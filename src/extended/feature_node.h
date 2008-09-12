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

#ifndef FEATURE_NODE_H
#define FEATURE_NODE_H

/* implements the ``genome node'' interface */
typedef struct GtFeatureNode GtFeatureNode;

#include "core/range.h"
#include "core/phase.h"
#include "core/strand.h"
#include "core/strarray.h"
#include "extended/feature_type.h"
#include "extended/genome_node.h"
#include "extended/transcript_feature_type.h"

typedef void (*AttributeIterFunc)(const char *attr_name, const char *attr_value,
                                  void *data);

const GtGenomeNodeClass* gt_feature_node_class(void);
/* Create an new <GtFeatureNode*> on sequence with ID <seqid> and type <type>
   which lies from <start> to <end> on strand <strand>.
   <start> and <end> always refer to the forward strand, therefore <start> has
   to be smaller or equal than <end>. */
GtGenomeNode*        gt_feature_node_new(GtStr *seqid, const char *type,
                                            unsigned long start,
                                            unsigned long end,
                                            GtStrand strand);
GtGenomeNode*        gt_feature_node_new_pseudo(GtFeatureNode*);
/* Return the ``standard gene'' (mainly for testing purposes). */
GtGenomeNode*        gt_feature_node_new_standard_gene(void);
const char*           gt_feature_node_get_source(GtFeatureNode*);
const char*           gt_feature_node_get_attribute(GtGenomeNode *gn,
                                                      const char *attr_name);
/* Return a GtStrArray containing the used attribute names. */
GtStrArray*          gt_feature_node_get_attribute_list(GtFeatureNode*);
const char*           gt_feature_node_get_type(GtFeatureNode*);
bool                  gt_feature_node_has_type(GtFeatureNode*,
                                                 const char*);
bool                  gt_feature_node_score_is_defined(const
                                                         GtFeatureNode*);
bool                  gt_feature_node_is_multi(const GtFeatureNode*);
bool                  gt_feature_node_is_pseudo(const GtFeatureNode*);
void                  gt_feature_node_make_multi_representative(const
                                                             GtFeatureNode*);
void                  gt_feature_node_set_multi_representative(
                                                             GtFeatureNode*,
                                                             GtFeatureNode*);
GtFeatureNode*     gt_feature_node_get_multi_representative(
                                                             GtFeatureNode*);
float                 gt_feature_node_get_score(GtFeatureNode*);
GtStrand             gt_feature_node_get_strand(GtFeatureNode*);
Phase                 gt_feature_node_get_phase(GtFeatureNode*);
void                  gt_feature_node_get_exons(GtFeatureNode*,
                                                  GtArray *exon_features);
void                  gt_feature_node_determine_transcripttypes(
                                                             GtFeatureNode*);
TranscriptFeatureType gt_feature_node_get_transcriptfeaturetype(
                                                             GtFeatureNode*);
void                  gt_feature_node_set_source(GtGenomeNode*,
                                                   GtStr *source);
void                  gt_feature_node_set_phase(GtGenomeNode*, Phase);
void                  gt_feature_node_set_end(GtFeatureNode*,
                                                unsigned long);
void                  gt_feature_node_set_score(GtFeatureNode*, float);
void                  gt_feature_node_unset_score(GtFeatureNode*);
void                  gt_feature_node_add_attribute(GtFeatureNode*,
                                                      const char *attr_name,
                                                      const char *attr_value);
void                  gt_feature_node_foreach_attribute(GtFeatureNode*,
                                                          AttributeIterFunc,
                                                          void *data);
bool                  gt_feature_node_has_CDS(const GtFeatureNode*);
bool                  gt_feature_node_has_splice_site(const
                                                        GtFeatureNode*);
double                gt_feature_node_average_splice_site_prob(const
                                                             GtFeatureNode*);
/* Returns true, if the given features have the same seqid, feature type, range,
   strand, and phase. */
bool                  gt_genome_features_are_similar(GtFeatureNode*,
                                                     GtFeatureNode*);
int                   gt_feature_node_unit_test(GtError*);

#define gt_feature_node_cast(genome_node) \
        gt_genome_node_cast(gt_feature_node_class(), genome_node)

#define gt_feature_node_try_cast(genome_node) \
        gt_genome_node_try_cast(gt_feature_node_class(), genome_node)

#endif
