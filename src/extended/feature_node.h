/*
  Copyright (c) 2006-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/bittab.h"
#include "core/range.h"
#include "core/strand_api.h"
#include "core/str_array.h"
#include "extended/feature_node_api.h"
#include "extended/feature_type.h"
#include "extended/genome_node.h"
#include "extended/transcript_feature_type.h"

/* Besides the ``mere'' feature nodes two ``special'' feature nodes exist:
   multi-features and pseudo-features.

   Multi-features represent features which span multiple lines (it is indicated
   in GFF3 files by the fact, that each line has the same ID attribute).

   To check if a feature is a multi-feature use the method
   <gt_feature_node_is_multi()>.
   Multi-features are connected via a ``representative''. That is, two features
   are part of the same multi-feature if they have the same representative
   (which can be retrieved via  <gt_feature_node_get_multi_representative()>).

   Pseudo-features became a technical necessity to be able to pass related
   top-level features as a single entity through the streaming machinery.
   There are two cases in which a pseudo-feature has to be introduced.

   First, if a multi-feature has no parent. In this case all features which
   comprise the multi-feature become the children of a pseudo-feature.

   Second, if two or more top-level features have the same childen (and are
   thereby connected). In this case all these top-level features become the
   children of a pseudo-feature.

   It should be clear from the explanation above that pseudo-features make only
   sense as top-level features (a fact which is enforced in the code).

   Pseudo-features are typically ignored during a traversal to give the illusion
   that they do not exist.
*/

typedef void (*AttributeIterFunc)(const char *attr_name, const char *attr_value,
                                  void *data);
typedef int (*GtGenomeNodeTraverseFunc)(GtGenomeNode*, void*, GtError*);

const GtGenomeNodeClass* gt_feature_node_class(void);

/* Create a new pseudo-<GtFeatureNode*> on sequence with ID <seqid> which lies
   from <start> to <end> on strand <strand>. Pseudo-features do not have a type.
   The <GtFeatureNode*> stores a new reference to <seqid>, so make sure you do
   not modify the original <seqid> afterwards.
   <start> and <end> always refer to the forward strand, therefore <start> has
   to be smaller or equal than <end>. */
GtGenomeNode*  gt_feature_node_new_pseudo(GtStr *seqid, unsigned long start,
                                          unsigned long end, GtStrand strand);
/* Create a new pseudo-<GtFeatureNode*> node which uses <feature_node> as
   template.  That is, the sequence ID, range, strand, and source are taken from
   <feature_node>. */
GtGenomeNode*  gt_feature_node_new_pseudo_template(GtFeatureNode *feature_node);
/* Return <true> if <feature_node> is a multi-feature, <false> otherwise. */
bool           gt_feature_node_is_multi(const GtFeatureNode *feature_node);
/* Return <true> if <feature_node> is a pseudo-feature, <false> otherwise. */
bool           gt_feature_node_is_pseudo(const GtFeatureNode *feature_node);
/* Make <feature_node> the representative of a multi-feature.
   Thereby <feature_node> becomes a multi-feature. */
void           gt_feature_node_make_multi_representative(GtFeatureNode
                                                         *feature_node);
/* Set the multi-feature representative of <feature_node> to <representative>.
   Thereby <feature_node> becomes a multi-feature. */
void           gt_feature_node_set_multi_representative(GtFeatureNode
                                                        *feature_node,
                                                        GtFeatureNode
                                                        *representative);
/* Unset the multi-feature status of <feature_node> and remove its multi-feature
   representative. */
void           gt_feature_node_unset_multi(GtFeatureNode *feature_node);
/* Return the representative of the multi-feature <feature_node>. */
GtFeatureNode* gt_feature_node_get_multi_representative(GtFeatureNode
                                                        *feature_node);
void           gt_feature_node_get_exons(GtFeatureNode*,
                                         GtArray *exon_features);
void           gt_feature_node_determine_transcripttypes(GtFeatureNode*);
GtTranscriptFeatureType
               gt_feature_node_get_transcriptfeaturetype(GtFeatureNode*);
void           gt_feature_node_set_end(GtFeatureNode*, unsigned long);
void           gt_feature_node_foreach_attribute(GtFeatureNode*,
                                                 AttributeIterFunc, void *data);
bool           gt_feature_node_has_CDS(const GtFeatureNode*);
bool           gt_feature_node_has_splice_site(const GtFeatureNode*);
double         gt_feature_node_average_splice_site_prob(const GtFeatureNode*);
/* Returns true, if the given features have the same seqid, feature type, range,
   strand, and phase. */
bool           gt_feature_nodes_are_similar(GtFeatureNode*, GtFeatureNode*);
int            gt_feature_node_unit_test(GtError*);

/* perform depth first traversal of the given genome node */
int            gt_genome_node_traverse_children(GtGenomeNode*, void*,
                                                GtGenomeNodeTraverseFunc,
                                                bool traverse_only_once,
                                                GtError*);
/* perform breadth first traversal of the given genome node  */
int            gt_genome_node_traverse_children_breadth(GtGenomeNode*, void*,
                                                       GtGenomeNodeTraverseFunc,
                                                        bool traverse_only_once,
                                                        GtError*);
int            gt_genome_node_traverse_direct_children(GtGenomeNode*, void*,
                                                       GtGenomeNodeTraverseFunc,
                                                       GtError*);
unsigned long  gt_genome_node_number_of_children(const GtGenomeNode*);
unsigned long  gt_genome_node_number_of_children_of_type(const GtGenomeNode
                                                         *parent,
                                                         const GtGenomeNode
                                                         *node);
/* does not free the leaf, do not use during traversal! */
void           gt_genome_node_remove_leaf(GtGenomeNode *tree,
                                          GtGenomeNode *leafn);
void           gt_genome_node_mark(GtGenomeNode*);
/* returns true if the (top-level) node is marked */
bool           gt_genome_node_is_marked(const GtGenomeNode*);

/* returns true if the given node graph contains a marked node */
bool           gt_genome_node_contains_marked(GtGenomeNode*);
bool           gt_genome_node_has_children(GtGenomeNode*);
bool           gt_genome_node_direct_children_do_not_overlap(GtGenomeNode*);
/* returns true if all direct childred of <parent> with the same type (s.t.) as
   <child> do not overlap */
bool           gt_genome_node_direct_children_do_not_overlap_st(GtGenomeNode
                                                                *parent,
                                                                GtGenomeNode
                                                                *child);
bool           gt_genome_node_is_tree(GtGenomeNode*);
/* returns true if the genome node overlaps at least one of the nodes given in
   the array. O(gt_array_size) */
bool           gt_genome_node_overlaps_nodes(GtGenomeNode*, GtArray*);
/* similar interface to gt_genome_node_overlaps_nodes(). Aditionally, if a
   bittab is given (which must have the same size as the array), the bits
   corresponding to overlapped nodes are marked (i.e., set) */
bool           gt_genome_node_overlaps_nodes_mark(GtGenomeNode*, GtArray*,
                                                  GtBittab*);

#define gt_feature_node_cast(genome_node) \
        gt_genome_node_cast(gt_feature_node_class(), genome_node)

#define gt_feature_node_try_cast(genome_node) \
        gt_genome_node_try_cast(gt_feature_node_class(), genome_node)

#endif
