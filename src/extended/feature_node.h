/*
  Copyright (c) 2006-2012 Gordon Gremme <gordon@gremme.org>
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
#include "extended/feature_node_observer.h"
#include "extended/feature_type.h"
#include "extended/genome_node.h"
#include "extended/transcript_feature_type.h"

typedef int (*GtFeatureNodeTraverseFunc)(GtFeatureNode*, void*, GtError*);

const GtGenomeNodeClass* gt_feature_node_class(void);

GtFeatureNode* gt_feature_node_clone(const GtFeatureNode*);
void           gt_feature_node_get_exons(GtFeatureNode*,
                                         GtArray *exon_features);
void           gt_feature_node_determine_transcripttypes(GtFeatureNode*);
GtTranscriptFeatureType
               gt_feature_node_get_transcriptfeaturetype(GtFeatureNode*);
void           gt_feature_node_set_end(GtFeatureNode*, GtUword);
bool           gt_feature_node_has_CDS(const GtFeatureNode*);
bool           gt_feature_node_has_splice_site(const GtFeatureNode*);
double         gt_feature_node_average_splice_site_prob(const GtFeatureNode*,
                                                        GtUword *num_ss);
void           gt_feature_node_set_observer(GtFeatureNode*,
                                            GtFeatureNodeObserver*);
void           gt_feature_node_unset_observer(GtFeatureNode*);
int            gt_feature_node_unit_test(GtError*);

/* Perform depth first traversal of the given <feature_node>. */
int            gt_feature_node_traverse_children(GtFeatureNode *feature_node,
                                                 void*,
                                                 GtFeatureNodeTraverseFunc,
                                                 bool traverse_only_once,
                                                 GtError*);
/* Perform topological-sorted depth first traversal of the given <feature_node>.
   Currently, this method can only be used once (for performance reasons)! */
int            gt_feature_node_traverse_children_top(GtFeatureNode
                                                     *feature_node, void*,
                                                     GtFeatureNodeTraverseFunc,
                                                     GtError*);
int            gt_feature_node_traverse_direct_children(GtFeatureNode*, void*,
                                                      GtFeatureNodeTraverseFunc,
                                                        GtError*);

/* Returns <true> if the given <feature_node> graph contains a marked node. */
bool           gt_feature_node_contains_marked(GtFeatureNode *feature_node);
bool           gt_feature_node_has_children(const GtFeatureNode*);
bool           gt_feature_node_direct_children_do_not_overlap(GtFeatureNode*);
/* Returns <true> if all direct childred of <parent> with the same type (s.t.)
   as <child> do not overlap. */
bool           gt_feature_node_direct_children_do_not_overlap_st(GtFeatureNode
                                                                 *parent,
                                                                 GtFeatureNode
                                                                 *child);
bool           gt_feature_node_is_tree(GtFeatureNode*);
/* Returns <true> if the <feature_node> overlaps at least one of the nodes given
   in the <array>. O(<gt_array_size(array)>). */
bool           gt_feature_node_overlaps_nodes(GtFeatureNode *feature_node,
                                              GtArray *array);
/* Similar interface to <gt_genome_node_overlaps_nodes()>. Aditionally, if a
   <bittab> is given (which must have the same size as the <array>), the bits
   corresponding to overlapped nodes are marked (i.e., set). */
bool           gt_feature_node_overlaps_nodes_mark(GtFeatureNode *feature_node,
                                                   GtArray *array,
                                                   GtBittab *bttab);

#endif
