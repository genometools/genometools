/*
  Copyright (c) 2006-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef FEATURE_NODE_API_H
#define FEATURE_NODE_API_H

#include "core/phase_api.h"
#include "core/str_api.h"
#include "core/strand_api.h"
#include "core/str_array_api.h"
#include "extended/genome_node_api.h"

/* Implements the <GtGenomeNode> interface. A single feature node corresponds
   to a regular GFF3 line (i.e., a line which does not start with <#>).
   Part-of relationships (which are realised in GFF3 with the <Parent> and <ID>
   attributes) are realised in the C API with the
   <gt_feature_node_add_child()> method. */
typedef struct GtFeatureNode GtFeatureNode;

/* Create an new <GtFeatureNode*> on sequence with ID <seqid> and type <type>
   which lies from <start> to <end> on strand <strand>.
   The <GtFeatureNode*> stores a new reference to <seqid>, so make sure you do
   not modify the original <seqid> afterwards.
   <start> and <end> always refer to the forward strand, therefore <start> has
   to be smaller or equal than <end>. */
GtGenomeNode* gt_feature_node_new(GtStr *seqid, const char *type,
                                  unsigned long start, unsigned long end,
                                  GtStrand strand);
/* Add <child> node to <parent> node. <parent> takes ownership of <child>.*/
void          gt_feature_node_add_child(GtFeatureNode *parent,
                                        GtFeatureNode *child);

/* Return the source of <feature_node>. If no source has been set, "." is
   returned. Corresponds to column 2 of regular GFF3 lines. */
const char*   gt_feature_node_get_source(const GtFeatureNode *feature_node);

/* Set the <source> of <feature_node>. Stores a new reference to <source>.
   Corresponds to column 2 of regular GFF3 lines. */
void          gt_feature_node_set_source(GtFeatureNode *feature_node,
                                         GtStr *source);

/* Return <true> if <feature_node> has a defined source (i.e., on different
   from "."). <false> otherwise. */
bool          gt_feature_node_has_source(const GtFeatureNode *feature_node);

/* Return the type of <feature_node>.
   Corresponds to column 3 of regular GFF3 lines. */
const char*   gt_feature_node_get_type(const GtFeatureNode *feature_node);

/* Set the type of <feature_node> to <type>. */
void          gt_feature_node_set_type(GtFeatureNode *feature_node,
                                       const char *type);

/* Return <true> if <feature_node> has given <type>, <false> otherwise. */
bool          gt_feature_node_has_type(GtFeatureNode *feature_node,
                                       const char *type);

/* Return <true> if the score of <feature_node> is defined, <false>
   otherwise. */
bool          gt_feature_node_score_is_defined(const GtFeatureNode
                                               *feature_node);
/* Return the score of <feature_node>. The score has to be defined.
   Corresponds to column 6 of regular GFF3 lines. */
float         gt_feature_node_get_score(const GtFeatureNode *feature_node);

/* Set the score of <feature_node> to <score>. */
void          gt_feature_node_set_score(GtFeatureNode *feature_node,
                                        float score);

/* Unset the score of <feature_node>. */
void          gt_feature_node_unset_score(GtFeatureNode *feature_node);

/* Return the strand of <feature_node>.
   Corresponds to column 7 of regular GFF3 lines. */
GtStrand      gt_feature_node_get_strand(const GtFeatureNode *feature_node);

/* Set the strand of <feature_node> to <strand>. */
void          gt_feature_node_set_strand(GtFeatureNode *feature_node,
                                         GtStrand strand);

/* Return the phase of <feature_node>.
   Corresponds to column 8 of regular GFF3 lines. */
GtPhase       gt_feature_node_get_phase(const GtFeatureNode *feature_node);

/* Set the phase of <feature_node> to <phase>. */
void          gt_feature_node_set_phase(GtFeatureNode *feature_node,
                                        GtPhase phase);

/* Return the attribute of <feature_node> with the given <name>.
   If no such attribute has been added, <NULL> is returned.
   The attributes are stored in column 9 of regular GFF3 lines. */
const char*   gt_feature_node_get_attribute(const GtFeatureNode *feature_node,
                                            const char *name);

/* Return a string array containing the used attribute names of <feature_node>.
   The caller is responsible to free the returned <GtStrArray*>. */
GtStrArray*   gt_feature_node_get_attribute_list(const GtFeatureNode
                                                 *feature_node);

/* Add attribute <tag>=<value> to <feature_node>. <tag> and <value> must at
   least have length 1. <feature_node> must not contain an attribute with the
   given <tag> already. You should not add Parent and ID attributes, use
   <gt_feature_node_add_child()> to denote part-of relationships. */
void          gt_feature_node_add_attribute(GtFeatureNode *feature_node,
                                            const char *tag, const char *value);

#endif
