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

#ifndef GENOME_NODE_H
#define GENOME_NODE_H

/* the ``genome node'' interface */
typedef struct GtGenomeNodeClass GtGenomeNodeClass;
typedef struct GtGenomeNode GtGenomeNode;

#include "core/bittab.h"
#include "core/phase.h"
#include "core/range.h"
#include "core/str.h"
#include "extended/genome_visitor.h"

typedef int (*GtGenomeNodeTraverseFunc)(GtGenomeNode*, void*, GtError*);

void           gt_genome_node_set_origin(GtGenomeNode*,
                                         GtStr *filename,
                                         unsigned int line_number);
GtGenomeNode* gt_genome_node_ref(GtGenomeNode*);
GtGenomeNode* gt_genome_node_rec_ref(GtGenomeNode*);
void*          gt_genome_node_cast(const GtGenomeNodeClass*, GtGenomeNode*);
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
const char*    gt_genome_node_get_filename(const GtGenomeNode*);
unsigned int   gt_genome_node_get_line_number(const GtGenomeNode*);
unsigned long  gt_genome_node_number_of_children(const GtGenomeNode*);
GtStr*        gt_genome_node_get_seqid(GtGenomeNode*);
/* used to sort nodes */
GtStr*        gt_genome_node_get_idstr(GtGenomeNode*);
unsigned long  gt_genome_node_get_start(GtGenomeNode*);
unsigned long  gt_genome_node_get_end(GtGenomeNode*);
GtRange       gt_genome_node_get_range(GtGenomeNode*);
void           gt_genome_node_set_range(GtGenomeNode*, GtRange);
void           gt_genome_node_change_seqid(GtGenomeNode*, GtStr*);
int            gt_genome_node_accept(GtGenomeNode*, GenomeVisitor*, GtError*);
/* Add <child> node to <parent> node. <parent> takes ownership of <child>.*/
void           gt_genome_node_add_child(GtGenomeNode *parent,
                                        GtGenomeNode *child);
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
int            gt_genome_node_cmp(GtGenomeNode*, GtGenomeNode*);
int            gt_genome_node_compare(GtGenomeNode**, GtGenomeNode**);
int            gt_genome_node_compare_with_data(GtGenomeNode**,
                                                GtGenomeNode**, void *unused);
/* <delta> has to point to a variable of type unsigned long */
int            gt_genome_node_compare_delta(GtGenomeNode**, GtGenomeNode**,
                                            void *delta);
void           gt_genome_node_delete(GtGenomeNode*);
void           gt_genome_node_rec_delete(GtGenomeNode*);

void           gt_genome_nodes_sort(GtArray*);
void           gt_genome_nodes_sort_stable(GtArray*);
bool           gt_genome_nodes_are_equal_sequence_regions(GtGenomeNode*,
                                                          GtGenomeNode*);
bool           gt_genome_nodes_are_sorted(const GtArray*);

#endif
