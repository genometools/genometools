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
typedef struct GT_GenomeNodeClass GT_GenomeNodeClass;
typedef struct GT_GenomeNode GT_GenomeNode;

#include "core/bittab.h"
#include "core/phase.h"
#include "core/range.h"
#include "core/str.h"
#include "extended/genome_visitor.h"

typedef int (*GT_GenomeNodeTraverseFunc)(GT_GenomeNode*, void*, GT_Error*);

void           gt_genome_node_set_origin(GT_GenomeNode*,
                                         GT_Str *filename,
                                         unsigned int line_number);
GT_GenomeNode* gt_genome_node_ref(GT_GenomeNode*);
GT_GenomeNode* gt_genome_node_rec_ref(GT_GenomeNode*);
void*          gt_genome_node_cast(const GT_GenomeNodeClass*, GT_GenomeNode*);
/* perform depth first traversal of the given genome node */
int            gt_genome_node_traverse_children(GT_GenomeNode*, void*,
                                                GT_GenomeNodeTraverseFunc,
                                                bool traverse_only_once,
                                                GT_Error*);
/* perform breadth first traversal of the given genome node  */
int            gt_genome_node_traverse_children_breadth(GT_GenomeNode*, void*,
                                                      GT_GenomeNodeTraverseFunc,
                                                        bool traverse_only_once,
                                                        GT_Error*);
int            gt_genome_node_traverse_direct_children(GT_GenomeNode*, void*,
                                                      GT_GenomeNodeTraverseFunc,
                                                       GT_Error*);
const char*    gt_genome_node_get_filename(const GT_GenomeNode*);
unsigned int   gt_genome_node_get_line_number(const GT_GenomeNode*);
unsigned long  gt_genome_node_number_of_children(const GT_GenomeNode*);
GT_Str*        gt_genome_node_get_seqid(GT_GenomeNode*);
/* used to sort nodes */
GT_Str*        gt_genome_node_get_idstr(GT_GenomeNode*);
unsigned long  gt_genome_node_get_start(GT_GenomeNode*);
unsigned long  gt_genome_node_get_end(GT_GenomeNode*);
GT_Range       gt_genome_node_get_range(GT_GenomeNode*);
void           gt_genome_node_set_range(GT_GenomeNode*, GT_Range);
void           gt_genome_node_change_seqid(GT_GenomeNode*, GT_Str*);
int            gt_genome_node_accept(GT_GenomeNode*, GenomeVisitor*, GT_Error*);
/* <parent> takes ownership of <child> */
void           gt_genome_node_is_part_of_genome_node(GT_GenomeNode *parent,
                                                     GT_GenomeNode *child);
/* does not free the leaf, do not use during traversal! */
void           gt_genome_node_remove_leaf(GT_GenomeNode *tree,
                                          GT_GenomeNode *leafn);
void           gt_genome_node_mark(GT_GenomeNode*);
/* returns true if the (top-level) node is marked */
bool           gt_genome_node_is_marked(const GT_GenomeNode*);
/* returns true if the given node graph contains a marked node */
bool           gt_genome_node_contains_marked(GT_GenomeNode*);
bool           gt_genome_node_has_children(GT_GenomeNode*);
bool           gt_genome_node_direct_children_do_not_overlap(GT_GenomeNode*);
/* returns true if all direct childred of <parent> with the same type (s.t.) as
   <child> do not overlap */
bool           gt_genome_node_direct_children_do_not_overlap_st(GT_GenomeNode
                                                                *parent,
                                                                GT_GenomeNode
                                                                *child);
bool           gt_genome_node_is_tree(GT_GenomeNode*);
/* returns true if the genome node overlaps at least one of the nodes given in
   the array. O(gt_array_size) */
bool           gt_genome_node_overlaps_nodes(GT_GenomeNode*, GT_Array*);
/* similar interface to gt_genome_node_overlaps_nodes(). Aditionally, if a
   bittab is given (which must have the same size as the array), the bits
   corresponding to overlapped nodes are marked (i.e., set) */
bool           gt_genome_node_overlaps_nodes_mark(GT_GenomeNode*, GT_Array*,
                                                  GT_Bittab*);
int            gt_genome_node_cmp(GT_GenomeNode*, GT_GenomeNode*);
int            gt_genome_node_compare(GT_GenomeNode**, GT_GenomeNode**);
int            gt_genome_node_compare_with_data(GT_GenomeNode**,
                                                GT_GenomeNode**, void *unused);
/* <delta> has to point to a variable of type unsigned long */
int            gt_genome_node_compare_delta(GT_GenomeNode**, GT_GenomeNode**,
                                            void *delta);
void           gt_genome_node_delete(GT_GenomeNode*);
void           gt_genome_node_rec_delete(GT_GenomeNode*);

void           gt_genome_nodes_sort(GT_Array*);
void           gt_genome_nodes_sort_stable(GT_Array*);
bool           gt_genome_nodes_are_equal_sequence_regions(GT_GenomeNode*,
                                                          GT_GenomeNode*);
bool           gt_genome_nodes_are_sorted(const GT_Array*);

#endif
