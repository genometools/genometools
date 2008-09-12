/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef INTERVAL_TREE_H
#define INTERVAL_TREE_H

#include "core/array.h"
#include "core/error.h"
#include "core/fptr_api.h"

/* Interval tree data structure, implemented according to
   Cormen et al., Introduction to Algorithms, 2nd edition, MIT Press,
   Cambridge, MA, USA, 2001 */

typedef struct GT_IntervalTree GT_IntervalTree;
typedef struct GT_IntervalTreeNode GT_IntervalTreeNode;

typedef int (*GT_IntervalTreeIteratorFunc)(GT_IntervalTreeNode*, void*);

/* transfers ownership of <data> to interval tree
   if GT_IntervalTreeDataFreeFunc is given */
GT_IntervalTreeNode* gt_interval_tree_node_new(void *data,
                                               unsigned long low,
                                               unsigned long high);
void*                gt_interval_tree_node_get_data(GT_IntervalTreeNode* n);

GT_IntervalTree*     gt_interval_tree_new(GT_FreeFunc);
unsigned long        gt_interval_tree_size(GT_IntervalTree*);
GT_IntervalTreeNode* gt_interval_tree_find_first_overlapping(GT_IntervalTree*,
                                                            unsigned long start,
                                                            unsigned long end);
void                 gt_interval_tree_insert(GT_IntervalTree*,
                                             GT_IntervalTreeNode*);
/* collects data pointers of all overlapping nodes in array */
void              gt_interval_tree_find_all_overlapping(GT_IntervalTree*,
                                                        unsigned long start,
                                                        unsigned long end,
                                                        GtArray*);
int               gt_interval_tree_traverse(GT_IntervalTree *it,
                                            GT_IntervalTreeIteratorFunc func,
                                            void *data);
void              gt_interval_tree_delete(GT_IntervalTree*);
int               gt_interval_tree_unit_test(GtError*);

#endif
