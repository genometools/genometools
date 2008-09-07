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
#include "core/fptr.h"

/* Interval tree data structure, implemented according to
   Cormen et al., Introduction to Algorithms, 2nd edition, MIT Press,
   Cambridge, MA, USA, 2001 */

typedef struct IntervalTree IntervalTree;
typedef struct IntervalTreeNode IntervalTreeNode;

typedef int (*IntervalTreeIteratorFunc)(IntervalTreeNode*, void*);

/* transfers ownership of <data> to interval tree
   if IntervalTreeDataFreeFunc is given */
IntervalTreeNode* interval_tree_node_new(void *data,
                                         unsigned long low,
                                         unsigned long high);
void*             interval_tree_node_get_data(IntervalTreeNode* n);

IntervalTree*     interval_tree_new(FreeFunc);
unsigned long     interval_tree_size(IntervalTree*);
IntervalTreeNode* interval_tree_find_first_overlapping(IntervalTree*,
                                                       unsigned long start,
                                                       unsigned long end);
void              interval_tree_insert(IntervalTree*, IntervalTreeNode*);
/* collects data pointers of all overlapping nodes in array */
void              interval_tree_find_all_overlapping(IntervalTree*,
                                                     unsigned long start,
                                                     unsigned long end,
                                                     GT_Array*);
int               interval_tree_traverse(IntervalTree *it,
                                         IntervalTreeIteratorFunc func,
                                         void *data);
void              interval_tree_delete(IntervalTree*);
int               interval_tree_unit_test(GT_Error*);

#endif
