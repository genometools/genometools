/*
  Copyright (c)      2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008-2012 Center for Bioinformatics, University of Hamburg

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

#ifndef RBTREE_H
#define RBTREE_H

#include "core/error.h"
#include "core/fptr_api.h"

typedef enum
{
  GT_RBTREE_PREORDER,
  GT_RBTREE_POSTORDER,
  GT_RBTREE_ENDORDER,
  GT_RBTREE_LEAF
} GtRBTreeContext;

typedef struct GtRBTree GtRBTree;
typedef struct GtRBTreeIter GtRBTreeIter;

typedef int   (*GtRBTreeAction)(void *key, GtRBTreeContext, unsigned long,
                                void*);
typedef void  (*GtRBTreeFreeFunc)(void *p);

GtRBTree*      gt_rbtree_new(GtCompareWithData cmp, GtRBTreeFreeFunc free,
                             void *info);
void           gt_rbtree_delete(GtRBTree *tree);
void           gt_rbtree_clear(GtRBTree *tree);
void*          gt_rbtree_find(GtRBTree *tree, void *key);
void*          gt_rbtree_find_with_cmp(GtRBTree *tree, void *key,
                                       GtCompareWithData cmpfunc, void *info);
int            gt_rbtree_insert(GtRBTree *tree, void *key);
int            gt_rbtree_insert_with_cmp(GtRBTree *tree, void *key,
                                         GtCompareWithData cmpfunc,
                                         void *info);
void*          gt_rbtree_search(GtRBTree *tree, void *key, bool *nodecreated);
void*          gt_rbtree_search_with_cmp(GtRBTree *tree, void *key,
                                         GtCompareWithData cmpfunc,
                                         void *info, bool *nodecreated);
int            gt_rbtree_erase(GtRBTree *tree, void *key);
size_t         gt_rbtree_size(GtRBTree *tree);
int            gt_rbtree_walk(GtRBTree *tree, GtRBTreeAction action,
                              void *actinfo);
int            gt_rbtree_walk_stop(GtRBTree *tree, GtRBTreeAction action,
                                   void *actinfo);
int            gt_rbtree_walk_reverse(GtRBTree *tree, GtRBTreeAction action,
                                    void *actinfo);
void*          gt_rbtree_minimum_key(GtRBTree *tree);
void*          gt_rbtree_maximum_key(GtRBTree *tree);
void*          gt_rbtree_root_key(GtRBTree *tree);
void*          gt_rbtree_next_key(GtRBTree *tree, void *key,
                                  GtCompareWithData cmpfun,
                                  void *cmpinfo);
void*          gt_rbtree_next_equal_key(GtRBTree *tree, void *key,
                                        GtCompareWithData cmpfun,
                                        void *cmpinfo);
void*          gt_rbtree_previous_key(GtRBTree *tree, void *key,
                                      GtCompareWithData cmpfun,
                                      void *cmpinfo);
void*          gt_rbtree_previous_equal_key(GtRBTree *tree, void *key,
                                            GtCompareWithData cmpfun,
                                            void *cmpinfo);
int            gt_rbtree_unit_test(GtError *err);

GtRBTreeIter*  gt_rbtree_iter_new_from_first(GtRBTree *tree);
GtRBTreeIter*  gt_rbtree_iter_new_from_last(GtRBTree *tree);
void*          gt_rbtree_iter_next(GtRBTreeIter *trav);
void*          gt_rbtree_iter_prev(GtRBTreeIter *trav);
void           gt_rbtree_iter_delete(GtRBTreeIter *trav);

#endif
