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
/* GtRBTree.c - red-black tree
 * Copyright (c) 2009 Charles Wilson
 * This file is part of the C Data Structures Library.
 *
 * Derived from code released into the public domain by Julienne Walker:
 * http://eternallyconfuzzled.com/tuts/datastructures/jsw_tut_rbtree.aspx
 *   > Created (Julienne Walker): August 23, 2003
 *   > Modified (Julienne Walker): March 14, 2008
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/ensure.h"
#include "core/error_api.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/rbtree.h"

#ifndef HEIGHT_LIMIT
#define HEIGHT_LIMIT 64
#endif

typedef struct GtRBTreeNode
{
  int red;
  void *key;
  struct GtRBTreeNode *link[2];  /* Left (0) and right (1) links */
} GtRBTreeNode;

struct GtRBTree
{
  GtRBTreeNode *root;
  GtCompareWithData cmp;
  GtRBTreeFreeFunc free;
  size_t size;
  void *info;
};

struct GtRBTreeIter
{
  GtRBTree *tree;
  GtRBTreeNode *it;
  GtRBTreeNode *path[HEIGHT_LIMIT];
  size_t top;
};

#define GT_RBTREE_NODE_IS_RED(node)\
        (node != NULL && node->red == 1)

static inline GtRBTreeNode *gt_rbtree_single(GtRBTreeNode *root, int dir)
{
  GtRBTreeNode *save = root->link[(int) !dir];
  root->link[(int) !dir] = save->link[dir];
  save->link[dir] = root;

  root->red = 1;
  save->red = 0;

  return save;
}

static inline GtRBTreeNode *gt_rbtree_double(GtRBTreeNode *root, int dir)
{
  root->link[(int) !dir] = gt_rbtree_single(root->link[(int) !dir], (int) !dir);
  return gt_rbtree_single (root, dir);
}

static GtRBTreeNode *gt_rbtree_new_node(void *key)
{
  GtRBTreeNode *rn = (GtRBTreeNode *) gt_malloc(sizeof *rn);

  if (rn == NULL)
    return NULL;

  rn->red = 1;
  rn->key = key;
  rn->link[0] = rn->link[1] = NULL;

  return rn;
}

GtRBTree *gt_rbtree_new(GtCompareWithData cmp, GtRBTreeFreeFunc free,
                        void *info)
{
  GtRBTree *rt = (GtRBTree *) gt_malloc(sizeof *rt);
  gt_assert(cmp);

  if (rt == NULL)
    return NULL;

  rt->root = NULL;
  rt->cmp = cmp;
  rt->free = free;
  rt->info = info;
  rt->size = 0;

  return rt;
}

void gt_rbtree_delete(GtRBTree *tree)
{
  gt_rbtree_clear(tree);
  gt_free(tree);
}

void gt_rbtree_clear(GtRBTree *tree)
{
  GtRBTreeNode *it;
  GtRBTreeNode *save;
  if (tree != NULL) {
    it = tree->root;
    if (it != NULL) {

      /* Rotate away the left links so that we can treat this like the
         destruction of a linked list */
      while (it != NULL) {
        if (it->link[0] == NULL) {
          /* No left links, just kill the node and move on */
          save = it->link[1];
          if (tree->free != NULL)
            tree->free(it->key);
          gt_free(it);
        } else {
          /* Rotate away the left link and check again */
          save = it->link[0];
          it->link[0] = save->link[1];
          save->link[1] = it;
        }
        it = save;
      }
      tree->size = 0;
      tree->root = NULL;
    }
  }
}

static inline void *gt_rbtree_find_with_cmp_g(GtRBTree *tree, void *key,
                                              GtCompareWithData cmpfunc,
                                              void *info)
{
  GtRBTreeNode *it = tree->root;

  while (it != NULL) {
    int cmp = cmpfunc(it->key, key, info);

    if (cmp == 0)
      break;

    /* If the tree supports duplicates, they should be
       chained to the right subtree for this to work */
    it = it->link[(int) (cmp < 0)];
  }

  return it == NULL ? NULL : it->key;
}

void *gt_rbtree_find_with_cmp(GtRBTree *tree, void *key,
                              GtCompareWithData cmpfunc, void *info)
{
  gt_assert(tree);
  gt_assert(cmpfunc);
  gt_assert(key);
  return gt_rbtree_find_with_cmp_g(tree, key, cmpfunc, info);
}

void *gt_rbtree_find(GtRBTree *tree, void *key)
{
  gt_assert(tree);
  gt_assert(key);
  return gt_rbtree_find_with_cmp_g(tree, key, tree->cmp, tree->info);
}

static inline int gt_rbtree_insert_g(GtRBTree *tree, void *key,
                                     bool *nodecreated,
                                     GtCompareWithData cmpfunc,
                                     void *info)
{
  *nodecreated = false;
  if (tree->root == NULL) {
    /* We have an empty tree; attach the new node directly to the root */
    tree->root = gt_rbtree_new_node(key);
    *nodecreated = true;

    if (tree->root == NULL)
      return -1;
  } else {
    GtRBTreeNode head = {0,0,{NULL, NULL}}; /* False tree root */
    GtRBTreeNode *g, *t;     /* Grandparent & parent */
    GtRBTreeNode *p, *q;     /* Iterator & parent */
    int dir = 0, last = 0;

    /* Set up our helpers */
    t = &head;
    g = p = NULL;
    q = t->link[1] = tree->root;

    /* Search down the tree for a place to insert */
    for (;;) {
      if (q == NULL) {
        /* Insert a new node at the first null link */
        gt_assert(p != NULL);
        p->link[dir] = q = gt_rbtree_new_node(key);
        *nodecreated = true;

        if (q == NULL)
          return -1;
      } else if (GT_RBTREE_NODE_IS_RED(q->link[0] )
                 && GT_RBTREE_NODE_IS_RED(q->link[1])) {
        /* Simple red violation: color flip */
        q->red = 1;
        q->link[0]->red = 0;
        q->link[1]->red = 0;
      }

      if ( GT_RBTREE_NODE_IS_RED(q) && GT_RBTREE_NODE_IS_RED(p)) {
        /* Hard red violation: rotations necessary */
        int dir2 = (int) (t->link[1] == g);

        if (q == p->link[last])
          t->link[dir2] = gt_rbtree_single(g, (int) !last);
        else
          t->link[dir2] = gt_rbtree_double(g, (int) !last);
      }

      /*
        Stop working if we inserted a node. This
        check also disallows duplicates in the tree
      */
      if (cmpfunc(q->key, key, info) == 0)
        break;

      last = dir;
      dir = (int) (cmpfunc(q->key, key, info) < 0);

      /* Move the helpers down */
      if (g != NULL)
        t = g;

      g = p;
      p = q;
      q = q->link[dir];
    }

    /* Update the root (it may be different) */
    tree->root = head.link[1];
  }

  /* Make the root black for simplified logic */
  tree->root->red = 0;
  ++tree->size;

  return 0;
}

int gt_rbtree_insert(GtRBTree *tree, void *key)
{
  GT_UNUSED bool nodecreated;
  gt_assert(tree);
  gt_assert(key);
  return gt_rbtree_insert_g(tree, key, &nodecreated, tree->cmp, tree->info);
}

int gt_rbtree_insert_with_cmp(GtRBTree *tree, void *key,
                              GtCompareWithData cmpfunc, void *info)
{
  GT_UNUSED bool nodecreated;
  gt_assert(tree);
  gt_assert(key);
  gt_assert(cmpfunc);
  return gt_rbtree_insert_g(tree, key, &nodecreated, cmpfunc, info);
}

void* gt_rbtree_search(GtRBTree *tree, void *key, bool *nodecreated)
{
  gt_assert(tree);
  gt_assert(key);
  gt_assert(nodecreated);
  return gt_rbtree_insert_g(tree,
                            key,
                            nodecreated,
                            tree->cmp,
                            tree->info) == 0
    ? key
    : NULL;
}

void* gt_rbtree_search_with_cmp(GtRBTree *tree, void *key,
                                GtCompareWithData cmpfunc, void *info,
                                bool *nodecreated)
{
  gt_assert(tree);
  gt_assert(nodecreated);
  gt_assert(cmpfunc);
  gt_assert(key);
  return gt_rbtree_insert_g(tree, key, nodecreated, cmpfunc, info) == 0
    ? key
    : NULL;
}

int gt_rbtree_erase(GtRBTree *tree, void *key)
{
  int rv = -1;
  gt_assert(tree);
  gt_assert(key);
  if (tree->root != NULL) {
    GtRBTreeNode head = {0,0,{NULL, NULL}}; /* False tree root */
    GtRBTreeNode *q, *p, *g; /* Helpers */
    GtRBTreeNode *f = NULL;  /* Found item */
    int dir = 1;

    /* Set up our helpers */
    q = &head;
    g = p = NULL;
    q->link[1] = tree->root;

    /* Search and push a red node down to fix red violations as we go */
    while (q->link[dir] != NULL) {
      int last = dir;

      /* Move the helpers down */
      g = p, p = q;
      q = q->link[dir];
      dir = (int) (tree->cmp(q->key, key, tree->info) < 0);

      /* Save the node with matching key and keep going; we'll do removal
         tasks at the end */
      if (tree->cmp(q->key, key, tree->info) == 0)
        f = q;

      /* Push the red node down with rotations and color flips */
      if (!GT_RBTREE_NODE_IS_RED(q)
            && !GT_RBTREE_NODE_IS_RED(q->link[dir])) {
        if (GT_RBTREE_NODE_IS_RED(q->link[(int) !dir])) {
          p = p->link[last] = gt_rbtree_single(q, dir);
        } else if (!GT_RBTREE_NODE_IS_RED(q->link[(int) !dir])) {
          GtRBTreeNode *s = p->link[(int) !last];

          if (s != NULL) {
            if (!GT_RBTREE_NODE_IS_RED(s->link[(int) !last])
                  && !GT_RBTREE_NODE_IS_RED(s->link[last])) {
              /* Color flip */
              p->red = 0;
              s->red = 1;
              q->red = 1;
            } else {
              int dir2;
              gt_assert(g != NULL);
              dir2 = (int) (g->link[1] == p);

              if (GT_RBTREE_NODE_IS_RED(s->link[last]))
                g->link[dir2] = gt_rbtree_double(p, last);
              else if (GT_RBTREE_NODE_IS_RED(s->link[(int) !last]))
                g->link[dir2] = gt_rbtree_single(p, last);

              /* Ensure correct coloring */
              q->red = g->link[dir2]->red = 1;
              g->link[dir2]->link[0]->red = 0;
              g->link[dir2]->link[1]->red = 0;
            }
          }
        }
      }
    }

    /* Replace and remove the saved node */
    if (f != NULL) {
      if (tree->free != NULL)
        tree->free(f->key);
      f->key = q->key;
      gt_assert(p);
      p->link[(int) (p->link[1] == q)] = q->link[(int) (q->link[0] == NULL)];
      gt_free(q);
      rv = 0;
      --tree->size;
    }

    /* Update the root (it may be different) */
    tree->root = head.link[1];

    /* Make the root black for simplified logic */
    if (tree->root != NULL)
      tree->root->red = 0;
  }

  return rv;
}

size_t gt_rbtree_size(GtRBTree *tree)
{
  gt_assert(tree);
  return tree->size;
}

static int gt_rbtree_recurse(GtRBTreeNode *root, GtRBTreeAction action,
                             unsigned long level, void *actinfo)
{
  if (root->link[0] == NULL && root->link[1] == NULL) {
    if (action(root->key, GT_RBTREE_LEAF, level, actinfo) != 0) {
      return -1;
    }
  } else {
    if (action(root->key, GT_RBTREE_PREORDER, level, actinfo) != 0) {
      return -2;
    }
    if (root->link[0] != NULL) {
      if (gt_rbtree_recurse(root->link[0], action, level + 1, actinfo) != 0) {
        return -3;
      }
    }
    if (action(root->key, GT_RBTREE_POSTORDER, level, actinfo) != 0) {
      return -4;
    }
    if (root->link[1] != NULL) {
      if (gt_rbtree_recurse(root->link[1], action, level + 1, actinfo) != 0) {
        return -5;
      }
    }
    if (action(root->key, GT_RBTREE_ENDORDER, level, actinfo) != 0) {
      return -6;
    }
  }
  return 0;
}

#define RBT_CHECK_RETURN_CODE\
        if (retcode < 0 || retcode == 1) {\
          return retcode;\
        }

static int gt_rbtree_recursewithstop(GtRBTreeNode *root, GtRBTreeAction action,
                                     unsigned long level, void *actinfo)
{
  int retcode;

  if (root->link[0] == NULL && root->link[1] == NULL) {
    retcode = action(root->key, GT_RBTREE_LEAF, level, actinfo);
    RBT_CHECK_RETURN_CODE;
  } else {
    retcode = action(root->key, GT_RBTREE_PREORDER, level, actinfo);
    RBT_CHECK_RETURN_CODE;
    if (root->link[0] != NULL) {
      retcode = gt_rbtree_recursewithstop(root->link[0], action,
                                      level + 1, actinfo);
      RBT_CHECK_RETURN_CODE;
    }
    retcode = action(root->key, GT_RBTREE_POSTORDER, level, actinfo);
    RBT_CHECK_RETURN_CODE;
    if (root->link[1] != NULL) {
      retcode = gt_rbtree_recursewithstop(root->link[1], action,
                                      level + 1, actinfo);
      RBT_CHECK_RETURN_CODE;
    }
    retcode = action(root->key, GT_RBTREE_ENDORDER, level, actinfo);
    RBT_CHECK_RETURN_CODE;
  }
  return 0;
}

static int gt_rbtree_recursereverseorder(GtRBTreeNode *root,
                                         GtRBTreeAction action,
                                         unsigned long level, void *actinfo)
{
  if (root->link[0] == NULL && root->link[1] == NULL) {
    if (action(root->key, GT_RBTREE_LEAF, level, actinfo) != 0)
    {
      return -1;
    }
  } else {
    if (action(root->key, GT_RBTREE_PREORDER, level, actinfo) != 0) {
      return -2;
    }
    if (root->link[1] != NULL) {
      if (gt_rbtree_recursereverseorder(root->link[1], action, level + 1,
                                        actinfo) != 0) {
        return -3;
      }
    }
    if (action(root->key, GT_RBTREE_POSTORDER, level, actinfo) != 0) {
      return -4;
    }
    if (root->link[0] != NULL) {
      if (gt_rbtree_recursereverseorder(root->link[0], action, level + 1,
                                        actinfo) != 0) {
        return -5;
      }
    }
    if (action (root->key, GT_RBTREE_ENDORDER, level, actinfo) != 0) {
      return -6;
    }
  }
  return 0;
}

int gt_rbtree_walk(GtRBTree *tree, GtRBTreeAction action, void *actinfo)
{
  gt_assert(tree);
  gt_assert(action);
  if (tree->root != NULL) {
    if (gt_rbtree_recurse(tree->root, action, 0, actinfo) != 0) {
      return -1;
    }
  }
  return 0;
}

int gt_rbtree_walk_stop(GtRBTree *tree, GtRBTreeAction action, void *actinfo)
{
  gt_assert(tree);
  gt_assert(action);
  if (tree->root != NULL) {
    int retcode = gt_rbtree_recursewithstop(tree->root, action, 0, actinfo);
    RBT_CHECK_RETURN_CODE;
  }
  return 0;
}

int gt_rbtree_walk_reverse(GtRBTree *tree, GtRBTreeAction action, void *actinfo)
{
  gt_assert(tree);
  gt_assert(action);
  if (tree->root != NULL) {
    if (gt_rbtree_recursereverseorder(tree->root, action, 0, actinfo) != 0) {
      return -1;
    }
  }
  return 0;
}

static inline void* gt_rbtree_minimum_key_for_node(GtRBTreeNode *root)
{
  if (root == NULL) {
    return NULL;
  }
  while (root->link[0] != NULL) {
    root = root->link[0];
  }
  return root->key;
}

static inline void* gt_rbtree_maximum_key_for_node(GtRBTreeNode *root)
{
  if (root == NULL) {
    return NULL;
  }
  while (root->link[1] != NULL) {
    root = root->link[1];
  }
  return root->key;
}

void* gt_rbtree_minimum_key(GtRBTree *tree)
{
  gt_assert(tree);
  return gt_rbtree_minimum_key_for_node(tree->root);
}

void* gt_rbtree_maximum_key(GtRBTree *tree)
{
  gt_assert(tree);
  return gt_rbtree_maximum_key_for_node(tree->root);
}

void* gt_rbtree_root_key(GtRBTree *tree)
{
  gt_assert(tree);
  if (tree->size == 0 || tree->root == NULL)
    return NULL;
  else
    return tree->root->key;
}

void* gt_rbtree_previous_key(GtRBTree *tree, void *key,
                             GtCompareWithData cmpfun, void *cmpinfo)
{
  int cmp;
  const GtRBTreeNode *current,
                     *found = NULL;
  gt_assert(tree);
  gt_assert(key);
  gt_assert(cmpfun);
  current = tree->root;

  while (current != NULL) {
    cmp = cmpfun(key, current->key, cmpinfo);
    if (cmp == 0) {
      if (current->link[0] == NULL) {
        if (found == NULL) {
          return NULL;
        }
        return found->key;
      }
      return gt_rbtree_maximum_key_for_node(current->link[0]);
    } else {
      if (cmp < 0) {
        current = current->link[0];
      } else {
        found = current;
        current = current->link[1];
      }
    }
  }
  if (found == NULL) {
    return NULL;
  }
  return found->key;
}

void* gt_rbtree_previous_equal_key(GtRBTree *tree, void *key,
                                   GtCompareWithData cmpfun, void *cmpinfo)
{
  int cmp;
  const GtRBTreeNode *current,
                     *found = NULL;
  gt_assert(tree);
  gt_assert(key);
  gt_assert(cmpfun);
  current = tree->root;

  while (current != NULL) {
    cmp = cmpfun(key, current->key, cmpinfo);
    if (cmp == 0) {
      return current->key;
    } else {
      if (cmp < 0) {
        current = current->link[0];
      } else {
        found = current;
        current = current->link[1];
      }
    }
  }
  if (found == NULL) {
    return NULL;
  }
  return found->key;
}

void* gt_rbtree_next_key(GtRBTree *tree, void *key, GtCompareWithData cmpfun,
                         void *cmpinfo)
{
  int cmp;
  const GtRBTreeNode *current,
                     *found = NULL;
  gt_assert(tree);
  gt_assert(key);
  gt_assert(cmpfun);
  current = tree->root;

  while (current != NULL) {
    cmp = cmpfun(key, current->key, cmpinfo);
    if (cmp == 0) {
      if (current->link[1] == NULL) {
        if (found == NULL) {
          return NULL;
        }
        return found->key;
      }
      return gt_rbtree_minimum_key_for_node(current->link[1]);
    } else {
      if (cmp < 0) {
        found = current;
        current = current->link[0];
      } else {
        current = current->link[1];
      }
    }
  }
  if (found == NULL) {
    return NULL;
  }
  return found->key;
}

void* gt_rbtree_next_equal_key(GtRBTree *tree, void *key,
                               GtCompareWithData cmpfun, void *cmpinfo)
{
  int cmp;
  const GtRBTreeNode *current,
                     *found = NULL;
  gt_assert(tree);
  gt_assert(key);
  gt_assert(cmpfun);
  current = tree->root;

  while (current != NULL) {
    cmp = cmpfun(key, current->key, cmpinfo);
    if (cmp == 0) {
      return current->key;
    } else {
      if (cmp < 0) {
        found = current;
        current = current->link[0];
      } else {
        current = current->link[1];
      }
    }
  }
  if (found == NULL) {
    return NULL;
  }
  return found->key;
}

static inline void *start(GtRBTreeIter *trav, GtRBTree *tree, int dir)
{
  trav->tree = tree;
  trav->it = tree->root;
  trav->top = 0;

  /* Save the path for later traversal */
  if (trav->it != NULL)
    {
      while (trav->it->link[dir] != NULL)
        {
          trav->path[trav->top++] = trav->it;
          trav->it = trav->it->link[dir];
        }
    }

  return trav->it == NULL ? NULL : trav->it->key;
}

static inline void *move(GtRBTreeIter *trav, int dir, GT_UNUSED bool strict)
{
  if (trav->it->link[dir] != NULL) {
    /* Continue down this branch */
    trav->path[trav->top++] = trav->it;
    trav->it = trav->it->link[dir];

    while (trav->it->link[(int) !dir] != NULL) {
      trav->path[trav->top++] = trav->it;
      trav->it = trav->it->link[(int) !dir];
    }
  } else {
    /* Move to the next branch */
    GtRBTreeNode *last;

    do {
      if (trav->top == 0) {
        trav->it = NULL;
        break;
      }

      last = trav->it;
      trav->it = trav->path[--trav->top];
    } while (last == trav->it->link[dir]);
  }

  return trav->it == NULL ? NULL : trav->it->key;
}

GtRBTreeIter* gt_rbtree_iter_new_from_first(GtRBTree *tree)
{
  GtRBTreeIter *trav = gt_malloc(sizeof (GtRBTreeIter));
  gt_assert(tree);
  (void) start(trav, tree, 0);
  return trav;
}

GtRBTreeIter* gt_rbtree_iter_new_from_last(GtRBTree *tree)
{
  GtRBTreeIter *trav = gt_malloc(sizeof (GtRBTreeIter));
  gt_assert(tree);
  (void) start(trav, tree, 1);
  return trav;
}

void *gt_rbtree_iter_next(GtRBTreeIter *trav)
{
  gt_assert(trav);
  return move(trav, 1, false); /* Toward larger items */
}

void *gt_rbtree_iter_prev(GtRBTreeIter *trav)
{
  gt_assert(trav);
  return move(trav, 0, false); /* Toward smaller items */
}

void *gt_rbtree_iter_data (GtRBTreeIter *trav)
{
  GtRBTreeNode *node;
  gt_assert(trav);
  node = trav->it;
  return (node != NULL ? node->key : NULL);
}

void gt_rbtree_iter_delete (GtRBTreeIter *trav)
{
  gt_free(trav);
}

/* Unit test for the red black tree datastructure */

#define GT_RBTREE_PASSES 100
#define GT_RBTREE_SIZE 100

/* Shuffle an array of integers.  */
static void gt_rbtree_permuteintarray (unsigned long *arr)
{
  unsigned long i, j, c;

  for (i = 0; i < (unsigned long) GT_RBTREE_SIZE; ++i) {
    j = (unsigned long) (random () % GT_RBTREE_SIZE);
    c = arr[i];
    arr[i] = arr[j];
    arr[j] = c;
  }
}

typedef enum
{
  GT_RBTREE_ASCENDING,
  GT_RBTREE_DESCENDING,
  GT_RBTREE_RANDOMORDER
} GtRBTreeSearchOrder;

typedef enum
{
 GtRBTreeBuild,
  GtRBTreeBuildDelete,
  GtRBTreeDelete,
  GtRBTreeFind
} GtRBTreeDoAction;

/* The keys we add to the tree.  */
static unsigned long *gt_rbtree_xtab;

/*
 * Pointers into the key array, possibly permutated, to define an order for
 * insertion/removal.
 */
static unsigned long *gt_rbtree_ytab;

/* Flags set for each element visited during a tree walk.  */
static unsigned long *gt_rbtree_ztab;

/*
 * Depths for all the elements, to check that the depth is constant for all
 * three visits.
 */
static unsigned long *gt_rbtree_depths;

/* Maximum depth during a tree walk.  */
static unsigned long gt_rbtree_max_depth;

static int nrbt_walk_action(void* nodekey, const GtRBTreeContext which,
                            const unsigned long depth, GT_UNUSED void *actinfo)
{
  unsigned long key = *(unsigned long *) nodekey;

  if (depth > gt_rbtree_max_depth) {
    gt_rbtree_max_depth = depth;
  }
  if (which == GT_RBTREE_LEAF || which == GT_RBTREE_PREORDER) {
    ++gt_rbtree_ztab[key];
    gt_rbtree_depths[key] = depth;
  } else {
    if (gt_rbtree_depths[key] != depth) {
      gt_xfputs("Depth for one element is not constant during tree walk.\n",
                stdout);
      return -1;
    }
  }
  return 0;
}

static int nrbt_cmp_fn(const void* a, const void* b, GT_UNUSED void *info)
{
  unsigned long va, vb;

  va = *(unsigned long *) a;
  vb = *(unsigned long *) b;

  if (va < vb) {
    return -1;
  }
  if (va > vb) {
    return 1;
  }
  return 0;
}

static int nrbt_walk_tree(GtRBTree *tree,unsigned long expected_count)
{
  unsigned long i;
  int error = 0;

  memset(gt_rbtree_ztab, 0, (GT_RBTREE_SIZE * sizeof (*gt_rbtree_ztab)));
  gt_rbtree_max_depth = 0;

  if (gt_rbtree_walk(tree, nrbt_walk_action, NULL) != 0) {
    gt_xfputs("walk failed\n", stdout);
    error = 1;
  }
  for (i = 0; i < expected_count; ++i) {
    if (gt_rbtree_ztab[i] != 1UL) {
      gt_xfputs("Node was not visited.\n", stdout);
      error = 1;
    }
  }
  if (gt_rbtree_max_depth >
      (unsigned long) (log ((double) expected_count) * 2.0 +
                       2.0)) {
    printf("Depth too large during tree walk (%lu vs %lu).\n",
                         gt_rbtree_max_depth,
                         (unsigned long) (log ((double) expected_count) * 2.0 +
                         2.0) );
    error = 1;
  }
  return error;
}

static bool nrbt_mangle_tree(GtRBTreeSearchOrder how, GtRBTreeDoAction what,
                             GtRBTree *tree, unsigned long lag, GtError *err)
{
  unsigned long i;
  bool nodecreated, haserr = false;

  if (how == GT_RBTREE_RANDOMORDER) {
    for (i = 0; i < (unsigned long) GT_RBTREE_SIZE; ++i) {
      gt_rbtree_ytab[i] = i;
    }
    gt_rbtree_permuteintarray (gt_rbtree_ytab);
  }
  for (i = 0; i < GT_RBTREE_SIZE + lag; ++i) {
    void *elem;
    unsigned long j, k;

    switch (how) {
      case GT_RBTREE_RANDOMORDER:
        if (i >= lag) {
          k = gt_rbtree_ytab[i - lag];
        } else {
          k = gt_rbtree_ytab[GT_RBTREE_SIZE - i - 1 + lag];
        }
        j = gt_rbtree_ytab[i];
        break;

      case GT_RBTREE_ASCENDING:
        k = i - lag;
        j = i;
        break;

      case GT_RBTREE_DESCENDING:
        k = GT_RBTREE_SIZE - i - 1 + lag;
        j = GT_RBTREE_SIZE - i - 1;
        break;

      default:
        abort ();
    }

    switch (what) {
      case GtRBTreeBuildDelete:
      case GtRBTreeBuild:
        if (i < (unsigned long) GT_RBTREE_SIZE) {
          if (gt_rbtree_find (tree, gt_rbtree_xtab + j) != NULL) {
            gt_error_set(err,"Found element which is not in tree yet");
            haserr = true;
          }
          if (!haserr) {
            elem = gt_rbtree_search(tree, gt_rbtree_xtab + j, &nodecreated);
            if (elem == NULL || gt_rbtree_find(tree,
                                               gt_rbtree_xtab + j) == NULL) {
              gt_error_set(err,"Couldn't find element after it was added");
              haserr = true;
            }
          }
        }
        if (haserr || what ==GtRBTreeBuild || i < lag)
          break;
        j = k;
        /*@fallthrough@*/
      case GtRBTreeDelete:
        elem = gt_rbtree_find (tree, gt_rbtree_xtab + j);
        if (elem == NULL || gt_rbtree_erase (tree, gt_rbtree_xtab + j) != 0) {
          gt_error_set(err,"Error deleting element");
          haserr = true;
        }
        break;

      case GtRBTreeFind:
        if (gt_rbtree_find (tree, gt_rbtree_xtab + j) == NULL) {
          gt_error_set(err,"Couldn't find element after it was added");
          haserr = true;
        }
        break;
    }
  }
  return haserr;
}

#define NRBT_MANGLECHECK(ORDER,MODE,LAG)\
        gt_ensure(had_err,!nrbt_mangle_tree(ORDER, MODE, tree, LAG,err))

#define NRBT_WALKCHECK\
        gt_ensure(had_err,!nrbt_walk_tree(tree, (unsigned long) GT_RBTREE_SIZE))

int gt_rbtree_unit_test(GtError *err)
{
  int had_err = 0;
  GtRBTree *tree = NULL;
  unsigned long i, j, k, *v;
  gt_error_check (err);

  gt_rbtree_xtab = gt_malloc(GT_RBTREE_SIZE * sizeof (*gt_rbtree_xtab));
  gt_rbtree_ytab = gt_malloc(2*GT_RBTREE_SIZE * sizeof (*gt_rbtree_ytab));
  gt_rbtree_ztab = gt_malloc(GT_RBTREE_SIZE * sizeof (*gt_rbtree_ztab));
  gt_rbtree_depths = gt_malloc(GT_RBTREE_SIZE * sizeof (*gt_rbtree_depths));

  tree = gt_rbtree_new(nrbt_cmp_fn, NULL, NULL);
  gt_ensure(had_err, tree != NULL);
  for (i = 0; i < (unsigned long) GT_RBTREE_SIZE; ++i) {
    gt_rbtree_xtab[i] = i;
  }
  /* Do this loop several times to get different permutations for the random
     case. */
  for (i = 0; i < (unsigned long) GT_RBTREE_PASSES; ++i) {
    NRBT_MANGLECHECK(GT_RBTREE_ASCENDING, GtRBTreeBuild, 0);
    NRBT_MANGLECHECK(GT_RBTREE_ASCENDING, GtRBTreeFind, 0);
    NRBT_MANGLECHECK(GT_RBTREE_DESCENDING, GtRBTreeFind, 0);
    NRBT_MANGLECHECK(GT_RBTREE_RANDOMORDER, GtRBTreeFind, 0);
    NRBT_WALKCHECK;
    NRBT_MANGLECHECK (GT_RBTREE_ASCENDING, GtRBTreeDelete, 0);

    NRBT_MANGLECHECK (GT_RBTREE_ASCENDING, GtRBTreeBuild, 0);
    NRBT_WALKCHECK;
    NRBT_MANGLECHECK (GT_RBTREE_DESCENDING, GtRBTreeDelete, 0);

    NRBT_MANGLECHECK (GT_RBTREE_ASCENDING, GtRBTreeBuild, 0);
    NRBT_WALKCHECK;
    NRBT_MANGLECHECK (GT_RBTREE_RANDOMORDER, GtRBTreeDelete, 0);

    NRBT_MANGLECHECK (GT_RBTREE_DESCENDING, GtRBTreeBuild, 0);
    NRBT_MANGLECHECK (GT_RBTREE_ASCENDING, GtRBTreeFind, 0);
    NRBT_MANGLECHECK (GT_RBTREE_DESCENDING, GtRBTreeFind, 0);
    NRBT_MANGLECHECK (GT_RBTREE_RANDOMORDER, GtRBTreeFind, 0);
    NRBT_WALKCHECK;
    NRBT_MANGLECHECK (GT_RBTREE_DESCENDING, GtRBTreeDelete, 0);

    NRBT_MANGLECHECK (GT_RBTREE_DESCENDING,GtRBTreeBuild, 0);
    NRBT_WALKCHECK;
    NRBT_MANGLECHECK (GT_RBTREE_DESCENDING, GtRBTreeDelete, 0);

    NRBT_MANGLECHECK (GT_RBTREE_DESCENDING, GtRBTreeBuild, 0);
    NRBT_WALKCHECK;
    NRBT_MANGLECHECK (GT_RBTREE_RANDOMORDER, GtRBTreeDelete, 0);

    NRBT_MANGLECHECK (GT_RBTREE_RANDOMORDER, GtRBTreeBuild, 0);
    NRBT_MANGLECHECK (GT_RBTREE_ASCENDING, GtRBTreeFind, 0);
    NRBT_MANGLECHECK (GT_RBTREE_DESCENDING, GtRBTreeFind, 0);
    NRBT_MANGLECHECK (GT_RBTREE_RANDOMORDER, GtRBTreeFind, 0);
    NRBT_WALKCHECK;
    NRBT_MANGLECHECK (GT_RBTREE_RANDOMORDER, GtRBTreeDelete, 0);

    for (j = 1UL; j < (unsigned long) GT_RBTREE_SIZE; j *= 2) {
      NRBT_MANGLECHECK (GT_RBTREE_RANDOMORDER, GtRBTreeBuildDelete, j);
    }
  }
  for (i = 1UL; i < (unsigned long) GT_RBTREE_SIZE; i *= 2) {
    NRBT_MANGLECHECK (GT_RBTREE_ASCENDING, GtRBTreeBuildDelete, i);
    NRBT_MANGLECHECK (GT_RBTREE_DESCENDING, GtRBTreeBuildDelete, i);
    NRBT_MANGLECHECK (GT_RBTREE_ASCENDING, GtRBTreeBuildDelete, i);
    NRBT_MANGLECHECK (GT_RBTREE_DESCENDING, GtRBTreeBuildDelete, i);
    NRBT_MANGLECHECK (GT_RBTREE_ASCENDING, GtRBTreeBuildDelete, i);
    NRBT_MANGLECHECK (GT_RBTREE_DESCENDING, GtRBTreeBuildDelete, i);
    NRBT_MANGLECHECK (GT_RBTREE_ASCENDING, GtRBTreeBuildDelete, i);
    NRBT_MANGLECHECK (GT_RBTREE_DESCENDING, GtRBTreeBuildDelete, i);
  }
  gt_rbtree_delete(tree);
  gt_free(gt_rbtree_xtab);
  gt_free(gt_rbtree_ytab);
  gt_free(gt_rbtree_ztab);
  gt_free(gt_rbtree_depths);

  i = 0; j = 1UL, k = 2UL;
  tree = gt_rbtree_new(nrbt_cmp_fn, NULL, NULL);
  v = gt_rbtree_root_key(tree);
  gt_ensure(had_err, !v);
  (void) gt_rbtree_insert(tree, &i);
  v = gt_rbtree_root_key(tree);
  gt_ensure(had_err, v == &i);
  gt_ensure(had_err, *v == i);
  (void) gt_rbtree_insert(tree, &j);
  (void) gt_rbtree_insert(tree, &k);
  v = gt_rbtree_root_key(tree);
  gt_ensure(had_err, v == &j);
  gt_ensure(had_err, *v == j);
  gt_rbtree_delete(tree);

  return had_err;
}
