/*
  Copyright (c) 2008-2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#include <limits.h>
#include <string.h>
#include "core/ensure.h"
#include "core/interval_tree.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/range.h"
#include "core/str.h"
#include "core/unused_api.h"

typedef enum GtIntervalTreeNodeColor {
  BLACK,
  RED
} GtIntervalTreeNodeColor ;

struct GtIntervalTreeNode {
  GtIntervalTreeNode *parent, *left, *right;
  void *data;
  GtIntervalTreeNodeColor color;
  unsigned long low, high, max;
};

struct GtIntervalTree {
  GtIntervalTreeNode *root, sentinel, *nil;
  unsigned long size;
  GtFree free_func;
};

GtIntervalTreeNode* gt_interval_tree_node_new(void *data,
                                         unsigned long low,
                                         unsigned long high)
{
  GtIntervalTreeNode* n;
  n = gt_calloc(1, sizeof (GtIntervalTreeNode));
  n->low = low;
  n->high = high;
  n->data = data;
  return n;
}

GtIntervalTree* gt_interval_tree_new(GtFree func)
{
  GtIntervalTree *it;
  it = gt_calloc(1, sizeof (GtIntervalTree));
  it->free_func = func;
  it->nil = &it->sentinel;
  it->root = it->nil;
  return it;
}

unsigned long gt_interval_tree_size(GtIntervalTree *it)
{
  gt_assert(it);
  return it->size;
}

void* gt_interval_tree_node_get_data(GtIntervalTreeNode *n)
{
  gt_assert(n);
  return n->data;
}

void gt_interval_tree_node_delete(GtIntervalTree *it, GtIntervalTreeNode *n)
{
  if (n == it->nil) return;
  if (n->data && it->free_func)
    it->free_func(n->data);
  gt_free(n);
}

static void interval_tree_node_rec_delete(GtIntervalTree *it,
                                          GtIntervalTreeNode *n)
{
  if (n == it->nil) return;
  interval_tree_node_rec_delete(it, n->left);
  interval_tree_node_rec_delete(it, n->right);
  gt_interval_tree_node_delete(it, n);
}

static GtIntervalTreeNode* interval_tree_search_internal(GtIntervalTree *it,
                                                         GtIntervalTreeNode
                                                          *node,
                                                         unsigned long low,
                                                         unsigned long high)
{
  GtIntervalTreeNode *x;
  x = node;

  while (x != it->nil && !(low <= x->high && x->low <= high)) {
    if (x->left != it->nil && x->left->max >= low)
      x = x->left;
    else
      x = x->right;
  }
  return (x == it->nil) ? NULL : x;
}

GtIntervalTreeNode* gt_interval_tree_find_first_overlapping(GtIntervalTree *it,
                                                            unsigned long low,
                                                            unsigned long high)
{
  gt_assert(it);
  if (it->root == it->nil)
    return NULL;
  return interval_tree_search_internal(it, it->root, low, high);
}

static int interval_tree_traverse_internal(GtIntervalTree *it,
                                           GtIntervalTreeNode *node,
                                           GtIntervalTreeIteratorFunc func,
                                           void *data)
{
  int had_err = 0;
  if (node == it->nil) return 0;
  if (!had_err)
    had_err = interval_tree_traverse_internal(it, node->left, func, data);
  if (!had_err)
    had_err = interval_tree_traverse_internal(it, node->right, func, data);
  if (!had_err)
    had_err = func(node, data);
  return had_err;
}

int gt_interval_tree_traverse(GtIntervalTree *it,
                              GtIntervalTreeIteratorFunc func, void *data)
{
  if (it->root == it->nil)
    return 0;
  return interval_tree_traverse_internal(it, it->root, func, data);
}

static int store_interval_node_in_array(GtIntervalTreeNode *x, void *data)
{
  GtArray *a = (GtArray*) data;
  gt_array_add(a, x->data);
  return 0;
}

static void interval_tree_find_all_internal(GtIntervalTree *it,
                                            GtIntervalTreeNode *node,
                                            GtIntervalTreeIteratorFunc func,
                                            unsigned long low,
                                            unsigned long high,
                                            void *data)
{
  GtIntervalTreeNode* x;
  if (node == it->nil) return;
  x = node;
  if (low <= x->high && x->low <= high)
    func(node, data);
  /* recursively search left and right subtrees */
  if (x->left != it->nil && low <= x->left->max)
    interval_tree_find_all_internal(it, x->left, func, low, high, data);
  if (x->right != it->nil && low <= x->right->max)
    interval_tree_find_all_internal(it, x->right, func, low, high, data);
}

void gt_interval_tree_find_all_overlapping(GtIntervalTree *it,
                                           unsigned long start,
                                           unsigned long end, GtArray* a)
{
  gt_assert(it && a && start <= end);
  if (it->root == it->nil) return;
  interval_tree_find_all_internal(it, it->root, store_interval_node_in_array,
                                  start, end, a);
}

void gt_interval_tree_iterate_overlapping(GtIntervalTree *it,
                                          GtIntervalTreeIteratorFunc func,
                                          unsigned long start,
                                          unsigned long end,
                                          void *data)
{
  gt_assert(it && func && start <= end);
  interval_tree_find_all_internal(it, it->root, func, start, end, data);
}

static void interval_tree_left_rotate(GtIntervalTree *it,
                                      GtIntervalTreeNode **root,
                                      GtIntervalTreeNode *x)
{
  GtIntervalTreeNode *y;
  y = x->right;
  x->right = y->left;
  if (y->left != it->nil)
    y->left->parent = x;
  y->parent = x->parent;
  if (x->parent == it->nil)
    *root = y;
  else {
    if (x == x->parent->left)
      x->parent->left = y;
    else
      x->parent->right = y;
  }
  y->left = x;
  x->parent = y;
  /* interval tree augmentation */
  x->max = x->high;
  if (x->left != it->nil && x->left->max > x->max)
    x->max = x->left->max;
  if (x->right != it->nil && x->right->max > x->max)
    x->max = x->right->max;
  y->max = y->high;
  if (y->left != it->nil && y->left->max > y->max)
    y->max = y->left->max;
  if (y->right != it->nil && y->right->max > y->max)
    y->max = y->right->max;
}

static void interval_tree_right_rotate(GtIntervalTree *it,
                                       GtIntervalTreeNode **root,
                                       GtIntervalTreeNode *y)
{
  GtIntervalTreeNode *x;
  x = y->left;
  y->left = x->right;
  if (x->right != it->nil)
    x->right->parent = y;
  x->parent = y->parent;
  if (y->parent == it->nil)
    *root = x;
  else {
    if (y == y->parent->left)
      y->parent->left = x;
    else
      y->parent->right = x;
  }
  x->right = y;
  y->parent = x;
  /* interval tree augmentation */
  x->max = x->high;
  if (x->left != it->nil && x->left->max > x->max)
    x->max = x->left->max;
  if (x->right != it->nil && x->right->max > x->max)
    x->max = x->right->max;
  y->max = y->high;
  if (y->left != it->nil && y->left->max > y->max)
    y->max = y->left->max;
  if (y->right != it->nil && y->right->max > y->max)
    y->max = y->right->max;
}

static void interval_tree_max_fixup(GtIntervalTree *it,
                                    GtIntervalTreeNode *root,
                                    GtIntervalTreeNode *x)
{
  while (x != it->nil && x != root) {
    if (x->left == it->nil && x->right != it->nil)
      x->max = x->right->max;
    else if (x->right == it->nil && x->left != it->nil)
      x->max = x->left->max;
    else if (x->right != it->nil && x->left != it->nil)
      x->max = MAX(x->left->max, x->right->max);
    x = x->parent;
  }
}

/* this is the insert routine from Cormen et al, p. 280*/
static void interval_tree_insert(GtIntervalTree *it,
                                 GtIntervalTreeNode **root,
                                 GtIntervalTreeNode *z)
{
  GtIntervalTreeNode *x, *y;
  y = it->nil;
  x = *root;
  z->max = z->high;
  while (x != it->nil)
  {
    y = x;
    /* interval tree augmentation */
    if (x->max < z->max)
      x->max = z->max;
    if (z->low < x->low)
      x = x->left;
    else
      x = x->right;
  }
  z->parent = y;
  if (y == it->nil)
    *root = z;
  else
  {
    if (z->low < y->low)
      y->left = z;
    else
      y->right = z;
  }
}

/* this is the fixup routine from Cormen et al, p. 281*/
static void interval_tree_insert_internal(GtIntervalTree *it,
                                          GtIntervalTreeNode **root,
                                          GtIntervalTreeNode *z)
{
  GtIntervalTreeNode* y;
  interval_tree_insert(it, root, z);
  z->color = RED;
  while (z != *root && z->parent->color == RED)
  {
    if (z->parent == z->parent->parent->left)
    {
      y = z->parent->parent->right;
      if (y != it->nil && y->color == RED)
      {
        z->parent->color = BLACK;
        y->color = BLACK;
        z->parent->parent->color = RED;
        z = z->parent->parent;
      }
      else
      {
        if (z == z->parent->right)
        {
          z = z->parent;
          interval_tree_left_rotate(it, root, z);
        }
        z->parent->color = BLACK;
        z->parent->parent->color = RED;
        interval_tree_right_rotate(it, root, z->parent->parent);
      }
    }
    else
    {
      y = z->parent->parent->left;
      if (y != it->nil && y->color == RED)
      {
        z->parent->color = BLACK;
        y->color = BLACK;
        z->parent->parent->color = RED;
        z = z->parent->parent;
      }
      else
      {
        if (z == z->parent->left)
        {
          z = z->parent;
          interval_tree_right_rotate(it, root, z);
        }
        z->parent->color = BLACK;
        z->parent->parent->color = RED;
        interval_tree_left_rotate(it, root, z->parent->parent);
      }
    }
  }
  (*root)->color = BLACK;
}

void gt_interval_tree_insert(GtIntervalTree *it, GtIntervalTreeNode *n)
{
  gt_assert(it && n);
  n->parent = it->nil;
  n->left = it->nil;
  n->right = it->nil;
  if (it->root == it->nil)
  {
    it->root = n;
  } else interval_tree_insert_internal(it, &(it->root), n);
  it->size++;
}

void gt_interval_tree_delete(GtIntervalTree *it)
{
  if (!it) return;
  interval_tree_node_rec_delete(it, it->root);
  gt_free(it);
}

GtIntervalTreeNode* gt_interval_tree_get_successor(GtIntervalTree *it,
                                                   GtIntervalTreeNode *x)
{
  GtIntervalTreeNode *y;

  if ((y = x->right)) {
    while (y->left != it->nil) {
      y = y->left;
    }
    return y;
  } else {
    y = x->parent;
    while (y != it->nil && x == y->right) {
      x = y;
      y = y->parent;
    }
    return y;
  }
}

static inline void interval_tree_delete_fixup(GtIntervalTree *it,
                                              GtIntervalTreeNode *x)
{
  GtIntervalTreeNode *w;

  while ((x->color == BLACK) && (it->root != x)) {
    if (x == x->parent->left) {
      w = x->parent->right;
      if (w->color == RED) {
        w->color = BLACK;
        x->parent->color = RED;
        interval_tree_left_rotate(it, &it->root, x->parent);
        w = x->parent->right;
      }
      if ( (w->right->color == BLACK) && (w->left->color == BLACK) ) {
        w->color = RED;
        x = x->parent;
      } else {
        if (w->right->color == BLACK) {
          w->left->color = BLACK;
          w->color = RED;
          interval_tree_right_rotate(it, &it->root, w);
          w = x->parent->right;
        }
        w->color = x->parent->color;
        x->parent->color = BLACK;
        w->right->color = BLACK;
        interval_tree_left_rotate(it, &it->root, x->parent);
        x = it->root;
      }
    } else {
      w = x->parent->left;
      if (w->color == RED) {
        w->color = BLACK;
        x->parent->color = RED;
        interval_tree_right_rotate(it, &it->root, x->parent);
        w=x->parent->left;
      }
      if ( (w->right->color == BLACK) && (w->left->color == BLACK) ) {
        w->color = RED;
        x = x->parent;
      } else {
        if (w->left->color == BLACK) {
          w->right->color = BLACK;
          w->color = RED;
          interval_tree_left_rotate(it, &it->root, w);
          w=x->parent->left;
        }
        w->color = x->parent->color;
        x->parent->color = BLACK;
        w->left->color = BLACK;
        interval_tree_right_rotate(it, &it->root, x->parent);
        x = it->root;
      }
    }
  }
  x->color = BLACK;
}

void gt_interval_tree_remove(GtIntervalTree *it, GtIntervalTreeNode *z)
{
  GtIntervalTreeNode *y, *x;
  gt_assert(it && it->size > 0);
  y = (z->left == it->nil || z->right == it->nil)
    ? z
    : gt_interval_tree_get_successor(it, z);
  x = (y->left != it->nil) ? y->left : y->right;
  gt_assert(y);

  x->parent = y->parent;

  if (y->parent == it->nil) {
    it->root = x;
  } else {
    if (y == y->parent->left) {
      y->parent->left = x;
    } else {
      y->parent->right = x;
    }
  }

  if (y != z) {
    z->max = y->max;
    z->low = y->low;
    z->high = y->high;
    z->data = y->data;
  }
  interval_tree_max_fixup(it, it->root, z->parent);
  if (y->color == BLACK) {
    y->color = z->color;
    interval_tree_delete_fixup(it, x);
  }
  if (y != it->nil)
    gt_interval_tree_node_delete(it, y);
  it->size--;
}

static void gt_interval_tree_print_rec(GtIntervalTree *it,
                                       GtIntervalTreeNode *n)
{
  if (n == it->nil) return;
  printf("(");
  gt_interval_tree_print_rec(it, n->left);
  printf("[%lu,%lu]", n->low, n->high);
  gt_interval_tree_print_rec(it, n->right);
  printf(")");
}

void gt_interval_tree_print(GtIntervalTree *it)
{
  gt_assert(it);
  gt_interval_tree_print_rec(it, it->root);
}

static int range_ptr_compare(const void *r1p, const void *r2p)
{
  int ret;
  gt_assert(r1p && r2p);
  ret = gt_range_compare(*(GtRange**) r1p, *(GtRange**) r2p);
  /* It could be that two identical ranges with different pointers are
     present. If so, compare pointers instead to get a canonical ordering. */
  if (ret == 0 && *(GtRange**) r1p != *(GtRange**) r2p)
  {
    if (*(GtRange**) r1p < *(GtRange**) r2p)
      ret = -1;
    else
      ret = 1;
  }
  return ret;
}

static int itree_test_get_node(GtIntervalTreeNode *x, void *data)
{
  GtArray *a = (GtArray*) data;
  gt_array_add(a, x);
  return 0;
}

int gt_interval_tree_unit_test(GT_UNUSED GtError *err)
{
  GtIntervalTree *it = NULL;
  GtIntervalTreeNode *res = NULL;
  unsigned long i = 0;
  int had_err = 0, num_testranges = 3000,
      num_samples = 300000, num_find_all_samples = 10000,
      gt_range_max_basepos = 90000, width = 700,
      query_width = 5000;
  GtRange *res_rng = NULL, qrange;
  GtArray *arr = NULL, *narr = NULL;

  arr = gt_array_new(sizeof (GtRange*));

  /* generate test ranges */
  for (i = 0;i<num_testranges;i++)
  {
    unsigned long start;
    GtRange *rng;
    rng  = gt_calloc(1, sizeof (GtRange));
    start = gt_rand_max(gt_range_max_basepos);
    rng->start = start;
    rng->end = start + gt_rand_max(width);
    gt_array_add(arr, rng);
  }

  it = gt_interval_tree_new(gt_free_func);

  /* insert ranges */
  for (i = 0; i < num_testranges && !had_err; i++)
  {
    GtIntervalTreeNode *new_node;
    GtRange *rng;
    rng = *(GtRange**) gt_array_get(arr, i);
    new_node = gt_interval_tree_node_new(rng, rng->start, rng->end);
    gt_interval_tree_insert(it, new_node);
  }
  gt_ensure(had_err, gt_interval_tree_size(it) == num_testranges);

  /* perform test queries */
  for (i = 0; i < num_samples && !had_err; i++)
  {
    unsigned long start = gt_rand_max(gt_range_max_basepos);
    qrange.start = start;
    qrange.end = start + gt_rand_max(width);
    res = gt_interval_tree_find_first_overlapping(it, qrange.start, qrange.end);
    if (res)
    {
      /* we have a hit, check if really overlapping */
      res_rng = (GtRange*) gt_interval_tree_node_get_data(res);
      gt_ensure(had_err, gt_range_overlap(&qrange, res_rng));
    } else {
      /* no hit, check whether there really is no overlapping
         interval in tree */
      GtRange *this_rng;
      unsigned long j;
      bool found = false;
      for (j = 0; j < gt_array_size(arr); j++)
      {
        this_rng = *(GtRange**) gt_array_get(arr, j);
        if (gt_range_overlap(this_rng, &qrange))
        {
          found = true;
          break;
        }
      }
      gt_ensure(had_err, !found);
    }
  }

  /* test searching for all overlapping intervals */
  for (i = 0; i < num_find_all_samples && !had_err; i++)
  {
    unsigned long start = gt_rand_max(gt_range_max_basepos);
    qrange.start = start;
    qrange.end = start + gt_rand_max(query_width);
    GtArray *res = gt_array_new(sizeof (GtRange*));
    gt_interval_tree_find_all_overlapping(it, qrange.start, qrange.end, res);
    if (res)
    {
      /* generate reference overlapping interval list by linear search */
      GtArray *ref;
      unsigned long j;
      ref = gt_array_new(sizeof (GtRange*));
      for (j = 0; j < gt_array_size(arr); j++)
      {
        GtRange *this_rng;
        this_rng = *(GtRange**) gt_array_get(arr, j);
        if (gt_range_overlap(this_rng, &qrange))
        {
          gt_array_add(ref, this_rng);
        }
      }
      /* compare reference with interval tree query result */
      gt_array_sort_stable(ref, range_ptr_compare);
      gt_array_sort_stable(res, range_ptr_compare);
      /* must be equal */
      gt_ensure(had_err, gt_array_cmp(ref, res)==0);
      gt_array_delete(ref);
    }
    gt_array_delete(res);
  }
  gt_interval_tree_delete(it);

  it = gt_interval_tree_new(NULL);
  gt_array_reset(arr);

  /* generate test ranges */
  for (i = 0;i<num_testranges && !had_err;i++)
  {
    unsigned long start;
    GtIntervalTreeNode *new_node;
    start = gt_rand_max(gt_range_max_basepos);
    new_node = gt_interval_tree_node_new((void*) i, start,
                                          start + gt_rand_max(width));
    gt_interval_tree_insert(it, new_node);
  }
  gt_ensure(had_err, gt_interval_tree_size(it) == num_testranges);

  narr = gt_array_new(sizeof (GtIntervalTreeNode*));
  for (i = 0; i < num_testranges && !had_err; i++) {
    unsigned long idx, n, val;
    GtIntervalTreeNode *node = NULL;

    /* get all nodes referenced by the interval tree */
    interval_tree_find_all_internal(it, it->root, itree_test_get_node, 0,
                                    gt_range_max_basepos+width, narr);

    /* remove a random node */
    idx = gt_rand_max(gt_array_size(narr)-1);
    node = *(GtIntervalTreeNode**) gt_array_get(narr, idx);
    gt_ensure(had_err, node != NULL);
    val = (unsigned long) gt_interval_tree_node_get_data(node);
    gt_interval_tree_remove(it, node);
    gt_array_reset(narr);

    /* make sure that the node has disappeared */
    gt_ensure(had_err, gt_interval_tree_size(it) == num_testranges - (i+1));
    interval_tree_find_all_internal(it, it->root, itree_test_get_node, 0,
                                    gt_range_max_basepos+width, narr);
    gt_ensure(had_err, gt_array_size(narr) == num_testranges - (i+1));
    for (n = 0; !had_err && n < gt_array_size(narr); n++) {
      GtIntervalTreeNode *onode = *(GtIntervalTreeNode**) gt_array_get(narr, n);
      gt_ensure(had_err, (unsigned long) gt_interval_tree_node_get_data(onode)
                           != val);
    }
  }

  gt_array_delete(arr);
  gt_array_delete(narr);
  gt_interval_tree_delete(it);
  return had_err;
}
