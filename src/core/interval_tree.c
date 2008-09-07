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

#include "core/ensure.h"
#include "core/interval_tree.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/range.h"
#include "core/str.h"
#include "core/unused.h"
#include <string.h>

typedef enum IntervalTreeNodeColor {
  BLACK,
  RED
} IntervalTreeNodeColor ;

struct IntervalTree {
  IntervalTreeNode *root;
  unsigned long size;
  GT_FreeFunc free_func;
};

struct IntervalTreeNode {
  IntervalTreeNode *parent, *left, *right;
  void *data;
  IntervalTreeNodeColor color;
  unsigned long low, high, max;
};

IntervalTreeNode* interval_tree_node_new(void *data,
                                         unsigned long low,
                                         unsigned long high)
{
  IntervalTreeNode* n;
  n = ma_calloc(1, sizeof (IntervalTreeNode));
  n->low = low;
  n->high = high;
  n->data = data;
  return n;
}

IntervalTree* interval_tree_new(GT_FreeFunc func)
{
  IntervalTree *it;
  it = ma_calloc(1, sizeof (IntervalTree));
  it->free_func = func;
  return it;
}

unsigned long interval_tree_size(IntervalTree *it)
{
  assert(it);
  return it->size;
}

void* interval_tree_node_get_data(IntervalTreeNode *n)
{
  assert(n);
  return n->data;
}

void interval_tree_node_delete(IntervalTree *it, IntervalTreeNode *n)
{
  if (!n) return;
  if (n->data && it->free_func)
    it->free_func(n->data);
  ma_free(n);
}

static void interval_tree_node_rec_delete(IntervalTree *it, IntervalTreeNode *n)
{
  if (!n) return;
  interval_tree_node_rec_delete(it, n->left);
  interval_tree_node_rec_delete(it, n->right);
  interval_tree_node_delete(it, n);
}

static IntervalTreeNode* interval_tree_search_internal(IntervalTreeNode *node,
                                                       unsigned long low,
                                                       unsigned long high)
{
  IntervalTreeNode *x;
  x = node;
  while (x && !(low <= x->high && x->low <= high)) {
    if (x->left && x->left->max >= low)
      x = x->left;
    else
      x = x->right;
  }
  return x;
}

IntervalTreeNode* interval_tree_find_first_overlapping(IntervalTree *it,
                                                       unsigned long low,
                                                       unsigned long high)
{
  assert(it);
  if (!it->root)
    return NULL;
  return interval_tree_search_internal(it->root, low, high);
}

static int interval_tree_traverse_internal(IntervalTreeNode *node,
                                           IntervalTreeIteratorFunc func,
                                           void *data)
{
  int had_err = 0;
  if (!node) return 0;
  if (!had_err)
    had_err = interval_tree_traverse_internal(node->left, func, data);
  if (!had_err)
    had_err = interval_tree_traverse_internal(node->right, func, data);
  if (!had_err)
    had_err = func(node, data);
  return had_err;
}

int interval_tree_traverse(IntervalTree *it, IntervalTreeIteratorFunc func,
                           void *data)
{
  if (!it->root)
    return 0;
  return interval_tree_traverse_internal(it->root, func, data);
}

static void interval_tree_find_all_internal(IntervalTreeNode *node,
                                            unsigned long low,
                                            unsigned long high, GT_Array *a)
{
  IntervalTreeNode* x;
  if (!node) return;
  x = node;
  if (low <= x->high && x->low <= high)
    gt_array_add(a, x->data);
  /* recursively search left and right subtrees */
  if (x->left && low <= x->left->max)
    interval_tree_find_all_internal(x->left, low, high, a);
  if (x->right && low <= x->right->max)
    interval_tree_find_all_internal(x->right, low, high, a);
}

void interval_tree_find_all_overlapping(IntervalTree *it, unsigned long start,
                                        unsigned long end, GT_Array* a)
{
  assert(it && a && start <= end);
  if (!it->root) return;
  interval_tree_find_all_internal(it->root, start, end, a);
}

static void interval_tree_left_rotate(IntervalTreeNode **root,
                                      IntervalTreeNode *x)
{
  IntervalTreeNode *y;
  y = x->right;
  x->right = y->left;
  if (y->left)
    y->left->parent = x;
  y->parent = x->parent;
  if (!x->parent)
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
  if (x->left && x->left->max > x->max)
    x->max = x->left->max;
  if (x->right && x->right->max > x->max)
    x->max = x->right->max;
  y->max = y->high;
  if (y->left && y->left->max > y->max)
    y->max = y->left->max;
  if (y->right && y->right->max > y->max)
    y->max = y->right->max;
}

static void interval_tree_right_rotate(IntervalTreeNode **root,
                                       IntervalTreeNode *y)
{
  IntervalTreeNode *x;
  x = y->left;
  y->left = x->right;
  if (x->right)
    x->right->parent = y;
  x->parent = y->parent;
  if (!y->parent)
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
  if (x->left && x->left->max > x->max)
    x->max = x->left->max;
  if (x->right && x->right->max > x->max)
    x->max = x->right->max;
  y->max = y->high;
  if (y->left && y->left->max > y->max)
    y->max = y->left->max;
  if (y->right && y->right->max > y->max)
    y->max = y->right->max;
}

/* this is the insert routine from Cormen et al, p. 280*/
static void tree_insert(IntervalTreeNode **root, IntervalTreeNode *z)
{
  IntervalTreeNode *x, *y;
  y = NULL;
  x = *root;
  z->max = z->high;
  while (x)
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
  if (!y)
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
static void interval_tree_insert_internal(IntervalTreeNode **root,
                                          IntervalTreeNode *z)
{
  IntervalTreeNode* y;
  tree_insert(root, z);
  z->color = RED;
  while (z != *root && z->parent->color == RED)
  {
    if (z->parent == z->parent->parent->left)
    {
      y = z->parent->parent->right;
      if (y && y->color == RED)
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
          interval_tree_left_rotate(root, z);
        }
        z->parent->color = BLACK;
        z->parent->parent->color = RED;
        interval_tree_right_rotate(root, z->parent->parent);
      }
    }
    else
    {
      y = z->parent->parent->left;
      if (y && y->color == RED)
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
          interval_tree_right_rotate(root, z);
        }
        z->parent->color = BLACK;
        z->parent->parent->color = RED;
        interval_tree_left_rotate(root, z->parent->parent);
      }
    }
  }
  (*root)->color = BLACK;
}

void interval_tree_insert(IntervalTree *it, IntervalTreeNode *n)
{
  assert(it && n);
  if (!it->root)
  {
    it->root = n;
  } else interval_tree_insert_internal(&(it->root), n);
  it->size++;
}

void interval_tree_delete(IntervalTree *it)
{
  if (!it) return;
  interval_tree_node_rec_delete(it, it->root);
  ma_free(it);
}

static int gt_range_ptr_compare(const void *r1p, const void *r2p)
{
  int ret;
  assert(r1p && r2p);
  ret = gt_range_compare(**(GT_Range**) r1p,**(GT_Range**) r2p);
  /* It could be that two identical ranges with different pointers are
     present. If so, compare pointers instead to get a canonical ordering. */
  if (ret == 0 && *(GT_Range**) r1p != *(GT_Range**) r2p)
  {
    if (*(GT_Range**) r1p < *(GT_Range**) r2p)
      ret = -1;
    else
      ret = 1;
  }
  return ret;
}

int interval_tree_unit_test(UNUSED GT_Error *err)
{
  IntervalTree *it = NULL;
  IntervalTreeNode *res = NULL;
  int had_err = 0, i = 0;
  int num_testranges = 3000;
  int num_samples = 300000;
  int num_find_all_samples = 10000;
  int gt_range_max_basepos = 90000;
  int width = 700;
  int query_width = 5000;

  GT_Range *res_rng = NULL, qrange;
  GT_Array *arr;

  arr = gt_array_new(sizeof (GT_Range*));

  /* generate test ranges */
  for (i = 0;i<num_testranges;i++)
  {
    unsigned long start;
    GT_Range *rng;
    rng  = ma_calloc(1, sizeof (GT_Range));
    start = rand_max(gt_range_max_basepos);
    rng->start = start;
    rng->end = start + rand_max(width);
    gt_array_add(arr, rng);
  }

  it = interval_tree_new(ma_free_func);

  /* insert ranges */
  for (i = 0; i < num_testranges && !had_err; i++)
  {
    IntervalTreeNode *new_node;
    GT_Range *rng;
    rng = *(GT_Range**) gt_array_get(arr, i);
    new_node = interval_tree_node_new(rng, rng->start, rng->end);
    interval_tree_insert(it, new_node);
  }

  ensure(had_err, interval_tree_size(it) == num_testranges);

  /* perform test queries */
  for (i = 0; i < num_samples && !had_err; i++)
  {
    unsigned long start = rand_max(gt_range_max_basepos);
    qrange.start = start;
    qrange.end = start + rand_max(width);
    res = interval_tree_find_first_overlapping(it, qrange.start, qrange.end);
    if (res)
    {
      /* we have a hit, check if really overlapping */
      res_rng = (GT_Range*) interval_tree_node_get_data(res);
      ensure(had_err, gt_range_overlap(qrange, *res_rng));
    } else {
      /* no hit, check whether there really is no overlapping
         interval in tree */
      GT_Range *this_rng;
      unsigned long j;
      bool found = false;
      for (j = 0; j < gt_array_size(arr); j++)
      {
        this_rng = *(GT_Range**) gt_array_get(arr, j);
        if (gt_range_overlap(*this_rng, qrange))
        {
          found = true;
          break;
        }
      }
      ensure(had_err, !found);
    }
  }

  /* test searching for all overlapping intervals */
  for (i = 0; i < num_find_all_samples && !had_err; i++)
  {
    unsigned long start = rand_max(gt_range_max_basepos);
    qrange.start = start;
    qrange.end = start + rand_max(query_width);
    GT_Array *res = gt_array_new(sizeof (GT_Range*));
    interval_tree_find_all_overlapping(it, qrange.start, qrange.end, res);
    if (res)
    {
      /* generate reference overlapping interval list by linear search */
      GT_Array *ref;
      unsigned long j;
      ref = gt_array_new(sizeof (GT_Range*));
      for (j = 0; j < gt_array_size(arr); j++)
      {
        GT_Range *this_rng;
        this_rng = *(GT_Range**) gt_array_get(arr, j);
        if (gt_range_overlap(*this_rng, qrange))
        {
          gt_array_add(ref, this_rng);
        }
      }
      /* compare reference with interval tree query result */
      gt_array_sort(ref, gt_range_ptr_compare);
      gt_array_sort(res, gt_range_ptr_compare);
      /* must be equal */
      ensure(had_err, gt_array_cmp(ref, res)==0);
      gt_array_delete(ref);
    }
    gt_array_delete(res);
  }
  gt_array_delete(arr);
  interval_tree_delete(it);
  return had_err;
}
