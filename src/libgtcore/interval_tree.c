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

#include "libgtcore/ensure.h"
#include "libgtcore/interval_tree.h"
#include "libgtcore/ma.h"
#include "libgtcore/mathsupport.h"
#include "libgtcore/range.h"
#include "libgtcore/str.h"
#include "libgtcore/unused.h"
#include <string.h>

typedef enum IntervalTreeNodeColor {
  BLACK,
  RED
} IntervalTreeNodeColor ;

struct IntervalTree {
  IntervalTreeNode *root;
  unsigned long size;
};

struct IntervalTreeNode {
  IntervalTreeNode *parent;
  IntervalTreeNode *left;
  IntervalTreeNode *right;
  FreeFunc free_func;
  void *data;
  IntervalTreeNodeColor color;
  unsigned long low;
  unsigned long high;
  unsigned long max;
};

IntervalTreeNode* interval_tree_node_new(void *data,
                                         unsigned long low,
                                         unsigned long high,
                                         FreeFunc free_func)
{
  IntervalTreeNode* n;
  n = ma_calloc(1, sizeof (IntervalTreeNode));
  n->low = low;
  n->high = high;
  n->data = data;
  n->free_func = free_func;
  return n;
}

IntervalTree* interval_tree_new(void)
{
  IntervalTree *it;
  it = ma_calloc(1, sizeof (IntervalTree));
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

void interval_tree_node_delete(IntervalTreeNode *n)
{
  if (!n) return;
  if (n->data && n->free_func)
    n->free_func(n->data);
  ma_free(n);
}

static void interval_tree_node_rec_delete(IntervalTreeNode *n)
{
  if (!n) return;
  interval_tree_node_rec_delete(n->left);
  interval_tree_node_rec_delete(n->right);
  interval_tree_node_delete(n);
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
    had_err = (int) func(node, data);
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
                                            unsigned long high, Array *a)
{
  IntervalTreeNode* x;
  if (!node) return;
  x = node;
  if (low <= x->high && x->low <= high)
    array_add(a, x->data);
  /* recursively search left and right subtrees */
  if (x->left && low <= x->left->max)
    interval_tree_find_all_internal(x->left, low, high, a);
  if (x->right && low <= x->right->max)
    interval_tree_find_all_internal(x->right, low, high, a);
}

void interval_tree_find_all_overlapping(IntervalTree *it, unsigned long start,
                                        unsigned long end, Array* a)
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
  interval_tree_node_rec_delete(it->root);
  ma_free(it);
}

static int range_ptr_compare(const void *r1p, const void *r2p)
{
  int ret;
  assert(r1p && r2p);
  ret = range_compare(**(Range**) r1p,**(Range**) r2p);
  /* It could be that two identical ranges with different pointers are
     present. If so, compare pointers instead to get a canonical ordering. */
  if (ret == 0 && *(Range**) r1p != *(Range**) r2p)
  {
    if (*(Range**) r1p < *(Range**) r2p)
      ret = -1;
    else
      ret = 1;
  }
  return ret;
}

int interval_tree_unit_test(UNUSED Error *err)
{
  IntervalTree *it = NULL;
  IntervalTreeNode *res = NULL;
  int had_err = 0, i = 0;
  int num_testranges = 3000;
  int num_samples = 300000;
  int num_find_all_samples = 10000;
  int range_max_basepos = 90000;
  int width = 700;
  int query_width = 5000;

  Range *res_rng = NULL, qrange;
  Array *arr;

  arr = array_new(sizeof (Range*));

  /* generate test ranges */
  for (i = 0;i<num_testranges;i++)
  {
    unsigned long start;
    Range *rng;
    rng  = ma_calloc(1, sizeof (Range));
    start = rand_max(range_max_basepos);
    rng->start = start;
    rng->end = start + rand_max(width);
    array_add(arr, rng);
  }

  it = interval_tree_new();

  /* insert ranges */
  for (i = 0; i < num_testranges && !had_err; i++)
  {
    IntervalTreeNode *new_node;
    Range *rng;
    rng = *(Range**) array_get(arr, i);
    new_node = interval_tree_node_new(rng, rng->start, rng->end, ma_free_func);
    interval_tree_insert(it, new_node);
  }

  ensure(had_err, interval_tree_size(it) == num_testranges);

  /* perform test queries */
  for (i = 0; i < num_samples && !had_err; i++)
  {
    unsigned long start = rand_max(range_max_basepos);
    qrange.start = start;
    qrange.end = start + rand_max(width);
    res = interval_tree_find_first_overlapping(it, qrange.start, qrange.end);
    if (res)
    {
      /* we have a hit, check if really overlapping */
      res_rng = (Range*) interval_tree_node_get_data(res);
      ensure(had_err, range_overlap(qrange, *res_rng));
    } else {
      /* no hit, check whether there really is no overlapping
         interval in tree */
      Range *this_rng;
      unsigned long j;
      bool found = false;
      for (j=0;j<array_size(arr);j++)
      {
        this_rng = *(Range**) array_get(arr, j);
        if (range_overlap(*this_rng, qrange))
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
    unsigned long start = rand_max(range_max_basepos);
    qrange.start = start;
    qrange.end = start + rand_max(query_width);
    Array *res = array_new(sizeof (Range*));
    interval_tree_find_all_overlapping(it, qrange.start, qrange.end, res);
    if (res)
    {
      /* generate reference overlapping interval list by linear search */
      Array *ref;
      unsigned long j;
      ref = array_new(sizeof (Range*));
      for (j=0;j<array_size(arr);j++)
      {
        Range *this_rng;
        this_rng = *(Range**) array_get(arr, j);
        if (range_overlap(*this_rng, qrange))
        {
          array_add(ref, this_rng);
        }
      }
      /* compare reference with interval tree query result */
      array_sort(ref, range_ptr_compare);
      array_sort(res, range_ptr_compare);
      /* must be equal */
      ensure(had_err, array_cmp(ref, res)==0);
      array_delete(ref);
    }
    array_delete(res);
  }
  array_delete(arr);
  interval_tree_delete(it);
  return had_err;
}
