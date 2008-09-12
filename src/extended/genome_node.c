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

#include <assert.h>
#include <stdarg.h>
#include "core/hashtable.h"
#include "core/ma.h"
#include "core/msort.h"
#include "core/queue.h"
#include "core/unused_api.h"
#include "extended/genome_node_rep.h"

typedef enum {
  NO_PARENT,
  ONE_PARENT,
  MULTIPLE_PARENTS
} ParentStatus;

typedef enum {
  TREE_STATUS_UNDETERMINED,
  IS_TREE,
  IS_NOT_A_TREE
} TreeStatus;

typedef struct {
  GtGenomeNodeTraverseFunc func;
  void *data;
} Traverse_children_info;

static int compare_gt_genome_node_type(GtGenomeNode *gn_a, GtGenomeNode *gn_b)
{
  void *sr_a, *sr_b, *sn_a, *sn_b;

  /* sequence regions first */
  sr_a = gt_genome_node_cast(gt_sequence_region_class(), gn_a);
  sr_b = gt_genome_node_cast(gt_sequence_region_class(), gn_b);

  if (sr_a && !sr_b)
    return -1;
  if (!sr_a && sr_b)
    return 1;

  /* sequence nodes last */
  sn_a = gt_genome_node_cast(gt_sequence_node_class(), gn_a);
  sn_b = gt_genome_node_cast(gt_sequence_node_class(), gn_b);

  if (sn_a && !sn_b)
    return 1;
  if (!sn_a && sn_b)
    return -1;

  return 0;
}

int gt_genome_node_cmp(GtGenomeNode *gn_a, GtGenomeNode *gn_b)
{
  int rval;
  assert(gn_a && gn_b);
  /* ensure that sequence regions come first and sequence nodes come last,
     otherwise we don't get a valid GFF3 stream */
  if ((rval = compare_gt_genome_node_type(gn_a, gn_b)))
    return rval;

  if ((rval = gt_str_cmp(gt_genome_node_get_idstr(gn_a),
                      gt_genome_node_get_idstr(gn_b)))) {
    return rval;
  }
  return gt_range_compare(gt_genome_node_get_range(gn_a),
                       gt_genome_node_get_range(gn_b));
}

static int compare_genome_nodes_with_delta(GtGenomeNode *gn_a,
                                           GtGenomeNode *gn_b,
                                           unsigned long delta)
{
  int rval;
  assert(gn_a && gn_b);
  /* ensure that sequence regions come first, otherwise we don't get a valid
     gff3 stream */
  if ((rval = compare_gt_genome_node_type(gn_a, gn_b)))
    return rval;

  if ((rval = gt_str_cmp(gt_genome_node_get_idstr(gn_a),
                      gt_genome_node_get_idstr(gn_b)))) {
    return rval;
  }
  return gt_range_compare_with_delta(gt_genome_node_get_range(gn_a),
                                  gt_genome_node_get_range(gn_b), delta);
}

static void set_tree_status(unsigned int *bit_field, TreeStatus tree_status)
{
  *bit_field &= ~(TREE_STATUS_MASK << TREE_STATUS_OFFSET);
  *bit_field |= tree_status << TREE_STATUS_OFFSET;
}

static TreeStatus get_tree_status(unsigned int bit_field)
{
  return (bit_field >> TREE_STATUS_OFFSET) & TREE_STATUS_MASK;
}

static void set_parent_status(unsigned int *bit_field,
                              ParentStatus parent_status)
{
  *bit_field &= ~(PARENT_STATUS_MASK << PARENT_STATUS_OFFSET);
  *bit_field |= parent_status << PARENT_STATUS_OFFSET;
}

static ParentStatus get_parent_status(unsigned int bit_field)
{
  return (bit_field >> PARENT_STATUS_OFFSET) & PARENT_STATUS_MASK;
}

static bool multiple_parents(unsigned int bit_field)
{
  if (get_parent_status(bit_field) == MULTIPLE_PARENTS)
    return true;
  return false;
}

static void add_parent(unsigned int *bit_field)
{
  switch (get_parent_status(*bit_field)) {
    case NO_PARENT:
      set_parent_status(bit_field, ONE_PARENT);
      break;
    case ONE_PARENT:
      set_parent_status(bit_field, MULTIPLE_PARENTS);
      break;
    case MULTIPLE_PARENTS:
      break;
  }
}

GtGenomeNode* gt_genome_node_create(const GtGenomeNodeClass *gnc)
{
  GtGenomeNode *gn;
  assert(gnc && gnc->size);
  gn                  = gt_malloc(gnc->size);
  gn->c_class         = gnc;
  gn->filename        = NULL; /* means the node is generated */
  gn->line_number     = 0;
  gn->children        = NULL; /* the children list is created on demand */
  gn->reference_count = 0;
  gn->bit_field       = 0;
  set_tree_status(&gn->bit_field, IS_TREE);
  return gn;
}

void gt_genome_node_set_origin(GtGenomeNode *gn,
                            GtStr *filename, unsigned int line_number)
{
  assert(gn && filename && line_number);
  gt_str_delete(gn->filename);
  gn->filename = gt_str_ref(filename);
  gn->line_number =line_number;
}

void* gt_genome_node_cast(const GtGenomeNodeClass *gnc, GtGenomeNode *gn)
{
  assert(gnc && gn);
  if (gn->c_class == gnc)
    return gn;
  return NULL;
}

static int increase_reference_count(GtGenomeNode *gn, GT_UNUSED void *data,
                                    GT_UNUSED GtError *err)
{
  gt_error_check(err);
  assert(gn);
  gn->reference_count++;
  return 0;
}

GtGenomeNode* gt_genome_node_ref(GtGenomeNode *gn)
{
  int had_err;
  had_err = increase_reference_count(gn, NULL, NULL);
  assert(!had_err); /* cannot happen, increase_reference_count() is sane */
  return gn;
}

int gt_genome_node_traverse_children_generic(GtGenomeNode *genome_node,
                                          void *data,
                                          GtGenomeNodeTraverseFunc traverse,
                                          bool traverse_only_once,
                                          bool depth_first, bool with_pseudo,
                                          GtError *err)
{
  GtArray *node_stack = NULL, *list_of_children;
  GT_Queue *node_queue = NULL;
  GtGenomeNode *gn, *gn_ref, *child_feature;
  GT_Dlistelem *dlistelem;
  unsigned long i;
  Hashtable *traversed_nodes = NULL;
  bool has_node_with_multiple_parents = false;
  int had_err = 0;

  if (!genome_node)
    return 0;

  /* create additional reference to <genome_node> (necessary if genome_node is
     freed by <traverse>) */
  gn_ref = gt_genome_node_ref(genome_node);

  if (depth_first) {
    node_stack = gt_array_new(sizeof (GtGenomeNode*));
    if (!with_pseudo &&
        gt_genome_node_cast(gt_genome_feature_class(), genome_node) &&
        gt_genome_feature_is_pseudo((GT_GenomeFeature*) genome_node)) {
      /* add the children backwards to traverse in order */
      for (dlistelem = gt_dlist_last(genome_node->children); dlistelem != NULL;
           dlistelem = gt_dlistelem_previous(dlistelem)) {
        child_feature = (GtGenomeNode*) gt_dlistelem_get_data(dlistelem);
        gt_array_add(node_stack, child_feature);
      }
    }
    else
      gt_array_add(node_stack, genome_node);
    assert(gt_array_size(node_stack));
  }
  else {
    node_queue = gt_queue_new();
    if (!with_pseudo &&
        gt_genome_node_cast(gt_genome_feature_class(), genome_node) &&
        gt_genome_feature_is_pseudo((GT_GenomeFeature*) genome_node)) {
      for (dlistelem = gt_dlist_first(genome_node->children); dlistelem != NULL;
           dlistelem = gt_dlistelem_next(dlistelem)) {
        child_feature = (GtGenomeNode*) gt_dlistelem_get_data(dlistelem);
        gt_queue_add(node_queue, child_feature);
      }
    }
    else
      gt_queue_add(node_queue, genome_node);
    assert(gt_queue_size(node_queue));
  }
  list_of_children = gt_array_new(sizeof (GtGenomeNode*));

  if (traverse_only_once)
  {
    static const HashElemInfo node_hashtype
      = { ht_ptr_elem_hash, { NULL }, sizeof (GtGenomeNode *),
          ht_ptr_elem_cmp, NULL, NULL };
    traversed_nodes = hashtable_new(node_hashtype);
  }

  while ((depth_first ? gt_array_size(node_stack)
                      : gt_queue_size(node_queue))) {
    if (depth_first)
      gn = *(GtGenomeNode**) gt_array_pop(node_stack);
    else
      gn = gt_queue_get(node_queue);
    gt_array_reset(list_of_children);
    if (gn->children) {
      /* a backup of the children array is necessary if traverse() frees the
         node */
      for (dlistelem = gt_dlist_first(gn->children); dlistelem != NULL;
           dlistelem = gt_dlistelem_next(dlistelem)) {
        child_feature = (GtGenomeNode*) gt_dlistelem_get_data(dlistelem);
        gt_array_add(list_of_children, child_feature);
      }
    }
    /* store the implications of <gn> to the tree status of <genome_node> */
    if (multiple_parents(gn->bit_field))
      has_node_with_multiple_parents = true;
    /* call traverse function */
    if (traverse) {
      had_err = traverse(gn, data, err);
      if (had_err)
        break;
    }
    for (i = 0; i < gt_array_size(list_of_children); i++) {
      if (depth_first) {
        /* we go backwards to traverse in order */
        child_feature = *(GtGenomeNode**) gt_array_get(list_of_children,
                                       gt_array_size(list_of_children) - i - 1);
      }
      else {
        child_feature = *(GtGenomeNode**) gt_array_get(list_of_children, i);
      }
      if (!traverse_only_once ||
          !hashtable_get(traversed_nodes, &child_feature)) {
        /* feature has not been traversed or has to be traversed multiple
           times */
        if (depth_first)
          gt_array_add(node_stack, child_feature);
        else
          gt_queue_add(node_queue, child_feature);
        if (traverse_only_once)
          hashtable_add(traversed_nodes, &child_feature);
      }
    }
  }

  /* save the tree status of the genome node */
  if (!had_err) {
    if (has_node_with_multiple_parents) {
      set_tree_status(&gn_ref->bit_field, IS_NOT_A_TREE);
      assert(get_tree_status(gn_ref->bit_field) == IS_NOT_A_TREE);
    }
    else {
      set_tree_status(&gn_ref->bit_field, IS_TREE);
      assert(get_tree_status(gn_ref->bit_field) ==IS_TREE);
    }
  }

  /* free */
  gt_genome_node_delete(gn_ref);
  if (traverse_only_once)
    hashtable_delete(traversed_nodes);
  gt_array_delete(list_of_children);
  gt_array_delete(node_stack);
  gt_queue_delete(node_queue);

  return had_err;
}

static int gt_genome_node_traverse_children_with_pseudo(GtGenomeNode
                                                        *genome_node,
                                                        void *data,
                                                       GtGenomeNodeTraverseFunc
                                                        traverse,
                                                        bool traverse_only_once,
                                                        GtError *err)
{
  return gt_genome_node_traverse_children_generic(genome_node, data, traverse,
                                               traverse_only_once, true, true,
                                               err);
}

GtGenomeNode* gt_genome_node_rec_ref(GtGenomeNode *gn)
{
  int had_err;
  assert(gn);
  had_err = gt_genome_node_traverse_children_with_pseudo(gn, NULL,
                                                      increase_reference_count,
                                                      true, NULL);
  assert(!had_err); /* cannot happen, increase_reference_count() is sane */
  return gn;
}

int gt_genome_node_traverse_children(GtGenomeNode *genome_node, void *data,
                                  GtGenomeNodeTraverseFunc traverse,
                                  bool traverse_only_once, GtError *err)
{
  return gt_genome_node_traverse_children_generic(genome_node, data, traverse,
                                               traverse_only_once, true, false,
                                               err);
}

int gt_genome_node_traverse_children_breadth(GtGenomeNode *genome_node,
                                             void *data,
                                             GtGenomeNodeTraverseFunc traverse,
                                             bool traverse_only_once,
                                             GtError *err)
{
  return gt_genome_node_traverse_children_generic(genome_node, data, traverse,
                                               traverse_only_once, false, false,
                                               err);
}

int gt_genome_node_traverse_direct_children(GtGenomeNode *gn,
                                         void *traverse_func_data,
                                         GtGenomeNodeTraverseFunc traverse,
                                         GtError *err)
{
  GT_Dlistelem *dlistelem;
  int had_err = 0;
  gt_error_check(err);
  if (!gn || !traverse)
    return 0;
  if (gn->children) {
    for (dlistelem = gt_dlist_first(gn->children); dlistelem != NULL;
         dlistelem = gt_dlistelem_next(dlistelem)) {
      had_err = traverse((GtGenomeNode*) gt_dlistelem_get_data(dlistelem),
                          traverse_func_data, err);
      if (had_err)
        break;
    }
  }
  return had_err;
}

const char* gt_genome_node_get_filename(const GtGenomeNode *gn)
{
  assert(gn);
  if (gn->filename)
    return gt_str_get(gn->filename);
  return "generated";
}

unsigned int gt_genome_node_get_line_number(const GtGenomeNode *gn)
{
  assert(gn);
  return gn->line_number;
}

unsigned long gt_genome_node_number_of_children(const GtGenomeNode *gn)
{
  assert(gn);
  return gt_dlist_size(gn->children);
}

GtStr* gt_genome_node_get_seqid(GtGenomeNode *gn)
{
  assert(gn && gn->c_class);
  if (gn->c_class->get_seqid)
    return gn->c_class->get_seqid(gn);
  return NULL;
}

GtStr* gt_genome_node_get_idstr(GtGenomeNode *gn)
{
  assert(gn && gn->c_class && gn->c_class->get_idstr);
  return gn->c_class->get_idstr(gn);
}

unsigned long gt_genome_node_get_start(GtGenomeNode *gn)
{
  GtRange range = gt_genome_node_get_range(gn);
  return range.start;
}

unsigned long gt_genome_node_get_end(GtGenomeNode *gn)
{
  GtRange range = gt_genome_node_get_range(gn);
  return range.end;
}

GtRange gt_genome_node_get_range(GtGenomeNode *gn)
{
  assert(gn && gn->c_class && gn->c_class->get_range);
  return gn->c_class->get_range(gn);
}

void gt_genome_node_set_range(GtGenomeNode *gn, GtRange range)
{
  assert(gn && gn->c_class && gn->c_class->set_range);
  gn->c_class->set_range(gn, range);
}

void gt_genome_node_change_seqid(GtGenomeNode *gn, GtStr *seqid)
{
  assert(gn && gn->c_class && gn->c_class->change_seqid && seqid);
  gn->c_class->change_seqid(gn, seqid);
}

int gt_genome_node_accept(GtGenomeNode *gn, GenomeVisitor *gv, GtError *err)
{
  gt_error_check(err);
  assert(gn && gv && gn->c_class && gn->c_class->accept);
  return gn->c_class->accept(gn, gv, err);
}

void gt_genome_node_add_child(GtGenomeNode *parent, GtGenomeNode *child)
{
  assert(parent && child);
  /* <parent> and <child> have the same seqid */
  assert(!gt_str_cmp(gt_genome_node_get_seqid(parent),
                     gt_genome_node_get_seqid(child)));
#ifndef NDEBUG
  if (gt_genome_node_cast(gt_genome_feature_class(), child)) {
    /* pseudo-features have to be top-level */
    assert(!gt_genome_feature_is_pseudo((GT_GenomeFeature*) child));
  }
#endif
  /* create children list on demand */
  if (!parent->children)
    parent->children = gt_dlist_new((GT_Compare) gt_genome_node_cmp);
  gt_dlist_add(parent->children, child); /* XXX: check for circles */
  /* update tree status of <parent> */
  set_tree_status(&parent->bit_field, TREE_STATUS_UNDETERMINED);
  /* update parent info of <child> */
  add_parent(&child->bit_field);
}

static int remove_leaf(GtGenomeNode *node, void *data, GT_UNUSED GtError *err)
{
  GT_Dlistelem *dlistelem;
  GtGenomeNode *child, *leaf = (GtGenomeNode*) data;
  gt_error_check(err);
  if (node != leaf && node->children) {
    for (dlistelem = gt_dlist_first(node->children); dlistelem != NULL;
         dlistelem = gt_dlistelem_next(dlistelem)) {
      child = (GtGenomeNode*) gt_dlistelem_get_data(dlistelem);
      if (child == leaf) {
        gt_dlist_remove(node->children, dlistelem);
        break;
      }
    }
  }
  return 0;
}

void gt_genome_node_remove_leaf(GtGenomeNode *tree, GtGenomeNode *leafn)
{
  int had_err;
  assert(tree && leafn);
  assert(!gt_genome_node_number_of_children(leafn));
  had_err = gt_genome_node_traverse_children(tree, leafn, remove_leaf, true,
                                             NULL);
  assert(!had_err); /* cannot happen, remove_leaf() is sane */
}

void gt_genome_node_mark(GtGenomeNode *gn)
{
  assert(gn);
  gn->bit_field |= 1;
}

bool gt_genome_node_is_marked(const GtGenomeNode *gn)
{
  assert(gn);
  return gn->bit_field & 1 ? true : false;
}

static int check_marked_status(GtGenomeNode *gn, void *data,
                               GT_UNUSED GtError *err)
{
  bool *marked = data;
  if (gt_genome_node_is_marked(gn))
    *marked = true;
  return 0;
}

bool gt_genome_node_contains_marked(GtGenomeNode *gn)
{
  bool contains_marked = false;
  int rval;
  assert(gn);
  rval = gt_genome_node_traverse_children(gn, &contains_marked,
                                       check_marked_status, true, NULL);
  assert(!rval); /* check_marked_status() is sane */
  return contains_marked;
}

bool gt_genome_node_has_children(GtGenomeNode *gn)
{
  assert(gn);
  if (!gn->children || gt_dlist_size(gn->children) == 0)
    return false;
  return true;
}

bool gt_genome_node_direct_children_do_not_overlap_generic(GtGenomeNode
                                                           *parent,
                                                           GtGenomeNode
                                                           *child)
{
  GtArray *children_ranges;
  GT_Dlistelem *dlistelem;
  GT_GenomeFeature *gf = NULL, *child_gf;
  GtRange range;
  bool rval;

  assert(parent);

  if (child)
    gf = gt_genome_node_cast(gt_genome_feature_class(), child);

  if (!parent->children)
    return true;

  /* get children ranges */
  children_ranges = gt_array_new(sizeof (GtRange));
  assert(parent->children);
  for (dlistelem = gt_dlist_first(parent->children); dlistelem != NULL;
       dlistelem = gt_dlistelem_next(dlistelem)) {
    if (!gf ||
        ((child_gf = gt_genome_node_cast(gt_genome_feature_class(),
                                      gt_dlistelem_get_data(dlistelem))) &&
         gt_genome_feature_get_type(gf) ==
         gt_genome_feature_get_type(child_gf))) {
      range = gt_genome_node_get_range((GtGenomeNode*)
                                    gt_dlistelem_get_data(dlistelem));
      gt_array_add(children_ranges, range);
    }
  }

  ranges_sort(children_ranges);
  assert(ranges_are_sorted(children_ranges));
  rval = ranges_do_not_overlap(children_ranges);

  gt_array_delete(children_ranges);

  return rval;
}

bool gt_genome_node_direct_children_do_not_overlap(GtGenomeNode *gn)
{
  return gt_genome_node_direct_children_do_not_overlap_generic(gn, NULL);
}

bool gt_genome_node_direct_children_do_not_overlap_st(GtGenomeNode *parent,
                                                   GtGenomeNode *child)
{
  return gt_genome_node_direct_children_do_not_overlap_generic(parent, child);
}

bool gt_genome_node_is_tree(GtGenomeNode *gn)
{
  bool status = false;
  assert(gn);
  switch (get_tree_status(gn->bit_field)) {
    case IS_TREE:
      status = true;
      break;
    case IS_NOT_A_TREE:
      status = false;
      break;
    case TREE_STATUS_UNDETERMINED:
      /* not implemented, the tree status must have been determined by a
         previous gt_genome_node_traverse_children() invocation */
    default: assert(0);
  }
  return status;
}

bool gt_genome_node_overlaps_nodes(GtGenomeNode *gn, GtArray *nodes)
{
  return gt_genome_node_overlaps_nodes_mark(gn, nodes, NULL);
}

bool gt_genome_node_overlaps_nodes_mark(GtGenomeNode *gn, GtArray *nodes,
                                             GtBittab *b)
{
  unsigned long i;
  GtGenomeNode *node;
  GtRange gn_range;
  bool rval = false;
#ifndef NDEBUG
  GtStr *gn_id;
  assert(gn && nodes);
  assert(!b || gt_bittab_size(b) == gt_array_size(nodes));
  gn_id = gt_genome_node_get_idstr(gn);
#endif
  gn_range = gt_genome_node_get_range(gn);

  for (i = 0; i < gt_array_size(nodes); i++) {
    node = *(GtGenomeNode**) gt_array_get(nodes, i);
    assert(!gt_str_cmp(gn_id, gt_genome_node_get_idstr(node)));
    if (gt_range_overlap(gn_range, gt_genome_node_get_range(node))) {
      rval = true;
      if (b)
        gt_bittab_set_bit(b, i);
      else
        break;
    }
  }
  return rval;
}

int gt_genome_node_compare(GtGenomeNode **gn_a, GtGenomeNode **gn_b)
{
  return gt_genome_node_cmp(*gn_a, *gn_b);
}

int gt_genome_node_compare_with_data(GtGenomeNode **gn_a, GtGenomeNode **gn_b,
                                  GT_UNUSED void *unused)
{
  return gt_genome_node_cmp(*gn_a, *gn_b);
}

int gt_genome_node_compare_delta(GtGenomeNode **gn_a, GtGenomeNode **gn_b,
                              void *delta)
{
  unsigned long *deltaptr = delta;
  assert(delta);
  return compare_genome_nodes_with_delta(*gn_a, *gn_b, *deltaptr);
}

void gt_genome_node_delete(GtGenomeNode *gn)
{
  if (!gn) return;
  if (gn->reference_count) { gn->reference_count--; return; }
  assert(gn->c_class);
  if (gn->c_class->free) gn->c_class->free(gn);
  gt_str_delete(gn->filename);
  gt_dlist_delete(gn->children);
  gt_free(gn);
}

static int free_genome_node(GtGenomeNode *gn, GT_UNUSED void *data,
                            GT_UNUSED GtError *err)
{
  gt_genome_node_delete(gn);
  return 0;
}

void gt_genome_node_rec_delete(GtGenomeNode *gn)
{
  int had_err;
  if (!gn) return;
  had_err = gt_genome_node_traverse_children_with_pseudo(gn, NULL,
                                                      free_genome_node, true,
                                                      NULL);
  assert(!had_err); /* cannot happen, free_genome_node() is sane */
}

void gt_genome_nodes_sort(GtArray *nodes)
{
  qsort(gt_array_get_space(nodes), gt_array_size(nodes),
        sizeof (GtGenomeNode*), (GT_Compare) gt_genome_node_compare);
}

void gt_genome_nodes_sort_stable(GtArray *nodes)
{
  gt_msort(gt_array_get_space(nodes), gt_array_size(nodes),
           sizeof (GtGenomeNode*), (GT_Compare) gt_genome_node_compare);
}

bool gt_genome_nodes_are_equal_sequence_regions(GtGenomeNode *gn_a,
                                                GtGenomeNode *gn_b)
{
  void *sr_a, *sr_b;

  sr_a = gn_a ? gt_genome_node_cast(gt_sequence_region_class(), gn_a) : NULL;
  sr_b = gn_b ? gt_genome_node_cast(gt_sequence_region_class(), gn_b) : NULL;

  if (sr_a && sr_b && !gt_str_cmp(gt_genome_node_get_seqid(gn_a),
                               gt_genome_node_get_seqid(gn_b))) {
    return true;
  }
  return false;
}

bool gt_genome_nodes_are_sorted(const GtArray *nodes)
{
  unsigned long i;
  assert(nodes);
  for (i = 1; i < gt_array_size(nodes); i++) {
    if (gt_genome_node_compare(gt_array_get(nodes, i-1),
                               gt_array_get(nodes, i)) > 0) {
      return false;
    }
  }
  return true;
}
