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
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"
#include "libgtcore/msort.h"
#include "libgtcore/queue.h"
#include "libgtcore/unused.h"
#include "libgtext/genome_node_rep.h"

typedef struct {
  GenomeNodeTraverseFunc func;
  void *data;
} Traverse_children_info;

static int compare_genome_node_type(GenomeNode *gn_a, GenomeNode *gn_b)
{
  void *sr_a, *sr_b;

  sr_a = genome_node_cast(sequence_region_class(), gn_a);
  sr_b = genome_node_cast(sequence_region_class(), gn_b);

  if (sr_a && !sr_b)
    return -1;
  if (!sr_a && sr_b)
    return 1;
  return 0;
}

int genome_node_cmp(GenomeNode *gn_a, GenomeNode *gn_b)
{
  int rval;
  assert(gn_a && gn_b);
  /* ensure that sequence regions come first, otherwise we don't get a valid
     gff3 stream */
  if ((rval = compare_genome_node_type(gn_a, gn_b)))
    return rval;

  if ((rval = str_cmp(genome_node_get_idstr(gn_a),
                      genome_node_get_idstr(gn_b)))) {
    return rval;
  }
  return range_compare(genome_node_get_range(gn_a),
                       genome_node_get_range(gn_b));
}

static int compare_genome_nodes_with_delta(GenomeNode *gn_a, GenomeNode *gn_b,
                                           unsigned long delta)
{
  int rval;
  assert(gn_a && gn_b);
  /* ensure that sequence regions come first, otherwise we don't get a valid
     gff3 stream */
  if ((rval = compare_genome_node_type(gn_a, gn_b)))
    return rval;

  if ((rval = str_cmp(genome_node_get_idstr(gn_a),
                      genome_node_get_idstr(gn_b)))) {
    return rval;
  }
  return range_compare_with_delta(genome_node_get_range(gn_a),
                                  genome_node_get_range(gn_b), delta);
}

GenomeNode* genome_node_create(const GenomeNodeClass *gnc, Str *filename,
                               unsigned long line_number)
{
  GenomeNode *gn;
  assert(gnc && gnc->size);
  gn                  = ma_malloc(gnc->size);
  gn->c_class         = gnc;
  gn->filename        = filename ? str_ref(filename)
                                 : str_new_cstr("generated");
  gn->line_number     = line_number;
  gn->children        = NULL; /* the children list is created on demand */
  gn->reference_count = 0;
  gn->info            = 0;
  genome_node_info_set_tree_status(&gn->info, GENOME_NODE_IS_TREE);
  gn->mark            = false;
  return gn;
}

void* genome_node_cast(const GenomeNodeClass *gnc, GenomeNode *gn)
{
  assert(gnc && gn);
  if (gn->c_class == gnc)
    return gn;
  return NULL;
}

static int increase_reference_count(GenomeNode *gn, UNUSED void *data,
                                    UNUSED Error *err)
{
  error_check(err);
  assert(gn);
  gn->reference_count++;
  return 0;
}

GenomeNode* genome_node_ref(GenomeNode *gn)
{
  int had_err;
  had_err = increase_reference_count(gn, NULL, NULL);
  assert(!had_err); /* cannot happen, increase_reference_count() is sane */
  return gn;
}

GenomeNode* genome_node_rec_ref(GenomeNode *gn)
{
  int had_err;
  assert(gn);
  had_err = genome_node_traverse_children(gn, NULL, increase_reference_count,
                                          true, NULL);
  assert(!had_err); /* cannot happen, increase_reference_count() is sane */
  return gn;
}

int genome_node_traverse_children_generic(GenomeNode *genome_node,
                                          void *data,
                                          GenomeNodeTraverseFunc traverse,
                                          bool traverse_only_once,
                                          bool depth_first, Error *e)
{
  Array *node_stack = NULL, *list_of_children;
  Queue *node_queue = NULL;
  GenomeNode *gn, *gn_ref, *child_feature;
  Dlistelem *dlistelem;
  unsigned long i;
  Hashtable *traversed_nodes = NULL;
  bool has_node_with_multiple_parents = false;
  int had_err = 0;

  if (!genome_node)
    return 0;

  /* create additional reference to <genome_node> (necessary if genome_node is
     freed by <traverse>) */
  gn_ref = genome_node_ref(genome_node);

  if (depth_first) {
    node_stack = array_new(sizeof (GenomeNode*));
    array_add(node_stack, genome_node);
  }
  else {
    node_queue = queue_new();
    queue_add(node_queue, genome_node);
  }
  list_of_children = array_new(sizeof (GenomeNode*));

  if (traverse_only_once)
    traversed_nodes = hashtable_new(HASH_DIRECT, NULL, NULL);

  while ((depth_first ? array_size(node_stack): queue_size(node_queue))) {
    if (depth_first)
      gn = *(GenomeNode**) array_pop(node_stack);
    else
      gn = queue_get(node_queue);
    array_reset(list_of_children);
    if (gn->children) {
      /* a backup of the children array is necessary if traverse() frees the
         node */
      for (dlistelem = dlist_first(gn->children); dlistelem != NULL;
           dlistelem = dlistelem_next(dlistelem)) {
        child_feature = (GenomeNode*) dlistelem_get_data(dlistelem);
        array_add(list_of_children, child_feature);
      }
    }
    /* store the implications of <gn> to the tree status of <genome_node> */
    if (genome_node_info_multiple_parents(&gn->info))
      has_node_with_multiple_parents = true;
    /* call traverse function */
    if (traverse) {
      had_err = traverse(gn, data, e);
      if (had_err)
        break;
    }
    for (i = 0; i < array_size(list_of_children); i++) {
      if (depth_first) {
        /* we go backwards to traverse in order */
        child_feature = *(GenomeNode**) array_get(list_of_children,
                                                  array_size(list_of_children)
                                                  - i - 1);
      }
      else {
        child_feature = *(GenomeNode**) array_get(list_of_children, i);
      }
      if (!traverse_only_once ||
          !hashtable_get(traversed_nodes, child_feature)) {
        /* feature has not been traversed or has to be traversed multiple
           times */
        if (depth_first)
          array_add(node_stack, child_feature);
        else
          queue_add(node_queue, child_feature);
        if (traverse_only_once)
          hashtable_add(traversed_nodes, child_feature, child_feature);
      }
    }
  }

  /* save the tree status of the genome node */
  if (!had_err) {
    if (has_node_with_multiple_parents) {
      genome_node_info_set_tree_status(&gn_ref->info,
                                       GENOME_NODE_IS_NOT_A_TREE);
      assert(genome_node_info_get_tree_status(&gn_ref->info) ==
             GENOME_NODE_IS_NOT_A_TREE);
    }
    else {
      genome_node_info_set_tree_status(&gn_ref->info, GENOME_NODE_IS_TREE);
      assert(genome_node_info_get_tree_status(&gn_ref->info) ==
             GENOME_NODE_IS_TREE);
    }
  }

  /* free */
  genome_node_delete(gn_ref);
  if (traverse_only_once)
    hashtable_delete(traversed_nodes);
  array_delete(list_of_children);
  array_delete(node_stack);
  queue_delete(node_queue);

  return had_err;
}

int genome_node_traverse_children(GenomeNode *genome_node, void *data,
                                  GenomeNodeTraverseFunc traverse,
                                  bool traverse_only_once, Error *e)
{
  return genome_node_traverse_children_generic(genome_node, data, traverse,
                                               traverse_only_once, true, e);
}

int genome_node_traverse_children_breadth(GenomeNode *genome_node, void *data,
                                          GenomeNodeTraverseFunc traverse,
                                          bool traverse_only_once, Error *e)
{
  return genome_node_traverse_children_generic(genome_node, data, traverse,
                                               traverse_only_once, false, e);
}

int genome_node_traverse_direct_children(GenomeNode *gn,
                                         void *traverse_func_data,
                                         GenomeNodeTraverseFunc traverse,
                                         Error *e)
{
  Dlistelem *dlistelem;
  int had_err = 0;
  error_check(e);
  if (!gn || !traverse)
    return 0;
  if (gn->children) {
    for (dlistelem = dlist_first(gn->children); dlistelem != NULL;
         dlistelem = dlistelem_next(dlistelem)) {
      had_err = traverse((GenomeNode*) dlistelem_get_data(dlistelem),
                          traverse_func_data, e);
      if (had_err)
        break;
    }
  }
  return had_err;
}

const char* genome_node_get_filename(const GenomeNode *gn)
{
  assert(gn);
  return str_get(gn->filename);
}

unsigned long genome_node_get_line_number(const GenomeNode *gn)
{
  assert(gn);
  return gn->line_number;
}

unsigned long genome_node_number_of_children(const GenomeNode *gn)
{
  assert(gn);
  return dlist_size(gn->children);
}

Str* genome_node_get_seqid(GenomeNode *gn)
{
  assert(gn && gn->c_class);
  if (gn->c_class->get_seqid)
    return gn->c_class->get_seqid(gn);
  return NULL;
}

Str* genome_node_get_idstr(GenomeNode *gn)
{
  assert(gn && gn->c_class && gn->c_class->get_idstr);
  return gn->c_class->get_idstr(gn);
}

unsigned long genome_node_get_start(GenomeNode *gn)
{
  Range range = genome_node_get_range(gn);
  return range.start;
}

unsigned long genome_node_get_end(GenomeNode *gn)
{
  Range range = genome_node_get_range(gn);
  return range.end;
}

Range genome_node_get_range(GenomeNode *gn)
{
  assert(gn && gn->c_class && gn->c_class->get_range);
  return gn->c_class->get_range(gn);
}

void genome_node_set_range(GenomeNode *gn, Range range)
{
  assert(gn && gn->c_class && gn->c_class->set_range);
  gn->c_class->set_range(gn, range);
}

void genome_node_set_seqid(GenomeNode *gn, Str *seqid)
{
  assert(gn && gn->c_class && gn->c_class->set_seqid && seqid);
  gn->c_class->set_seqid(gn, seqid);
}

int genome_node_accept(GenomeNode *gn, GenomeVisitor *gv, Error *e)
{
  error_check(e);
  assert(gn && gv && gn->c_class && gn->c_class->accept);
  return gn->c_class->accept(gn, gv, e);
}

void genome_node_is_part_of_genome_node(GenomeNode *parent, GenomeNode *child)
{
  assert(parent && child);
  /* create children list on demand */
  if (!parent->children)
    parent->children = dlist_new((Compare) genome_node_cmp);
  dlist_add(parent->children, child); /* XXX: check for circles */
  /* update tree status of <parent> */
  genome_node_info_set_tree_status(&parent->info,
                                   GENOME_NODE_TREE_STATUS_UNDETERMINED);
  /* update parent info of <child> */
  genome_node_info_add_parent(&child->info);
}

static int remove_leaf(GenomeNode *node, void *data, UNUSED Error *err)
{
  Dlistelem *dlistelem;
  GenomeNode *child, *leaf = (GenomeNode*) data;
  error_check(err);
  if (node != leaf && node->children) {
    for (dlistelem = dlist_first(node->children); dlistelem != NULL;
         dlistelem = dlistelem_next(dlistelem)) {
      child = (GenomeNode*) dlistelem_get_data(dlistelem);
      if (child == leaf) {
        dlist_remove(node->children, dlistelem);
        break;
      }
    }
  }
  return 0;
}

void genome_node_remove_leaf(GenomeNode *tree, GenomeNode *leafn)
{
  int had_err;
  assert(tree && leafn);
  assert(!genome_node_number_of_children(leafn));
  had_err = genome_node_traverse_children(tree, leafn, remove_leaf, true, NULL);
  assert(!had_err); /* cannot happen, remove_leaf() is sane */
}

void genome_node_mark(GenomeNode *gn)
{
  assert(gn);
  gn->mark = true;
}

bool genome_node_is_marked(const GenomeNode *gn)
{
  assert(gn);
  return gn->mark;
}

static int check_marked_status(GenomeNode *gn, void *data, UNUSED Error *err)
{
  bool *marked = data;
  if (gn->mark)
    *marked = true;
  return 0;
}

bool genome_node_contains_marked(GenomeNode *gn)
{
  bool contains_marked = false;
  int rval;
  assert(gn);
  rval = genome_node_traverse_children(gn, &contains_marked,
                                       check_marked_status, true, NULL);
  assert(!rval); /* check_marked_status() is sane */
  return contains_marked;
}

bool genome_node_has_children(GenomeNode *gn)
{
  assert(gn);
  if (!gn->children || dlist_size(gn->children) == 0)
    return false;
  return true;
}

bool genome_node_direct_children_do_not_overlap_generic(GenomeNode *parent,
                                                        GenomeNode *child)
{
  Array *children_ranges;
  Dlistelem *dlistelem;
  GenomeFeature *gf = NULL, *child_gf;
  Range range;
  bool rval;

  assert(parent);

  if (child)
    gf = genome_node_cast(genome_feature_class(), child);

  if (!parent->children)
    return true;

  /* get children ranges */
  children_ranges = array_new(sizeof (Range));
  assert(parent->children);
  for (dlistelem = dlist_first(parent->children); dlistelem != NULL;
       dlistelem = dlistelem_next(dlistelem)) {
    if (!gf ||
        ((child_gf = genome_node_cast(genome_feature_class(),
                                      dlistelem_get_data(dlistelem))) &&
         genome_feature_get_type(gf) == genome_feature_get_type(child_gf))) {
      range = genome_node_get_range((GenomeNode*)
                                    dlistelem_get_data(dlistelem));
      array_add(children_ranges, range);
    }
  }

  ranges_sort(children_ranges);
  assert(ranges_are_sorted(children_ranges));
  rval = ranges_do_not_overlap(children_ranges);

  array_delete(children_ranges);

  return rval;
}

bool genome_node_direct_children_do_not_overlap(GenomeNode *gn)
{
  return genome_node_direct_children_do_not_overlap_generic(gn, NULL);
}

bool genome_node_direct_children_do_not_overlap_st(GenomeNode *parent,
                                                   GenomeNode *child)
{
  return genome_node_direct_children_do_not_overlap_generic(parent, child);
}

bool genome_node_is_tree(GenomeNode *gn)
{
  bool status = false;
  assert(gn);
  switch (genome_node_info_get_tree_status(&gn->info)) {
    case GENOME_NODE_IS_TREE:
      status = true;
      break;
    case GENOME_NODE_IS_NOT_A_TREE:
      status = false;
      break;
    case GENOME_NODE_TREE_STATUS_UNDETERMINED:
      /* not implemented, the tree status must have been determined by a
         previous genome_node_traverse_children() invocation */
    default: assert(0);
  }
  return status;
}

bool genome_node_overlaps_nodes(GenomeNode *gn, Array *nodes)
{
  return genome_node_overlaps_nodes_mark(gn, nodes, NULL);
}

bool genome_node_overlaps_nodes_mark(GenomeNode *gn, Array *nodes,
                                             Bittab *b)
{
  unsigned long i;
  GenomeNode *node;
  Range gn_range;
  bool rval = false;
#ifndef NDEBUG
  Str *gn_id;
  assert(gn && nodes);
  assert(!b || bittab_size(b) == array_size(nodes));
  gn_id = genome_node_get_idstr(gn);
#endif
  gn_range = genome_node_get_range(gn);

  for (i = 0; i < array_size(nodes); i++) {
    node = *(GenomeNode**) array_get(nodes, i);
    assert(!str_cmp(gn_id, genome_node_get_idstr(node)));
    if (range_overlap(gn_range, genome_node_get_range(node))) {
      rval = true;
      if (b)
        bittab_set_bit(b, i);
      else
        break;
    }
  }
  return rval;
}

int genome_node_compare(GenomeNode **gn_a, GenomeNode **gn_b)
{
  return genome_node_cmp(*gn_a, *gn_b);
}

int genome_node_compare_with_data(GenomeNode **gn_a, GenomeNode **gn_b,
                                  UNUSED void *unused)
{
  return genome_node_cmp(*gn_a, *gn_b);
}

int genome_node_compare_delta(GenomeNode **gn_a, GenomeNode **gn_b,
                              void *delta)
{
  unsigned long *deltaptr = delta;
  assert(delta);
  return compare_genome_nodes_with_delta(*gn_a, *gn_b, *deltaptr);
}

void genome_node_delete(GenomeNode *gn)
{
  if (!gn) return;
  if (gn->reference_count) { gn->reference_count--; return; }
  assert(gn->c_class);
  if (gn->c_class->free) gn->c_class->free(gn);
  str_delete(gn->filename);
  dlist_delete(gn->children);
  ma_free(gn);
}

static int free_genome_node(GenomeNode *gn, UNUSED void *data,
                            UNUSED Error *err)
{
  genome_node_delete(gn);
  return 0;
}

void genome_node_rec_delete(GenomeNode *gn)
{
  int had_err;
  if (!gn) return;
  had_err = genome_node_traverse_children(gn, NULL, free_genome_node, true,
                                          NULL);
  assert(!had_err); /* cannot happen, free_genome_node() is sane */
}

void genome_nodes_sort(Array *nodes)
{
  qsort(array_get_space(nodes), array_size(nodes), sizeof (GenomeNode*),
        (Compare) genome_node_compare);
}

void genome_nodes_sort_stable(Array *nodes)
{
  msort(array_get_space(nodes), array_size(nodes), sizeof (GenomeNode*),
        (Compare) genome_node_compare);

}

bool genome_nodes_are_equal_sequence_regions(GenomeNode *gn_a, GenomeNode *gn_b)
{
  void *sr_a, *sr_b;

  sr_a = gn_a ? genome_node_cast(sequence_region_class(), gn_a) : NULL;
  sr_b = gn_b ? genome_node_cast(sequence_region_class(), gn_b) : NULL;

  if (sr_a && sr_b && !str_cmp(genome_node_get_seqid(gn_a),
                               genome_node_get_seqid(gn_b))) {
    return true;
  }
  return false;
}

bool genome_nodes_are_sorted(const Array *nodes)
{
  unsigned long i;
  assert(nodes);
  for (i = 1; i < array_size(nodes); i++) {
    if (genome_node_compare(array_get(nodes, i-1), array_get(nodes, i)) > 0)
      return false;
  }
  return true;
}
