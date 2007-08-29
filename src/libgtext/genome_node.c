/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdarg.h>
#include "libgtcore/hashtable.h"
#include "libgtcore/msort.h"
#include "libgtcore/queue.h"
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

static int compare_genome_nodes(GenomeNode *gn_a, GenomeNode *gn_b)
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

GenomeNode* genome_node_create(const GenomeNodeClass *gnc,
                               const char *filename,
                               unsigned long line_number, Env *env)
{
  GenomeNode *gn;
  assert(gnc && gnc->size && filename);
  gn                  = env_ma_malloc(env, gnc->size);
  gn->c_class         = gnc;
  gn->filename        = filename;
  gn->line_number     = line_number;
  gn->children        = NULL; /* the children list is created on demand */
  gn->reference_count = 0;
  gn->info            = 0;
  genome_node_info_set_tree_status(&gn->info, GENOME_NODE_IS_TREE);
  return gn;
}

void* genome_node_cast(const GenomeNodeClass *gnc, GenomeNode *gn)
{
  assert(gnc && gn);
  if (gn->c_class == gnc)
    return gn;
  return NULL;
}

static int increase_reference_count(GenomeNode *gn, /*@unused@*/ void *data,
                                    Env *env)
{
  env_error_check(env);
  assert(gn);
  gn->reference_count++;
  return 0;
}

static GenomeNode* genome_node_ref(GenomeNode *gn)
{
  int had_err;
  had_err = increase_reference_count(gn, NULL, NULL);
  assert(!had_err); /* cannot happen, increase_reference_count() is sane */
  return gn;
}

GenomeNode* genome_node_rec_ref(GenomeNode *gn, Env *env)
{
  int had_err;
  assert(gn);
  had_err = genome_node_traverse_children(gn, NULL, increase_reference_count,
                                          true, env);
  assert(!had_err); /* cannot happen, increase_reference_count() is sane */
  return gn;
}

int genome_node_traverse_children_generic(GenomeNode *genome_node,
                                          void *data,
                                          GenomeNodeTraverseFunc traverse,
                                          bool traverse_only_once,
                                          bool depth_first, Env *env)
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
    node_stack = array_new(sizeof (GenomeNode*), env);
    array_add(node_stack, genome_node, env);
  }
  else {
    node_queue = queue_new(sizeof (GenomeNode*), env);
    queue_add(node_queue, genome_node, env);
  }
  list_of_children = array_new(sizeof (GenomeNode*), env);

  if (traverse_only_once)
    traversed_nodes = hashtable_new(HASH_DIRECT, NULL, NULL, env);

  while ((depth_first ? array_size(node_stack): queue_size(node_queue))) {
    if (depth_first)
      gn = *(GenomeNode**) array_pop(node_stack);
    else
      gn = *(GenomeNode**) queue_get(node_queue);
    array_reset(list_of_children);
    if (gn->children) {
      /* a backup of the children array is necessary if traverse() frees the
         node */
      for (dlistelem = dlist_first(gn->children); dlistelem != NULL;
           dlistelem = dlistelem_next(dlistelem)) {
        child_feature = (GenomeNode*) dlistelem_get_data(dlistelem);
        array_add(list_of_children, child_feature, env);
      }
    }
    /* store the implications of <gn> to the tree status of <genome_node> */
    if (genome_node_info_multiple_parents(&gn->info))
      has_node_with_multiple_parents = true;
    /* call traverse function */
    if (traverse) {
      had_err = traverse(gn, data, env);
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
          array_add(node_stack, child_feature, env);
        else
          queue_add(node_queue, child_feature, env);
        if (traverse_only_once)
          hashtable_add(traversed_nodes, child_feature, child_feature, env);
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
  genome_node_delete(gn_ref, env);
  if (traverse_only_once)
    hashtable_delete(traversed_nodes, env);
  array_delete(list_of_children, env);
  array_delete(node_stack, env);
  queue_delete(node_queue, env);

  return had_err;
}

int genome_node_traverse_children(GenomeNode *genome_node, void *data,
                                  GenomeNodeTraverseFunc traverse,
                                  bool traverse_only_once, Env *env)
{
  return genome_node_traverse_children_generic(genome_node, data, traverse,
                                               traverse_only_once, true, env);
}

int genome_node_traverse_children_breadth(GenomeNode *genome_node, void *data,
                                          GenomeNodeTraverseFunc traverse,
                                          bool traverse_only_once,
                                          Env *env)
{
  return genome_node_traverse_children_generic(genome_node, data, traverse,
                                               traverse_only_once, false, env);
}

int genome_node_traverse_direct_children(GenomeNode *gn,
                                         void *traverse_func_data,
                                         GenomeNodeTraverseFunc traverse,
                                         Env *env)
{
  Dlistelem *dlistelem;
  int had_err = 0;
  env_error_check(env);
  if (!gn || !traverse)
    return 0;
  if (gn->children) {
    for (dlistelem = dlist_first(gn->children); dlistelem != NULL;
         dlistelem = dlistelem_next(dlistelem)) {
      had_err = traverse((GenomeNode*) dlistelem_get_data(dlistelem),
                          traverse_func_data, env);
      if (had_err)
        break;
    }
  }
  return had_err;
}

const char* genome_node_get_filename(const GenomeNode *gn)
{
  assert(gn);
  return gn->filename;
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
  assert(gn && gn->c_class && gn->c_class->get_seqid);
  return gn->c_class->get_seqid(gn);
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

void genome_node_set_source(GenomeNode *gn, Str *source)
{
  assert(gn && gn->c_class && gn->c_class->set_source && source);
  gn->c_class->set_source(gn, source);
}

void genome_node_set_phase(GenomeNode *gn, Phase p)
{
  assert(gn && gn->c_class && gn->c_class->set_phase);
  gn->c_class->set_phase(gn, p);
}

int genome_node_accept(GenomeNode *gn, GenomeVisitor *gv, Env *env)
{
  env_error_check(env);
  assert(gn && gv && gn->c_class && gn->c_class->accept);
  return gn->c_class->accept(gn, gv, env);
}

void genome_node_is_part_of_genome_node(GenomeNode *parent, GenomeNode *child,
                                        Env *env)
{
  assert(parent && child);
  /* create children list  on demand */
  if (!parent->children)
    parent->children = dlist_new((Compare) compare_genome_nodes, env);
  dlist_add(parent->children, child, env); /* XXX: check for circles */
  /* update tree status of <parent> */
  genome_node_info_set_tree_status(&parent->info,
                                   GENOME_NODE_TREE_STATUS_UNDETERMINED);
  /* update parent info of <child> */
  genome_node_info_add_parent(&child->info);
}

static int remove_leaf(GenomeNode *node, void *data, Env *env)
{
  Dlistelem *dlistelem;
  GenomeNode *child, *leaf = (GenomeNode*) data;
  env_error_check(env);
  if (node != leaf && node->children) {
    for (dlistelem = dlist_first(node->children); dlistelem != NULL;
         dlistelem = dlistelem_next(dlistelem)) {
      child = (GenomeNode*) dlistelem_get_data(dlistelem);
      if (child == leaf) {
        dlist_remove(node->children, dlistelem, env);
        break;
      }
    }
  }
  return 0;
}

void genome_node_remove_leaf(GenomeNode *tree, GenomeNode *leafn, Env *env)
{
  int had_err;
  assert(tree && leafn);
  assert(!genome_node_number_of_children(leafn));
  had_err = genome_node_traverse_children(tree, leafn, remove_leaf, true, env);
  assert(!had_err); /* cannot happen, remove_leaf() is sane */
}

bool genome_node_has_children(GenomeNode *gn)
{
  assert(gn);
  if (!gn->children || dlist_size(gn->children) == 0)
    return false;
  return true;
}

bool genome_node_direct_children_do_not_overlap_generic(GenomeNode *parent,
                                                        GenomeNode *child,
                                                        Env *env)
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
  children_ranges = array_new(sizeof (Range), env);
  assert(parent->children);
  for (dlistelem = dlist_first(parent->children); dlistelem != NULL;
       dlistelem = dlistelem_next(dlistelem)) {
    if (!gf ||
        ((child_gf = genome_node_cast(genome_feature_class(),
                                      dlistelem_get_data(dlistelem))) &&
         genome_feature_get_type(gf) == genome_feature_get_type(child_gf))) {
      range = genome_node_get_range((GenomeNode*)
                                    dlistelem_get_data(dlistelem));
      array_add(children_ranges, range, env);
    }
  }

  ranges_sort(children_ranges);
  assert(ranges_are_sorted(children_ranges));
  rval = ranges_do_not_overlap(children_ranges);

  array_delete(children_ranges, env);

  return rval;
}

bool genome_node_direct_children_do_not_overlap(GenomeNode *gn, Env *env)
{
  return genome_node_direct_children_do_not_overlap_generic(gn, NULL, env);
}

bool genome_node_direct_children_do_not_overlap_st(GenomeNode *parent,
                                                   GenomeNode *child, Env *env)
{
  return genome_node_direct_children_do_not_overlap_generic(parent, child, env);
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

bool genome_node_tree_is_sorted(GenomeNode **buffer, GenomeNode *current_node,
                                Env *env)
{
  assert(buffer && current_node);

  if (*buffer) {
    /* the last node is not larger than the current one */
    if (genome_node_compare(buffer, &current_node) == 1)
      return false;
    genome_node_delete(*buffer, env);
  }
  *buffer = genome_node_ref(current_node);
  return true;
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
  return compare_genome_nodes(*gn_a, *gn_b);
}

int genome_node_compare_with_data(GenomeNode **gn_a, GenomeNode **gn_b,
                                  void *unused)
{
  return compare_genome_nodes(*gn_a, *gn_b);
}

int genome_node_compare_delta(GenomeNode **gn_a, GenomeNode **gn_b,
                              void *delta)
{
  unsigned long *deltaptr = delta;
  assert(delta);
  return compare_genome_nodes_with_delta(*gn_a, *gn_b, *deltaptr);
}

void genome_node_delete(GenomeNode *gn, Env *env)
{
  if (!gn) return;
  if (gn->reference_count) { gn->reference_count--; return; }
  assert(gn->c_class);
  if (gn->c_class->free) gn->c_class->free(gn, env);
  dlist_delete(gn->children, env);
  env_ma_free(gn, env);
}

static int free_genome_node(GenomeNode *gn, /*@unused@*/ void *data, Env *env)
{
  genome_node_delete(gn, env);
  return 0;
}

void genome_node_rec_delete(GenomeNode *gn, Env *env)
{
  int had_err;
  if (!gn) return;
  had_err = genome_node_traverse_children(gn, NULL, free_genome_node, true,
                                          env);
  assert(!had_err); /* cannot happen, free_genome_node() is sane */
}

void genome_nodes_sort(Array *nodes)
{
  qsort(array_get_space(nodes), array_size(nodes), sizeof (GenomeNode*),
        (Compare) genome_node_compare);
}

void genome_nodes_sort_stable(Array *nodes, Env *env)
{
  msort(array_get_space(nodes), array_size(nodes), sizeof (GenomeNode*),
        (Compare) genome_node_compare, env);

}

bool genome_nodes_are_sorted(const Array *nodes)
{
  unsigned long i;
  assert(nodes);
  for (i = 1; i < array_size(nodes); i++) {
    if (genome_node_compare(array_get(nodes, i-1), array_get(nodes, i)) == 1)
      return false;
  }
  return true;
}
