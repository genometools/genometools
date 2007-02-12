/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdarg.h>
#include "fptr.h"
#include "genome_node_rep.h"
#include "hashtable.h"
#include "msort.h"
#include "xansi.h"

typedef struct {
  Genome_node_traverse_func func;
  void *data;
} Traverse_children_info;

static int compare_genome_node_type(Genome_node *gn_a, Genome_node *gn_b)
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

static int compare_genome_nodes(Genome_node *gn_a, Genome_node *gn_b)
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

void genome_node_class_init(Genome_node_class *gnc, size_t size, ...)
{
  va_list ap;
  Fptr func, meth, *mm;

  gnc->size = size;
  va_start(ap, size);
  while ((func = va_arg(ap, Fptr))) {
    meth = va_arg(ap, Fptr);
    assert(meth);
    if (func == (Fptr) genome_node_free) {
      mm  = (Fptr*) &gnc->free; *mm = meth;
    }
    else assert(0);
  }
  va_end(ap);
}

Genome_node* genome_node_create(const Genome_node_class *gnc,
                                const char *filename,
                                unsigned long line_number)
{
  Genome_node *gn;

  assert(gnc && gnc->size);

  gn                  = xmalloc(gnc->size);
  gn->c_class         = gnc;
  gn->filename        = filename;
  gn->line_number     = line_number;
  gn->children        = NULL; /* the children list is created on demand */
  gn->reference_count = 0;

  return gn;
}

void* genome_node_cast(const Genome_node_class *gnc, Genome_node *gn)
{
  assert(gnc && gn);
  if (gn->c_class == gnc)
    return gn;
  return NULL;
}

static void increase_reference_count(Genome_node *gn, /*@unused@*/ void *data)
{
  assert(gn);
  gn->reference_count++;
}

static Genome_node* genome_node_ref(Genome_node *gn)
{
  increase_reference_count(gn, NULL);
  return gn;
}

Genome_node* genome_node_rec_ref(Genome_node *gn)
{
  assert(gn);
  genome_node_traverse_children(gn, NULL, increase_reference_count, 1);
  return gn;
}

void genome_node_traverse_children(Genome_node *genome_node,
                                   void *data,
                                   Genome_node_traverse_func traverse,
                                   bool traverse_only_once)
{
  Array *node_stack, *list_of_children;
  Genome_node *gn, *child_feature;
  Dlistelem *dlistelem;
  unsigned long i;
  Hashtable *traversed_nodes = NULL;

  if (!genome_node) return;

  node_stack = array_new(sizeof(Genome_node*));
  list_of_children = array_new(sizeof(Genome_node*));
  array_add(node_stack, genome_node);

  if (traverse_only_once)
    traversed_nodes = hashtable_new(HASH_DIRECT, NULL, NULL);

  while (array_size(node_stack)) {
    gn = *(Genome_node**) array_pop(node_stack);
    array_set_size(list_of_children, 0);
    if (gn->children) {
      /* a backup of the children array is necessary if traverse() frees the
         node */
      for (dlistelem = dlist_first(gn->children); dlistelem != NULL;
           dlistelem = dlistelem_next(dlistelem)) {
        child_feature = (Genome_node*) dlistelem_get_data(dlistelem);
        array_add(list_of_children, child_feature);
      }
    }
    if (traverse) traverse(gn, data);
    /* we go backwards to traverse in order */
    for (i = array_size(list_of_children); i > 0; i--) {
      child_feature = *((Genome_node**) array_get(list_of_children, i-1));
      if (!traverse_only_once ||
          !hashtable_get(traversed_nodes, child_feature)) {
        /* feature has not been traversed or has to be traversed multiple
           times */
        array_add(node_stack, child_feature);
        if (traverse_only_once)
          hashtable_add(traversed_nodes, child_feature, child_feature);
      }
    }
  }

  if (traverse_only_once)
    hashtable_free(traversed_nodes);
  array_free(list_of_children);
  array_free(node_stack);
}

void genome_node_traverse_direct_children(Genome_node *gn,
                                          void *traverse_func_data,
                                          Genome_node_traverse_func traverse)
{
  Dlistelem *dlistelem;
  if (!gn || !traverse) return;
  if (gn->children) {
    for (dlistelem = dlist_first(gn->children); dlistelem != NULL;
         dlistelem = dlistelem_next(dlistelem)) {
      traverse((Genome_node*) dlistelem_get_data(dlistelem),
               traverse_func_data);
    }
  }
}

const char* genome_node_get_filename(const Genome_node *gn)
{
  assert(gn);
  return gn->filename;
}

unsigned long genome_node_get_line_number(const Genome_node *gn)
{
  assert(gn);
  return gn->line_number;
}

unsigned long genome_node_number_of_children(const Genome_node *gn)
{
  assert(gn);
  return dlist_size(gn->children);
}

Str* genome_node_get_seqid(Genome_node *gn)
{
  assert(gn && gn->c_class && gn->c_class->get_seqid);
  return gn->c_class->get_seqid(gn);
}

Str* genome_node_get_idstr(Genome_node *gn)
{
  assert(gn && gn->c_class && gn->c_class->get_idstr);
  return gn->c_class->get_idstr(gn);
}

unsigned long genome_node_get_start(Genome_node *gn)
{
  Range range = genome_node_get_range(gn);
  return range.start;
}

unsigned long genome_node_get_end(Genome_node *gn)
{
  Range range = genome_node_get_range(gn);
  return range.end;
}

Range genome_node_get_range(Genome_node *gn)
{
  assert(gn && gn->c_class && gn->c_class->get_range);
  return gn->c_class->get_range(gn);
}

void genome_node_set_range(Genome_node *gn, Range range)
{
  assert(gn && gn->c_class && gn->c_class->set_range);
  gn->c_class->set_range(gn, range);
}

void genome_node_set_seqid(Genome_node *gn, Str *seqid)
{
  assert(gn && gn->c_class && gn->c_class->set_seqid && seqid);
  gn->c_class->set_seqid(gn, seqid);
}

void genome_node_set_source(Genome_node *gn, Str *source)
{
  assert(gn && gn->c_class && gn->c_class->set_source && source);
  gn->c_class->set_source(gn, source);
}

void genome_node_set_phase(Genome_node *gn, Phase p)
{
  assert(gn && gn->c_class && gn->c_class->set_phase);
  gn->c_class->set_phase(gn, p);
}

void genome_node_accept(Genome_node *gn, Genome_visitor *gv, Log *l)
{
  assert(gn && gv && gn->c_class && gn->c_class->accept);
  gn->c_class->accept(gn, gv, l);
}

void genome_node_is_part_of_genome_node(Genome_node *parent,
                                        Genome_node *child)
{
  assert(parent && child);
  /* create children list  on demand */
  if (!parent->children)
    parent->children = dlist_new((Compare) compare_genome_nodes);
  dlist_add(parent->children, child); /* XXX: check for circles */
}

static void remove_leaf(Genome_node *node, void *data)
{
  Dlistelem *dlistelem;
  Genome_node *child, *leaf = (Genome_node*) data;
  if (node != leaf && node->children) {
    for (dlistelem = dlist_first(node->children); dlistelem != NULL;
         dlistelem = dlistelem_next(dlistelem)) {
      child = (Genome_node*) dlistelem_get_data(dlistelem);
      if (child == leaf) {
        dlist_remove(node->children, dlistelem);
        break;
      }
    }
  }
}

void genome_node_remove_leaf(Genome_node *tree, Genome_node *leafn)
{
  assert(tree && leafn);
  assert(!genome_node_number_of_children(leafn));
  genome_node_traverse_children(tree, leafn, remove_leaf, 1);
}

bool genome_node_has_children(Genome_node *gn)
{
  assert(gn);
  if (!gn->children || dlist_size(gn->children) == 0)
    return false;
  return true;
}

bool genome_node_direct_children_do_not_overlap(Genome_node *gn)
{
  Array *children_ranges = array_new(sizeof(Range));
  Dlistelem *dlistelem;
  Range range;
  bool rval;

  assert(gn);

  if (!gn->children)
    return 1;

  /* get children ranges */
  if (gn->children) {
    for (dlistelem = dlist_first(gn->children); dlistelem != NULL;
         dlistelem = dlistelem_next(dlistelem)) {
      range = genome_node_get_range((Genome_node*)
                                    dlistelem_get_data(dlistelem));
      array_add(children_ranges, range);
    }
  }

  ranges_sort(children_ranges);
  assert(ranges_are_sorted(children_ranges));
  rval = ranges_do_not_overlap(children_ranges);

  array_free(children_ranges);

  return rval;
}

bool genome_node_tree_is_sorted(Genome_node **buffer,
                                        Genome_node *current_node)
{
  assert(buffer && current_node);

  if (*buffer) {
    /* the last node is not larger than the current one */
    if (genome_node_compare(buffer, &current_node) == 1)
      return false;
    genome_node_free(*buffer);
  }
  *buffer = genome_node_ref(current_node);
  return true;
}

bool genome_node_overlaps_nodes(Genome_node *gn, Array *nodes)
{
  return genome_node_overlaps_nodes_mark(gn, nodes, NULL);
}

bool genome_node_overlaps_nodes_mark(Genome_node *gn, Array *nodes,
                                             Bittab *b)
{
  unsigned long i;
  Genome_node *node;
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
    node = *(Genome_node**) array_get(nodes, i);
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

int genome_node_compare(Genome_node **gn_a, Genome_node **gn_b)
{
  return compare_genome_nodes(*gn_a, *gn_b);
}

void genome_node_free(Genome_node *gn)
{
  if (!gn) return;
  if (gn->reference_count) { gn->reference_count--; return; }
  assert(gn->c_class);
  if (gn->c_class->free) gn->c_class->free(gn);
  dlist_free(gn->children);
  free(gn);
}

static void free_genome_node(Genome_node *gn, /*@unused@*/ void *data)
{
  genome_node_free(gn);
}

void genome_node_rec_free(Genome_node *gn)
{
  if (!gn) return;
  genome_node_traverse_children(gn, NULL, free_genome_node, 1);
}

void genome_nodes_sort(Array *nodes)
{
  qsort(array_get_space(nodes), array_size(nodes), sizeof(Genome_node*),
        (Compare) genome_node_compare);
}

void genome_nodes_sort_stable(Array *nodes)
{
  msort(array_get_space(nodes), array_size(nodes), sizeof(Genome_node*),
        (Compare) genome_node_compare);

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
