/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/genome_node_iterator.h"
#include "extended/genome_node_rep.h"
#include "extended/feature_type_factory_builtin.h"

struct GT_GenomeNodeIterator {
  GT_GenomeNode *gn;
  GT_Array *node_stack;
  bool direct;
};

static GT_GenomeNodeIterator* gt_genome_node_iterator_new_base(GT_GenomeNode *gn)
{
  GT_GenomeNodeIterator *gni;
  assert(gn);
  gni = ma_malloc(sizeof *gni);
  gni->gn = gt_genome_node_rec_ref(gn);
  gni->node_stack = gt_array_new(sizeof (GT_GenomeNode*));
  return gni;
}

GT_GenomeNodeIterator* gt_genome_node_iterator_new(GT_GenomeNode *gn)
{
  GT_GenomeNodeIterator *gni;
  GT_GenomeNode *child_feature;
  Dlistelem *dlistelem;
  assert(gn);
  gni = gt_genome_node_iterator_new_base(gn);
  if (gt_genome_node_cast(gt_genome_feature_class(), gn) &&
      gt_genome_feature_is_pseudo((GT_GenomeFeature*) gn)) {
    /* add the children backwards to traverse in order */
    for (dlistelem = dlist_last(gn->children); dlistelem != NULL;
         dlistelem = dlistelem_previous(dlistelem)) {
      child_feature = (GT_GenomeNode*) dlistelem_get_data(dlistelem);
      gt_array_add(gni->node_stack, child_feature);
    }
  }
  else
    gt_array_add(gni->node_stack, gni->gn);
  assert(gt_array_size(gni->node_stack));
  gni->direct = false;
  return gni;
}

static void add_children_to_stack(GT_Array *node_stack, GT_GenomeNode *gn)
{
  GT_GenomeNode *child;
  Dlistelem *dlistelem;
  assert(node_stack && gn && gn->children);
  /* add the children backwards to traverse in order */
  for (dlistelem = dlist_last(gn->children); dlistelem != NULL;
       dlistelem = dlistelem_previous(dlistelem)) {
    child = dlistelem_get_data(dlistelem);
    gt_array_add(node_stack, child);
  }
}

GT_GenomeNodeIterator* gt_genome_node_iterator_new_direct(GT_GenomeNode *gn)
{
  GT_GenomeNodeIterator *gni;
  assert(gn);
  gni = gt_genome_node_iterator_new_base(gn);
  if (gn->children)
    add_children_to_stack(gni->node_stack, gn);
  gni->direct = true;
  return gni;
}

GT_GenomeNode* gt_genome_node_iterator_next(GT_GenomeNodeIterator *gni)
{
  GT_GenomeNode *gn;
  assert(gni);
  if (!gt_array_size(gni->node_stack))
    return NULL;
  /* pop */
  gn = *(GT_GenomeNode**) gt_array_pop(gni->node_stack);
  /* push children on stack */
  if (!gni->direct && gn->children)
    add_children_to_stack(gni->node_stack, gn);
  return gn;
}

int gt_genome_node_iterator_example(GT_UNUSED GT_Error *err)
{
  GT_FeatureTypeFactory *feature_type_factory;
  GT_GenomeNodeIterator *gni;
  GT_GenomeNode *gn, *node;
  feature_type_factory = gt_feature_type_factory_builtin_new();
  gn = gt_genome_feature_new_standard_gene(feature_type_factory);

  /* an example genome node iterator use case */
  gni = gt_genome_node_iterator_new(gn);
  while ((node = gt_genome_node_iterator_next(gni))) {
    /* do something with <node> */
  }
  gt_genome_node_iterator_delete(gni);

  gt_genome_node_rec_delete(gn);
  gt_feature_type_factory_delete(feature_type_factory);
  return 0;
}

void gt_genome_node_iterator_delete(GT_GenomeNodeIterator *gni)
{
  if (!gni) return;
  gt_genome_node_rec_delete(gni->gn);
  gt_array_delete(gni->node_stack);
  ma_free(gni);
}
