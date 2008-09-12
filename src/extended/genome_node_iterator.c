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

#include "core/assert.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/genome_node_iterator.h"
#include "extended/feature_node_rep.h"

struct GtGenomeNodeIterator {
  GtGenomeNode *gn;
  GtArray *feature_stack;
  bool direct;
};

static GtGenomeNodeIterator* genome_node_iterator_new_base(GtGenomeNode *gn)
{
  GtGenomeNodeIterator *gni;
  gt_assert(gn);
  gni = gt_malloc(sizeof *gni);
  gni->gn = gt_genome_node_rec_ref(gn);
  gni->feature_stack = gt_array_new(sizeof (GtFeatureNode*));
  return gni;
}

GtGenomeNodeIterator* gt_genome_node_iterator_new(GtGenomeNode *gn)
{
  GtFeatureNode *fn;
  GtGenomeNodeIterator *gni;
  GtGenomeNode *child_feature;
  GtDlistelem *dlistelem;
  gt_assert(gn);
  fn = gt_feature_node_cast(gn); /* XXX */
  gni = genome_node_iterator_new_base(gn);
  if (gt_genome_node_cast(gt_feature_node_class(), gn) &&
      gt_feature_node_is_pseudo((GtFeatureNode*) gn)) {
    /* add the children backwards to traverse in order */
    for (dlistelem = gt_dlist_last(fn->children); dlistelem != NULL;
         dlistelem = gt_dlistelem_previous(dlistelem)) {
      child_feature = (GtGenomeNode*) gt_dlistelem_get_data(dlistelem);
      gt_array_add(gni->feature_stack, child_feature);
    }
  }
  else
    gt_array_add(gni->feature_stack, gni->gn);
  gt_assert(gt_array_size(gni->feature_stack));
  gni->direct = false;
  return gni;
}

static void add_children_to_stack(GtArray *feature_stack, GtFeatureNode *fn)
{
  GtGenomeNode *child;
  GtDlistelem *dlistelem;
  gt_assert(feature_stack && fn && fn->children);
  /* add the children backwards to traverse in order */
  for (dlistelem = gt_dlist_last(fn->children); dlistelem != NULL;
       dlistelem = gt_dlistelem_previous(dlistelem)) {
    child = gt_dlistelem_get_data(dlistelem);
    gt_array_add(feature_stack, child);
  }
}

GtGenomeNodeIterator* gt_genome_node_iterator_new_direct(GtGenomeNode *gn)
{
  GtFeatureNode *fn;
  GtGenomeNodeIterator *gni;
  gt_assert(gn);
  gni = genome_node_iterator_new_base(gn);
  fn = gt_feature_node_cast(gn); /* XXX */
  if (fn->children)
    add_children_to_stack(gni->feature_stack, fn);
  gni->direct = true;
  return gni;
}

GtGenomeNode* gt_genome_node_iterator_next(GtGenomeNodeIterator *gni)
{
  GtFeatureNode *fn;
  gt_assert(gni);
  if (!gt_array_size(gni->feature_stack))
    return NULL;
  /* pop */
  fn = *(GtFeatureNode**) gt_array_pop(gni->feature_stack);
  /* push children on stack */
  if (!gni->direct && fn->children)
    add_children_to_stack(gni->feature_stack, fn);
  return (GtGenomeNode*) fn; /* XXX */
}

int gt_genome_node_iterator_example(GT_UNUSED GtError *err)
{
  GtGenomeNodeIterator *gni;
  GtGenomeNode *gn, *node;
  gn = gt_feature_node_new_standard_gene();

  /* an example genome node iterator use case */
  gni = gt_genome_node_iterator_new(gn);
  while ((node = gt_genome_node_iterator_next(gni))) {
    /* do something with <node> */
  }
  gt_genome_node_iterator_delete(gni);

  gt_genome_node_rec_delete(gn);
  return 0;
}

void gt_genome_node_iterator_delete(GtGenomeNodeIterator *gni)
{
  if (!gni) return;
  gt_genome_node_rec_delete(gni->gn);
  gt_array_delete(gni->feature_stack);
  gt_free(gni);
}
