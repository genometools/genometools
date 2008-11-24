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

#include "core/assert_api.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/feature_node_rep.h"

struct GtFeatureNodeIterator {
  GtFeatureNode *fn;
  GtArray *feature_stack;
  bool direct;
};

static GtFeatureNodeIterator* feature_node_iterator_new_base(const GtFeatureNode
                                                             *fn)
{
  GtFeatureNodeIterator *fni;
  gt_assert(fn);
  fni = gt_malloc(sizeof *fni);
  fni->fn = (GtFeatureNode*) gt_genome_node_ref((GtGenomeNode*) fn);
  fni->feature_stack = gt_array_new(sizeof (GtFeatureNode*));
  return fni;
}

GtFeatureNodeIterator* gt_feature_node_iterator_new(const GtFeatureNode *fn)
{
  GtFeatureNodeIterator *fni;
  GtFeatureNode *child_feature;
  GtDlistelem *dlistelem;
  gt_assert(fn);
  fni = feature_node_iterator_new_base(fn);
  if (gt_feature_node_is_pseudo((GtFeatureNode*) fn)) {
    /* add the children backwards to traverse in order */
    for (dlistelem = gt_dlist_last(fn->children); dlistelem != NULL;
         dlistelem = gt_dlistelem_previous(dlistelem)) {
      child_feature = (GtFeatureNode*) gt_dlistelem_get_data(dlistelem);
      gt_array_add(fni->feature_stack, child_feature);
    }
  }
  else
    gt_array_add(fni->feature_stack, fni->fn);
  gt_assert(gt_array_size(fni->feature_stack));
  fni->direct = false;
  return fni;
}

static void add_children_to_stack(GtArray *feature_stack,
                                  const GtFeatureNode *fn)
{
  GtFeatureNode *child;
  GtDlistelem *dlistelem;
  gt_assert(feature_stack && fn && fn->children);
  /* add the children backwards to traverse in order */
  for (dlistelem = gt_dlist_last(fn->children); dlistelem != NULL;
       dlistelem = gt_dlistelem_previous(dlistelem)) {
    child = gt_dlistelem_get_data(dlistelem);
    gt_array_add(feature_stack, child);
  }
}

GtFeatureNodeIterator* gt_feature_node_iterator_new_direct(const GtFeatureNode
                                                           *fn)
{
  GtFeatureNodeIterator *fni;
  gt_assert(fn);
  fni = feature_node_iterator_new_base(fn);
  if (fn->children)
    add_children_to_stack(fni->feature_stack, fn);
  fni->direct = true;
  return fni;
}

GtFeatureNode* gt_feature_node_iterator_next(GtFeatureNodeIterator *fni)
{
  GtFeatureNode *fn;
  gt_assert(fni);
  if (!gt_array_size(fni->feature_stack))
    return NULL;
  /* pop */
  fn = *(GtFeatureNode**) gt_array_pop(fni->feature_stack);
  /* push children on stack */
  if (!fni->direct && fn->children)
    add_children_to_stack(fni->feature_stack, fn);
  return fn;
}

int gt_feature_node_iterator_example(GT_UNUSED GtError *err)
{
  GtFeatureNodeIterator *fni;
  GtFeatureNode *fn, *node;
  fn = (GtFeatureNode*) gt_feature_node_new_standard_gene();

  /* an example genome node iterator use case */
  fni = gt_feature_node_iterator_new(fn);
  while ((node = gt_feature_node_iterator_next(fni))) {
    /* do something with <node> */
  }
  gt_feature_node_iterator_delete(fni);

  gt_genome_node_delete((GtGenomeNode*) fn);
  return 0;
}

void gt_feature_node_iterator_delete(GtFeatureNodeIterator *fni)
{
  if (!fni) return;
  gt_genome_node_delete((GtGenomeNode*) fni->fn);
  gt_array_delete(fni->feature_stack);
  gt_free(fni);
}
