/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/ma.h"
#include "libgtext/genome_node_iterator.h"
#include "libgtext/genome_node_rep.h"

struct GenomeNodeIterator {
  Array *node_stack;
};

GenomeNodeIterator* genome_node_iterator_new(GenomeNode *gn, Env *env)
{
  GenomeNodeIterator *gni;
  env_error_check(env);
  assert(gn);
  gni = ma_malloc(sizeof *gni);
  gni->node_stack = array_new(sizeof (GenomeNode*), env);
  array_add(gni->node_stack, gn, env);
  return gni;
}

GenomeNode* genome_node_iterator_next(GenomeNodeIterator *gni, Env *env)
{
  GenomeNode *gn, *child;
  Dlistelem *dlistelem;
  assert(gni);
  if (!array_size(gni->node_stack))
    return NULL;
  /* pop */
  gn = *(GenomeNode**) array_pop(gni->node_stack);
  /* push children on stack */
  if (gn->children) {
    /* add the children backwards to traverse in order */
    for (dlistelem = dlist_last(gn->children); dlistelem != NULL;
         dlistelem = dlistelem_previous(dlistelem)) {
      child = dlistelem_get_data(dlistelem);
      array_add(gni->node_stack, child, env);
    }
  }
  return gn;
}

int genome_node_iterator_example(Env *env)
{
  GenomeNodeIterator *gni;
  GenomeNode *gn, *node;
  env_error_check(env);
  gn = genome_feature_new_standard_gene(env);

  /* an example genome node iterator use case */
  gni = genome_node_iterator_new(gn, env);
  while ((node = genome_node_iterator_next(gni, env))) {
    /* do something with <node> */
  }
  genome_node_iterator_delete(gni, env);

  genome_node_rec_delete(gn, env);
  return 0;
}

void genome_node_iterator_delete(GenomeNodeIterator *gni, Env *env)
{
  if (!gni) return;
  array_delete(gni->node_stack, env);
  ma_free(gni);
}
