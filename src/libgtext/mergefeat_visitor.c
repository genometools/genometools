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
#include "libgtcore/hashtable.h"
#include "libgtcore/undef.h"
#include "libgtcore/unused.h"
#include "libgtext/mergefeat_visitor.h"
#include "libgtext/genome_visitor_rep.h"

struct MergefeatVisitor {
  const GenomeVisitor parent_instance;
  GenomeNode *current_tree;
  Hashtable *ht; /* type -> previous node */
  Array *nodes_to_remove;
};

#define mergefeat_visitor_cast(GV)\
        genome_visitor_cast(mergefeat_visitor_class(), GV)

static void mergefeat_visitor_free(GenomeVisitor *gv)
{
  MergefeatVisitor *mergefeat_visitor = mergefeat_visitor_cast(gv);
  assert(mergefeat_visitor);
  hashtable_delete(mergefeat_visitor->ht);
  array_delete(mergefeat_visitor->nodes_to_remove);
}

static int mergefeat_in_children(GenomeNode *gn, void *data, UNUSED Error *err)
{
  MergefeatVisitor *v = (MergefeatVisitor*) data;
  GenomeFeature *previous_feature, *current_feature;
  Range previous_range, current_range;
  error_check(err);
  current_feature = genome_node_cast(genome_feature_class(), gn);
  assert(current_feature);
  if ((previous_feature = hashtable_get(v->ht, genome_feature_type_get_cstr(
                                  genome_feature_get_type(current_feature))))) {
    /* previous feature found -> check if merging is necessary */
    assert(genome_feature_get_type(previous_feature) ==
           genome_feature_get_type(current_feature));
    previous_range = genome_node_get_range((GenomeNode*) previous_feature);
    current_range = genome_node_get_range((GenomeNode*) current_feature);
    assert(range_compare(previous_range, current_range) <= 0); /* sorted */
    if (previous_range.end + 1 == current_range.start) {
      /* merge nodes */
      genome_feature_set_end(previous_feature, current_range.end);
      /* XXX: compute average score ? */
      genome_feature_set_score(previous_feature, UNDEF_DOUBLE);
      assert(!genome_node_number_of_children((GenomeNode*) current_feature));
      array_add(v->nodes_to_remove, current_feature);
    }
    /* remove previous feature */
    hashtable_remove(v->ht, (char*) genome_feature_type_get_cstr(
                     genome_feature_get_type(previous_feature)));
  }
  /* add current feature */
  hashtable_add(v->ht, (char*) genome_feature_type_get_cstr(
                genome_feature_get_type(current_feature)), current_feature);
  return 0;
}

static int mergefeat_if_necessary(GenomeNode *gn, void *data, Error *e)
{
  MergefeatVisitor *v = (MergefeatVisitor*) data;
  GenomeFeature *gf;
  error_check(e);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf);
  v->current_tree = gn;
  hashtable_reset(v->ht);
  return genome_node_traverse_direct_children(gn, v, mergefeat_in_children, e);
}

static int mergefeat_visitor_genome_feature(GenomeVisitor *gv,
                                            GenomeFeature *gf, Error *e)
{
  MergefeatVisitor *v;
  GenomeNode *leaf;
  unsigned long i;
  int had_err = 0;
  error_check(e);
  v = mergefeat_visitor_cast(gv);
  array_reset(v->nodes_to_remove);
  had_err = genome_node_traverse_children((GenomeNode*) gf, v,
                                          mergefeat_if_necessary, false, e);
  if (!had_err) {
    for (i = 0; i < array_size(v->nodes_to_remove); i++) {
      leaf = *(GenomeNode**) array_get(v->nodes_to_remove, i);
      genome_node_remove_leaf((GenomeNode*) gf, leaf);
      genome_node_delete(leaf);
    }
  }
  return had_err;
}

const GenomeVisitorClass* mergefeat_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (MergefeatVisitor),
                                          mergefeat_visitor_free,
                                          NULL,
                                          mergefeat_visitor_genome_feature,
                                          NULL };
  return &gvc;
}

GenomeVisitor* mergefeat_visitor_new(void)
{
  GenomeVisitor *gv = genome_visitor_create(mergefeat_visitor_class());
  MergefeatVisitor *mergefeat_visitor = mergefeat_visitor_cast(gv);
  mergefeat_visitor->ht = hashtable_new(HASH_STRING, NULL, NULL);
  mergefeat_visitor->nodes_to_remove = array_new(sizeof (GenomeNode*));
  return gv;
}
