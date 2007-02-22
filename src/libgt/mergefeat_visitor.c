/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "mergefeat_visitor.h"
#include "genome_visitor_rep.h"
#include "hashtable.h"
#include "undef.h"

struct MergefeatVisitor {
  const GenomeVisitor parent_instance;
  GenomeNode *current_tree;
  Hashtable *ht; /* type -> previous node */
};

#define mergefeat_visitor_cast(GV)\
        genome_visitor_cast(mergefeat_visitor_class(), GV)

static void mergefeat_visitor_free(GenomeVisitor *gv, Env *env)
{
  MergefeatVisitor *mergefeat_visitor = mergefeat_visitor_cast(gv);
  assert(mergefeat_visitor);
  hashtable_delete(mergefeat_visitor->ht, env);
}

static int mergefeat_in_children(GenomeNode *gn, void *data, Env *env)
{
  MergefeatVisitor *v = (MergefeatVisitor*) data;
  GenomeFeature *previous_feature, *current_feature;
  Range previous_range, current_range;
  env_error_check(env);
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
      genome_feature_set_score(previous_feature, UNDEFDOUBLE);
      assert(!genome_node_number_of_children((GenomeNode*) current_feature));
      /* XXX: */
      genome_node_remove_leaf(v->current_tree, (GenomeNode*) current_feature,
                              env);
#if 0
      genome_node_delete((GenomeNode*) current_feature);
#endif
    }
    /* remove previous feature */
    hashtable_remove(v->ht, (char*) genome_feature_type_get_cstr(
                     genome_feature_get_type(previous_feature)), env);
  }
  /* add current feature */
  hashtable_add(v->ht, (char*) genome_feature_type_get_cstr(
                genome_feature_get_type(current_feature)), current_feature,
                env);
  return 0;
}

static int mergefeat_if_necessary(GenomeNode *gn, void *data, Env *env)
{
  MergefeatVisitor *v = (MergefeatVisitor*) data;
  GenomeFeature *gf;
  env_error_check(env);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf);
  v->current_tree = gn;
  hashtable_reset(v->ht, env);
  return genome_node_traverse_direct_children(gn, v, mergefeat_in_children,
                                              env);
}

static int mergefeat_visitor_genome_feature(GenomeVisitor *gv,
                                            GenomeFeature *gf, Env *env)
{
  MergefeatVisitor *v;
  env_error_check(env);
  v = mergefeat_visitor_cast(gv);
  return genome_node_traverse_children((GenomeNode*) gf, v,
                                       mergefeat_if_necessary, false, env);
}

const GenomeVisitorClass* mergefeat_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (MergefeatVisitor),
                                          mergefeat_visitor_free,
                                          NULL,
                                          mergefeat_visitor_genome_feature,
                                          NULL,
                                          NULL };
  return &gvc;
}

GenomeVisitor* mergefeat_visitor_new(Env *env)
{
  GenomeVisitor *gv = genome_visitor_create(mergefeat_visitor_class());
  MergefeatVisitor *mergefeat_visitor = mergefeat_visitor_cast(gv);
  mergefeat_visitor->ht = hashtable_new(HASH_STRING, NULL, NULL, env);
  return gv;
}
