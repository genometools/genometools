/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "addintrons_visitor.h"
#include "genome_visitor_rep.h"

struct AddIntronsVisitor {
  const GenomeVisitor parent_instance;
  GenomeFeature *previous_feature;
};

#define addintrons_visitor_cast(GV)\
        genome_visitor_cast(addintrons_visitor_class(), GV)

static void addintrons_visitor_free(GenomeVisitor *gv, Env *env)
{
  AddIntronsVisitor *addintrons_visitor = addintrons_visitor_cast(gv);
  assert(addintrons_visitor);
}

static int addintrons_in_children(GenomeNode *gn, void *data, Env *env)
{
  AddIntronsVisitor *v = (AddIntronsVisitor*) data;
  GenomeFeature *current_feature;
  env_error_check(env);
  current_feature = genome_node_cast(genome_feature_class(), gn);
  assert(current_feature);
  if (v->previous_feature) {
    /* XXX: add introns here  (previous_feature -> previous_gene_feature) */
  }
  v->previous_feature = current_feature;
  return 0;
}

static int addintrons_if_necessary(GenomeNode *gn, void *data, Env *env)
{
  AddIntronsVisitor *v = (AddIntronsVisitor*) data;
  GenomeFeature *gf;
  env_error_check(env);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf);
  return genome_node_traverse_direct_children(gn, v, addintrons_in_children,
                                              env);
}

static int addintrons_visitor_genome_feature(GenomeVisitor *gv,
                                            GenomeFeature *gf, Env *env)
{
  AddIntronsVisitor *v;
  env_error_check(env);
  v = addintrons_visitor_cast(gv);
  v->previous_feature = NULL;
  return genome_node_traverse_children((GenomeNode*) gf, v,
                                       addintrons_if_necessary, false, env);
}

const GenomeVisitorClass* addintrons_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (AddIntronsVisitor),
                                          addintrons_visitor_free,
                                          NULL,
                                          addintrons_visitor_genome_feature,
                                          NULL,
                                          NULL };
  return &gvc;
}

GenomeVisitor* addintrons_visitor_new(Env *env)
{
  return genome_visitor_create(addintrons_visitor_class(), env);
}
