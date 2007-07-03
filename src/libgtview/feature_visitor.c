/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <gtcore.h>
#include <libgtview/feature_visitor.h>
#include <libgtview/feature_index.h>
#include <libgtext/genome_visitor_rep.h>
#include <libgtext/sequence_region.h>

struct FeatureVisitor {
  const GenomeVisitor parent_instance;
        FeatureIndex *features;
};

#define feature_visitor_cast(GV)\
        genome_visitor_cast(feature_visitor_class(), GV)

static void feature_visitor_free(GenomeVisitor *gv,
                                 Env *env)
{
  FeatureVisitor *feature_visitor = feature_visitor_cast(gv);
  assert(feature_visitor);
}

static int feature_visitor_genome_feature(GenomeVisitor *gv,
                                          GenomeFeature *gf,
                                          Env *env)
{
  FeatureVisitor *v = feature_visitor_cast(gv);
  env_error_check(env);

  feature_index_add_genome_feature(v->features, gf, env);

  return 0;
}

static int feature_visitor_sequence_region(GenomeVisitor *gv,
                                           SequenceRegion *sr,
                                           Env *env)
{
  FeatureVisitor *v = feature_visitor_cast(gv);
  env_error_check(env);
  return feature_index_add_sequence_region(v->features, sr, env);
}

const GenomeVisitorClass* feature_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (FeatureVisitor),
                                          feature_visitor_free,
                                          NULL,
                                          feature_visitor_genome_feature,
                                          feature_visitor_sequence_region,
                                          NULL };
  return &gvc;
}

GenomeVisitor* feature_visitor_new(FeatureIndex *fi,
                                   Env *env)
{
  GenomeVisitor *gv;
  FeatureVisitor *feature_visitor;
  env_error_check(env);
  assert(fi);
  gv = genome_visitor_create(feature_visitor_class(), env);
  feature_visitor = feature_visitor_cast(gv);
  feature_visitor->features = fi;
  return gv;
}
