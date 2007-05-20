/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <gtcore.h>
#include <libgtext/genome_visitor_rep.h>
#include <libgtext/regioncov_visitor.h>

struct RegionCovVisitor {
  const GenomeVisitor parent_instance;
  unsigned long max_feature_dist;
};

#define regioncov_visitor_cast(GV)\
        genome_visitor_cast(regioncov_visitor_class(), GV)

static void regioncov_visitor_free(GenomeVisitor *gv, Env *env)
{
  /* RegionCovVisitor *regioncov_visitor = regioncov_visitor_cast(gv); */
}

static int regioncov_visitor_genome_feature(GenomeVisitor *gv,
                                            GenomeFeature *gf, Env *env)
{
  RegionCovVisitor *regioncov_visitor;
  env_error_check(env);
  regioncov_visitor = regioncov_visitor_cast(gv);
  return 0;
}

const GenomeVisitorClass* regioncov_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (RegionCovVisitor),
                                          regioncov_visitor_free,
                                          NULL,
                                          regioncov_visitor_genome_feature,
                                          NULL,
                                          NULL };
  return &gvc;
}

GenomeVisitor* regioncov_visitor_new(unsigned long max_feature_dist, Env *env)
{
  GenomeVisitor *gv = genome_visitor_create(regioncov_visitor_class(), env);
  RegionCovVisitor *regioncov_visitor = regioncov_visitor_cast(gv);
  regioncov_visitor->max_feature_dist = max_feature_dist;
  return gv;
}

void regioncov_visitor_show_coverage(GenomeVisitor *gv)
{
  /* RegionCovVisitor *regioncov_visitor = regioncov_visitor_cast(gv); */
  printf("region coverage:\n");
}
