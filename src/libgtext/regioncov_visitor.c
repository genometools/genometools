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
  Hashtable *region2rangelist;
};

#define regioncov_visitor_cast(GV)\
        genome_visitor_cast(regioncov_visitor_class(), GV)

static void regioncov_visitor_free(GenomeVisitor *gv, Env *env)
{
  RegionCovVisitor *regioncov_visitor = regioncov_visitor_cast(gv);
  hashtable_delete(regioncov_visitor->region2rangelist, env);
}

static int regioncov_visitor_genome_feature(GenomeVisitor *gv,
                                            GenomeFeature *gf, Env *env)
{
  RegionCovVisitor *regioncov_visitor;
  env_error_check(env);
  regioncov_visitor = regioncov_visitor_cast(gv);
  return 0;
}

static int regioncov_visitor_sequence_region(GenomeVisitor *gv,
                                             SequenceRegion *sr, Env *env)
{
  RegionCovVisitor *regioncov_visitor;
  Array *rangelist;
  env_error_check(env);
  regioncov_visitor = regioncov_visitor_cast(gv);
  rangelist = array_new(sizeof (Range), env);
  hashtable_add(regioncov_visitor->region2rangelist,
                cstr_dup(str_get(genome_node_get_seqid((GenomeNode*) sr)), env),
                rangelist, env);
  return 0;
}

const GenomeVisitorClass* regioncov_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (RegionCovVisitor),
                                          regioncov_visitor_free,
                                          NULL,
                                          regioncov_visitor_genome_feature,
                                          regioncov_visitor_sequence_region,
                                          NULL };
  return &gvc;
}

GenomeVisitor* regioncov_visitor_new(unsigned long max_feature_dist, Env *env)
{
  GenomeVisitor *gv = genome_visitor_create(regioncov_visitor_class(), env);
  RegionCovVisitor *regioncov_visitor = regioncov_visitor_cast(gv);
  regioncov_visitor->max_feature_dist = max_feature_dist;
  regioncov_visitor->region2rangelist = hashtable_new(HASH_STRING,
                                                      env_ma_free_func,
                                                      (FreeFunc) array_delete,
                                                      env);
  return gv;
}

void regioncov_visitor_show_coverage(GenomeVisitor *gv)
{
  /*RegionCovVisitor *regioncov_visitor = regioncov_visitor_cast(gv); */
  printf("region coverage:\n");
  /*
  hashtable_foreach(regioncov_visitor->region2rangelist, show_rangelist, NULL,
                    env);
  */

  /* XXX:
    - save all sequence region & range list pairs
    - sort them
    - show them
  */
}
