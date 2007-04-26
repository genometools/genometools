/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <libgtext/genome_node.h>
#include <libgtext/genome_visitor_rep.h>
#include <libgtext/gtf_visitor.h>

struct GTFVisitor {
  const GenomeVisitor parent_instance;
  GenFile *outfp;
};

#define gtf_visitor_cast(GV)\
        genome_visitor_cast(gtf_visitor_class(), GV)

static void gtf_visitor_free(GenomeVisitor *gv, Env *env)
{
  GTFVisitor *gtf_visitor = gtf_visitor_cast(gv);
  assert(gtf_visitor);
}

static int gtf_visitor_comment(GenomeVisitor *gv, Comment *c, Env *env)
{
  GTFVisitor *gtf_visitor;
  env_error_check(env);
  gtf_visitor = gtf_visitor_cast(gv);
  return 0;
}

static int gtf_visitor_genome_feature(GenomeVisitor *gv, GenomeFeature *gf,
                                       Env *env)
{
  GTFVisitor *gtf_visitor;
  env_error_check(env);
  gtf_visitor = gtf_visitor_cast(gv);
  return 0;
}

static int gtf_visitor_sequence_region(GenomeVisitor *gv, SequenceRegion *sr,
                                        Env *env)
{
  GTFVisitor *gtf_visitor;
  env_error_check(env);
  gtf_visitor = gtf_visitor_cast(gv);
  assert(gv && sr);
  return 0;
}

const GenomeVisitorClass* gtf_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (GTFVisitor),
                                          gtf_visitor_free,
                                          gtf_visitor_comment,
                                          gtf_visitor_genome_feature,
                                          gtf_visitor_sequence_region,
                                          NULL };
  return &gvc;
}

GenomeVisitor* gtf_visitor_new(GenFile *outfp, Env *env)
{
  GenomeVisitor *gv = genome_visitor_create(gtf_visitor_class(), env);
  GTFVisitor *gtf_visitor = gtf_visitor_cast(gv);
  gtf_visitor->outfp = outfp;
  return gv;
}
