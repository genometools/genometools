/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdlib.h>
#include <libgt/genome_visitor_rep.h>
#include <libgt/xansi.h>

GenomeVisitor* genome_visitor_create(const GenomeVisitorClass *gvc, Env *env)
{
  GenomeVisitor *gv;
  assert(gvc && gvc->size);
  gv = env_ma_calloc(env, 1, gvc->size);
  gv->c_class = gvc;
  return gv;
}

void* genome_visitor_cast(const GenomeVisitorClass *gvc, GenomeVisitor *gv)
{
  assert(gvc && gv && gv->c_class == gvc);
  return gv;
}

int genome_visitor_visit_comment(GenomeVisitor *gv, Comment *c, Env *env)
{
  env_error_check(env);
  assert(gv && c && gv->c_class);
  if (gv->c_class->comment)
    return gv->c_class->comment(gv, c, env);
  else if (gv->c_class->default_func)
    return gv->c_class->default_func(gv, (GenomeNode*) c, env);
  return 0;
}

int genome_visitor_visit_genome_feature(GenomeVisitor *gv, GenomeFeature *gf,
                                        Env *env)
{
  env_error_check(env);
  assert(gv && gf && gv->c_class);
  if (gv->c_class->genome_feature)
    return gv->c_class->genome_feature(gv, gf, env);
  else if (gv->c_class->default_func)
    return gv->c_class->default_func(gv, (GenomeNode*) gf, env);
  return 0;
}

int genome_visitor_visit_sequence_region(GenomeVisitor *gv, SequenceRegion *sr,
                                         Env *env)
{
  env_error_check(env);
  assert(gv && sr && gv->c_class);
  if (gv->c_class->sequence_region)
    return gv->c_class->sequence_region(gv, sr, env);
  else if (gv->c_class->default_func)
    return gv->c_class->default_func(gv, (GenomeNode*) sr, env);
  return 0;
}

void genome_visitor_delete(GenomeVisitor *gv, Env *env)
{
  if (!gv) return;
  assert(gv->c_class);
  if (gv->c_class->free)
    gv->c_class->free(gv, env);
  env_ma_free(gv, env);
}
