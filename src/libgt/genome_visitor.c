/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdlib.h>
#include "genome_visitor_rep.h"
#include "xansi.h"

GenomeVisitor* genome_visitor_create(const GenomeVisitorClass *gvc)
{
  GenomeVisitor *gv;
  assert(gvc && gvc->size);
  gv = xcalloc(1, gvc->size);
  gv->c_class = gvc;
  return gv;
}

void* genome_visitor_cast(const GenomeVisitorClass *gvc, GenomeVisitor *gv)
{
  assert(gvc && gv && gv->c_class == gvc);
  return gv;
}

int genome_visitor_visit_comment(GenomeVisitor *gv, Comment *c, Log *l,
                                 Error *err)
{
  error_check(err);
  assert(gv && c && gv->c_class);
  if (gv->c_class->comment)
    return gv->c_class->comment(gv, c, l, err);
  else if (gv->c_class->default_func)
    return gv->c_class->default_func(gv, (GenomeNode*) c, l, err);
  return 0;
}

int genome_visitor_visit_genome_feature(GenomeVisitor *gv, Genome_feature *gf,
                                        Log *l, Error *err)
{
  error_check(err);
  assert(gv && gf && gv->c_class);
  if (gv->c_class->genome_feature)
    return gv->c_class->genome_feature(gv, gf, l, err);
  else if (gv->c_class->default_func)
    return gv->c_class->default_func(gv, (GenomeNode*) gf, l, err);
  return 0;
}

int genome_visitor_visit_sequence_region(GenomeVisitor *gv, SequenceRegion *sr,
                                         Log *l, Error *err)
{
  error_check(err);
  assert(gv && sr && gv->c_class);
  if (gv->c_class->sequence_region)
    return gv->c_class->sequence_region(gv, sr, l, err);
  else if (gv->c_class->default_func)
    return gv->c_class->default_func(gv, (GenomeNode*) sr, l, err);
  return 0;
}

void genome_visitor_free(GenomeVisitor *gv)
{
  if (!gv) return;
  assert(gv->c_class);
  if (gv->c_class->free)
    gv->c_class->free(gv);
  free(gv);
}
