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

void genome_visitor_visit_comment(GenomeVisitor *gv, Comment *c, Log *l)
{
  assert(gv && c && gv->c_class);
  if (gv->c_class->comment)
    gv->c_class->comment(gv, c, l);
  else if (gv->c_class->default_func)
    gv->c_class->default_func(gv, (GenomeNode*) c, l);
}

void genome_visitor_visit_genome_feature(GenomeVisitor *gv, Genome_feature *gf,
                                         Log *l)
{
  assert(gv && gf && gv->c_class);
  if (gv->c_class->genome_feature)
    gv->c_class->genome_feature(gv, gf, l);
  else if (gv->c_class->default_func)
    gv->c_class->default_func(gv, (GenomeNode*) gf, l);
}

void genome_visitor_visit_sequence_region(GenomeVisitor *gv, SequenceRegion *sr,
                                          Log *l)
{
  assert(gv && sr && gv->c_class);
  if (gv->c_class->sequence_region)
    gv->c_class->sequence_region(gv, sr, l);
  else if (gv->c_class->default_func)
    gv->c_class->default_func(gv, (GenomeNode*) sr, l);
}

void genome_visitor_free(GenomeVisitor *gv)
{
  if (!gv) return;
  assert(gv->c_class);
  if (gv->c_class->free)
    gv->c_class->free(gv);
  free(gv);
}
