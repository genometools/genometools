/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdlib.h>
#include "libgtext/sequence_region.h"
#include "libgtext/genome_node_rep.h"

struct SequenceRegion
{
  const GenomeNode parent_instance;
  Str *seqid;
  Range range;
};

#define sequence_region_cast(GN)\
        genome_node_cast(sequence_region_class(), GN)

static void sequence_region_free(GenomeNode *gn, Env *env)
{
  SequenceRegion *sr = sequence_region_cast(gn);
  assert(sr && sr->seqid);
  str_delete(sr->seqid, env);
}

static Str* sequence_region_get_seqid(GenomeNode *gn)
{
  SequenceRegion *sr = sequence_region_cast(gn);
  return sr->seqid;
}

static Range sequence_region_get_range(GenomeNode *gn)
{
  SequenceRegion *sr = sequence_region_cast(gn);
  return sr->range;
}

static void sequence_region_set_range(GenomeNode *gn, Range range)
{
  SequenceRegion *sr = sequence_region_cast(gn);
  sr->range = range;
}

static int sequence_region_accept(GenomeNode *gn, GenomeVisitor *gv, Env *env)
{
  SequenceRegion *sr;
  env_error_check(env);
  sr = sequence_region_cast(gn);
  return genome_visitor_visit_sequence_region(gv, sr, env);
}

const GenomeNodeClass* sequence_region_class()
{
  static const GenomeNodeClass gnc = { sizeof (SequenceRegion),
                                       sequence_region_free,
                                       sequence_region_get_seqid,
                                       sequence_region_get_seqid,
                                       sequence_region_get_range,
                                       sequence_region_set_range,
                                       NULL,
                                       NULL,
                                       NULL,
                                       sequence_region_accept };
  return &gnc;
}

GenomeNode* sequence_region_new(Str *seqid, Range range, Str *filename,
                                unsigned long line_number, Env *env)
{
  GenomeNode *gn = genome_node_create(sequence_region_class(), filename,
                                      line_number, env);
  SequenceRegion *sr = sequence_region_cast(gn);
  assert(seqid);
  sr->seqid = str_ref(seqid);
  sr->range = range;
  return gn;
}
