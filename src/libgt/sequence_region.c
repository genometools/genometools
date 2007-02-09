/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdlib.h>
#include "sequence_region.h"
#include "genome_node_rep.h"

struct Sequence_region
{
  const Genome_node parent_instance;
  Str *seqid;
  Range range;
};

#define sequence_region_cast(GN)\
        genome_node_cast(sequence_region_class(), GN)

static void sequence_region_free(Genome_node *gn)
{
  Sequence_region *sr = sequence_region_cast(gn);
  assert(sr && sr->seqid);
  str_free(sr->seqid);
}

static Str* sequence_region_get_seqid(Genome_node *gn)
{
  Sequence_region *sr = sequence_region_cast(gn);
  return sr->seqid;
}

static Range sequence_region_get_range(Genome_node *gn)
{
  Sequence_region *sr = sequence_region_cast(gn);
  return sr->range;
}

static void sequence_region_set_range(Genome_node *gn, Range range)
{
  Sequence_region *sr = sequence_region_cast(gn);
  sr->range = range;
}

static void sequence_region_accept(Genome_node *gn, Genome_visitor *gv, Log *l)
{
  Sequence_region *sr = sequence_region_cast(gn);
  genome_visitor_visit_sequence_region(gv, sr, l);
}

const Genome_node_class* sequence_region_class()
{
  static const Genome_node_class gnc = { sizeof(Sequence_region),
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

Genome_node* sequence_region_new(Str *seqid,
                                 Range range,
                                 const char *filename,
                                 unsigned long line_number)
{
  Genome_node *gn = genome_node_create(sequence_region_class(), filename,
                                       line_number);
  Sequence_region *sr = sequence_region_cast(gn);
  assert(seqid);
  sr->seqid = str_ref(seqid);
  sr->range = range;
  return gn;
}
