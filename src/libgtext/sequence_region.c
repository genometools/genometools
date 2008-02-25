/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
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

static void sequence_region_free(GenomeNode *gn)
{
  SequenceRegion *sr = sequence_region_cast(gn);
  assert(sr && sr->seqid);
  str_delete(sr->seqid);
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

static void sequence_region_set_seqid(GenomeNode *gn, Str *seqid)
{
  SequenceRegion *sr = sequence_region_cast(gn);
  assert(sr && seqid);
  str_delete(sr->seqid);
  sr->seqid = str_ref(seqid);
}

static int sequence_region_accept(GenomeNode *gn, GenomeVisitor *gv, Error *e)
{
  SequenceRegion *sr;
  error_check(e);
  sr = sequence_region_cast(gn);
  return genome_visitor_visit_sequence_region(gv, sr, e);
}

const GenomeNodeClass* sequence_region_class()
{
  static const GenomeNodeClass gnc = { sizeof (SequenceRegion),
                                       sequence_region_free,
                                       sequence_region_get_seqid,
                                       sequence_region_get_seqid,
                                       sequence_region_get_range,
                                       sequence_region_set_range,
                                       sequence_region_set_seqid,
                                       sequence_region_accept };
  return &gnc;
}

GenomeNode* sequence_region_new(Str *seqid, Range range, Str *filename,
                                unsigned long line_number)
{
  GenomeNode *gn = genome_node_create(sequence_region_class(), filename,
                                      line_number);
  SequenceRegion *sr = sequence_region_cast(gn);
  assert(seqid);
  sr->seqid = str_ref(seqid);
  sr->range = range;
  return gn;
}

void sequence_regions_consolidate(GenomeNode *gn_a, GenomeNode *gn_b)
{
  Range range_a, range_b;

  assert(gn_a);
  assert(gn_b);
  assert(genome_node_cast(sequence_region_class(), gn_a));
  assert(genome_node_cast(sequence_region_class(), gn_b));
  assert(!str_cmp(genome_node_get_seqid(gn_a), genome_node_get_seqid(gn_b)));

  range_a = genome_node_get_range(gn_a);
  range_b = genome_node_get_range(gn_b);
  range_a = range_join(range_a, range_b);
  genome_node_set_range(gn_a, range_a);
}
