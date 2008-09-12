/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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
#include "extended/sequence_region.h"
#include "extended/genome_node_rep.h"

struct GT_SequenceRegion
{
  const GtGenomeNode parent_instance;
  GtStr *seqid;
  GtRange range;
};

#define gt_sequence_region_cast(GN)\
        gt_genome_node_cast(gt_sequence_region_class(), GN)

static void gt_sequence_region_free(GtGenomeNode *gn)
{
  GT_SequenceRegion *sr = gt_sequence_region_cast(gn);
  assert(sr && sr->seqid);
  gt_str_delete(sr->seqid);
}

static GtStr* gt_sequence_region_get_seqid(GtGenomeNode *gn)
{
  GT_SequenceRegion *sr = gt_sequence_region_cast(gn);
  return sr->seqid;
}

static GtRange gt_sequence_region_get_range(GtGenomeNode *gn)
{
  GT_SequenceRegion *sr = gt_sequence_region_cast(gn);
  return sr->range;
}

static void gt_sequence_region_set_range(GtGenomeNode *gn, GtRange range)
{
  GT_SequenceRegion *sr = gt_sequence_region_cast(gn);
  sr->range = range;
}

static void gt_sequence_region_change_seqid(GtGenomeNode *gn, GtStr *seqid)
{
  GT_SequenceRegion *sr = gt_sequence_region_cast(gn);
  assert(sr && seqid);
  gt_str_delete(sr->seqid);
  sr->seqid = gt_str_ref(seqid);
}

static int gt_sequence_region_accept(GtGenomeNode *gn, GenomeVisitor *gv,
                                     GtError *err)
{
  GT_SequenceRegion *sr;
  gt_error_check(err);
  sr = gt_sequence_region_cast(gn);
  return genome_visitor_visit_sequence_region(gv, sr, err);
}

const GtGenomeNodeClass* gt_sequence_region_class()
{
  static const GtGenomeNodeClass gnc = { sizeof (GT_SequenceRegion),
                                       gt_sequence_region_free,
                                       gt_sequence_region_get_seqid,
                                       gt_sequence_region_get_seqid,
                                       gt_sequence_region_get_range,
                                       gt_sequence_region_set_range,
                                       gt_sequence_region_change_seqid,
                                       gt_sequence_region_accept };
  return &gnc;
}

GtGenomeNode* gt_sequence_region_new(GtStr *seqid, GtRange range)
{
  GtGenomeNode *gn = gt_genome_node_create(gt_sequence_region_class());
  GT_SequenceRegion *sr = gt_sequence_region_cast(gn);
  assert(seqid);
  sr->seqid = gt_str_ref(seqid);
  sr->range = range;
  return gn;
}

void gt_sequence_regions_consolidate(GtGenomeNode *gn_a, GtGenomeNode *gn_b)
{
  GtRange gt_range_a, gt_range_b;

  assert(gn_a);
  assert(gn_b);
  assert(gt_genome_node_cast(gt_sequence_region_class(), gn_a));
  assert(gt_genome_node_cast(gt_sequence_region_class(), gn_b));
  assert(!gt_str_cmp(gt_genome_node_get_seqid(gn_a),
          gt_genome_node_get_seqid(gn_b)));

  gt_range_a = gt_genome_node_get_range(gn_a);
  gt_range_b = gt_genome_node_get_range(gn_b);
  gt_range_a = gt_range_join(gt_range_a, gt_range_b);
  gt_genome_node_set_range(gn_a, gt_range_a);
}
