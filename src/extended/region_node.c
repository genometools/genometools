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

#include <stdlib.h>
#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "extended/genome_node_rep.h"
#include "extended/region_node.h"

struct GtRegionNode
{
  const GtGenomeNode parent_instance;
  GtStr *seqid;
  GtRange range;
};

static void region_node_free(GtGenomeNode *gn)
{
  GtRegionNode *rn = gt_region_node_cast(gn);
  gt_assert(rn && rn->seqid);
  gt_str_delete(rn->seqid);
}

static GtStr* region_node_get_seqid(GtGenomeNode *gn)
{
  GtRegionNode *rn = gt_region_node_cast(gn);
  return rn->seqid;
}

static GtRange region_node_get_range(GtGenomeNode *gn)
{
  GtRegionNode *rn = gt_region_node_cast(gn);
  return rn->range;
}

static void region_node_set_range(GtGenomeNode *gn, const GtRange *range)
{
  GtRegionNode *rn = gt_region_node_cast(gn);
  rn->range = *range;
}

static void region_node_change_seqid(GtGenomeNode *gn, GtStr *seqid)
{
  GtRegionNode *rn = gt_region_node_cast(gn);
  gt_assert(rn && seqid);
  gt_str_delete(rn->seqid);
  rn->seqid = gt_str_ref(seqid);
}

static int region_node_accept(GtGenomeNode *gn, GtNodeVisitor *nv, GtError *err)
{
  GtRegionNode *rn;
  gt_error_check(err);
  rn = gt_region_node_cast(gn);
  return gt_node_visitor_visit_region_node(nv, rn, err);
}

const GtGenomeNodeClass* gt_region_node_class()
{
  static const GtGenomeNodeClass *gnc = NULL;
  gt_class_alloc_lock_enter();
  if (!gnc) {
    gnc = gt_genome_node_class_new(sizeof (GtRegionNode),
                                   region_node_free,
                                   region_node_get_seqid,
                                   region_node_get_seqid,
                                   region_node_get_range,
                                   region_node_set_range,
                                   region_node_change_seqid,
                                   region_node_accept);
  }
  gt_class_alloc_lock_leave();
  return gnc;
}

GtGenomeNode* gt_region_node_new(GtStr *seqid, unsigned long start,
                                               unsigned long end)
{
  GtGenomeNode *gn = gt_genome_node_create(gt_region_node_class());
  GtRegionNode *rn = gt_region_node_cast(gn);
  gt_assert(seqid);
  gt_assert(start <= end);
  rn->seqid = gt_str_ref(seqid);
  rn->range.start = start;
  rn->range.end   = end;
  return gn;
}

void gt_region_node_consolidate(GtRegionNode *rn_a, GtRegionNode *rn_b)
{
  GtRange range_a, range_b;
  gt_assert(rn_a);
  gt_assert(rn_b);
  gt_assert(!gt_str_cmp(gt_genome_node_get_seqid((GtGenomeNode*) rn_a),
                        gt_genome_node_get_seqid((GtGenomeNode*) rn_b)));
  range_a = gt_genome_node_get_range((GtGenomeNode*) rn_a);
  range_b = gt_genome_node_get_range((GtGenomeNode*) rn_b);
  range_a = gt_range_join(&range_a, &range_b);
  gt_genome_node_set_range((GtGenomeNode*) rn_a, &range_a);
}
