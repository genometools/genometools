/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include "core/cstr_api.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/str_api.h"
#include "core/unused_api.h"
#include "extended/region_node.h"
#include "extended/region_node_builder.h"

struct GtRegionNodeBuilder {
  GtHashmap *sequence_region_to_range;
};

GtRegionNodeBuilder* gt_region_node_builder_new(void)
{
  GtRegionNodeBuilder *rnb = gt_malloc(sizeof *rnb);
  rnb->sequence_region_to_range = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                                 gt_free_func);
  return rnb;
}

void gt_region_node_builder_delete(GtRegionNodeBuilder *rnb)
{
  if (!rnb) return;
  gt_hashmap_delete(rnb->sequence_region_to_range);
  gt_free(rnb);
}

void gt_region_node_builder_add_region(GtRegionNodeBuilder *rnb,
                                       const char *region_name,
                                       GtRange region_range)
{
  GtRange *rangeptr;
  gt_assert(rnb && region_name);
  if ((rangeptr = gt_hashmap_get(rnb->sequence_region_to_range, region_name))) {
    /* sequence region is already defined -> update range */
    *rangeptr = gt_range_join(&region_range, rangeptr);
  }
  else {
    /* sequence region is not already defined -> define it */
    rangeptr = gt_malloc(sizeof (GtRange));
    *rangeptr = region_range;
    gt_hashmap_add(rnb->sequence_region_to_range, gt_cstr_dup(region_name),
                   rangeptr);
  }
}

static int build_region_nodes(void *key, void *value, void *data,
                              GT_UNUSED GtError *err)
{
  GtStr *seqid;
  GtRange range;
  GtGenomeNode *gn;
  GtQueue *genome_nodes = (GtQueue*) data;
  gt_error_check(err);
  gt_assert(key && value && data);
  seqid = gt_str_new_cstr(key);
  range = *(GtRange*) value;
  gn = gt_region_node_new(seqid, range.start, range.end);
  gt_queue_add(genome_nodes, gn);
  gt_str_delete(seqid);
  return 0;
}

void gt_region_node_builder_build(const GtRegionNodeBuilder *rnb,
                                  GtQueue *genome_nodes)
{
  GT_UNUSED int had_err;
  gt_assert(rnb && genome_nodes);
  had_err = gt_hashmap_foreach(rnb->sequence_region_to_range,
                               build_region_nodes, genome_nodes, NULL);
  gt_assert(!had_err); /* build_region_nodes() is sane */
}

void gt_region_node_builder_reset(GtRegionNodeBuilder *rnb)
{
  gt_assert(rnb);
  gt_hashmap_reset(rnb->sequence_region_to_range);
}
