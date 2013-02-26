/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/cstr_api.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/range.h"
#include "core/unused_api.h"
#include "extended/node_visitor_api.h"
#include "extended/regioncov_visitor.h"

struct GtRegionCovVisitor {
  const GtNodeVisitor parent_instance;
  unsigned long max_feature_dist;
  GtHashmap *region2rangelist;
};

#define gt_regioncov_visitor_cast(GV)\
        gt_node_visitor_cast(gt_regioncov_visitor_class(), GV)

static void gt_regioncov_visitor_free(GtNodeVisitor *nv)
{
  GtRegionCovVisitor *regioncov_visitor = gt_regioncov_visitor_cast(nv);
  gt_hashmap_delete(regioncov_visitor->region2rangelist);
}

static int gt_regioncov_visitor_feature_node(GtNodeVisitor *nv,
                                             GtFeatureNode *fn,
                                             GT_UNUSED GtError *err)
{
  GtRange *old_range_ptr, old_range, new_range;
  GtArray *ranges;
  GtRegionCovVisitor *regioncov_visitor;
  gt_error_check(err);
  regioncov_visitor = gt_regioncov_visitor_cast(nv);
  ranges = gt_hashmap_get(regioncov_visitor->region2rangelist,
                       gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*)
                                                           fn)));
  gt_assert(ranges);
  new_range = gt_genome_node_get_range((GtGenomeNode*) fn);
  if (!gt_array_size(ranges))
    gt_array_add(ranges, new_range);
  else {
    old_range_ptr = gt_array_get_last(ranges);
    old_range = *old_range_ptr;
    old_range.end += regioncov_visitor->max_feature_dist;
    if (gt_range_overlap(&old_range, &new_range)) {
      old_range_ptr->end = MAX(old_range_ptr->end, new_range.end);
    }
    else
      gt_array_add(ranges, new_range);
  }
  return 0;
}

static int gt_regioncov_visitor_region_node(GtNodeVisitor *nv, GtRegionNode *rn,
                                         GT_UNUSED GtError *err)
{
  GtRegionCovVisitor *regioncov_visitor;
  GtArray *rangelist;
  gt_error_check(err);
  regioncov_visitor = gt_regioncov_visitor_cast(nv);
  rangelist = gt_array_new(sizeof (GtRange));
  gt_hashmap_add(regioncov_visitor->region2rangelist,
              gt_cstr_dup(gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*)
                                                              rn))),
              rangelist);
  return 0;
}

const GtNodeVisitorClass* gt_regioncov_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtRegionCovVisitor),
                                    gt_regioncov_visitor_free,
                                    NULL,
                                    gt_regioncov_visitor_feature_node,
                                    gt_regioncov_visitor_region_node,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_regioncov_visitor_new(unsigned long max_feature_dist)
{
  GtNodeVisitor *nv = gt_node_visitor_create(gt_regioncov_visitor_class());
  GtRegionCovVisitor *regioncov_visitor = gt_regioncov_visitor_cast(nv);
  regioncov_visitor->max_feature_dist = max_feature_dist;
  regioncov_visitor->region2rangelist =
    gt_hashmap_new(GT_HASH_STRING, gt_free_func, (GtFree) gt_array_delete);
  return nv;
}

static int show_rangelist(void *key, void *value, GT_UNUSED void *data,
                          GT_UNUSED GtError *err)
{
  unsigned long i;
  GtArray *rangelist;
  GtRange *rangeptr;
  gt_error_check(err);
  gt_assert(key && value);
  rangelist = (GtArray*) value;
  if (gt_array_size(rangelist)) {
    gt_assert(gt_ranges_are_sorted_and_do_not_overlap(rangelist));
    printf("%s:\n", (char*) key);
    for (i = 0; i < gt_array_size(rangelist); i++) {
      rangeptr = gt_array_get(rangelist, i);
      printf("%lu, %lu\n", rangeptr->start, rangeptr->end);
    }
  }
  return 0;
}

void gt_regioncov_visitor_show_coverage(GtNodeVisitor *nv)
{
  GtRegionCovVisitor *regioncov_visitor = gt_regioncov_visitor_cast(nv);
  GT_UNUSED int had_err;
  had_err = gt_hashmap_foreach_in_key_order(regioncov_visitor->region2rangelist,
                                         show_rangelist, NULL, NULL);
  gt_assert(!had_err); /* show_rangelist() is sane */
}
