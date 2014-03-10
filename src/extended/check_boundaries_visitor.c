/*
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>

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

#include <string.h>
#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/cstr_api.h"
#include "core/hashmap.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#include "core/queue.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/check_boundaries_visitor_api.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_visitor_api.h"

struct GtCheckBoundariesVisitor {
  const GtNodeVisitor parent_instance;
};

static int check_boundaries_visitor_check_rec(GtFeatureNode *parent,
                                              GtFeatureNode *child,
                                              GtError *err)
{
  GtFeatureNodeIterator *fni;
  GtFeatureNode *node;
  GtRange range,
          p_range;
  int had_err = 0;

  range = gt_genome_node_get_range((GtGenomeNode*) child);
  p_range = gt_genome_node_get_range((GtGenomeNode*) parent);

  if (range.start < p_range.start || range.end > p_range.end) {
    gt_warning("%s child range " GT_WU "-" GT_WU " (file %s, line %u) not "
               "contained in %s parent range " GT_WU "-" GT_WU " (file %s, "
               "line %u)",
               gt_feature_node_get_type(child),
               range.start, range.end,
               gt_genome_node_get_filename((GtGenomeNode*) child),
               gt_genome_node_get_line_number((GtGenomeNode*) child),
               gt_feature_node_get_type(parent),
               p_range.start, p_range.end,
               gt_genome_node_get_filename((GtGenomeNode*) parent),
               gt_genome_node_get_line_number((GtGenomeNode*) parent));
  }

  fni = gt_feature_node_iterator_new_direct(child);
  while ((node = gt_feature_node_iterator_next(fni))) {
    had_err = check_boundaries_visitor_check_rec(child, node, err);
  }
  gt_feature_node_iterator_delete(fni);

  return had_err;
}

static int check_boundaries_visitor_feature_node(GT_UNUSED GtNodeVisitor *nv,
                                                 GtFeatureNode *fn,
                                                 GT_UNUSED GtError *err)
{
  GtFeatureNodeIterator *fni;
  GtFeatureNode *node;
  int had_err = 0;

  fni = gt_feature_node_iterator_new_direct(fn);
  while (!had_err && (node = gt_feature_node_iterator_next(fni))) {
    had_err = check_boundaries_visitor_check_rec(fn, node, err);
  }
  gt_feature_node_iterator_delete(fni);

  return 0;
}

const GtNodeVisitorClass* gt_check_boundaries_visitor_class()
{
  static GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtCheckBoundariesVisitor),
                                    NULL,
                                    NULL,
                                    check_boundaries_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_check_boundaries_visitor_new(void)
{
  GtNodeVisitor *nv;
  nv = gt_node_visitor_create(gt_check_boundaries_visitor_class());
  return nv;
}
