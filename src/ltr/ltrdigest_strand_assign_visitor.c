/*
  Copyright (c) 2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/strand_api.h"
#include "core/unused_api.h"
#include "extended/node_visitor_api.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/feature_type.h"
#include "ltr/ltrdigest_def.h"
#include "ltr/ltrdigest_strand_assign_visitor.h"

struct GtLTRdigestStrandAssignVisitor {
  const GtNodeVisitor parent_instance;
  GtStrand strand;
};

const GtNodeVisitorClass* gt_ltrdigest_strand_assign_visitor_class(void);

#define gt_ltrdigest_strand_assign_visitor_cast(GV)\
        gt_node_visitor_cast(gt_ltrdigest_strand_assign_visitor_class(), GV)

static int gt_ltrdigest_strand_assign_visitor_feature_node(GtNodeVisitor *nv,
                                                 GtFeatureNode *fn,
                                                 GtError *err)
{
  GtLTRdigestStrandAssignVisitor *lv;
  GtFeatureNodeIterator *fni;
  GtFeatureNode *curnode = NULL;
  int had_err = 0;
  lv = gt_ltrdigest_strand_assign_visitor_cast(nv);
  gt_assert(lv);
  gt_error_check(err);

  lv->strand = GT_STRAND_UNKNOWN;
  fni = gt_feature_node_iterator_new(fn);
  while (!had_err && (curnode = gt_feature_node_iterator_next(fni))) {
    GtStrand node_strand = gt_feature_node_get_strand(curnode);
    if (lv->strand == GT_STRAND_UNKNOWN && node_strand != lv->strand)
      lv->strand = node_strand;
    else {
      if (node_strand != GT_STRAND_UNKNOWN && node_strand != lv->strand) {
        gt_error_set(err, "inconsistent strands encountered in `%s' feature "
                        "in file %s, line %u: found %c, expected %c",
                        gt_feature_node_get_type(curnode),
                        gt_genome_node_get_filename((GtGenomeNode*) curnode),
                        gt_genome_node_get_line_number((GtGenomeNode*) curnode),
                        GT_STRAND_CHARS[node_strand],
                        GT_STRAND_CHARS[lv->strand]);
        had_err = -1;
      }
    }
  }
  gt_feature_node_iterator_delete(fni);

  if (!had_err && lv->strand != GT_STRAND_UNKNOWN) {
    gt_feature_node_set_strand(fn, lv->strand);
    fni = gt_feature_node_iterator_new(fn);
    while (!had_err && (curnode = gt_feature_node_iterator_next(fni))) {
      gt_feature_node_set_strand(curnode, lv->strand);
    }
    gt_feature_node_iterator_delete(fni);
  }
  return had_err;
}

GT_UNUSED void gt_ltrdigest_strand_assign_visitor_free(GtNodeVisitor *nv)
{
  GT_UNUSED GtLTRdigestStrandAssignVisitor *lv;
  if (!nv) return;
  lv = gt_ltrdigest_strand_assign_visitor_cast(nv);
}

const GtNodeVisitorClass* gt_ltrdigest_strand_assign_visitor_class(void)
{
  static const GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtLTRdigestStrandAssignVisitor),
                                NULL,
                                NULL,
                                gt_ltrdigest_strand_assign_visitor_feature_node,
                                NULL,
                                NULL,
                                NULL);
  }
  return nvc;
}

GtNodeVisitor* gt_ltrdigest_strand_assign_visitor_new()
{
  GtNodeVisitor *nv = NULL;
  GtLTRdigestStrandAssignVisitor *sav;
  nv = gt_node_visitor_create(gt_ltrdigest_strand_assign_visitor_class());
  sav = gt_ltrdigest_strand_assign_visitor_cast(nv);
  sav->strand = GT_STRAND_UNKNOWN;
  return nv;
}
