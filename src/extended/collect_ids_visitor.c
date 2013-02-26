/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "extended/collect_ids_visitor.h"
#include "extended/node_visitor_api.h"

struct GtCollectIDsVisitor {
  const GtNodeVisitor parent_instance;
  GtCstrTable *seqids;
};

const GtNodeVisitorClass* gt_collect_ids_visitor_class(void);

#define collect_ids_visitor_cast(GV)\
        gt_node_visitor_cast(gt_collect_ids_visitor_class(), GV)

static int collect_ids_visitor_feature_node(GtNodeVisitor *nv,
                                            GtFeatureNode *fn,
                                            GT_UNUSED GtError *err)
{
  GtCollectIDsVisitor *civ;
  const char *seqid;
  civ = collect_ids_visitor_cast(nv);
  seqid = gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) fn));
  if (!gt_cstr_table_get(civ->seqids, seqid))
    gt_cstr_table_add(civ->seqids, seqid);
  return 0;
}

static int collect_ids_visitor_region_node(GtNodeVisitor *nv, GtRegionNode *rn,
                                           GT_UNUSED GtError *err)
{
  GtCollectIDsVisitor *civ;
  const char *seqid;
  int had_err = 0;
  gt_error_check(err);
  civ = collect_ids_visitor_cast(nv);
  seqid = gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*) rn));
  if (!gt_cstr_table_get(civ->seqids, seqid))
    gt_cstr_table_add(civ->seqids, seqid);
  return had_err;
}

const GtNodeVisitorClass* gt_collect_ids_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtCollectIDsVisitor),
                                    NULL,
                                    NULL,
                                    collect_ids_visitor_feature_node,
                                    collect_ids_visitor_region_node,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_collect_ids_visitor_new(GtCstrTable *cst)
{
  GtNodeVisitor *nv = gt_node_visitor_create(gt_collect_ids_visitor_class());
  GtCollectIDsVisitor *collect_ids_visitor = collect_ids_visitor_cast(nv);
  collect_ids_visitor->seqids = cst;
  return nv;
}
