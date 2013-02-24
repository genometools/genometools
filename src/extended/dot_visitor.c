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

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/genome_node.h"
#include "extended/gff3_defines.h"
#include "extended/dot_visitor.h"

struct GtDotVisitor {
  const GtNodeVisitor parent_instance;
  GtHashmap *ids;
  unsigned long idctr;
};

#define dv_cast(GV)\
        gt_node_visitor_cast(gt_dot_visitor_class(), GV)

static void dv_free(GtNodeVisitor *nv)
{
  GtDotVisitor *dv = dv_cast(nv);
  gt_assert(dv);
  gt_hashmap_delete(dv->ids);
}

static int dv_output_edge(GtFeatureNode *fn, void *data,
                          GT_UNUSED GtError *err)
{
  GtFeatureNodeIterator *fni;
  GtFeatureNode *curnode = NULL;
  unsigned long pid, cid;
  int had_err = 0;
  GtDotVisitor *dv = (GtDotVisitor*) data;

  if (!(pid = (unsigned long) gt_hashmap_get(dv->ids, fn))) {
    pid = dv->idctr++;
    gt_hashmap_add(dv->ids, fn, (void*) pid);
  }
  printf("%lu [shape=record, label=\"{ %lu | %s | %p }\"]\n",
             pid, pid, gt_feature_node_get_type(fn), fn);

  fni = gt_feature_node_iterator_new_direct(fn);
  while ((curnode = gt_feature_node_iterator_next(fni))) {
    if (!(cid = (unsigned long) gt_hashmap_get(dv->ids, curnode))) {
      cid = dv->idctr++;
      gt_hashmap_add(dv->ids, curnode, (void*) cid);
    }
    printf("%lu [shape=record, label=\"{ %lu | %s | %p }\"]\n",
             cid, cid, gt_feature_node_get_type(curnode), curnode);
    printf("%lu -> %lu\n", pid, cid);
  }
  gt_feature_node_iterator_delete(fni);

  return had_err;
}

static int dv_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                            GtError *err)
{
  GT_UNUSED GtDotVisitor *dv;
  int had_err = 0;
  unsigned long pid;
  gt_error_check(err);
  dv = dv_cast(nv);

  if (!(pid = (unsigned long) gt_hashmap_get(dv->ids, fn))) {
    pid = dv->idctr++;
    printf("subgraph %lu {\n", pid);
    gt_hashmap_add(dv->ids, fn, (void*) pid);
  }
  had_err = gt_feature_node_traverse_children(fn, dv,
                                              dv_output_edge,
                                              true, err);
  printf("}\n");
  gt_assert(!had_err);

  return had_err;
}

const GtNodeVisitorClass* gt_dot_visitor_class()
{
  static GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtDotVisitor),
                                    dv_free,
                                    NULL,
                                    dv_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  return nvc;
}

void gt_dot_visitor_finalize(GT_UNUSED GtNodeVisitor *v)
{
  printf("}\n");
}

GtNodeVisitor* gt_dot_visitor_new()
{
  GtNodeVisitor *nv = gt_node_visitor_create(gt_dot_visitor_class());
  GtDotVisitor *dv = dv_cast(nv);

  dv->ids = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
  dv->idctr = 1;
  printf("\ndigraph nodes {\n");

  return nv;
}
