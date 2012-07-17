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

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/feature_node.h"
#include "extended/genome_node.h"
#include "extended/gff3_defines.h"
#include "extended/multi_sanitizer_visitor.h"

struct GtMultiSanitizerVisitor {
  const GtNodeVisitor parent_instance;
  GtHashmap *first_elems;
};

#define msv_cast(GV)\
        gt_node_visitor_cast(gt_multi_sanitizer_visitor_class(), GV)

static void msv_free(GtNodeVisitor *nv)
{
  GtMultiSanitizerVisitor *msv = msv_cast(nv);
  gt_assert(msv);
  gt_hashmap_delete(msv->first_elems);
}

static int msv_fixup_representatives(GtFeatureNode *fn, void *data,
                                     GT_UNUSED GtError *err)
{
  int had_err = 0;
  GtMultiSanitizerVisitor *msv = (GtMultiSanitizerVisitor*) data;

  if (gt_feature_node_is_multi(fn)) {
    GtFeatureNode *rep = NULL;
    if (!(rep = gt_hashmap_get(msv->first_elems,
                               gt_feature_node_get_multi_representative(fn)))) {
      gt_hashmap_add(msv->first_elems,
                     gt_feature_node_get_multi_representative(fn),
                     fn);
      rep = fn;
      gt_feature_node_unset_multi(fn);
      gt_feature_node_make_multi_representative(fn);
    } else {
      gt_assert(rep);
      gt_feature_node_unset_multi(fn);
      gt_feature_node_set_multi_representative(fn, rep);
    }
  }
  return had_err;
}

static int msv_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                            GtError *err)
{
  GtMultiSanitizerVisitor *msv;
  int had_err = 0;
  gt_error_check(err);
  msv = msv_cast(nv);

  gt_hashmap_reset(msv->first_elems);
  had_err = gt_feature_node_traverse_children(fn, msv,
                                              msv_fixup_representatives,
                                              true, err);
  gt_assert(!had_err);   /* msv_fixup_representatives() is sane */

  return had_err;
}

const GtNodeVisitorClass* gt_multi_sanitizer_visitor_class()
{
  static GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtMultiSanitizerVisitor),
                                    msv_free,
                                    NULL,
                                    msv_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  return nvc;
}

GtNodeVisitor* gt_multi_sanitizer_visitor_new()
{
  GtNodeVisitor *nv =
                     gt_node_visitor_create(gt_multi_sanitizer_visitor_class());
  GtMultiSanitizerVisitor *msv = msv_cast(nv);
  msv->first_elems = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);

  return nv;
}
