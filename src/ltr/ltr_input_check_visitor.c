/*
  Copyright (c) 2014 Sascha Steinbiss <sascha@steinbiss.name>

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
#include "core/range_api.h"
#include "core/str_api.h"
#include "core/strand_api.h"
#include "core/symbol_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/node_visitor_api.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/feature_type.h"
#include "ltr/ltr_input_check_visitor.h"

struct GtLTRInputCheckVisitor {
  const GtNodeVisitor parent_instance;
  bool only_ltrs;
};

const GtNodeVisitorClass* gt_ltr_input_check_visitor_class(void);

#define gt_ltr_input_check_visitor_cast(GV)\
        gt_node_visitor_cast(gt_ltr_input_check_visitor_class(), GV)

static int gt_ltr_input_check_visitor_feature_node(GtNodeVisitor *nv,
                                                   GtFeatureNode *fn,
                                                   GtError *err)
{
  GT_UNUSED GtLTRInputCheckVisitor *lv;
  GtFeatureNodeIterator *fni;
  bool seen_left = false;
  GtFeatureNode *curnode = NULL,
                *ltr_retrotrans = NULL,
                *lltr = NULL,
                *rltr = NULL;
  int had_err = 0;
  lv = gt_ltr_input_check_visitor_cast(nv);
  gt_assert(lv);
  gt_error_check(err);

  /* traverse annotation subgraph and find LTR components */
  fni = gt_feature_node_iterator_new(fn);
  while (!had_err && (curnode = gt_feature_node_iterator_next(fni))) {
    if (strcmp(gt_feature_node_get_type(curnode),
               gt_ft_LTR_retrotransposon) == 0) {
      ltr_retrotrans = curnode;
    }
    if (strcmp(gt_feature_node_get_type(curnode),
               gt_ft_long_terminal_repeat) == 0) {
      if (seen_left)
        rltr = curnode;
      else {
        lltr = curnode;
        seen_left = true;
      }
    }
  }
  gt_feature_node_iterator_delete(fni);

  if (lv->only_ltrs) {
    if (!had_err && !ltr_retrotrans) {
      gt_error_set(err, "connected component with %s entry node (%s, line %u) "
                        "does not contain a '%s' node, which is required",
                   gt_feature_node_get_type(fn),
                   gt_genome_node_get_filename((GtGenomeNode*) fn),
                   gt_genome_node_get_line_number((GtGenomeNode*) fn),
                   gt_ft_LTR_retrotransposon);
      had_err = -1;
    }
  }

  if (!had_err && ltr_retrotrans && (!lltr || !rltr)) {
    gt_error_set(err, "LTR_retrotransposon feature (%s, line %u) "
                      "does not contain two %s child features, both of which "
                      "are required",
                 gt_genome_node_get_filename((GtGenomeNode*) ltr_retrotrans),
                 gt_genome_node_get_line_number((GtGenomeNode*) ltr_retrotrans),
                 gt_ft_long_terminal_repeat);
    had_err = -1;
  }

  return had_err;
}

const GtNodeVisitorClass* gt_ltr_input_check_visitor_class(void)
{
  static const GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtLTRInputCheckVisitor),
                                    NULL,
                                    NULL,
                                    gt_ltr_input_check_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  return nvc;
}

void gt_ltr_input_check_visitor_disallow_non_ltr(GtLTRInputCheckVisitor *liv)
{
  gt_assert(liv);
  liv->only_ltrs = true;
}

GtNodeVisitor* gt_ltr_input_check_visitor_new(void)
{
  GtNodeVisitor *nv = NULL;
  GtLTRInputCheckVisitor *lv;
  nv = gt_node_visitor_create(gt_ltr_input_check_visitor_class());
  gt_assert(nv);
  lv = gt_ltr_input_check_visitor_cast(nv);
  lv->only_ltrs = false;
  return nv;
}
