/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "core/cstr_api.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/range.h"
#include "core/unused_api.h"
#include "extended/node_visitor_rep.h"
#include "extended/region_node.h"
#include "ltr/ltr_visitor.h"

struct GtLTRVisitor {
  const GtNodeVisitor parent_instance;
  GtLTRElement *element;
};

#define gt_ltr_visitor_cast(GV)\
        gt_node_visitor_cast(gt_ltr_visitor_class(), GV)

static int gt_ltr_visitor_feature_node(GtNodeVisitor *gv, GtFeatureNode *gf,
                                    GT_UNUSED GtError *err)
{
  GtLTRVisitor *lv;
  GtRange node_range;
  GtArray *pdomarr = NULL;
  const char *pfamname;
  const char *gft;
  lv = gt_ltr_visitor_cast(gv);
  gt_assert(lv);
  gt_error_check(err);

  gft = gt_feature_node_get_type(gf);

  if (strcmp(gft, "LTR_retrotransposon") == 0)
  {
    lv->element->mainnode = gf;
  } else if (strcmp(gft, "long_terminal_repeat") == 0)
  {
    if (lv->element->leftLTR == NULL)
    {
      node_range = gt_genome_node_get_range((GtGenomeNode*) gf);
      lv->element->leftLTR = gf;
      /* compensate for 1-based node coords */
      lv->element->leftLTR_5 = node_range.start - 1;
      lv->element->leftLTR_3 = node_range.end - 1;
    }
    else
    {
      node_range = gt_genome_node_get_range((GtGenomeNode*) gf);
      lv->element->rightLTR = gf;
      /* compensate for 1-based node coords */
      lv->element->rightLTR_5 = node_range.start - 1;
      lv->element->rightLTR_3 = node_range.end - 1;
    }
  } else if (strcmp(gft, "target_site_duplication") == 0)
  {
    if (lv->element->leftTSD == NULL)
    {
      lv->element->leftTSD = gf;
    }
    else
    {
      lv->element->rightTSD = gf;
    }
  } else if (strcmp(gft, "RR_tract") == 0)
  {
    if (lv->element->ppt == NULL)
    {
      lv->element->ppt = gf;
    }
  } else if (strcmp(gft, "primer_binding_site") == 0)
  {
    if (lv->element->pbs == NULL)
    {
      lv->element->pbs = gf;
    }
  } else if (strcmp(gft,"protein_match") == 0)
  {
    if (!lv->element->pdoms)
    {
      lv->element->pdoms = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                          (GtFree) gt_array_delete);
    }
    pfamname = gt_feature_node_get_attribute(gf, "name");
    if (!(pdomarr = (GtArray*) gt_hashmap_get(lv->element->pdoms, pfamname)))
    {
      char *pfamcpy = gt_cstr_dup(pfamname);
      pdomarr = gt_array_new(sizeof (GtFeatureNode*));
      gt_hashmap_add(lv->element->pdoms, pfamcpy, pdomarr);
      if (lv->element->pdomorder)
        gt_array_add(lv->element->pdomorder, pfamcpy);
    }
    gt_array_add(pdomarr, gf);
  }
  return 0;
}

const GtNodeVisitorClass* gt_ltr_visitor_class(void)
{
  static const GtNodeVisitorClass *gvc = NULL;
  if (!gvc) {
    gvc = gt_node_visitor_class_new(sizeof (GtLTRVisitor),
                                    NULL,
                                    NULL,
                                    gt_ltr_visitor_feature_node,
                                    NULL,
                                    NULL);
  }
  return gvc;
}

GtNodeVisitor* gt_ltr_visitor_new(GtLTRElement *element)
{
  GtNodeVisitor *gv;
  GtLTRVisitor *lv;
  gt_assert(element);
  gv = gt_node_visitor_create(gt_ltr_visitor_class());
  lv = gt_ltr_visitor_cast(gv);
  lv->element = element;
  gt_assert(lv);
  return gv;
}
