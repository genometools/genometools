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
#include "core/class_alloc_lock.h"
#include "core/cstr_api.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/range.h"
#include "core/unused_api.h"
#include "extended/feature_type.h"
#include "extended/node_visitor_api.h"
#include "extended/region_node.h"
#include "ltr/ltr_visitor.h"

struct GtLTRVisitor {
  const GtNodeVisitor parent_instance;
  GtLTRElement *element;
};

#define gt_ltr_visitor_cast(GV)\
        gt_node_visitor_cast(gt_ltr_visitor_class(), GV)

static int gt_ltr_visitor_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                       GT_UNUSED GtError *err)
{
  GtLTRVisitor *lv;
  GtRange node_range;
  GtArray *pdomarr = NULL;
  const char *pfamname;
  const char *fnt;
  lv = gt_ltr_visitor_cast(nv);
  gt_assert(lv);
  gt_error_check(err);

  fnt = gt_feature_node_get_type(fn);

  if (strcmp(fnt, gt_ft_LTR_retrotransposon) == 0)
  {
    lv->element->mainnode = fn;
  } else if (strcmp(fnt, gt_ft_long_terminal_repeat) == 0)
  {
    if (lv->element->leftLTR == NULL)
    {
      node_range = gt_genome_node_get_range((GtGenomeNode*) fn);
      lv->element->leftLTR = fn;
      /* compensate for 1-based node coords */
      lv->element->leftLTR_5 = node_range.start - 1;
      lv->element->leftLTR_3 = node_range.end - 1;
    }
    else
    {
      node_range = gt_genome_node_get_range((GtGenomeNode*) fn);
      lv->element->rightLTR = fn;
      /* compensate for 1-based node coords */
      lv->element->rightLTR_5 = node_range.start - 1;
      lv->element->rightLTR_3 = node_range.end - 1;
    }
  } else if (strcmp(fnt, gt_ft_target_site_duplication) == 0)
  {
    if (lv->element->leftTSD == NULL)
    {
      lv->element->leftTSD = fn;
    }
    else
    {
      lv->element->rightTSD = fn;
    }
  } else if (strcmp(fnt, gt_ft_RR_tract) == 0)
  {
    if (lv->element->ppt == NULL)
    {
      lv->element->ppt = fn;
    }
  } else if (strcmp(fnt, gt_ft_primer_binding_site) == 0)
  {
    if (lv->element->pbs == NULL)
    {
      lv->element->pbs = fn;
    }
  } else if (strcmp(fnt, gt_ft_protein_match) == 0)
  {
    if (!lv->element->pdoms)
    {
      lv->element->pdoms = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                          (GtFree) gt_array_delete);
    }
    pfamname = gt_feature_node_get_attribute(fn, "name");
    if (!(pdomarr = (GtArray*) gt_hashmap_get(lv->element->pdoms, pfamname)))
    {
      char *pfamcpy = gt_cstr_dup(pfamname);
      pdomarr = gt_array_new(sizeof (GtFeatureNode*));
      gt_hashmap_add(lv->element->pdoms, pfamcpy, pdomarr);
      if (lv->element->pdomorder != NULL)
        gt_array_add(lv->element->pdomorder, pfamcpy);
    }
    gt_array_add(pdomarr, fn);
  }
  return 0;
}

const GtNodeVisitorClass* gt_ltr_visitor_class(void)
{
  static const GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtLTRVisitor),
                                    NULL,
                                    NULL,
                                    gt_ltr_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_ltr_visitor_new(GtLTRElement *element)
{
  GtNodeVisitor *nv;
  GtLTRVisitor *lv;
  gt_assert(element);
  nv = gt_node_visitor_create(gt_ltr_visitor_class());
  lv = gt_ltr_visitor_cast(nv);
  lv->element = element;
  gt_assert(lv);
  return nv;
}
