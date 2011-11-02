/*
  Copyright (c) 2007      Christin Schaerfer <schaerfer@zbh.uni-hamburg.de>
  Copyright (c)      2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#include <string.h>
#include "annotationsketch/element.h"
#include "annotationsketch/style.h"
#include "core/array.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/strand_api.h"
#include "extended/feature_node.h"

struct GtElement {
  const char *type;
  unsigned long refcount;
  GtStrand strand;
  GtFeatureNode *gn;
  GtRange range;
  GtDrawingRange drange;
  bool mark;
};

GtElement* gt_element_new(GtFeatureNode *node)
{
  GtElement *element;
  gt_assert(node);
  element = gt_element_new_empty();
  gt_element_set_type(element, gt_feature_node_get_type(node));
  gt_element_set_range(element, gt_genome_node_get_range((GtGenomeNode*) node));
  element->strand = gt_feature_node_get_strand(node);
  element->mark = gt_feature_node_is_marked(node);
  element->gn = (GtFeatureNode*) gt_genome_node_ref((GtGenomeNode*) node);
  return element;
}

GtElement* gt_element_ref(GtElement *elem)
{
  gt_assert(elem);
  elem->refcount++;
  return elem;
}

GtElement* gt_element_new_empty(void)
{
  return gt_calloc(1, sizeof (GtElement));
}

GtRange gt_element_get_range(const GtElement *element)
{
  gt_assert(element);
  return element->range;
}

void gt_element_set_range(GtElement *element, GtRange r)
{
  gt_assert(element);
  element->range = r;
}

const char* gt_element_get_type(const GtElement *element)
{
  gt_assert(element);
  return element->type;
}

void gt_element_set_type(GtElement *element, const char *type)
{
  gt_assert(element);
  element->type = type;
}

GtStrand gt_element_get_strand(const GtElement *element)
{
  gt_assert(element);
  return element->strand;
}
bool gt_element_is_marked(const GtElement *element)
{
  gt_assert(element);
  return element->mark;
}

static bool elements_are_equal(const GtElement *e1, const GtElement *e2)
{
  gt_assert(e1 && e2);
  if (e1->type == e2->type && !gt_range_compare(&e1->range, &e2->range))
    return true;
  return false;
}

int gt_element_sketch(GtElement *elem, GtCanvas *canvas, GtError *err)
{
  int had_err = 0;
  gt_assert(elem && canvas);
  had_err = gt_canvas_visit_element(canvas, elem, err);
  return had_err;
}

GtFeatureNode* gt_element_get_node_ref(const GtElement *elem)
{
  gt_assert(elem);
  return elem->gn;
}

int gt_element_unit_test(GtError *err)
{
  GtRange r1, r2, r_temp;
  GtGenomeNode *gn, *gn2;
  GtElement *e, *e2, *e3;
  GtStr *seqid;
  int had_err = 0;
  gt_error_check(err);

  r1.start = 10UL;
  r1.end = 50UL;

  r2.start = 20UL;
  r2.end = 50UL;

  seqid = gt_str_new_cstr("seqid");
  gn = gt_feature_node_new(seqid, gt_ft_exon, r1.start, r1.end, GT_STRAND_BOTH);
  gn2 = gt_feature_node_new(seqid, gt_ft_exon, r2.start, r2.end,
                             GT_STRAND_BOTH);

  e = gt_element_new((GtFeatureNode*) gn);
  e2 = gt_element_new((GtFeatureNode*)gn);
  e3 = gt_element_new((GtFeatureNode*)gn2);

  /* tests gt_element_get_range */
  r_temp = gt_element_get_range(e);
  gt_ensure(had_err, (0 == gt_range_compare(&r1, &r_temp)));
  gt_ensure(had_err, (1 == gt_range_compare(&r2, &r_temp)));

  /* tests gt_element_get_type and gt_element_set_type*/
  gt_ensure(had_err, !strcmp(gt_ft_exon, gt_element_get_type(e)));
  gt_ensure(had_err, strcmp(gt_ft_intron, gt_element_get_type(e)));
  gt_element_set_type(e, gt_ft_intron);
  gt_ensure(had_err, !strcmp(gt_ft_intron, gt_element_get_type(e)));
  gt_element_set_type(e2, gt_ft_intron);

  /* tests elements_are_equal */
  gt_ensure(had_err, elements_are_equal(e, e2));
  gt_ensure(had_err, !elements_are_equal(e, e3));
  gt_ensure(had_err, !elements_are_equal(e2, e3));

  gt_element_delete(e);
  gt_element_delete(e2);
  gt_element_delete(e3);
  gt_genome_node_delete(gn);
  gt_genome_node_delete(gn2);
  gt_str_delete(seqid);

  return had_err;

}

void gt_element_delete(GtElement *element)
{
  if (!element) return;
  if (element->refcount) {
    element->refcount--;
    return;
  }
  gt_genome_node_delete((GtGenomeNode*) element->gn);
  gt_free(element);
}
