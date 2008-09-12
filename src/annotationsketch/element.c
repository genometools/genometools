/*
  Copyright (c) 2007      Christin Schaerfer <cschaerfer@zbh.uni-hamburg.de>
  Copyright (c)      2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
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
#include "core/strand.h"
#include "extended/feature_type.h"
#include "extended/genome_feature.h"

struct GtElement {
  const char *type;
  GtStrand strand;
  GtGenomeNode *gn;
  GtRange range;
  GT_DrawingRange drange;
  bool mark;
};

GtElement* gt_element_new(GtGenomeNode *gn)
{
  GtElement *element;
  GT_GenomeFeature *gf = (GT_GenomeFeature*) gn;
  assert(gn);
  element = gt_element_new_empty();
  gt_element_set_type(element, gt_genome_feature_get_type(gf));
  gt_element_set_range(element, gt_genome_node_get_range(gn));
  element->strand = gt_genome_feature_get_strand(gf);
  element->mark = gt_genome_node_is_marked(gn);
  element->gn = gt_genome_node_ref(gn);
  return element;
}

GtElement* gt_element_new_empty(void)
{
  return gt_calloc(1, sizeof (GtElement));
}

GT_DrawingRange gt_element_calculate_drawing_range(GtElement *element,
                                             GtCanvas *canvas)
{
  assert(element && canvas);
  element->drange = gt_canvas_convert_coords(canvas, element->range);
  return element->drange;
}

GtRange gt_element_get_range(const GtElement *element)
{
  assert(element);
  return element->range;
}

void gt_element_set_range(GtElement *element, GtRange r)
{
  assert(element);
  element->range = r;
}

const char* gt_element_get_type(const GtElement *element)
{
  assert(element);
  return element->type;
}

void gt_element_set_type(GtElement *element, const char *type)
{
  assert(element);
  element->type = type;
}

GtStrand gt_element_get_strand(const GtElement *element)
{
  assert(element);
  return element->strand;
}
bool gt_element_is_marked(const GtElement *element)
{
  assert(element);
  return element->mark;
}

bool elements_are_equal(const GtElement *e1, const GtElement *e2)
{
  assert(e1 && e2);
  if (e1->type == e2->type && !gt_range_compare(e1->range, e2->range))
    return true;
  return false;
}

int gt_element_sketch(GtElement *elem, GtCanvas *canvas)
{
  int had_err = 0;
  assert(elem && canvas);
  gt_canvas_visit_element(canvas, elem);
  return had_err;
}

GtGenomeNode* gt_element_get_node_ref(const GtElement *elem)
{
  assert(elem);
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
  gn = gt_genome_feature_new(seqid, gft_exon, r1.start, r1.end, GT_STRAND_BOTH);
  gn2 = gt_genome_feature_new(seqid, gft_exon, r2.start, r2.end,
                             GT_STRAND_BOTH);

  e = gt_element_new(gn);
  e2 = gt_element_new(gn);
  e3 = gt_element_new(gn2);

  /* tests gt_element_get_range */
  r_temp = gt_element_get_range(e);
  ensure(had_err, (0 == gt_range_compare(r1, r_temp)));
  ensure(had_err, (1 == gt_range_compare(r2, r_temp)));

  /* tests gt_element_get_type and gt_element_set_type*/
  ensure(had_err, !strcmp(gft_exon, gt_element_get_type(e)));
  ensure(had_err, strcmp(gft_intron, gt_element_get_type(e)));
  gt_element_set_type(e, gft_intron);
  ensure(had_err, !strcmp(gft_intron, gt_element_get_type(e)));
  gt_element_set_type(e2, gft_intron);

  /* tests elements_are_equal */
  ensure(had_err, elements_are_equal(e, e2));
  ensure(had_err, !elements_are_equal(e, e3));
  ensure(had_err, !elements_are_equal(e2, e3));

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
  if (element->gn)
    gt_genome_node_delete(element->gn);
  gt_free(element);
}
