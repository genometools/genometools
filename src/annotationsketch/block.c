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
#include "core/ensure.h"
#include "core/log.h"
#include "core/ma.h"
#include "extended/feature_type_factory_builtin.h"
#include "annotationsketch/block.h"
#include "annotationsketch/element.h"

struct GT_Block {
  Dlist *elements;
  GT_Range range;
  Str *caption;
  bool show_caption;
  Strand strand;
  GT_GenomeFeatureType *type;
  GT_GenomeNode *top_level_feature;
  unsigned long reference_count;
};

/* Compare function used to insert Elements into dlist, order by type */
static int elemcmp(const void *a, const void *b)
{
  Element *elem_a = (Element*) a;
  Element *elem_b = (Element*) b;

  GT_GenomeFeatureType *ta = element_get_type(elem_a);
  GT_GenomeFeatureType *tb = element_get_type(elem_b);

  if (ta == tb)
    return 0;
  else if (strcmp(genome_feature_type_get_cstr(ta),
                  genome_feature_type_get_cstr(tb)) < 0) {
    return 1;
  }
  return -1;
}

int gt_block_compare(const GT_Block *block1, const GT_Block *block2)
{
  int ret;
  assert(block1 && block2);
  ret = gt_range_compare(gt_block_get_range(block1), gt_block_get_range(block2));
  if (ret == 0 && block1 != block2)
    ret = (block1 < block2 ? -1 : 1);
  return ret;
}

GT_Block* gt_block_ref(GT_Block *block)
{
  assert(block);
  block->reference_count++;
  return block;
}

GT_Block* gt_block_new(void)
{
  GT_Block *block = ma_calloc(1, sizeof (GT_Block));
  block->elements = dlist_new(elemcmp);
  block->caption = NULL;
  block->show_caption = true;
  block->strand = STRAND_UNKNOWN;
  block->top_level_feature = NULL;
  return block;
}

GT_Block* gt_block_new_from_node(GT_GenomeNode *node)
{
  GT_Block *block;
  assert(node);
  block = gt_block_new();
  block->range = gt_genome_node_get_range(node);
  block->strand = genome_feature_get_strand((GT_GenomeFeature*) node);
  block->type = genome_feature_get_type((GT_GenomeFeature*) node);
  block->top_level_feature = gt_genome_node_ref(node);
  return block;
}

void gt_block_insert_element(GT_Block *block, GT_GenomeNode *gn)
{
  Element *element;
  assert(block && gn);
  if (!block->top_level_feature)
    block->top_level_feature = gt_genome_node_ref(gn);
  element = element_new(gn);
  dlist_add(block->elements, element);
}

GT_GenomeNode* gt_block_get_top_level_feature(const GT_Block *block)
{
  assert(block);
  return block->top_level_feature;
}

GT_Range gt_block_get_range(const GT_Block *block)
{
   assert(block);
   return block->range;
}

GT_Range* gt_block_get_range_ptr(const GT_Block *block)
{
   assert(block);
   return (GT_Range*) &(block->range);
}

void gt_block_set_range(GT_Block *block, GT_Range r)
{
  assert(block && r.start <= r.end);
  block->range = r;
}

bool gt_block_has_only_one_fullsize_element(const GT_Block *block)
{
  bool ret = false;
  assert(block);
  if (dlist_size(block->elements) == 1UL) {
    GT_Range elem_range, block_range;
    assert(dlist_first(block->elements) == dlist_last(block->elements));
    elem_range = element_get_range(
                   dlistelem_get_data(dlist_first(block->elements)));
    block_range = gt_block_get_range(block);
    ret = (gt_range_compare(block_range, elem_range) == 0);
  }
  return ret;
}

void gt_block_set_caption_visibility(GT_Block *block, bool val)
{
  assert(block);
  block->show_caption = val;
}

bool gt_block_caption_is_visible(const GT_Block *block)
{
  assert(block);
  return (block->caption && block->show_caption);
}

void gt_block_set_caption(GT_Block *block, Str *caption)
{
  assert(block);
  block->caption = caption;
}

Str* gt_block_get_caption(const GT_Block *block)
{
  assert(block);
  return block->caption;
}

void gt_block_set_strand(GT_Block *block, Strand strand)
{
  assert(block);
  block->strand = strand;
}

Strand gt_block_get_strand(const GT_Block *block)
{
  assert(block);
  return block->strand;
}

void gt_block_set_type(GT_Block *block, GT_GenomeFeatureType *type)
{
  assert(block);
  block->type = type;
}

GT_GenomeFeatureType* gt_block_get_type(const GT_Block *block)
{
  assert(block);
  return block->type;
}

unsigned long gt_block_get_size(const GT_Block *block)
{
  assert(block && block->elements);
  return dlist_size(block->elements);
}

int gt_block_sketch(GT_Block *block, GT_Canvas *canvas)
{
 int had_err = 0;
 Dlistelem *delem;
 assert(block && canvas);
 /* if resulting block was too short,
    do not traverse this feature tree further */
 if (-1 == gt_canvas_visit_block(canvas, block))
   return had_err;
 for (delem = dlist_first(block->elements); delem;
      delem = dlistelem_next(delem)) {
    Element* elem = (Element*) dlistelem_get_data(delem);
    element_sketch(elem, canvas);
  }
  return had_err;
}

int gt_block_unit_test(GT_Error *err)
{
  FeatureTypeFactory *feature_type_factory;
  GT_GenomeFeatureType *gft;
  GT_Range r1, r2, r_temp, b_range;
  int had_err = 0;
  Strand s;
  GT_GenomeNode *gn1, *gn2;
  Element *e1, *e2;
  GT_Block * b;
  Str *seqid, *caption1, *caption2;
  gt_error_check(err);

  feature_type_factory = feature_type_factory_builtin_new();
  seqid = str_new_cstr("seqid");
  caption1 = str_new_cstr("foo");
  caption2 = str_new_cstr("bar");

  r1.start = 10UL;
  r1.end = 50UL;

  r2.start = 40UL;
  r2.end = 50UL;

  gft = feature_type_factory_create_gft(feature_type_factory, gft_gene);
  gn1 = genome_feature_new(seqid, gft, r1, STRAND_FORWARD);
  gft = feature_type_factory_create_gft(feature_type_factory, gft_exon);
  gn2 = genome_feature_new(seqid, gft, r2, STRAND_FORWARD);

  e1 = element_new(gn1);
  e2 = element_new(gn2);

  b = gt_block_new();

  /* test gt_block_insert_elements */
  ensure(had_err, (0UL == gt_block_get_size(b)));
  gt_block_insert_element(b, gn1);
  ensure(had_err, (1UL == gt_block_get_size(b)));
  gt_block_insert_element(b, gn2);
  ensure(had_err, (2UL == gt_block_get_size(b)));

  /* test gt_block_set_range & gt_block_get_range */
  r_temp = gt_range_join(r1, r2);
  gt_block_set_range(b, r_temp);
  b_range = gt_block_get_range(b);
  ensure(had_err, (0 == gt_range_compare(b_range, r_temp)));
  ensure(had_err, (1 == gt_range_compare(r2, r_temp)));

  /* tests gt_block_set_caption & gt_block_get_caption */
  gt_block_set_caption(b, caption1);
  ensure(had_err, (0 == str_cmp(gt_block_get_caption(b), caption1)));
  ensure(had_err, (0 != str_cmp(gt_block_get_caption(b), caption2)));

  /* tests gt_block_set_strand & gt_block_get_range */
  s = gt_block_get_strand(b);
  ensure(had_err, (STRAND_UNKNOWN == s));
  gt_block_set_strand(b, STRAND_FORWARD);
  s = gt_block_get_strand(b);
  ensure(had_err, (STRAND_FORWARD == s));

  str_delete(caption2);
  str_delete(seqid);
  element_delete(e1);
  element_delete(e2);
  gt_block_delete(b);
  gt_genome_node_delete(gn1);
  gt_genome_node_delete(gn2);
  feature_type_factory_delete(feature_type_factory);

  return had_err;
}

void gt_block_delete(GT_Block *block)
{
  Dlistelem *delem;
  if (!block) return;
  if (block->reference_count) {
    block->reference_count--;
    return;
  }
  for (delem = dlist_first(block->elements); delem;
       delem = dlistelem_next(delem)) {
    Element* elem = (Element*) dlistelem_get_data(delem);
    element_delete(elem);
  }
  if (block->caption)
    str_delete(block->caption);
  dlist_delete(block->elements);
  if (block->top_level_feature)
    gt_genome_node_delete(block->top_level_feature);
  ma_free(block);
}
