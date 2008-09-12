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
#include "annotationsketch/block.h"
#include "annotationsketch/element.h"
#include "core/ensure.h"
#include "core/log.h"
#include "core/ma.h"

struct GtBlock {
  GtDlist *elements;
  GtRange range;
  GtStr *caption;
  bool show_caption;
  GtStrand strand;
  const char *type;
  GtGenomeNode *top_level_feature;
  unsigned long reference_count;
};

/* GT_Compare function used to insert GtElements into dlist, order by type */
static int elemcmp(const void *a, const void *b)
{
  const char *type_a, *type_b;
  GtElement *elem_a = (GtElement*) a;
  GtElement *elem_b = (GtElement*) b;

  type_a = gt_element_get_type(elem_a);
  type_b = gt_element_get_type(elem_b);

  if (type_a == type_b)
    return 0;
  if (strcmp(type_a, type_b) < 0)
    return 1;
  return -1;
}

int gt_block_compare(const GtBlock *block1, const GtBlock *block2)
{
  int ret;
  assert(block1 && block2);
  ret = gt_range_compare(gt_block_get_range(block1),
                         gt_block_get_range(block2));
  if (ret == 0 && block1 != block2)
    ret = (block1 < block2 ? -1 : 1);
  return ret;
}

GtBlock* gt_block_ref(GtBlock *block)
{
  assert(block);
  block->reference_count++;
  return block;
}

GtBlock* gt_block_new(void)
{
  GtBlock *block = gt_calloc(1, sizeof (GtBlock));
  block->elements = gt_dlist_new(elemcmp);
  block->caption = NULL;
  block->show_caption = true;
  block->strand = GT_STRAND_UNKNOWN;
  block->top_level_feature = NULL;
  return block;
}

GtBlock* gt_block_new_from_node(GtGenomeNode *node)
{
  GtBlock *block;
  assert(node);
  block = gt_block_new();
  block->range = gt_genome_node_get_range(node);
  block->strand = gt_genome_feature_get_strand((GtGenomeFeature*) node);
  block->type = gt_genome_feature_get_type((GtGenomeFeature*) node);
  block->top_level_feature = gt_genome_node_ref(node);
  return block;
}

void gt_block_insert_element(GtBlock *block, GtGenomeNode *gn)
{
  GtElement *element;
  assert(block && gn);
  if (!block->top_level_feature)
    block->top_level_feature = gt_genome_node_ref(gn);
  element = gt_element_new(gn);
  gt_dlist_add(block->elements, element);
}

GtGenomeNode* gt_block_get_top_level_feature(const GtBlock *block)
{
  assert(block);
  return block->top_level_feature;
}

GtRange gt_block_get_range(const GtBlock *block)
{
   assert(block);
   return block->range;
}

GtRange* gt_block_get_range_ptr(const GtBlock *block)
{
   assert(block);
   return (GtRange*) &(block->range);
}

void gt_block_set_range(GtBlock *block, GtRange r)
{
  assert(block && r.start <= r.end);
  block->range = r;
}

bool gt_block_has_only_one_fullsize_element(const GtBlock *block)
{
  bool ret = false;
  assert(block);
  if (gt_dlist_size(block->elements) == 1UL) {
    GtRange elem_range, block_range;
    assert(gt_dlist_first(block->elements) == gt_dlist_last(block->elements));
    elem_range = gt_element_get_range(
                   gt_dlistelem_get_data(gt_dlist_first(block->elements)));
    block_range = gt_block_get_range(block);
    ret = (gt_range_compare(block_range, elem_range) == 0);
  }
  return ret;
}

void gt_block_set_caption_visibility(GtBlock *block, bool val)
{
  assert(block);
  block->show_caption = val;
}

bool gt_block_caption_is_visible(const GtBlock *block)
{
  assert(block);
  return (block->caption && block->show_caption);
}

void gt_block_set_caption(GtBlock *block, GtStr *caption)
{
  assert(block);
  block->caption = caption;
}

GtStr* gt_block_get_caption(const GtBlock *block)
{
  assert(block);
  return block->caption;
}

void gt_block_set_strand(GtBlock *block, GtStrand strand)
{
  assert(block);
  block->strand = strand;
}

GtStrand gt_block_get_strand(const GtBlock *block)
{
  assert(block);
  return block->strand;
}

void gt_block_set_type(GtBlock *block, const char *type)
{
  assert(block);
  block->type = type;
}

const char* gt_block_get_type(const GtBlock *block)
{
  assert(block);
  return block->type;
}

unsigned long gt_block_get_size(const GtBlock *block)
{
  assert(block && block->elements);
  return gt_dlist_size(block->elements);
}

int gt_block_sketch(GtBlock *block, GtCanvas *canvas)
{
 int had_err = 0;
 GtDlistelem *delem;
 assert(block && canvas);
 /* if resulting block was too short,
    do not traverse this feature tree further */
 if (-1 == gt_canvas_visit_block(canvas, block))
   return had_err;
 for (delem = gt_dlist_first(block->elements); delem;
      delem = gt_dlistelem_next(delem)) {
    GtElement* elem = (GtElement*) gt_dlistelem_get_data(delem);
    gt_element_sketch(elem, canvas);
  }
  return had_err;
}

int gt_block_unit_test(GtError *err)
{
  GtRange r1, r2, r_temp, b_range;
  GtStrand s;
  GtGenomeNode *gn1, *gn2;
  GtElement *e1, *e2;
  GtBlock * b;
  GtStr *seqid, *caption1, *caption2;
  int had_err = 0;
  gt_error_check(err);

  seqid = gt_str_new_cstr("seqid");
  caption1 = gt_str_new_cstr("foo");
  caption2 = gt_str_new_cstr("bar");

  r1.start = 10UL;
  r1.end = 50UL;

  r2.start = 40UL;
  r2.end = 50UL;

  gn1 = gt_genome_feature_new(seqid, gft_gene, r1.start, r1.end,
                              GT_STRAND_FORWARD);
  gn2 = gt_genome_feature_new(seqid, gft_exon, r2.start, r2.end,
                              GT_STRAND_FORWARD);

  e1 = gt_element_new(gn1);
  e2 = gt_element_new(gn2);

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
  ensure(had_err, (0 == gt_str_cmp(gt_block_get_caption(b), caption1)));
  ensure(had_err, (0 != gt_str_cmp(gt_block_get_caption(b), caption2)));

  /* tests gt_block_set_strand & gt_block_get_range */
  s = gt_block_get_strand(b);
  ensure(had_err, (GT_STRAND_UNKNOWN == s));
  gt_block_set_strand(b, GT_STRAND_FORWARD);
  s = gt_block_get_strand(b);
  ensure(had_err, (GT_STRAND_FORWARD == s));

  gt_str_delete(caption2);
  gt_str_delete(seqid);
  gt_element_delete(e1);
  gt_element_delete(e2);
  gt_block_delete(b);
  gt_genome_node_delete(gn1);
  gt_genome_node_delete(gn2);

  return had_err;
}

void gt_block_delete(GtBlock *block)
{
  GtDlistelem *delem;
  if (!block) return;
  if (block->reference_count) {
    block->reference_count--;
    return;
  }
  for (delem = gt_dlist_first(block->elements); delem;
       delem = gt_dlistelem_next(delem)) {
    GtElement* elem = (GtElement*) gt_dlistelem_get_data(delem);
    gt_element_delete(elem);
  }
  if (block->caption)
    gt_str_delete(block->caption);
  gt_dlist_delete(block->elements);
  if (block->top_level_feature)
    gt_genome_node_delete(block->top_level_feature);
  gt_free(block);
}
