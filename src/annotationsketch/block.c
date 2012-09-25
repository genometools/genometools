/*
  Copyright (c) 2007      Christin Schaerfer <schaerfer@zbh.uni-hamburg.de>
  Copyright (c) 2008-2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg

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
#include "annotationsketch/default_formats.h"
#include "annotationsketch/element.h"
#include "core/array.h"
#include "core/cstr_api.h"
#include "core/ensure.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/undef_api.h"
#include "core/unused_api.h"

struct GtBlock {
  GtArray *elements;
  GtRange range;
  GtStr *caption,
        *track_id;
  bool show_caption,
       sorted;
  GtStrand strand;
  const char *type;
  GtFeatureNode *top_level_feature;
  unsigned long reference_count;
};

/* This function orders GtElements by z-index or type. This enables the sketch
   function to paint elements in a specific order based on their type
   (painter's algorithm-like) */
static int elemcmp(const void *a, const void *b, void *data)
{
  const char *type_a, *type_b;
  double zindex_a = GT_UNDEF_DOUBLE, zindex_b= GT_UNDEF_DOUBLE;
  GtStyle *sty = (GtStyle*) data;
  GtElement *elem_a = *(GtElement**) a;
  GtElement *elem_b = *(GtElement**) b;

  type_a = gt_element_get_type(elem_a);
  type_b = gt_element_get_type(elem_b);

  /* same types are equal, no further action needed */
  if (type_a == type_b)
    return 0;
  /* if not, then get z-index from style */
  if (sty)
  {
    (void) gt_style_get_num(sty, type_a, "z_index", &zindex_a, NULL, NULL);
    (void) gt_style_get_num(sty, type_b, "z_index", &zindex_b, NULL, NULL);
  }
  /* only one is set -> put types with set z-indices always on top of others*/
  if (zindex_a == GT_UNDEF_DOUBLE && zindex_b != GT_UNDEF_DOUBLE)
    return -1;
  if (zindex_b == GT_UNDEF_DOUBLE && zindex_a != GT_UNDEF_DOUBLE)
    return 1;
  /* none is set, fall back to default alphabetic ordering */
  if (zindex_a == GT_UNDEF_DOUBLE && zindex_b == GT_UNDEF_DOUBLE)
  {
    if (strcmp(type_a, type_b) < 0)
      return 1;
    return -1;
  }
  /* both z-indices are set */
  if (gt_double_equals_double(zindex_a, zindex_b))
    return 0;
  if (gt_double_smaller_double(zindex_a, zindex_b))
    return -1;
  return 1;
}

int gt_block_compare(const GtBlock *block1, const GtBlock *block2,
                     GT_UNUSED void *data)
{
  GtRange range_a, range_b;
  int ret = 0;
  gt_assert(block1 && block2);
  range_a = gt_block_get_range(block1),
  range_b = gt_block_get_range(block2);
  ret = gt_range_compare(&range_a, &range_b);
  if (ret == 0 && block1 != block2) {
    GtStr *caption1, *caption2;
    caption1 = gt_block_get_caption(block1);
    caption2 = gt_block_get_caption(block2);
    /* blocks do not necessarily have captions. If both have a caption, we
       compare them. If only one block has a caption, this block comes first. */
    if (caption1 && caption2)
      ret = strcmp(gt_str_get(caption1), gt_str_get(caption2));
    else if (caption1)
      ret = -1;
    else if (caption2)
      ret = 1;
  }
  return ret;
}

GtBlock* gt_block_ref(GtBlock *block)
{
  gt_assert(block);
  block->reference_count++;
  return block;
}

void gt_block_print(const GtBlock* block)
{
  unsigned long i;
  gt_assert(block);
  for (i=0;i<gt_array_size(block->elements);i++)
  {
    GtElement *elem;
    elem = gt_element_ref(*(GtElement**) gt_array_get(block->elements, i));
    gt_assert(elem);
    GtRange r = gt_element_get_range(elem);
    printf("%s\t%lu-%lu\n", gt_element_get_type(elem),
                           r.start,
                           r.end);
  }
}

GtBlock* gt_block_new(void)
{
  GtBlock *block = gt_calloc(1, sizeof (GtBlock));
  block->elements = gt_array_new(sizeof (GtElement*));
  block->caption = NULL;
  block->show_caption = true;
  block->sorted = false;
  block->strand = GT_STRAND_UNKNOWN;
  block->top_level_feature = NULL;
  return block;
}

GtBlock* gt_block_new_from_node(GtFeatureNode *node)
{
  GtBlock *block;
  gt_assert(node);
  block = gt_block_new();
  block->range = gt_genome_node_get_range((GtGenomeNode*) node);
  block->strand = gt_feature_node_get_strand(node);
  block->type = gt_feature_node_get_type(node);
  if (!gt_feature_node_is_pseudo(node)) {
    block->top_level_feature = (GtFeatureNode*)
                               gt_genome_node_ref((GtGenomeNode*) node);
  }
  return block;
}

void gt_block_insert_element(GtBlock *block, GtFeatureNode *node)
{
  GtElement *element;
  gt_assert(block && node);
  if (!block->top_level_feature) {
    block->top_level_feature = (GtFeatureNode*)
                               gt_genome_node_ref((GtGenomeNode*) node);
  }
  element = gt_element_new(node);
  /* invalidate sortedness flag because insertion of element at the end
     may break ordering */
  block->sorted = false;
  gt_array_add(block->elements, element);
}

void gt_block_merge(GtBlock *b1, GtBlock *b2)
{
  unsigned int GT_UNUSED merged_size, i;
  gt_assert(b1 && b2);
  merged_size = gt_block_get_size(b1) + gt_block_get_size(b2);
  for (i=0;i<gt_array_size(b2->elements);i++)
  {
    GtElement *elem;
    elem = gt_element_ref(*(GtElement**) gt_array_get(b2->elements, i));
    gt_assert(elem);
    gt_array_add(b1->elements, elem);
  }
  gt_assert(gt_block_get_size(b1) == merged_size);
}

GtBlock* gt_block_clone(GtBlock *block)
{
  GtBlock* newblock;
  unsigned long i;
  gt_assert(block);
  newblock = gt_block_new();
  for (i=0;i<gt_array_size(block->elements);i++)
  {
    GtElement *elem;
    elem = gt_element_ref(*(GtElement**) gt_array_get(block->elements, i));
    gt_assert(elem);
    gt_array_add(newblock->elements, elem);
  }
  gt_assert(gt_block_get_size(newblock) == gt_block_get_size(block));
  newblock->caption = gt_str_ref(block->caption);
  newblock->type = block->type;
  newblock->range.start = block->range.start;
  newblock->range.end = block->range.end;
  newblock->show_caption = block->show_caption;
  newblock->strand = block->strand;
  newblock->top_level_feature = (GtFeatureNode*)
                                gt_genome_node_ref((GtGenomeNode*)
                                                   block->top_level_feature);
  return newblock;
}

GtFeatureNode* gt_block_get_top_level_feature(const GtBlock *block)
{
  gt_assert(block);
  return block->top_level_feature;
}

GtRange gt_block_get_range(const GtBlock *block)
{
   gt_assert(block);
   return block->range;
}

GtRange* gt_block_get_range_ptr(const GtBlock *block)
{
   gt_assert(block);
   return (GtRange*) &(block->range);
}

void gt_block_set_range(GtBlock *block, GtRange r)
{
  gt_assert(block && r.start <= r.end);
  block->range = r;
}

bool gt_block_has_only_one_fullsize_element(const GtBlock *block)
{
  bool ret = false;
  unsigned long bsize;
  gt_assert(block);
  bsize = gt_array_size(block->elements);
  if (bsize == 1) {
    GtRange elem_range, block_range;
    gt_assert(*(GtElement**) gt_array_get(block->elements, 0) ==
              *(GtElement**) gt_array_get(block->elements, bsize-1));
    elem_range = gt_element_get_range(*(GtElement**)
                                             gt_array_get(block->elements, 0));
    block_range = gt_block_get_range(block);
    ret = (gt_range_compare(&block_range, &elem_range) == 0);
  }
  return ret;
}

void gt_block_set_caption_visibility(GtBlock *block, bool val)
{
  gt_assert(block);
  block->show_caption = val;
}

bool gt_block_caption_is_visible(const GtBlock *block)
{
  gt_assert(block);
  return (block->caption && block->show_caption);
}

void gt_block_set_caption(GtBlock *block, GtStr *caption)
{
  gt_assert(block);
  block->caption = caption;
}

GtStr* gt_block_get_caption(const GtBlock *block)
{
  gt_assert(block);
  return block->caption;
}

void gt_block_set_strand(GtBlock *block, GtStrand strand)
{
  gt_assert(block);
  block->strand = strand;
}

GtStrand gt_block_get_strand(const GtBlock *block)
{
  gt_assert(block);
  return block->strand;
}

void gt_block_set_type(GtBlock *block, const char *type)
{
  gt_assert(block);
  block->type = type;
}

const char* gt_block_get_type(const GtBlock *block)
{
  gt_assert(block);
  return block->type;
}

int gt_block_get_max_height(const GtBlock *block, double *result,
                            const GtStyle *sty, GtError *err)
{
  unsigned long max_height = 0;
  GtStyleQueryStatus rval = GT_STYLE_QUERY_OK;
  unsigned long i;
  gt_assert(block && sty);
  for (i=0;i<gt_array_size(block->elements);i++) {
    GtElement *elem;
    double height = 0;
    elem = *(GtElement**) gt_array_get(block->elements, i);
    /* get default or image-wide bar height */
    rval = gt_style_get_num(sty, "format", "bar_height", &height, NULL, err);
    switch (rval) {
    case GT_STYLE_QUERY_ERROR:
      return -1;
      break; /* should never reach this */
    case GT_STYLE_QUERY_NOT_SET:
      height = BAR_HEIGHT_DEFAULT;
      break;
    default:
      /* do nothing */
      break;
    }
    /* try to get type-specific bar height */
    rval = gt_style_get_num(sty, gt_element_get_type(elem), "bar_height",
                            &height, gt_element_get_node_ref(elem),
                            err);
    if (rval == GT_STYLE_QUERY_ERROR) {
      return -1;
    }
    if (gt_double_smaller_double(max_height, height))
      max_height = height;
  }
  *result = max_height;
  return 0;
}

unsigned long gt_block_get_size(const GtBlock *block)
{
  gt_assert(block && block->elements);
  return gt_array_size(block->elements);
}

int gt_block_sketch(GtBlock *block, GtCanvas *canvas, GtError *err)
{
  int had_err = 0;
  unsigned long i;
  gt_assert(block && canvas && err);
  /* if resulting block was too short,
     do not traverse this feature tree further */
  had_err = gt_canvas_visit_block(canvas, block, err);
  if (had_err)
  {
    if (had_err == 1)  /* 'magic' value */
      return 0;
    else
      return had_err;
  } /* we have in any case returned if had_err was set */
  if (!block->sorted)
  {
    GtStyle *sty = gt_canvas_get_style(canvas);
    /* sort elements if they have changed since last sketch operation */
    gt_array_sort_with_data(block->elements, elemcmp, sty);
    block->sorted = true;
  }
  /* delegate sketch request to elements */
  for (i=0;i<gt_array_size(block->elements);i++) {
     GtElement *elem = *(GtElement**) gt_array_get(block->elements, i);
     had_err = gt_element_sketch(elem, canvas, err);
  }
  return had_err;
}

int gt_block_unit_test(GtError *err)
{
  GtRange r1, r2, r_temp, b_range;
  GtStrand s;
  GtGenomeNode *gn1, *gn2;
  GtElement *e1, *e2;
  double height;
  GtBlock *b;
  GtStr *seqid, *caption1, *caption2;
  int had_err = 0;
  GtStyle *sty;
  GtError *testerr;
  gt_error_check(err);

  seqid = gt_str_new_cstr("seqid");
  caption1 = gt_str_new_cstr("foo");
  caption2 = gt_str_new_cstr("bar");
  testerr = gt_error_new();

  r1.start = 10UL;
  r1.end = 50UL;

  r2.start = 40UL;
  r2.end = 50UL;

  gn1 = gt_feature_node_new(seqid, gt_ft_gene, r1.start, r1.end,
                            GT_STRAND_FORWARD);
  gn2 = gt_feature_node_new(seqid, gt_ft_exon, r2.start, r2.end,
                            GT_STRAND_FORWARD);

  e1 = gt_element_new((GtFeatureNode*) gn1);
  e2 = gt_element_new((GtFeatureNode*) gn2);

  b = gt_block_new();

  /* test gt_block_insert_elements */
  gt_ensure(had_err, (0UL == gt_block_get_size(b)));
  gt_block_insert_element(b, (GtFeatureNode*) gn1);
  gt_ensure(had_err, (1UL == gt_block_get_size(b)));
  gt_block_insert_element(b, (GtFeatureNode*) gn2);
  gt_ensure(had_err, (2UL == gt_block_get_size(b)));

  /* test gt_block_set_range & gt_block_get_range */
  r_temp = gt_range_join(&r1, &r2);
  gt_block_set_range(b, r_temp);
  b_range = gt_block_get_range(b);
  gt_ensure(had_err, (0 == gt_range_compare(&b_range, &r_temp)));
  gt_ensure(had_err, (1 == gt_range_compare(&r2, &r_temp)));

  /* tests gt_block_set_caption & gt_block_get_caption */
  gt_block_set_caption(b, caption1);
  gt_ensure(had_err, (0 == gt_str_cmp(gt_block_get_caption(b), caption1)));
  gt_ensure(had_err, (0 != gt_str_cmp(gt_block_get_caption(b), caption2)));

  /* tests gt_block_set_strand & gt_block_get_range */
  s = gt_block_get_strand(b);
  gt_ensure(had_err, (GT_STRAND_UNKNOWN == s));
  gt_block_set_strand(b, GT_STRAND_FORWARD);
  s = gt_block_get_strand(b);
  gt_ensure(had_err, (GT_STRAND_FORWARD == s));

  /* test gt_block_get_max_height() */
  sty = gt_style_new(err);
  gt_ensure(had_err, gt_block_get_max_height(b, &height, sty, err) == 0);
  gt_ensure(had_err, !gt_error_is_set(testerr));
  gt_ensure(had_err, height == BAR_HEIGHT_DEFAULT);
  gt_style_set_num(sty, "exon", "bar_height", 42);
  gt_ensure(had_err, gt_block_get_max_height(b, &height, sty, err) == 0);
  gt_ensure(had_err, !gt_error_is_set(testerr));
  gt_ensure(had_err, height == 42);
  gt_style_set_num(sty, "gene", "bar_height", 23);
  gt_ensure(had_err, gt_block_get_max_height(b, &height, sty, err) == 0);
  gt_ensure(had_err, !gt_error_is_set(testerr));
  gt_ensure(had_err, height == 42);
  gt_style_unset(sty, "exon", "bar_height");
  gt_ensure(had_err, gt_block_get_max_height(b, &height, sty, err) == 0);
  gt_ensure(had_err, !gt_error_is_set(testerr));
  gt_ensure(had_err, height == 23);

  gt_str_delete(caption2);
  gt_str_delete(seqid);
  gt_element_delete(e1);
  gt_element_delete(e2);
  gt_block_delete(b);
  gt_style_delete(sty);
  gt_error_delete(testerr);
  gt_genome_node_delete(gn1);
  gt_genome_node_delete(gn2);

  return had_err;
}

void gt_block_delete(GtBlock *block)
{
  unsigned long i;
  if (!block) return;
  if (block->reference_count) {
    block->reference_count--;
    return;
  }
  for (i=0;i<gt_array_size(block->elements);i++) {
    GtElement *elem = *(GtElement**) gt_array_get(block->elements, i);
    gt_element_delete(elem);
  }
  if (block->caption)
    gt_str_delete(block->caption);
  gt_array_delete(block->elements);
  if (block->top_level_feature)
    gt_genome_node_delete((GtGenomeNode*) block->top_level_feature);
  gt_free(block);
}
