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

#include "annotationsketch/block.h"
#include "annotationsketch/default_formats.h"
#include "annotationsketch/line.h"
#include "core/ensure.h"
#include "core/interval_tree.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/range.h"
#include "core/undef_api.h"
#include "extended/gff3_defines.h"

struct GtLine {
  bool has_captions;
  GtArray *blocks;
};

GtLine* gt_line_new(void)
{
  GtLine *line;
  line = gt_calloc(1, sizeof (GtLine));
  line->blocks = gt_array_new(sizeof (GtBlock*));
  line->has_captions = false;
  return line;
}

void gt_line_insert_block(GtLine *line, GtBlock *block)
{
  gt_assert(line && block);
  if (!line->has_captions && gt_block_get_caption(block) != NULL)
    line->has_captions = true;
  gt_array_add(line->blocks, block);
}

bool gt_line_has_captions(const GtLine *line)
{
  gt_assert(line);
  return line->has_captions;
}

GtArray* gt_line_get_blocks(GtLine* line)
{
  gt_assert(line);
  return line->blocks;
}

int gt_line_sketch(GtLine *line, GtCanvas *canvas, GtError *err)
{
  int i = 0, had_err = 0;
  gt_assert(line && canvas);
  had_err = gt_canvas_visit_line_pre(canvas, line, err);
  if (!had_err)
  {
    for (i = 0; !had_err && i < gt_array_size(line->blocks); i++) {
      GtBlock *block;
      block = *(GtBlock**) gt_array_get(line->blocks, i);
      had_err = gt_block_sketch(block, canvas, err);
    }
  }
  if (!had_err)
    had_err = gt_canvas_visit_line_post(canvas, line, err);
  return had_err;
}

int gt_line_get_height(GtLine *line, double *height, const GtStyle *sty,
                       GtError *err)
{
  double line_height = 0;
  unsigned long i;
  gt_assert(line && sty);
  for (i = 0; i < gt_array_size(line->blocks); i++) {
    GtBlock *block;
    double myheight = BAR_HEIGHT_DEFAULT;
    block = *(GtBlock**) gt_array_get(line->blocks, i);
    /* check again for caption presence, may have changed in the meantime */
    if (!line->has_captions && gt_block_get_caption(block) != NULL) {
      line->has_captions = true;
    }
    if (gt_block_get_max_height(block, &myheight, sty, err) < 0) {
      return -1;
    }
    if (gt_double_smaller_double(line_height, myheight)) {
      line_height = myheight;
    }
  }
  *height = line_height;
  return 0;
}

int gt_line_unit_test(GtError *err)
{
  GtRange r1, r2, r3;
  GtArray* blocks;
  GtStr *seqid1;
  int had_err = 0;
  GtGenomeNode *parent, *gn1, *gn2, *gn3;
  double height;
  GtLine *l1;
  GtBlock *b1, *b2;
  GtStyle *sty = NULL;
  const char* foo = "foo";
  const char* bar = "bar";
  const char* blub = "blub";
  gt_error_check(err);

  r1.start = 10;
  r1.end = 40;
  r2.start = 51;
  r2.end = 80;
  r3.start = 0;
  r3.end = 7;

  seqid1 = gt_str_new_cstr("test1");

  parent = gt_feature_node_new(seqid1, gt_ft_gene, r1.start, r2.end,
                               GT_STRAND_FORWARD);
  gn1 = gt_feature_node_new(seqid1, gt_ft_exon, r1.start, r1.end,
                            GT_STRAND_FORWARD);
  gn2 = gt_feature_node_new(seqid1, gt_ft_exon, r2.start, r2.end,
                            GT_STRAND_FORWARD);
  gn3 = gt_feature_node_new(seqid1, gt_ft_TF_binding_site, r3.start, r3.end,
                            GT_STRAND_FORWARD);

  gt_feature_node_add_child((GtFeatureNode*) parent, (GtFeatureNode*) gn1);
  gt_feature_node_add_child((GtFeatureNode*) parent, (GtFeatureNode*) gn2);

  l1 = gt_line_new();

  gt_feature_node_add_attribute((GtFeatureNode*) parent, GT_GFF_NAME, foo);
  gt_feature_node_add_attribute((GtFeatureNode*) gn1, GT_GFF_NAME, bar);
  gt_feature_node_add_attribute((GtFeatureNode*) gn2, GT_GFF_NAME, bar);
  gt_feature_node_add_attribute((GtFeatureNode*) gn3, GT_GFF_NAME, blub);

  b1 = gt_block_new();
  b2 = gt_block_new();

  gt_block_insert_element(b1, (GtFeatureNode*) parent);
  gt_block_insert_element(b1, (GtFeatureNode*) gn1);
  gt_block_insert_element(b1, (GtFeatureNode*) gn2);
  gt_block_insert_element(b2, (GtFeatureNode*) gn3);
  gt_block_set_range(b1, r1);
  gt_block_set_range(b2, r2);

  /* test gt_line_insert_block */
  gt_ensure(had_err,  (0 == gt_array_size(gt_line_get_blocks(l1))));
  gt_line_insert_block(l1, b1);
  gt_ensure(had_err,  (1 == gt_array_size(gt_line_get_blocks(l1))));
  gt_line_insert_block(l1, b2);
  gt_ensure(had_err,  (2 == gt_array_size(gt_line_get_blocks(l1))));

  /* test gt_line_get_blocks */
  blocks = gt_line_get_blocks(l1);
  gt_ensure(had_err, (2 == gt_array_size(blocks)));

  /* test gt_line_get_height() */
  if (!had_err)
  {
    sty = gt_style_new(err);
    gt_ensure(had_err, sty && !gt_error_is_set(err));
    gt_ensure(had_err, gt_line_get_height(l1, &height, sty, err) == 0);
    gt_ensure(had_err, height == BAR_HEIGHT_DEFAULT);
    gt_ensure(had_err, !gt_error_is_set(err));
    gt_style_set_num(sty, "exon", "bar_height", 42);
    gt_ensure(had_err, gt_line_get_height(l1, &height, sty, err) == 0);
    gt_ensure(had_err, height == 42);
    gt_ensure(had_err, !gt_error_is_set(err));
    gt_style_set_num(sty, "gene", "bar_height", 23);
    gt_ensure(had_err, gt_line_get_height(l1, &height, sty, err) == 0);
    gt_ensure(had_err, height == 42);
    gt_ensure(had_err, !gt_error_is_set(err));
    gt_style_unset(sty, "exon", "bar_height");
    gt_ensure(had_err, gt_line_get_height(l1, &height, sty, err) == 0);
    gt_ensure(had_err, height == 23);
    gt_ensure(had_err, !gt_error_is_set(err));
    gt_style_unset(sty, "gene", "bar_height");
    gt_style_set_num(sty, "format", "bar_height", 99);
    gt_ensure(had_err, gt_line_get_height(l1, &height, sty, err) == 0);
    gt_ensure(had_err, height == 99);
    gt_ensure(had_err, !gt_error_is_set(err));
  }

  gt_str_delete(seqid1);
  gt_line_delete(l1);
  gt_style_delete(sty);
  gt_genome_node_delete(parent);
  gt_genome_node_delete(gn3);
  return had_err;
}

void gt_line_delete(GtLine *line)
{
  unsigned long i;
  if (!line) return;
  for (i = 0; i < gt_array_size(line->blocks); i++)
    gt_block_delete(*(GtBlock**) gt_array_get(line->blocks, i));
  gt_array_delete(line->blocks);
  gt_free(line);
}
