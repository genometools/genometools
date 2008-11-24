/*
  Copyright (c) 2007      Christin Schaerfer <cschaerfer@zbh.uni-hamburg.de>
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

#include "core/ensure.h"
#include "core/interval_tree.h"
#include "core/ma.h"
#include "core/range.h"
#include "annotationsketch/block.h"
#include "annotationsketch/line.h"
#include "annotationsketch/style.h"

struct GtLine {
  bool has_captions;
  GtArray *blocks;
};

GtLine* gt_line_new(void)
{
  GtLine *line;
  line = gt_malloc(sizeof (GtLine));
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
    for (i = 0; i < gt_array_size(line->blocks); i++) {
      GtBlock *block;
      block = *(GtBlock**) gt_array_get(line->blocks, i);
      had_err = gt_block_sketch(block, canvas, err);
    }
  }
  if (!had_err)
    gt_canvas_visit_line_post(canvas, line, err);
  return had_err;
}

int gt_line_unit_test(GtError *err)
{
  GtRange r1, r2;
  GtArray* blocks;
  GtStr *seqid1, *seqid2, *seqid3;
  int had_err = 0;
  GtGenomeNode *parent, *gn1, *gn2, *gn3, *gn4;
  GtLine *l1, *l2;
  GtBlock *b1, *b2;
  gt_error_check(err);

  const char* foo = "foo";
  const char* bar = "bar";
  const char* blub = "blub";

  r1.start = 10;
  r1.end = 50;

  r2.start = 51;
  r2.end = 80;

  seqid1 = gt_str_new_cstr("test1");
  seqid2 = gt_str_new_cstr("test2");
  seqid3 = gt_str_new_cstr("foo");

  parent = gt_feature_node_new(seqid1, gft_gene, 10, 80, GT_STRAND_FORWARD);
  gn1 = gt_feature_node_new(seqid3, gft_exon, r1.start, r1.end,
                              GT_STRAND_FORWARD);
  gn2 = gt_feature_node_new(seqid3, gft_exon, r2.start, r2.end,
                             GT_STRAND_FORWARD);
  gn3 = gt_feature_node_new(seqid2, gft_exon, 70, 100, GT_STRAND_FORWARD);
  gn4 = gt_feature_node_new(seqid3, gft_TF_binding_site, 10, 20,
                              GT_STRAND_FORWARD);

  l1 = gt_line_new();
  l2 = gt_line_new();

  gt_feature_node_add_attribute((GtFeatureNode*) parent, "Name", foo);
  gt_feature_node_add_attribute((GtFeatureNode*) gn1, "Name", bar);
  gt_feature_node_add_attribute((GtFeatureNode*) gn2, "Name", bar);
  gt_feature_node_add_attribute((GtFeatureNode*) gn3, "Name", blub);
  gt_feature_node_add_attribute((GtFeatureNode*) gn4, "Name", bar);

  b1 = gt_block_new();
  b2 = gt_block_new();

  gt_block_insert_element(b1, (GtFeatureNode*) gn1);
  gt_block_insert_element(b2, (GtFeatureNode*) gn2);
  gt_block_set_range(b1, r1);
  gt_block_set_range(b2, r2);

  /* test gt_line_insert_block */
  ensure(had_err,  (0 == gt_array_size(gt_line_get_blocks(l1))));
  gt_line_insert_block(l1, b1);
  ensure(had_err,  (1 == gt_array_size(gt_line_get_blocks(l1))));
  gt_line_insert_block(l1, b2);
  ensure(had_err,  (2 == gt_array_size(gt_line_get_blocks(l1))));

  /* test gt_line_get_blocks */
  blocks = gt_line_get_blocks(l1);
  ensure(had_err, (2 == gt_array_size(blocks)));

  gt_str_delete(seqid1);
  gt_str_delete(seqid2);
  gt_str_delete(seqid3);
  gt_line_delete(l1);
  gt_line_delete(l2);
  gt_genome_node_delete(parent);
  gt_genome_node_delete(gn1);
  gt_genome_node_delete(gn2);
  gt_genome_node_delete(gn3);
  gt_genome_node_delete(gn4);

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
