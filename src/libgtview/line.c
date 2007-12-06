/*
  Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/ensure.h"
#include "libgtcore/ma.h"
#include "libgtcore/range.h"
#include "libgtview/block.h"
#include "libgtview/line.h"

struct Line {
  Array *blocks;
};

Line* line_new(void)
{
  Line *line;
  line = ma_malloc(sizeof (Line));
  line->blocks = array_new(sizeof (Block*));
  return line;
}

void line_insert_block(Line *line, Block *block)
{
  assert(line && block);
  array_add(line->blocks, block);
}

bool line_is_occupied(const Line *line, Range r)
{
  unsigned long i;
  Range r1;
  assert(line);
  for (i = 0; i < array_size(line->blocks); i++) {
    r1 = block_get_range(*(Block**) array_get(line->blocks, i));
    if (range_overlap(r1, r))
      return true;
  }
  return false;
}

Array* line_get_blocks(Line* line)
{
  assert(line);
  return line->blocks;
}

int line_unit_test(Error *err)
{
  Range r1, r2, r3, r4, r_parent;
  Array* blocks;
  Str *seqid1, *seqid2, *seqid3;
  int had_err = 0;
  Config *cfg;
  GenomeNode *parent, *gn1, *gn2, *gn3, *gn4;
  Line *l1, *l2;
  Block *b1, *b2;
  error_check(err);

  const char* foo = "foo";
  const char* bar = "bar";
  const char* blub = "blub";

  if (!(cfg = config_new(false, err)))
    had_err = -1;

  r_parent.start = 10UL;
  r_parent.end = 80UL;

  r1.start = 10UL;
  r1.end = 50UL;

  r2.start = 51UL;
  r2.end = 80UL;

  r3.start = 70UL;
  r3.end = 100UL;

  r4.start = 10UL;
  r4.end = 20UL;

  seqid1 = str_new_cstr("test1");
  seqid2 = str_new_cstr("test2");
  seqid3 = str_new_cstr("foo");

  parent = genome_feature_new(gft_gene, r_parent, STRAND_FORWARD, NULL, 0);
  gn1 = genome_feature_new(gft_exon, r1, STRAND_FORWARD, NULL, 0);
  gn2 = genome_feature_new(gft_exon, r2, STRAND_FORWARD, NULL, 0);
  gn3 = genome_feature_new(gft_exon, r3, STRAND_FORWARD, NULL, 0);
  gn4 = genome_feature_new(gft_TF_binding_site, r4, STRAND_FORWARD, NULL, 0);

  genome_node_set_seqid((GenomeNode*) parent, seqid1);
  genome_node_set_seqid((GenomeNode*) gn1, seqid3);
  genome_node_set_seqid((GenomeNode*) gn2, seqid3);
  genome_node_set_seqid((GenomeNode*) gn3, seqid2);
  genome_node_set_seqid((GenomeNode*) gn4, seqid3);

  l1 = line_new();
  l2 = line_new();

  genome_feature_add_attribute((GenomeFeature*) parent, "Name", foo);
  genome_feature_add_attribute((GenomeFeature*) gn1, "Name", bar);
  genome_feature_add_attribute((GenomeFeature*) gn2, "Name", bar);
  genome_feature_add_attribute((GenomeFeature*) gn3, "Name", blub);
  genome_feature_add_attribute((GenomeFeature*) gn4, "Name", bar);

  b1 = block_new();
  b2 = block_new();

  block_insert_element(b1, gn1, cfg);
  block_insert_element(b2, gn2, cfg);
  block_set_range(b1, r1);
  block_set_range(b2, r2);

  /* test line_insert_block */
  ensure(had_err,  (0 == array_size(line_get_blocks(l1))));
  line_insert_block(l1, b1);
  ensure(had_err,  (1 == array_size(line_get_blocks(l1))));
  line_insert_block(l1, b2);
  ensure(had_err,  (2 == array_size(line_get_blocks(l1))));

  /* test line_is_occupied */
  ensure(had_err, !line_is_occupied(l2, r3));
  ensure(had_err, line_is_occupied(l1, r3));

  /* test line_get_blocks */
  blocks = line_get_blocks(l1);
  ensure(had_err, (2 == array_size(blocks)));

  config_delete(cfg);
  str_delete(seqid1);
  str_delete(seqid2);
  str_delete(seqid3);
  line_delete(l1);
  line_delete(l2);
  genome_node_delete(parent);
  genome_node_delete(gn1);
  genome_node_delete(gn2);
  genome_node_delete(gn3);
  genome_node_delete(gn4);

  return had_err;
}

void line_delete(Line *line)
{
  unsigned long i;
  if (!line) return;
  for (i = 0; i < array_size(line->blocks); i++)
    block_delete(*(Block**) array_get(line->blocks, i));
  array_delete(line->blocks);
  ma_free(line);
}
