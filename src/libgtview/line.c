/*
   Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/
/**
 * \if INTERNAL \file line.c \endif
 * \author Christin Schaerfer <cschaerfer@zbh.uni-hamburg.de>
 */

#include "libgtcore/ensure.h"
#include "libgtcore/range.h"
#include "libgtview/block.h"
#include "libgtview/line.h"

struct Line
{
  Array *blocks;
};

Line* line_new(Env *env)
{
  Line *line;
  env_error_check(env);
  line = env_ma_malloc(env, sizeof (Line));
  line->blocks = array_new(sizeof (Block*), env);

  assert(line != NULL);
  return line;
}

void line_insert_block(Line *line, Block *block, Env *env)
{
  assert(line != NULL && block != NULL);

  array_add(line->blocks, block, env);
}

bool line_is_occupied(Line *line, Range r)
{
  int i;
  Range r1;

  assert(line != NULL);

  for (i=0; i<array_size(line->blocks); i++)
  {
    r1 = block_get_range(*(Block**)  array_get(line->blocks, i));
    if (range_overlap(r1, r))
    {
      return true;
    }
  }
  return false;
}

Array* line_get_blocks(Line* line)
{
  assert(line != NULL);
  return line->blocks;
}

int line_unit_test(Env* env)
{
  Range r1, r2, r3, r4, r_parent;
  Array* blocks;
  Str *seqid1, *seqid2, *seqid3;
  int had_err = 0;
  Config *cfg;
  GenomeNode *parent, *gn1, *gn2, *gn3, *gn4;
  Line *l1, *l2;
  Block *b1, *b2;

  const char* foo = "foo";
  const char* bar = "bar";
  const char* blub = "blub";

  cfg = config_new(env, false);

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

  seqid1 = str_new_cstr("test1", env);
  seqid2 = str_new_cstr("test2", env);
  seqid3 = str_new_cstr("foo", env);

  parent = genome_feature_new(gft_gene, r_parent,
                                          STRAND_FORWARD, NULL, 0, env);
  gn1 = genome_feature_new(gft_exon, r1,
                                       STRAND_FORWARD, NULL, 0, env);
  gn2 = genome_feature_new(gft_exon, r2,
                                       STRAND_FORWARD, NULL, 0, env);
  gn3 = genome_feature_new(gft_exon, r3,
                                       STRAND_FORWARD, NULL, 0, env);
  gn4 = genome_feature_new(gft_TF_binding_site, r4,
                                       STRAND_FORWARD, NULL, 0, env);

  genome_node_set_seqid((GenomeNode*) parent, seqid1);
  genome_node_set_seqid((GenomeNode*) gn1, seqid3);
  genome_node_set_seqid((GenomeNode*) gn2, seqid3);
  genome_node_set_seqid((GenomeNode*) gn3, seqid2);
  genome_node_set_seqid((GenomeNode*) gn4, seqid3);

  l1 = line_new(env);
  l2 = line_new(env);

  genome_feature_add_attribute((GenomeFeature*) parent, "Name", foo, env);
  genome_feature_add_attribute((GenomeFeature*) gn1, "Name", bar, env);
  genome_feature_add_attribute((GenomeFeature*) gn2, "Name", bar, env);
  genome_feature_add_attribute((GenomeFeature*) gn3, "Name", blub, env);
  genome_feature_add_attribute((GenomeFeature*) gn4, "Name", bar, env);

  b1 = block_new(env);
  b2 = block_new(env);

  block_insert_element(b1, gn1, cfg, env);
  block_insert_element(b2, gn2, cfg, env);
  block_set_range(b1, r1);
  block_set_range(b2, r2);

  /* test line_insert_block */
  ensure(had_err,  (0 == array_size(line_get_blocks(l1))));
  line_insert_block(l1, b1, env);
  ensure(had_err,  (1 == array_size(line_get_blocks(l1))));
  line_insert_block(l1, b2, env);
  ensure(had_err,  (2 == array_size(line_get_blocks(l1))));

  /* test line_is_occupied */
  ensure(had_err, !line_is_occupied(l2, r3));
  ensure(had_err, line_is_occupied(l1, r3));

  /* test line_get_blocks */
  blocks = line_get_blocks(l1);
  ensure(had_err, (2 == array_size(blocks)));

  config_delete(cfg, env);
  str_delete(seqid1, env);
  str_delete(seqid2, env);
  str_delete(seqid3, env);
  line_delete(l1, env);
  line_delete(l2, env);
  genome_node_delete(parent, env);
  genome_node_delete(gn1, env);
  genome_node_delete(gn2, env);
  genome_node_delete(gn3, env);
  genome_node_delete(gn4, env);

  return had_err;
}

void line_delete(Line *line, Env *env)
{
  int i;

  if (!line) return;

  for (i=0; i<array_size(line->blocks); i++)
  {
    block_delete(*(Block**) array_get(line->blocks, i), env);
  }

  array_delete(line->blocks, env);
  env_ma_free(line, env);
}

