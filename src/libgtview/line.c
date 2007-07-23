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

  assert(line);
  return line;
}

void line_insert_block(Line *line, Block *block, Env *env)
{
  assert(line && block);

  array_add(line->blocks, block, env);
}

bool line_is_occupied(Line *line, Range r)
{
   assert(line);

   int i;
   Range r1;

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
  return line->blocks;
}

int line_unit_test(Env* env)
{
  Range r1, r2, r3, r4, r_parent;
  Array* blocks;
  Str *seqid1, *seqid2, *seqid3;
  int had_err = 0;
  Config *cfg;

  cfg = config_new(env, false);

  r_parent.start = 10;
  r_parent.end = 80;

  r1.start = 10;
  r1.end = 50;

  r2.start = 51;
  r2.end = 80;

  r3.start = 70;
  r3.end = 100;

  r4.start = 10;
  r4.end = 20;

  seqid1 = str_new_cstr("test1", env);
  seqid2 = str_new_cstr("test2", env);
  seqid3 = str_new_cstr("foo", env);

  GenomeNode* parent = genome_feature_new(gft_gene, r_parent,
                                          STRAND_FORWARD, NULL, 0, env);
  GenomeNode* gn1 = genome_feature_new(gft_exon, r1,
                                       STRAND_FORWARD, NULL, 0, env);
  GenomeNode* gn2 = genome_feature_new(gft_exon, r2,
                                       STRAND_FORWARD, NULL, 0, env);
  GenomeNode* gn3 = genome_feature_new(gft_exon, r3,
                                       STRAND_FORWARD, NULL, 0, env);
  GenomeNode* gn4 = genome_feature_new(gft_TF_binding_site, r4,
                                       STRAND_FORWARD, NULL, 0, env);

  genome_node_set_seqid((GenomeNode*) parent, seqid1);
  genome_node_set_seqid((GenomeNode*) gn1, seqid3);
  genome_node_set_seqid((GenomeNode*) gn2, seqid3);
  genome_node_set_seqid((GenomeNode*) gn3, seqid2);
  genome_node_set_seqid((GenomeNode*) gn4, seqid3);

  Line* l1 = line_new(env);
  Line* l2 = line_new(env);

  const char* foo = "foo";
  const char* bar = "bar";
  const char* blub = "blub";

  genome_feature_add_attribute((GenomeFeature*) parent, "Name", foo, env);
  genome_feature_add_attribute((GenomeFeature*) gn1, "Name", bar, env);
  genome_feature_add_attribute((GenomeFeature*) gn2, "Name", bar, env);
  genome_feature_add_attribute((GenomeFeature*) gn3, "Name", blub, env);
  genome_feature_add_attribute((GenomeFeature*) gn4, "Name", bar, env);

  Block* b1 = block_new(env);
  Block* b2 = block_new(env);

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

