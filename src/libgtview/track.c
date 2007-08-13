/*
   Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/
/**
 * \if INTERNAL \file track.c \endif
 * \author Christin Schaerfer <cschaerfer@zbh.uni-hamburg.de>
 */

#include "libgtcore/ensure.h"
#include "libgtcore/hashtable.h"
#include "libgtview/track.h"

typedef struct
{
  void *ptr;
  bool collapse_target;
} CollapseCheckInfo;

struct Track
{
  Str *title;
  Array *lines;
};

Track* track_new(Str *title,
                 Env *env)
{
  assert(title && env);

  Track *track;
  env_error_check(env);
  track = env_ma_malloc(env, sizeof (Track));
  Str* t = title;
  track->title = t;
  track->lines = array_new(sizeof (Line*), env);
  assert(track);
  return track;
}

/*
Constructs Track title
*/
/*static int add_type_to_title(void *key, void *value, void* data, Env* env)
{
  Str *title = (Str*) data;
  const char *type = (const char*) key;
  str_append_cstr(title, type, env);
  str_append_cstr(title, ", ", env);
  return 0;
}*/

Str* track_get_title(Track *track)
{
  assert(track && track->title);
  return track->title;
}

Line* get_next_free_line(Track* track, Range r, Env* env)
{
  assert(track && env);

  int i;
  Line* line;

  for (i=0; i<array_size(track->lines); i++)
  {
    line = *(Line**) array_get(track->lines, i);
    if (!line_is_occupied(line, r))
    {
      return line;
    }
  }
  line = line_new(env);
  array_add(track->lines, line, env);

  assert(line);
  return line;
}

Array* track_get_lines(Track* track)
{
  return track->lines;
}

int track_get_number_of_lines(Track *track)
{
  assert(track);

  int nof_tracks;
  nof_tracks = array_size(track->lines);
  return nof_tracks;
}

/*
Sort block into free line
*/
void track_insert_block(Track *track, Block *block, Env *env)
{
  assert(track && block && env);

  Range r;
  Line *line;

  r = block_get_range(block);
  line = get_next_free_line(track, r, env);
  line_insert_block(line, block, env);
}

int track_unit_test(Env* env)
{
/*  Range r1, r2, r3, r_parent, r_mRNA1, r_mRNA2;
  Str *seqid1, *seqid2, *seqid3;
  Array* lines;*/
  int had_err = 0;
/*
  Config *cfg;

  cfg = config_new(env, false);

  r_parent.start = 10;
  r_parent.end = 80;

  r1.start = 10;
  r1.end = 50;

  r2.start = 60;
  r2.end = 80;

  r3.start = 70;
  r3.end = 100;

  r_mRNA1 = r_mRNA2 = r_parent;

  seqid1 = str_new_cstr("test1", env);
  seqid2 = str_new_cstr("test2", env);
  seqid3 = str_new_cstr("foo", env);

  GenomeNode* parent = genome_feature_new(gft_gene, r_parent,
                                          STRAND_FORWARD, NULL, 0, env);
  GenomeNode* gn1 = genome_feature_new(gft_exon, r1,
                                       STRAND_FORWARD, NULL, 0, env);
  GenomeNode* gn2 = genome_feature_new(gft_exon, r2,
                                       STRAND_FORWARD, NULL, 0, env);
  GenomeNode* gn3 = genome_feature_new(gft_intron, r3,
                                       STRAND_FORWARD, NULL, 0, env);
  GenomeNode* mRNA1 = genome_feature_new(gft_mRNA, r_mRNA1,
                                         STRAND_FORWARD, NULL, 0, env);
  GenomeNode* mRNA2 = genome_feature_new(gft_mRNA, r_mRNA2,
                                         STRAND_FORWARD, NULL, 0, env);

  genome_node_set_seqid((GenomeNode*) parent, seqid1);
  genome_node_set_seqid((GenomeNode*) gn1, seqid3);
  genome_node_set_seqid((GenomeNode*) gn2, seqid3);
  genome_node_set_seqid((GenomeNode*) gn3, seqid2);
  genome_node_set_seqid((GenomeNode*) mRNA1, seqid2);
  genome_node_set_seqid((GenomeNode*) mRNA2, seqid2);

  genome_node_is_part_of_genome_node(parent, mRNA1, env);
  genome_node_is_part_of_genome_node(parent, mRNA2, env);

  Line* l1 = line_new(env);
  Line* l2 = line_new(env);

  const char* foo = "foo";
  const char* bar = "bar";
  const char* blub = "blub";

  genome_feature_add_attribute((GenomeFeature*) parent, "Name", foo, env);
  genome_feature_add_attribute((GenomeFeature*) gn1, "Name", bar, env);
  genome_feature_add_attribute((GenomeFeature*) gn2, "Name", bar, env);
  genome_feature_add_attribute((GenomeFeature*) gn3, "Name", blub, env);
  genome_feature_add_attribute((GenomeFeature*) mRNA1, "Name", blub, env);
  genome_feature_add_attribute((GenomeFeature*) mRNA2, "Name", blub, env);

  Str* title = str_new(env);
  str_set(title, "exon", env);
  Str* s = str_new(env);
  str_set(s, "foo", env);

  Track* t = track_new(title, env);
  Track* t2 = track_new(s, env);
  Track* t3 = track_new(s, env);

  int nof_blocks;
  nof_blocks = track_get_number_of_blocks(t, env);
  ensure(had_err, (0 == nof_blocks));
  track_insert_element(t, gn1, cfg, mRNA1, env);
  nof_blocks = track_get_number_of_blocks(t, env);
  ensure(had_err, (1 == nof_blocks));
  track_insert_element(t, gn2, cfg, mRNA1, env);
  nof_blocks = track_get_number_of_blocks(t, env);
  ensure(had_err, (1 == nof_blocks));
  track_insert_element(t, gn3, cfg, gn1, env);
  nof_blocks = track_get_number_of_blocks(t, env);
  ensure(had_err, (2 == nof_blocks));
  nof_blocks = track_get_number_of_blocks(t3, env);
  ensure(had_err, (0 == nof_blocks));
  track_insert_element(t3, mRNA1, cfg, parent, env);
  nof_blocks = track_get_number_of_blocks(t3, env);
  ensure(had_err, (1 == nof_blocks));
  track_insert_element(t3, mRNA2, cfg, parent, env);
  nof_blocks = track_get_number_of_blocks(t3, env);
  ensure(had_err, (2 == nof_blocks));
  track_finish(t3, env);

  ensure(had_err, (0 == array_size(track_get_lines(t))));
  track_finish(t, env);
  ensure(had_err, (2 == array_size(track_get_lines(t))));

  ensure(had_err, (0 == str_cmp(title, track_get_title(t))));
  ensure(had_err, !(0 == str_cmp(s, track_get_title(t))));

  lines = track_get_lines(t2);
  ensure(had_err, (0 == array_size(lines)));
  ensure(had_err, (0 == track_get_number_of_lines(t2)));
  track_insert_element(t2, gn1, cfg, NULL, env);
  track_finish(t2, env);
  lines = track_get_lines(t2);
  ensure(had_err, (1 == array_size(lines)));
  ensure(had_err, (1 == track_get_number_of_lines(t2)));
  lines = track_get_lines(t);
  ensure(had_err, (2 == array_size(lines)));
  ensure(had_err, (2 == track_get_number_of_lines(t)));

  config_delete(cfg, env);
  line_delete(l1, env);
  line_delete(l2, env);
  track_delete(t, env);
  track_delete(t2, env);
  track_delete(t3, env);
  str_delete(title, env);
  str_delete(s, env);
  str_delete(seqid1, env);
  str_delete(seqid2, env);
  str_delete(seqid3, env);
  genome_node_delete(parent, env);
  genome_node_delete(gn1, env);
  genome_node_delete(gn2, env);
  genome_node_delete(gn3, env);
  genome_node_delete(mRNA1, env);
  genome_node_delete(mRNA2, env);
*/
  return had_err;
}

void track_delete(Track *track,
                  Env *env)
{
  int i;
  if (!track) return;
  for (i=0; i<array_size(track->lines); i++)
  {
    line_delete(*(Line**) array_get(track->lines, i), env);
  }
  array_delete(track->lines, env);
  str_delete(track->title, env);
  env_ma_free(track, env);
}

