/*
   Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtext/track.h>

GenomeNode* last_parent = NULL;

struct Track
{
  Str *title;
  Hashtable *types;
  Array *lines;
  Hashtable *blocks;
};

Track* track_new(Str *title,
                 Env *env)
{
  assert(title);

  Track *track;
  env_error_check(env);
  track = env_ma_malloc(env, sizeof (Track));
  Str* t = str_ref(title);
  track->title = t;
  track->lines = array_new(sizeof (Line*), env);
  track->blocks = hashtable_new(HASH_DIRECT, NULL, NULL, env);
  track->types = hashtable_new(HASH_STRING, NULL, NULL, env);
  assert(track);
  return track;
}

void track_insert_element(Track *track,
                          GenomeNode *gn,
                          Config *cfg,
                          GenomeNode *parent,
                          Env *env)
{
  assert(track && gn && cfg);
  Block *block;
  const char* caption;
  const char* parent_caption;
  Range gn_range;
  GenomeFeatureType gn_type;
  Strand gn_strand;
  Strand parent_strand;
  const char* feature_type;

  gn_type = genome_feature_get_type((GenomeFeature*) gn);
  feature_type = genome_feature_type_get_cstr(gn_type);
  gn_range = genome_node_get_range(gn);
  gn_strand = genome_feature_get_strand((GenomeFeature*) gn);

  if (!hashtable_get(track->types, feature_type))
    hashtable_add(track->types,
                  (char *) feature_type,
                  (char *) feature_type,
                  env);

  if (parent == NULL)
  {
    block = block_new(env);
    hashtable_add(track->blocks, gn, block, env);
    caption = genome_feature_get_attribute(gn, "Name");
    if (caption == NULL)
    {
      caption = genome_feature_get_attribute(gn, "ID");
    }
    block_set_caption(block, caption);
    block_set_range(block, genome_node_get_range(gn));
    block_set_type(block, gn_type);
    block_set_strand(block, gn_strand);
  }
  else
  {
    parent_strand = genome_feature_get_strand((GenomeFeature*) parent);
    block = (Block*) hashtable_get(track->blocks, parent);
    if (block == NULL)
    {
      block = block_new(env);
      if (gn_type == gft_mRNA)
      {
        hashtable_add(track->blocks, gn, block, env);
      }
      else
      {
        hashtable_add(track->blocks, parent, block, env);
      }
      caption = genome_feature_get_attribute(gn, "Name");
      if (caption == NULL)
      {
         caption = genome_feature_get_attribute(gn, "ID");
      }
      block_set_caption(block, caption);
      parent_caption =  genome_feature_get_attribute(parent, "Name");
      if (parent_caption == NULL)
      {
        parent_caption = genome_feature_get_attribute(parent, "ID");
      }
      block_set_parent_caption(block, parent_caption);
      block_set_range(block, genome_node_get_range(parent));
      block_set_type(block, gn_type);
      block_set_strand(block, parent_strand);
    }
  }
  block_insert_element(block, gn, cfg, env);
}

/*
Constructs Track title
*/
int add_type_to_title(void *key, void *value, void* data, Env* env)
{
  Str *title = (Str*) data;
  const char *type = (const char*) key;
  str_append_cstr(title, type, env);
  str_append_cstr(title, ", ", env);
  return 0;
}

Str* track_get_title(Track *track)
{
  assert(track && track->title);
  return track->title;
}

Line* get_next_free_line(Track* track, Range r, Env* env)
{
  assert(track);

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
foreach function to track_get_number_of_blocks
*/
int track_add_block(void *parent , void *block, void *nof_blocks, Env *env)
{
  int *add;

  add = nof_blocks;
  *add += 1;

  return 0;
}

/*
Returns number of Blocks of a Track object
*/
int track_get_number_of_blocks(Track *track, Env *env)
{
  assert(track);
  int nof_blocks;
  nof_blocks = 0;
  hashtable_foreach(track->blocks, track_add_block, &nof_blocks, env);

  return nof_blocks;
}

/*
Sort block into free line
*/
int process_block(void *parent, void *block, void *track, Env *env)
{
  Range r;
  Line *line;
  Block *cur_block = (Block*) block;
  Track *t = (Track*) track;

  r = block_get_range(cur_block);
  line = get_next_free_line(t, r, env);
  line_insert_block(line, cur_block, env);

  return 0;
}

void track_finish(Track *track, Env *env)
{
  hashtable_foreach(track->blocks, process_block, track, env);
  hashtable_foreach(track->types, add_type_to_title, track->title, env);
  str_set_length(track->title, str_length(track->title)-2);
  hashtable_delete(track->blocks, env);
}

int track_unit_test(Env* env)
{
  Range r1, r2, r3, r_parent, r_mRNA1, r_mRNA2;
  Str *seqid1, *seqid2, *seqid3;
  Array* lines;
  int has_err = 0;

  Config *cfg;
  Str *luafile = str_new_cstr("config.lua",env);

  cfg = config_new(env, false);
  config_load_file(cfg, luafile, env);

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

  /* test track_insert_elements
     (implicit test of get_next_free_line) */
  int nof_blocks;
  nof_blocks = track_get_number_of_blocks(t, env);
  ensure(has_err, (0 == nof_blocks));
  track_insert_element(t, gn1, cfg, mRNA1, env);
  nof_blocks = track_get_number_of_blocks(t, env);
  ensure(has_err, (1 == nof_blocks));
  track_insert_element(t, gn2, cfg, mRNA1, env);
  nof_blocks = track_get_number_of_blocks(t, env);
  ensure(has_err, (1 == nof_blocks));
  track_insert_element(t, gn3, cfg, gn1, env);
  nof_blocks = track_get_number_of_blocks(t, env);
  ensure(has_err, (2 == nof_blocks));
  nof_blocks = track_get_number_of_blocks(t3, env);
  ensure(has_err, (0 == nof_blocks));
  track_insert_element(t3, mRNA1, cfg, parent, env);
  nof_blocks = track_get_number_of_blocks(t3, env);
  ensure(has_err, (1 == nof_blocks));
  track_insert_element(t3, mRNA2, cfg, parent, env);
  nof_blocks = track_get_number_of_blocks(t3, env);
  ensure(has_err, (2 == nof_blocks));
  track_finish(t3, env);

  /* test track_finish*/
  ensure(has_err, (0 == array_size(track_get_lines(t))));
  track_finish(t, env);
  ensure(has_err, (2 == array_size(track_get_lines(t))));

  /* test track_get_title */
  ensure(has_err, (0 == str_cmp(title, track_get_title(t))));
  ensure(has_err, !(0 == str_cmp(s, track_get_title(t))));

  /* test track_get_lines and track_get_number_of_lines */
  lines = track_get_lines(t2);
  ensure(has_err, (0 == array_size(lines)));
  ensure(has_err, (0 == track_get_number_of_lines(t2)));
  track_insert_element(t2, gn1, cfg, NULL, env);
  track_finish(t2, env);
  lines = track_get_lines(t2);
  ensure(has_err, (1 == array_size(lines)));
  ensure(has_err, (1 == track_get_number_of_lines(t2)));
  lines = track_get_lines(t);
  ensure(has_err, (2 == array_size(lines)));
  ensure(has_err, (2 == track_get_number_of_lines(t)));

  config_delete(cfg, env);
  str_delete(luafile, env);
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

  return has_err;
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
  hashtable_delete(track->types, env);
  array_delete(track->lines, env);
  str_delete(track->title, env);
  env_ma_free(track, env);
}

