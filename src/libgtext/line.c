/*
   Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>   
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtcore/range.h>
#include <libgtext/block.h>
#include <libgtext/line.h>

GenomeNode* last_parent = NULL;

struct Line
{
  Array *blocks;
};

/*!
Creates a new Line object.
\param env Pointer to Environment object.
\return Pointer to a new Line object.
*/
Line* line_new(Env *env)
{
  Line *line;
  env_error_check(env);
  line = env_ma_malloc(env, sizeof (Line));
  line->blocks = array_new(sizeof (Block*), env);

  assert(line);
  return line;
}

/*!
Inserts an element into a Line object
\param line Line in which the element shall be insert
\param gn Pointer to GenomeNode object
\param cfg Pointer to Config file
\param env Pointer to Environment object
*/
void line_insert_element(Line *line,
                         GenomeNode *gn,
			 Config *cfg,
			 GenomeNode *parent,
			 Env *env)
{
  assert(line && gn && cfg);
  Block *block;
  Str* caption;

  if((last_parent != NULL)
     && (parent != NULL)
     && (0 == genome_node_compare(&parent, &last_parent)))
  {
    block = *(Block**) array_get(line->blocks, (array_size(line->blocks) -1));
    if(!range_overlap(genome_node_get_range(gn), block_get_range(block)))
    {
      caption = genome_node_get_idstr(parent);
    }
    else
    {
      block = block_new(env);
      array_add(line->blocks, block, env);
      caption = genome_node_get_idstr(gn);
    }
  }
  else
  {
    block = block_new(env);
    array_add(line->blocks, block, env);
    caption = genome_node_get_idstr(gn);
  }

  block_insert_element(block, gn, cfg, env);
  block_set_caption(block, caption);
  last_parent = parent;
}

/*!
Checks if Line is occupied
\param line Line object to check
\param gn Pointer to GenomeNode object
\return True or False
*/
bool line_is_occupied(Line *line, GenomeNode *gn)
{
   assert(line && gn);

   int i;
   Range r1, r2;

   r2 = genome_node_get_range(gn);

   for(i=0; i<array_size(line->blocks); i++)
   {
     r1 = block_get_range(*(Block**)  array_get(line->blocks, i));
     if(range_overlap(r1, r2))
     {
       return true;
     }
   }
   return false;
}

/*!
Returns Array with Pointer to Block objects
\param line Pointer to Line object
\return Pointer to Array
*/
Array* line_get_blocks(Line* line)
{
  return line->blocks;
}


/*!
Delets Line
\param line Pointer to Line object to delete
\param env Pointer to Environment object
*/
void line_delete(Line *line,
                 Env *env)
{
  int i;

  if(!line) return;

  for(i=0; i<array_size(line->blocks); i++)
  {
    block_delete(*(Block**) array_get(line->blocks, i), env);
  }

  array_delete(line->blocks, env);
  env_ma_free(line, env);
}

/*!
Prints all Blocks of a Line object
\param line Pointer to Line object
*/
void print_line(Line* line)
{
  assert(line);
  int i;

  for(i=0; i<array_size(line->blocks); i++)
  {
    printf("(");
    print_block(*(Block**) array_get(line->blocks, i));
    printf(")");
  }
}

/*!
Tests Line Class
\param env Pointer to Environment object
*/
int line_unit_test(Env* env)
{
  Range r1, r2, r3, r_parent;
  Array* blocks;
  Str *seqid1, *seqid2, *seqid3;
  int has_err = 0;
  Config *cfg;
  Str *luafile = str_new_cstr("config.lua",env);

  /* do not show warnings during the unit test */
  bool verbose = false;

  cfg = config_new(env, &verbose);
  config_load_file(cfg, luafile, env);

  r_parent.start = 10;
  r_parent.end = 80;
  
  r1.start = 10;
  r1.end = 50;

  r2.start = 51;
  r2.end = 80;

  r3.start = 70;
  r3.end = 100;

  seqid1 = str_new_cstr("test1", env);
  seqid2 = str_new_cstr("test2", env);
  seqid3 = str_new_cstr("foo", env);

  GenomeNode* parent = genome_feature_new(gft_gene, r_parent, STRAND_FORWARD, NULL, 0, env);
  GenomeNode* gn1 = genome_feature_new(gft_exon, r1, STRAND_FORWARD, NULL, 0, env);
  GenomeNode* gn2 = genome_feature_new(gft_exon, r2, STRAND_FORWARD, NULL, 0, env);
  GenomeNode* gn3 = genome_feature_new(gft_exon, r3, STRAND_FORWARD, NULL, 0, env);

  genome_node_set_seqid((GenomeNode*) parent, seqid1);
  genome_node_set_seqid((GenomeNode*) gn1, seqid3);
  genome_node_set_seqid((GenomeNode*) gn2, seqid3);
  genome_node_set_seqid((GenomeNode*) gn3, seqid2);

  Line* l1 = line_new(env);
  Line* l2 = line_new(env);
			    
  /* test line_insert_elements */
  ensure(has_err, (0 == array_size(line_get_blocks(l1))));
  line_insert_element(l1, gn1, cfg, parent, env);
  ensure(has_err, (1 == array_size(line_get_blocks(l1))));
  line_insert_element(l1, gn2, cfg, parent, env);
  blocks = line_get_blocks(l1);
  ensure(has_err, (1 == array_size(blocks)));
  Block* b = *(Block**) array_get(blocks, 0);
  ensure(has_err, (0 == str_cmp(block_get_caption(b), genome_node_get_idstr(parent))));
  line_insert_element(l1, gn3, cfg, gn1, env);
  blocks = line_get_blocks(l1);
  ensure(has_err, (2 == array_size(blocks)));
  b = *(Block**) array_get(blocks, 1);
  ensure(has_err, (0 == str_cmp(block_get_caption(b), genome_node_get_idstr(gn3))));

  /* test line_is_occupied */
  ensure(has_err, !line_is_occupied(l2, gn3));
  ensure(has_err, line_is_occupied(l1, gn3));

  /* test line_get_blocks */
  blocks = line_get_blocks(l1);
  ensure(has_err, (2 == array_size(blocks)));

  config_delete(cfg, env);
  str_delete(luafile, env);
  str_delete(seqid1, env);
  str_delete(seqid2, env);
  str_delete(seqid3, env);
  line_delete(l1, env);
  line_delete(l2, env);
  genome_node_delete(parent, env);
  genome_node_delete(gn1, env);
  genome_node_delete(gn2, env);
  genome_node_delete(gn3, env);

  return has_err;
}


