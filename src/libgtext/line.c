/*
   Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>   
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtcore/range.h>
#include <libgtext/block.h>
#include <libgtext/line.h>

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
			 Env *env)
{
  Block *block;

  block = block_new(env);
  array_add(line->blocks, block, env);

  block_insert_element(block, gn, cfg, env);
}

/*!
Checks if Line is occupied
\param line Line object to check
\param gn Pointer to GenomeNode object
\return True or False
*/
bool line_is_occupied(Line *line, GenomeNode *gn)
{
   int i;
   Range r1, r2;

   r1 = genome_node_get_range(gn);

   for(i=0; i<array_size(line->blocks); i++)
   {
     r2 = block_get_range(*(Block**)  array_get(line->blocks, i));
     if(!range_overlap(r1, r2))
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
  int i;

  for(i=0; i<array_size(line->blocks); i++)
  {
    printf("(");
    print_block(*(Block**) array_get(line->blocks, i));
    printf(")");
  }
}

/*!
Tests Block Class
\param env Pointer to Environment object
*/
int line_unit_test(Env* env)
{
  Range r1, r2, r3; 
  Array* blocks;
  int has_err = 0;

  r1.start = 10;
  r1.end = 50;

  r2.start = 60;
  r2.end = 80;

  r3.start = 70;
  r3.end = 100;

  GenomeNode* gn1 = genome_feature_new(gft_exon, r1, STRAND_FORWARD, NULL, 0, env);
  GenomeNode* gn2 = genome_feature_new(gft_intron, r2, STRAND_FORWARD, NULL, 0, env);
  GenomeNode* gn3 = genome_feature_new(gft_intron, r3, STRAND_FORWARD, NULL, 0, env);
		    
  Line* l1 = line_new(env);
  Line* l2 = line_new(env);
			    
  /* test line_insert_elements */
  ensure(has_err, (0 == array_size(line_get_blocks(l1))));
  line_insert_element(l1, gn1, NULL, env);
  ensure(has_err, (1 == array_size(line_get_blocks(l1))));
  line_insert_element(l1, gn2, NULL, env);
  ensure(has_err, (2 == array_size(line_get_blocks(l1))));

  /* test line_is_occupied */
  ensure(has_err, !line_is_occupied(l2, gn3));
  ensure(has_err, line_is_occupied(l1, gn3));

  /* test line_get_blocks */
  blocks = line_get_blocks(l1);
  ensure(has_err, (2 == array_size(blocks)));

  line_delete(l1, env);
  line_delete(l2, env);
  genome_node_delete(gn1, env);
  genome_node_delete(gn2, env);
  genome_node_delete(gn3, env);

  return has_err;
}


