/*
   Copyright (c) 2007 Sascha Steinbiss, Christin Schaerfer, Malte Mader
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
  array_add(lines->blocks, block, env);

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
     r2 = (Block*) array_get(line->blocks, i);
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
    block_delete((Block*) array_get(line->blocks, i));
  }

  array_delete(line->blocks, env);
  env_ma_free(line, env);
}

