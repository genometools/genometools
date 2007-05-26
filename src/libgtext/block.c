/*
   Copyright (c) 2007 Sascha Steinbiss, Christin Schaerfer, Malte Mader
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtext/block.h>
#include <libgtext/element.h>

struct Block
{
  Array *elements;
};

/*!
Creates a new Block object.
\param env Pointer to Environment object.
\return Pointer to a new Block object.
*/
Block* block_new(Env *env)
{
  Block *block;
  env_error_check(env);
  block = env_ma_malloc(env, sizeof (Block));
  block->elements = array_new(sizeof (Element*), env);
  return block;
}

/*!
Inserts an element into a Block object
\param block Block in which the element shall be insert
\param gn Pointer to GenomeNode object
\param cfg Pointer to Config file
\param env Pointer to Environment object
*/
void block_insert_element(Block *block,
                          GenomeNode *gn,
			  Config *cfg,
			  Env *env)
{
  Element *element;

  element = element_new(gn, cfg, env);
  array_add(block->elements, element, env);
}

/*!
Returns range of a Block object
\param block Pointer to Block object
\return Pointer to Range object
*/
Range block_get_range(Block *block)
{  

   Range r1, r2;
   int i;

   r1 = element_get_range(*(Element**) array_get(block->elements, 0));

   for(i=1; i<array_size(block->elements); i++)
   {
     r2 = element_get_range(*(Element**) array_get(block->elements, i));
     r1 = range_join(r1, r2);
   }

   return r1;
}

/*!
Returns Array with Pointer to Element objects
\param block Pointer to Block object
\return Pointer to Array
*/
Array* block_get_elements(Block* block)
{
  return block->elements;
}

/*!
Delets Block
\param block Pointer to Block object to delete
\param env Pointer to Environment object
*/
void block_delete(Block *block,
                  Env *env)
{
  int i;

  if(!block) return;

  for(i=0; i<array_size(block->elements); i++)
  {
    element_delete(*(Element**) array_get(block->elements, i), env);
  }

  array_delete(block->elements, env);
  env_ma_free(block, env);
}

