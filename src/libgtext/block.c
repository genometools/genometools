/*
   Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>   
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtext/block.h>
#include <libgtext/element.h>

struct Block
{
  Array *elements;
  Range range;
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
  Range r;
  r.start = 0;
  r.end = 0;
  block->range = r;
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
  Range r;

  element = element_new(gn, cfg, env);
  array_add(block->elements, element, env);
  r = element_get_range(element);
  if(1 == array_size(block->elements))
  {
    block->range = element_get_range(element);
  }
  else 
  {
    block->range = range_join(r, block->range);
  }

}

/*!
Returns range of a Block object
\param block Pointer to Block object
\return Pointer to Range object
*/
Range block_get_range(Block *block)
{  
   return block->range;
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

/*!
Prints all Elements of a Block object
\param block Pointer to Block object
*/
void print_block(Block* block)
{
  int i;

  for(i=0; i<array_size(block->elements); i++)
  {
    printf("[");
    print_element(*(Element**) array_get(block->elements, i));
    printf("]");
  }
}

/*!
 * Unit Test for Block Class
 * \param env Pointer to Environment object
 * */
int block_unit_test(Env* env)
{
  Range r1, r2, r_temp, b_range;
  Array* elements;
  int has_err = 0;

  r1.start = 10;
  r1.end = 50;

  r2.start = 50;
  r2.end = 80;


  GenomeNode* gn1 = genome_feature_new(gft_exon, r1, STRAND_FORWARD, NULL, 0, env);
  GenomeNode* gn2 = genome_feature_new(gft_intron, r2, STRAND_FORWARD, NULL, 0, env);

  Element* e1 = element_new(gn1, NULL, env);
  Element* e2 = element_new(gn2, NULL, env);

  Block* b = block_new(env);

  /* test block_insert_elements */
  ensure(has_err, (0 == array_size(block_get_elements(b))));
  block_insert_element(b, gn1, NULL, env);
  ensure(has_err, (1 == array_size(block_get_elements(b))));
  block_insert_element(b, gn2, NULL, env);
  ensure(has_err, (2 == array_size(block_get_elements(b))));

  /* test block_get_elements */
  elements = block_get_elements(b);
  ensure(has_err, elements_are_equal(e1, *(Element**) array_get(elements, 0)));
  ensure(has_err, !elements_are_equal(e1, *(Element**) array_get(elements, 1)));
  ensure(has_err, !elements_are_equal(e2, *(Element**) array_get(elements, 0)));
  ensure(has_err, elements_are_equal(e2, *(Element**) array_get(elements, 1)));

  /* test block_get_range */ 
  b_range = block_get_range(b);
  r_temp = range_join(r1, r2);
  ensure(has_err, (0 == range_compare(b_range, r_temp)));
  ensure(has_err, (1 == range_compare(r2, r_temp)));
	      
  element_delete(e1, env);
  element_delete(e2, env);
  block_delete(b, env);
  genome_node_delete(gn1, env);
  genome_node_delete(gn2, env);
	      
  return has_err;
}

