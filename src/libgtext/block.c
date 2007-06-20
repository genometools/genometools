/*
   Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>   
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtext/block.h>
#include <libgtext/element.h>

struct Block
{
  Dlist *elements;
  Range range;
  const char* caption;
};

/*!

*/
int elemcmp(const void *a, const void *b)
{
  Element *elem_a = (Element*) a;
  Element *elem_b = (Element*) b;

  Range ra = element_get_range(elem_a);
  Range rb = element_get_range(elem_b);

  return range_compare(ra, rb);
}

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
  block->elements = dlist_new(elemcmp, env);
  Range r;
  r.start = 0;
  r.end = 0;
  block->range = r;
  block->caption = NULL;

  assert(block);
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
  assert(block && gn);

  Element *element;

  element = element_new(gn, cfg, env);
  dlist_add(block->elements, element, env);
}

/*!
Returns range of a Block object
\param block Pointer to Block object
\return Pointer to Range object
*/
Range block_get_range(Block *block)
{  
   assert(block);

   return block->range;
}

/*!
Sets range of a Block object
\param block Pointer to Block object to set range
\param r Range to set
*/
void block_set_range(Block *block, Range r)
{
  assert(block && r.start && r.end);

  block->range = r;
}

/*!
Sets caption of a Block object
\param block Pointer to Block object to set caption
\param caption Pointer to String object
*/
void block_set_caption(Block *block,
                       const char* caption)
{
  assert(block);
  block->caption = caption;
}

/*!
Gets caption of a Block object
\param block Pointer to Block object 
\return caption Pointer to String object
*/
const char* block_get_caption(Block *block)
{
  assert(block);

  return block->caption;
}

/*!
Returns Array with Pointer to Element objects
\param block Pointer to Block object
\return Pointer to Array
*/
Dlist* block_get_elements(Block* block)
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
  Dlistelem *delem;

  if(!block) return;

  for(delem = dlist_first(block->elements); delem != NULL;
      delem = dlistelem_next(delem))
  {
    Element* elem = (Element*) dlistelem_get_data(delem);
    element_delete(elem, env);
  }
  dlist_delete(block->elements, env);
  env_ma_free(block, env);
}

/*!
Prints all Elements of a Block object
\param block Pointer to Block object
*/
void print_block(Block* block)
{
  assert(block);

  Dlistelem *elem;

  for(elem = dlist_first(block->elements); elem != NULL;
      elem = dlistelem_next(elem))
  {
    Element *element = (Element*) dlistelem_get_data(elem);
    printf("[");
    print_element(element);
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
  Dlist* elements;
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

  const char* caption1 = "foo";
  const char* caption2 = "bar";

  /* test block_insert_elements */
  ensure(has_err, (0 == dlist_size(block_get_elements(b))));
  block_insert_element(b, gn1, NULL, env);
  ensure(has_err, (1 == dlist_size(block_get_elements(b))));
  block_insert_element(b, gn2, NULL, env);
  ensure(has_err, (2 == dlist_size(block_get_elements(b))));

  /* test block_get_elements */
  elements = block_get_elements(b);
  Element *elem = (Element*) dlistelem_get_data(dlist_first(elements));
  ensure(has_err, elements_are_equal(e1, elem));
  ensure(has_err, !elements_are_equal(e2, (Element*) dlist_first(elements)));
  elem = (Element*) dlistelem_get_data(dlist_last(elements));
  ensure(has_err, !elements_are_equal(e1, elem));
  ensure(has_err, elements_are_equal(e2, elem));

  /* test block_set_range & block_get_range */
  r_temp = range_join(r1, r2);
  block_set_range(b, r_temp);
  b_range = block_get_range(b);
  ensure(has_err, (0 == range_compare(b_range, r_temp)));
  ensure(has_err, (1 == range_compare(r2, r_temp)));

  /* tests block_set_caption 
     & block_get_caption */
  block_set_caption(b, caption1);
  ensure(has_err, (0 == strcmp(block_get_caption(b), caption1)));
  ensure(has_err, (0 != strcmp(block_get_caption(b), caption2)));
	      
  element_delete(e1, env);
  element_delete(e2, env);
  block_delete(b, env);
  genome_node_delete(gn1, env);
  genome_node_delete(gn2, env);
	      
  return has_err;
}

