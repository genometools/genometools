/*
   Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>   
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtcore/array.h>
#include <libgtext/genome_feature.h>
#include <libgtext/genome_feature_type.h>
#include <libgtcore/ensure.h>
#include <libgtcore/strand.h>
#include <libgtext/element.h>

enum ArrowStatus{
  Left = 1,
  Right = 2,
  NoArrow = 3,
};

struct Element
{
  GenomeFeatureType type;
  Range range;
  int arrow_status;
  Config* cfg;
};

/*!
Creates a new Element object.
\param gn Pointer to GenomeNode object.
\param cfg Pointer to Config object.
\param env Pointer to Environment object.
\return Pointer to a new Element object.
*/
Element* element_new(GenomeNode *gn, Config *cfg, Env *env)
{ 
  assert(gn);

  Element *element;
  GenomeFeature *gf = (GenomeFeature*) gn;

  env_error_check(env);
  element = env_ma_malloc(env, sizeof (Element));
  element->type = genome_feature_get_type(gf);
  element->range = genome_node_get_range(gn);
  element->arrow_status = NoArrow;
  element->cfg = cfg;

  assert(element);
  return element;
}

/*!
Sets ArrowStatus of an Element object
\param element Element on which status to set.
\param status ArrowStatus.
*/
void element_set_arrow_status(Element *element, int status)
{
  assert(element && status);
  element->arrow_status = status;
}

/*!
Returns ArrowStatus of an Element object
\param element Pointer to Element object
\return ArrowStatus
*/
int element_get_arrow_status(Element *element)
{
  assert(element);

  assert(element->arrow_status);
  return element->arrow_status;
}

/*!
Returns range of an Element object
\param block Pointer to Element object
\return Pointer to Range object
*/
Range element_get_range(Element *element)
{
  assert(element);

  assert(element->range.start && element->range.end);
  return element->range;
}

/*!
Delets Element
\param element Pointer to Element object to delete
\param env Pointer to Environment object
*/
void element_delete(Element *element,
                    Env *env)
{
  if(!element) return;

  env_ma_free(element, env);
}

/*!
Prints Element type and range
\param element Pointer to Element object to print
*/
void print_element(Element* element)
{
  assert(element);
  printf("%s, %lu - %lu", genome_feature_type_get_cstr(element->type),
                            element->range.start, element->range.end);
}

/*!
Checks if two Element objects are equal
\param e1 Pointer to Element object
\param e2 Pointer to Element object
\returns true if e1 and e2 are equal
*/
bool elements_are_equal(Element* e1,
                        Element* e2)
{
  assert(e1 && e2);
  if((0 == strcmp(genome_feature_type_get_cstr(e1->type), genome_feature_type_get_cstr(e2->type)))
     && (0 == range_compare(e1->range, e2->range))
     && (e1->arrow_status == e2->arrow_status))
    return true;
  else 
    return false;
}

/*!
Unit Test for Element Class
\param env Pointer to Environment object
*/
int element_unit_test(Env* env)
{
  Range r1, r2, r_temp;
  int as1, as2, as_temp;
  int has_err = 0;

  r1.start = 10;
  r1.end = 50;

  r2.start = 20;
  r2.end = 50;

  as1 = NoArrow;
  as2 = Right;

  GenomeNode* gn = genome_feature_new(gft_exon, r1, STRAND_BOTH, NULL, 0, env);
  GenomeNode* gn2 = genome_feature_new(gft_exon, r2, STRAND_BOTH, NULL, 0, env);

  Element* e = element_new(gn, NULL, env);
  Element* e2 = element_new(gn, NULL, env);
  Element* e3 = element_new(gn2, NULL, env);

  /* tests element_get_range */
  r_temp = element_get_range(e);
  ensure(has_err, (0 == range_compare(r1, r_temp)));
  ensure(has_err, (1 == range_compare(r2, r_temp)));

  /* tests element_get_arrow_status */
  as_temp = element_get_arrow_status(e);
  ensure(has_err, (as1 == as_temp));
  ensure(has_err, !(as2 == as_temp));

  /* tests element_set_arrow_status */
  element_set_arrow_status(e, as2);
  element_set_arrow_status(e2, as2);
  element_set_arrow_status(e3, as2);
  as_temp = element_get_arrow_status(e);
  ensure(has_err, (as2 == as_temp));
  ensure(has_err, !(as1 == as_temp));
 
  /* tests elements_are_equal */
  ensure(has_err, elements_are_equal(e, e2));
  ensure(has_err, !elements_are_equal(e, e3));
  ensure(has_err, !elements_are_equal(e2, e3));

  element_delete(e, env);
  element_delete(e2, env);
  element_delete(e3, env);
  genome_node_delete(gn, env);
  genome_node_delete(gn2, env);

  return has_err;

}

