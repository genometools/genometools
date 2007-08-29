/*
   Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/
/**
 * \if INTERNAL \file element.c \endif
 * \author Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
 */

#include <string.h>
#include "libgtcore/array.h"
#include "libgtcore/ensure.h"
#include "libgtcore/strand.h"
#include "libgtext/genome_feature.h"
#include "libgtext/genome_feature_type.h"
#include "libgtview/element.h"

struct Element
{
  GenomeFeatureType type;
  Range range;
};

Element* element_new(GenomeNode *gn, Env *env)
{
  Element *element;
  GenomeFeature *gf = (GenomeFeature*) gn;

  assert(gn != NULL);

  env_error_check(env);
  element = element_new_empty(env);
  element_set_type(element,genome_feature_get_type(gf));
  element_set_range(element, genome_node_get_range(gn));
  return element;
}

Element* element_new_empty(Env *env)
{
  Element *element;

  env_error_check(env);
  element = env_ma_malloc(env, sizeof (Element));
  assert(element != NULL);
  return element;
}

GenomeFeatureType element_get_type(Element *element)
{
  assert(element != NULL);
  return element->type;
}

void element_set_type(Element *element,
                      GenomeFeatureType type)
{
  assert(element != NULL);
  element->type = type;
}

Range element_get_range(Element *element)
{
  assert(element != NULL);

  return element->range;
}

void element_set_range(Element *element,
                       Range r)
{
  assert(element != NULL);

  element->range = r;
}

bool elements_are_equal(Element* e1,
                        Element* e2)
{
  assert(e1 != NULL && e2 != NULL);
  if ((0 == strcmp(genome_feature_type_get_cstr(e1->type),
                   genome_feature_type_get_cstr(e2->type)))
     && (0 == range_compare(e1->range, e2->range)))
    return true;
  else
    return false;
}

int element_unit_test(Env* env)
{
  Range r1, r2, r_temp;
  int had_err = 0;
  GenomeNode *gn, *gn2;
  Element *e, *e2, *e3;

  r1.start = 10UL;
  r1.end = 50UL;

  r2.start = 20UL;
  r2.end = 50UL;

  gn = genome_feature_new(gft_exon, r1, STRAND_BOTH, "unit_test", 0, env);
  gn2 = genome_feature_new(gft_exon, r2, STRAND_BOTH, "unit_test", 0, env);

  e = element_new(gn, env);
  e2 = element_new(gn, env);
  e3 = element_new(gn2, env);

  /* tests element_get_range */
  r_temp = element_get_range(e);
  ensure(had_err, (0 == range_compare(r1, r_temp)));
  ensure(had_err, (1 == range_compare(r2, r_temp)));

  /* tests element_get_type and element_set_type*/
  ensure(had_err, (gft_exon == element_get_type(e)));
  ensure(had_err, (gft_intron != element_get_type(e)));
  element_set_type(e, gft_intron);
  ensure(had_err, (gft_intron == element_get_type(e)));
  element_set_type(e2, gft_intron);

  /* tests elements_are_equal */
  ensure(had_err, elements_are_equal(e, e2));
  ensure(had_err, !elements_are_equal(e, e3));
  ensure(had_err, !elements_are_equal(e2, e3));

  element_delete(e, env);
  element_delete(e2, env);
  element_delete(e3, env);
  genome_node_delete(gn, env);
  genome_node_delete(gn2, env);

  return had_err;

}

void element_delete(Element *element,
                    Env *env)
{
  if (!element) return;

  env_ma_free(element, env);
}

