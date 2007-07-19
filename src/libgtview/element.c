/*
   Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/
/**
 * \file element.c
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
  Config* cfg;
};

Element* element_new(GenomeNode *gn, Config *cfg, Env *env)
{
  assert(gn);

  Element *element;
  GenomeFeature *gf = (GenomeFeature*) gn;

  env_error_check(env);
  element = env_ma_malloc(env, sizeof (Element));
  element->type = genome_feature_get_type(gf);
  element->range = genome_node_get_range(gn);
  element->cfg = cfg;

  assert(element);
  return element;
}

Element* element_new_empty(Config *cfg,
                           Env *env)
{
  Element *element;

  env_error_check(env);
  element = env_ma_malloc(env, sizeof (Element));
  element->cfg = cfg;

  assert(element);
  return element;
}

GenomeFeatureType element_get_type(Element *element)
{
  assert(element);
  return element->type;
}

void element_set_type(Element *element,
                      GenomeFeatureType type)
{
  assert(element);
  element->type = type;
}

Range element_get_range(Element *element)
{
  assert(element);

  return element->range;
}

void element_set_range(Element *element,
                       Range r)
{
  assert(element);

  element->range = r;
}

bool elements_are_equal(Element* e1,
                        Element* e2)
{
  assert(e1 && e2);
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

  r1.start = 10;
  r1.end = 50;

  r2.start = 20;
  r2.end = 50;

  GenomeNode* gn = genome_feature_new(gft_exon, r1,
                                      STRAND_BOTH, NULL, 0, env);
  GenomeNode* gn2 = genome_feature_new(gft_exon, r2,
                                       STRAND_BOTH, NULL, 0, env);

  Element* e = element_new(gn, NULL, env);
  Element* e2 = element_new(gn, NULL, env);
  Element* e3 = element_new(gn2, NULL, env);

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

