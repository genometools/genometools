/*
   Copyright (c) 2007 Sascha Steinbiss, Christin Schaerfer, Malte Mader
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtcore/array.h>
#include <libgtext/genome_feature.h>
#include <libgtext/genome_feature_type.h>
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
  Element *element;
  GenomeFeature *gf = (GenomeFeature*) gn;

  env_error_check(env);
  element = env_ma_malloc(env, sizeof (Element));
  element->type = genome_feature_get_type(gf);
  element->range = genome_node_get_range(gn);
  element->arrow_status = NoArrow;
  element->cfg = cfg;
  return element;
}

/*!
Sets ArrowStatus of an Element object
\param element Element on which status to set.
\param status ArrowStatus.
*/
void element_set_arrow_status(Element *element, int status)
{
  element->arrow_status = status;
}

/*!
Returns ArrowStatus of an Element object
\param element Pointer to Element object
\return ArrowStatus
*/
int element_get_arrow_status(Element *element)
{
  return element->arrow_status;
}

/*!
Returns range of an Element object
\param block Pointer to Element object
\return Pointer to Range object
*/
Range element_get_range(Element *element)
{
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

