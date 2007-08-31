/*
  Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ELEMENT_H
#define ELEMENT_H

#include "libgtcore/range.h"
#include "libgtext/genome_node.h"
#include "libgtext/genome_feature_type.h"
#include "libgtview/config.h"

/* An element has a type, a range and a config object. */
typedef struct Element Element;

/* Creates a complete new Element object. */
Element*          element_new(GenomeNode*, Env*);
/* Creates an empty Element object. Range and type have to be set afterwards. */
Element*          element_new_empty(Env*);
Range             element_get_range(const Element*);
void              element_set_range(Element*, Range);
GenomeFeatureType element_get_type(const Element*);
void              element_set_type(Element*, GenomeFeatureType);
bool              element_is_marked(const Element*);
bool              elements_are_equal(const Element*, const Element*);
int               element_unit_test(Env* env);
void              element_delete(Element* element, Env* env);

#endif
