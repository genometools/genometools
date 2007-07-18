/*
   Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/
/**
 * \file element.h
 * \author Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
 */

#ifndef ELEMENT_H
#define ELEMENT_H

#include <libgtcore/range.h>
#include "libgtext/genome_node.h"
#include "libgtext/genome_feature_type.h"
#include <libgtview/config.h>

/*!
Element interface
An element has a type, a range and a config object.
*/
typedef struct Element Element;

/*!
Creates a new Element object.
\param gn Pointer to GenomeNode object.
\param cfg Pointer to Config object.
\param env Pointer to Environment object.
\return Pointer to a new Element object.
*/
Element* element_new(GenomeNode *gn,
                     Config *cfg,
                     Env* env);

/*!
Creates a new Element object.
range and type have to be set afterwards with
element_set_range and element_set_type
\param cfg Pointer to Config object.
\param env Pointer to Environment object.
\return Pointer to a new Element object.
*/
Element* element_new_empty(Config* cfg,
                           Env* env);

/*!
Returns range of an Element object
\param block Pointer to Element object
\return Pointer to Range object
*/
Range element_get_range(Element* element);

/*!
Sets Range of an Element object
\param element Element to set range
\param r Range to set
*/
void element_set_range(Element* element,
                       Range r);

/*!
Returns Type of an Element object
\param element Pointer to Element object
\return GenomeFeatureType
*/
GenomeFeatureType element_get_type(Element* element);

/*!
Sets Type of an Element object
\param element Element to set type
\param type GenomeFeatureType to set
*/
void element_set_type(Element *element,
                      GenomeFeatureType type);

/*!
Checks if two Element objects are equal
\param e1 Pointer to Element object
\param e2 Pointer to Element object
\returns true if e1 and e2 are equal
*/
bool elements_are_equal(Element* e1,
                        Element* e2);

/*!
Unit Test for Element Class
\param env Pointer to Environment object
*/
int element_unit_test(Env* env);

/*!
Delets Element
\param element Pointer to Element object to delete
\param env Pointer to Environment object
*/
void element_delete(Element* element,
                    Env* env);

#endif

