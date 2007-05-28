/*
   Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>   
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ELEMENT_H
#define ELEMENT_H

#include <libgtcore/range.h>
#include <libgtext/genome_node.h>
#include <libgtext/config.h>

typedef struct Element Element;

Element* element_new(GenomeNode *gn,
                     Config *cfg,
		     Env* env);
Range element_get_range(Element* element);
void element_set_arrow_status(Element* element,
                              int status);
int element_get_arrow_status(Element* element);
void element_delete(Element* element,
                    Env* env);
void print_element(Element* element);
bool elements_are_equal(Element* e1,
                        Element* e2);
int element_unit_test(Env* env);

#endif

