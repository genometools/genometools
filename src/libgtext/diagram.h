/*
  Copyright (c) 2007 Sascha Steinbiss, Christin Schaerfer, Malte Mader
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FEATUREINDEX_H
#define FEATUREINDEX_H

#include <libgtcore/array.h>
#include <libgtcore/str.h>
#include <libgtext/config.h>
#include <libgtcore/range.h>

typedef struct Diagram Diagram;

Diagram* diagram_new(Array*,Range,Config*,Env*);
void diagram_set_config(Diagram*,Config*,Env*);
void diagram_delete(Diagram*,Env*);

#endif
