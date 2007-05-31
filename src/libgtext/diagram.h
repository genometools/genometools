/*
  Copyright (c) 2007 Sascha Steinbiss, Christin Schaerfer, Malte Mader
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef DIAGRAM_H
#define DIAGRAM_H

#include <libgtcore/array.h>
#include <libgtext/config.h>
#include <libgtcore/range.h>
#include <libgtcore/hashtable.h>


typedef struct Diagram Diagram;

Diagram* diagram_new(Array*,Range,Config*,Env*);
void diagram_set_config(Diagram*,Config*,Env*);
Hashtable* diagram_get_tracks(Diagram*);
int* diagram_get_total_lines(Diagram*, Env*);
void diagram_delete(Diagram*,Env*);
int diagram_unit_test(Env*);

#endif
