/*
  Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef TRACK_H
#define TRACK_H

#include "libgtcore/str.h"
#include "libgtcore/array.h"
#include "libgtext/genome_node.h"
#include "libgtview/config.h"
#include "libgtview/line.h"

/* A track has a title and a type und contains line objects. */
typedef struct Track Track;

Track* track_new(Str* title, Env* env);
void   track_insert_block(Track*, Block*, Env*);
Str*   track_get_title(const Track*);
/* Returns Array containing Pointers to Line objects. */
Array* track_get_lines(const Track*);
int    track_get_number_of_lines(const Track*);
int    track_unit_test(Env*);
void   track_delete(Track*, Env*);

#endif

