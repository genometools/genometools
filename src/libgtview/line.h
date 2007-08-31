/*
  Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef LINE_H
#define LINE_H

#include "libgtcore/array.h"
#include "libgtext/genome_node.h"
#include "libgtview/config.h"
#include "libgtview/block.h"

/* A line contains block objects. */
typedef struct Line Line;

Line*  line_new(Env* env);
void   line_insert_block(Line*, Block*, Env*); /* takes ownership */
bool   line_is_occupied(const Line*, Range);
/* Returns Array containing Pointers to Block objects. */
Array* line_get_blocks(Line*);
int    line_unit_test(Env*);
void   line_delete(Line*, Env*);

#endif
