/*
   Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef LINE_H
#define LINE_H

#include <libgtcore/array.h>
#include <libgtext/genome_node.h>
#include <libgtext/config.h>
#include <libgtext/block.h>

typedef struct Line Line;

Line* line_new(Env* env);
void line_insert_block(Line* line,
                       Block* block,
		       Env* env);
bool line_is_occupied(Line* line,
                      Range r);
Array* line_get_blocks(Line* line);
void line_delete(Line* line,
                 Env* env);
void print_line(Line* line);
int line_unit_test(Env* env);

#endif

