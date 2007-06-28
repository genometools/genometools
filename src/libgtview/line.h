/*
   Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef LINE_H
#define LINE_H

#include <libgtcore/array.h>
#include <libgtext/genome_node.h>
#include <libgtview/config.h>
#include <libgtview/block.h>


/*!
Line interface
A line contains block objects
*/
typedef struct Line Line;

/*!
Creates a new Line object.
\param env Pointer to Environment object.
\return Pointer to a new Line object.
*/
Line* line_new(Env* env);

/*!
insert Block object in Line object
\param line Pointer to Line object
\param block Pointer to Block object to instert
\param env Pointer to environment object
*/
void line_insert_block(Line* line,
                       Block* block,
		       Env* env);

/*!
Checks if Line is occupied
\param line Line object to check
\param gn Pointer to GenomeNode object
\return True or False
*/
bool line_is_occupied(Line* line,
                      Range r);

/*!
Returns Array with Pointer to Block objects
\param line Pointer to Line object
\return Pointer to Array
*/
Array* line_get_blocks(Line* line);

/*!
Tests Line Class
\param env Pointer to Environment object
*/
int line_unit_test(Env* env);

/*!
Delets Line
\param line Pointer to Line object to delete
\param env Pointer to Environment object
*/
void line_delete(Line* line,
                 Env* env);

#endif

