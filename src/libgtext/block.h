/*
   Copyright (c) Sascha Steinbiss, Malte Mader, Christin Schaerfer
   Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef BLOCK_H
#define BLOCK_H

#include <libgtcore/range.h>
#include <libgtext/genome_node.h>
#include <libgtcore/array.h>
#include <libgtext/element.h>
#include <libgtext/config.h>

typedef struct Block Block;

Block* block_new(Env* env);
void block_insert_element(Block* block, 
                          GenomeNode* gn, 
			  Config* cfg, 
			  Env* env);
Range block_get_range(Block* block);
Array* block_get_elements(Block* block);
void block_delete(Block* block,
                  Env* env);

#endif

