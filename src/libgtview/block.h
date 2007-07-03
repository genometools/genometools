/*
   Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef BLOCK_H
#define BLOCK_H

#include <libgtcore/range.h>
#include <libgtext/genome_node.h>
#include <libgtcore/array.h>
#include <libgtview/element.h>
#include <libgtview/config.h>

/*!
Block interface
A block has a range, a caption, a parent caption,
a strand and a type and it contains element objects.
*/
typedef struct Block Block;

/*!
Creates a new Block object.
\param env Pointer to Environment object.
\return Pointer to a new Block object.
*/
Block* block_new(Env* env);

/*!
Inserts an element into a Block object
\param block Block in which the element shall be insert
\param gn Pointer to GenomeNode object
\param cfg Pointer to Config file
\param env Pointer to Environment object
*/
void block_insert_element(Block* block,
                          GenomeNode* gn,
			  Config* cfg,
			  Env* env);

/*!
Returns range of a Block object
\param block Pointer to Block object
\return Pointer to Range object
*/
Range block_get_range(Block* block);

/*!
Sets range of a Block object
\param block Pointer to Block object to set range
\param r Range to set
*/
void block_set_range(Block* block,
                     Range r);

/*!
Sets whether a block caption should be displayed or not.
\param block Pointer to Block object to set visibility for.
\param val TRUE to show caption, FALSE to disable
*/
void block_set_caption_visibility(Block* block, bool val);

/*!
Returns block caption visibility setting.
\param block Pointer to Block object to query.
\return TRUE if visible, FALSE otherwise.
*/
bool block_caption_is_visible(Block* block);

/*!
Sets caption of a Block object
\param block Pointer to Block object to set caption
\param caption Pointer to String object
*/
void block_set_caption(Block* block,
                       const char* caption);

/*!
Gets caption of a Block object
\param block Pointer to Block object
\return caption Pointer to String object
*/
const char* block_get_caption(Block* block);

/*!
Sets parent_caption of a Block object
\param block Pointer to Block object to set caption
\param caption Pointer to String object
*/
void block_set_parent_caption(Block* block,
                              const char* caption);

/*!
Gets parent_caption of a Block object
\param block Pointer to Block object
\return caption Pointer to String object
*/
const char* block_get_parent_caption(Block* block);

/*!
Sets strand of a Block object
\param block Block to set strand
\param strand Strand to set
*/
void block_set_strand(Block* block,
                      Strand strand);

/*!
Gets strand of a Block object
\param block Pointer to Block object
\return strand Strand
*/
Strand block_get_strand(Block* block);

/*!
Sets type of a Block object
\param block Block to set type
\param type GenomeFeatureType to set
*/
void block_set_type(Block* block,
                    GenomeFeatureType type);

/*!
Gets type of a Block object
\param block Pointer to Block object
\returen GenomeFeature Type
*/
GenomeFeatureType block_get_type(Block* block);

/*!
Returns Array with Pointer to Element objects
\param block Pointer to Block object
\return Pointer to Array
*/
Dlist* block_get_elements(Block* block);

/*!
Unit Test for Block Class
\param env Pointer to Environment object
*/
int block_unit_test(Env* env);

/*!
Delets Block
\param block Pointer to Block object to delete
\param env Pointer to Environment object
*/
void block_delete(Block* block,
                  Env* env);

#endif

