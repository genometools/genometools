/*
  Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef BLOCK_H
#define BLOCK_H

#include "libgtcore/dlist.h"
#include "libgtcore/range.h"
#include "libgtcore/array.h"
#include "libgtext/genome_node.h"
#include "libgtview/element.h"
#include "libgtview/config.h"

/* A block has a range, a caption, a parent caption, a strand, and a type and it
   contains element objects. */
typedef struct Block Block;

Block*            block_new(void);
/* Create a new Block object, setting block parameters (such as strand, range)
   from a given <node> template. */
Block*            block_new_from_node(GenomeNode *node);
/* Insert <node> into block. */
void              block_insert_element(Block*, GenomeNode *node, Config*);
Range             block_get_range(const Block*);
void              block_set_range(Block*, Range r);
/* Checks whether a Block is occupied completely by a single element. */
bool              block_has_only_one_fullsize_element(Block*);
/* Set whether a block caption should be displayed or not. */
void              block_set_caption_visibility(Block*, bool);
bool              block_caption_is_visible(const Block*);
void              block_set_caption(Block*, Str*);
Str*              block_get_caption(const Block*);
void              block_set_strand(Block*, Strand);
Strand            block_get_strand(const Block*);
void              block_set_type(Block*, GenomeFeatureType);
GenomeFeatureType block_get_type(const Block*);
/* Returns Dlist with Pointer to Element objects. */
Dlist*            block_get_elements(const Block*);
int               block_unit_test(Error*);
void              block_delete(Block*);

#endif
