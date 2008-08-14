/*
  Copyright (c) 2007      Christin Schaerfer <cschaerfer@zbh.uni-hamburg.de>
  Copyright (c)      2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

/* A block has a range, a caption, a parent caption, a strand, and a type and it
   contains element objects. */
typedef struct Block Block;

#include "libgtcore/dlist.h"
#include "libgtcore/range.h"
#include "libgtcore/array.h"
#include "libgtext/genome_node.h"
#include "libgtview/canvas.h"
#include "libgtview/diagram.h"
#include "libgtview/element.h"

Block*             block_new(void);
/* Create a new Block object, setting block parameters (such as strand, range)
   from a given <node> template. */
Block*             block_new_from_node(GenomeNode *node);
Block*             block_ref(Block*);
/* Insert <node> into block. */
void               block_insert_element(Block*, GenomeNode *node);
Range              block_get_range(const Block*);
Range*             block_get_range_ptr(const Block *block);
void               block_set_range(Block*, Range r);
/* Checks whether a Block is occupied completely by a single element. */
bool               block_has_only_one_fullsize_element(Block*);
/* Set whether a block caption should be displayed or not. */
void               block_set_caption_visibility(Block*, bool);
bool               block_caption_is_visible(const Block*);
void               block_set_caption(Block*, Str*);
Str*               block_get_caption(const Block*);
void               block_set_strand(Block*, Strand);
GenomeNode*        block_get_top_level_feature(Block*);
Strand             block_get_strand(const Block*);
void               block_set_type(Block*, GenomeFeatureType*);
GenomeFeatureType* block_get_type(const Block*);
unsigned long      block_get_size(Block*);
int                block_render(Block*, Canvas*);
int                block_unit_test(Error*);
void               block_delete(Block*);

#endif
