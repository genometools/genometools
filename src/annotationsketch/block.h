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
typedef struct GT_Block GT_Block;

#include "annotationsketch/canvas.h"
#include "annotationsketch/diagram.h"
#include "annotationsketch/element.h"
#include "core/dlist.h"
#include "core/range.h"
#include "core/array.h"
#include "extended/genome_node.h"

GT_Block*             gt_block_new(void);
/* Create a new GT_Block object, setting block parameters (such as strand,
   range) from a given <node> template. */
GT_Block*             gt_block_new_from_node(GT_GenomeNode *node);
GT_Block*             gt_block_ref(GT_Block*);
/* Insert <node> into block. */
void                  gt_block_insert_element(GT_Block*, GT_GenomeNode *node);
GT_Range              gt_block_get_range(const GT_Block*);
GT_Range*             gt_block_get_range_ptr(const GT_Block *block);
void                  gt_block_set_range(GT_Block*, GT_Range r);
/* Checks whether a GT_Block is occupied completely by a single element. */
bool                  gt_block_has_only_one_fullsize_element(const GT_Block*);
/* Set whether a block caption should be displayed or not. */
void                  gt_block_set_caption_visibility(GT_Block*, bool);
bool                  gt_block_caption_is_visible(const GT_Block*);
void                  gt_block_set_caption(GT_Block*, GT_Str*);
GT_Str*               gt_block_get_caption(const GT_Block*);
void                  gt_block_set_strand(GT_Block*, GT_Strand);
GT_GenomeNode*        gt_block_get_top_level_feature(const GT_Block*);
GT_Strand             gt_block_get_strand(const GT_Block*);
void                  gt_block_set_type(GT_Block*, GT_FeatureType*);
GT_FeatureType* gt_block_get_type(const GT_Block*);
unsigned long         gt_block_get_size(const GT_Block*);
int                   gt_block_sketch(GT_Block*, GT_Canvas*);
int                   gt_block_compare(const GT_Block*, const GT_Block*);
int                   gt_block_unit_test(GT_Error*);
void                  gt_block_delete(GT_Block*);

#endif
