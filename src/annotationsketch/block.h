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
typedef struct GtBlock GtBlock;

#include "annotationsketch/canvas.h"
#include "annotationsketch/diagram.h"
#include "annotationsketch/element.h"
#include "core/dlist.h"
#include "core/range.h"
#include "core/array.h"
#include "extended/genome_node.h"

GtBlock*              gt_block_new(void);
/* Create a new GtBlock object, setting block parameters (such as strand,
   range) from a given <node> template. */
GtBlock*              gt_block_new_from_node(GtGenomeFeature *node);
GtBlock*              gt_block_ref(GtBlock*);
/* Insert <node> into block. */
void                  gt_block_insert_element(GtBlock*, GtGenomeFeature *node);
GtRange               gt_block_get_range(const GtBlock*);
GtRange*              gt_block_get_range_ptr(const GtBlock *block);
void                  gt_block_set_range(GtBlock*, GtRange r);
/* Checks whether a GtBlock is occupied completely by a single element. */
bool                  gt_block_has_only_one_fullsize_element(const GtBlock*);
void                  gt_block_merge(GtBlock*, GtBlock*);
GtBlock*              gt_block_clone(GtBlock*);
/* Set whether a block caption should be displayed or not. */
void                  gt_block_set_caption_visibility(GtBlock*, bool);
bool                  gt_block_caption_is_visible(const GtBlock*);
void                  gt_block_set_caption(GtBlock*, GtStr*);
GtStr*                gt_block_get_caption(const GtBlock*);
void                  gt_block_set_strand(GtBlock*, GtStrand);
GtGenomeFeature*      gt_block_get_top_level_feature(const GtBlock*);
GtStrand              gt_block_get_strand(const GtBlock*);
void                  gt_block_set_type(GtBlock*, const char *type);
const char*           gt_block_get_type(const GtBlock*);
unsigned long         gt_block_get_size(const GtBlock*);
int                   gt_block_sketch(GtBlock*, GtCanvas*);
int                   gt_block_compare(const GtBlock*, const GtBlock*);
int                   gt_block_unit_test(GtError*);
void                  gt_block_delete(GtBlock*);

#endif
