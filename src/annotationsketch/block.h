/*
  Copyright (c) 2007      Christin Schaerfer <schaerfer@zbh.uni-hamburg.de>
  Copyright (c)      2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#include "annotationsketch/block_api.h"
#include "annotationsketch/canvas.h"
#include "annotationsketch/style.h"
#include "core/error.h"

void        gt_block_insert_element(GtBlock*, GtFeatureNode *node);
void        gt_block_set_range(GtBlock*, GtRange r);
void        gt_block_set_type(GtBlock*, const char *type);
int         gt_block_compare(const GtBlock *block1, const GtBlock *block2,
                             void *data);
int         gt_block_sketch(GtBlock*, GtCanvas*, GtError*);
int         gt_block_get_max_height(const GtBlock *block,
                                    double *result,
                                    const GtStyle *sty,
                                    GtError *err);
void        gt_block_print(const GtBlock* block);

int         gt_block_unit_test(GtError*);

#endif
