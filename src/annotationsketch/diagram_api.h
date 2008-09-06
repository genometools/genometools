/*
  Copyright (c) 2007 Malte Mader <mmader@zbh.uni-hamburg.de>
  Copyright (c) 2007 Sascha Steinbiss, Christin Schaerfer, Malte Mader
  Copyright (c) 2007 Christin Schaerfer <cschaerfer@zbh.uni-hamburg.de>
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

#ifndef DIAGRAM_API_H
#define DIAGRAM_API_H

#include "annotationsketch/feature_index.h"

typedef struct GT_Diagram GT_Diagram;

/* Create a new GT_Diagram object representing the genome nodes in
   <feature_index> in region <seqid> overlapping with <range>. */
GT_Diagram* gt_diagram_new(GT_FeatureIndex *feature_index, const char *seqid,
                           const Range*, GT_Style*);
GT_Diagram* gt_diagram_new_from_array(Array *features, const Range*, GT_Style*);
Range       gt_diagram_get_range(GT_Diagram*);
/* Render <diagram> on the given <canvas>. */
int         gt_diagram_sketch(GT_Diagram *diagram, GT_Canvas *canvas);
void        gt_diagram_delete(GT_Diagram*);

#endif
