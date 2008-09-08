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

#ifndef ELEMENT_H
#define ELEMENT_H

/* An element has a type, a range and a config object. */
typedef struct GT_Element GT_Element;

#include "annotationsketch/canvas.h"
#include "annotationsketch/drawing_range.h"
#include "core/range.h"
#include "core/strand.h"
#include "extended/genome_node.h"
#include "extended/genome_feature_type.h"

/* Creates a complete new GT_Element object. */
GT_Element*           gt_element_new(GT_GenomeNode*);
/* Creates an empty GT_Element object. GT_Range and type have to be set afterwards. */
GT_Element*           gt_element_new_empty(void);
GT_Range              gt_element_get_range(const GT_Element*);
void                  gt_element_set_range(GT_Element*, GT_Range);
GT_DrawingRange       gt_element_calculate_drawing_range(GT_Element*,
                                                         GT_Canvas*);
GT_GenomeFeatureType* gt_element_get_type(const GT_Element*);
void                  gt_element_set_type(GT_Element*, GT_GenomeFeatureType*);
GT_Strand             gt_element_get_strand(const GT_Element*);
GT_GenomeNode*        gt_element_get_node_ref(const GT_Element*);
bool                  gt_element_is_marked(const GT_Element*);
bool                  gt_elements_are_equal(const GT_Element*,
                                            const GT_Element*);
int                   gt_element_sketch(GT_Element*, GT_Canvas*);
int                   gt_element_unit_test(GT_Error*);
void                  gt_element_delete(GT_Element* element);

#endif
