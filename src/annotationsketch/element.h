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

#ifndef ELEMENT_H
#define ELEMENT_H

/* An element has a type, a range and a config object. */
typedef struct GtElement GtElement;

#include "annotationsketch/canvas.h"
#include "annotationsketch/drawing_range.h"
#include "core/range.h"
#include "core/strand_api.h"
#include "extended/feature_type.h"
#include "extended/genome_node.h"

/* Creates a complete new <GtElement> object. */
GtElement*     gt_element_new(GtFeatureNode*);
GtElement*     gt_element_ref(GtElement*);
/* Creates an empty <GtElement> object.
   Range and type have to be set afterwards. */
GtElement*     gt_element_new_empty(void);
GtRange        gt_element_get_range(const GtElement*);
void           gt_element_set_range(GtElement*, GtRange);
const char*    gt_element_get_type(const GtElement*);
void           gt_element_set_type(GtElement*, const char *type);
GtStrand       gt_element_get_strand(const GtElement*);
GtFeatureNode* gt_element_get_node_ref(const GtElement*);
bool           gt_element_is_marked(const GtElement*);
int            gt_element_sketch(GtElement*, GtCanvas*, GtError*);
int            gt_element_unit_test(GtError*);
void           gt_element_delete(GtElement* element);

#endif
