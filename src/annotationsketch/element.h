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
typedef struct Element Element;

#include "annotationsketch/canvas.h"
#include "annotationsketch/drawing_range.h"
#include "core/range.h"
#include "core/strand.h"
#include "extended/genome_node.h"
#include "extended/genome_feature_type.h"

/* Creates a complete new Element object. */
Element*           element_new(GenomeNode*);
/* Creates an empty Element object. Range and type have to be set afterwards. */
Element*           element_new_empty(void);
Range              element_get_range(const Element*);
void               element_set_range(Element*, Range);
DrawingRange       element_calculate_drawing_range(Element*, Canvas*);
GenomeFeatureType* element_get_type(const Element*);
void               element_set_type(Element*, GenomeFeatureType*);
Strand             element_get_strand(const Element*);
GenomeNode*        element_get_node_ref(const Element*);
bool               element_is_marked(const Element*);
bool               elements_are_equal(const Element*, const Element*);
int                element_sketch(Element*, Canvas*);
int                element_unit_test(Error*);
void               element_delete(Element* element);

#endif
