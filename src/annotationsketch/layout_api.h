/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef LAYOUT_API_H
#define LAYOUT_API_H

#include "annotationsketch/canvas_api.h"
#include "annotationsketch/diagram_api.h"
#include "annotationsketch/drawing_range.h"
#include "annotationsketch/style_api.h"
#include "annotationsketch/text_width_calculator_api.h"
#include "core/range.h"

/* The <GtLayout> class represents contents (tracks) of a <GtDiagram> broken up
   into lines such that a given horizontal space allotment given in pixels
   or points is used up most efficiently. This is done using the <GtLineBreaker>
   and <GtTextWidthCalculator> classes. As defaults, Cairo-based instances of
   these classes are used but can be specified separately.

   A <GtLayout> can be queried for the height of the laid out representation and
   finally be rendered to a <GtCanvas>. */
typedef struct GtLayout GtLayout;

/* Creates a new <GtLayout> object for the contents of <diagram>.
   The layout is done for a target image width of <width> and using the rules in
   <GtStyle> object <style>. */
GtLayout*              gt_layout_new(GtDiagram *diagram, unsigned int width,
                                     GtStyle*, GtError *err);
/* Like <gt_layout_new()>, but allows use of a different <GtTextWidthCalculator>
   implementation. */
GtLayout*              gt_layout_new_with_twc(GtDiagram*,
                                              unsigned int width,
                                              GtStyle*,
                                              GtTextWidthCalculator*,
                                              GtError *err);
/* Returns the height of the layout in pixels. */
unsigned long          gt_layout_get_height(const GtLayout*);
/* Renders the layout on the <target_canvas>. */
int                    gt_layout_sketch(GtLayout*, GtCanvas *target_canvas,
                                        GtError*);
/* Destroys a layout. */
void                   gt_layout_delete(GtLayout*);

#endif
