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

#ifndef LINE_BREAKER_CAPTIONS_H
#define LINE_BREAKER_CAPTIONS_H

#include "core/range.h"
#include "annotationsketch/canvas.h"
#include "annotationsketch/layout.h"
#include "annotationsketch/line_breaker.h"
#include "annotationsketch/style.h"

/* Implements the GtLineBreaker interface; returns new lines if captions
   overlap. */
typedef struct GtLineBreakerCaptions GtLineBreakerCaptions;

const GtLineBreakerClass* gt_line_breaker_captions_class(void);
GtLineBreaker*            gt_line_breaker_captions_new(GtLayout*,
                                                       unsigned long width,
                                                       GtStyle*);

#endif
