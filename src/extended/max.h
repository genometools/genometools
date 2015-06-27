/*
  Copyright (C) 2015 Annika Seidel, annika.seidel@studium.uni-hamburg.de

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
#ifndef MAX_H
#define MAX_H
#include "core/types_api.h"

typedef struct Gtmaxcoordvalue Gtmaxcoordvalue;

Gtmaxcoordvalue* gt_max_new(void);

void gt_max_delete(Gtmaxcoordvalue *max);

void gt_max_set_value(Gtmaxcoordvalue *max, GtWord value);

GtWord gt_max_get_value(Gtmaxcoordvalue *max);

void gt_max_set_start(Gtmaxcoordvalue*max, GtUwordPair start );

GtUwordPair gt_max_get_start(Gtmaxcoordvalue *max);

void gt_max_set_end_with_pair(Gtmaxcoordvalue *max, GtUwordPair end);

void gt_max_set_end(Gtmaxcoordvalue *max, GtUword a, GtUword b);

GtUwordPair gt_max_get_end(Gtmaxcoordvalue *max);

GtUword gt_max_get_row_length(Gtmaxcoordvalue *max);

GtUword gt_max_get_col_length(Gtmaxcoordvalue *max);
#endif
