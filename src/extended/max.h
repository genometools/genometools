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

typedef struct Gtmax Gtmax;

Gtmax* gt_max_new(void);

void gt_max_delete(Gtmax *max);

void gt_max_set_value(Gtmax *max, GtWord value);

GtWord gt_max_get_value(Gtmax *max);

void gt_max_set_start(Gtmax *max, GtUwordPair *start );

GtUwordPair* gt_max_get_start(Gtmax *max);

void gt_max_set_end_with_pair(Gtmax *max, GtUwordPair *end);

void gt_max_set_end(Gtmax *max, GtUword a, GtUword b);

GtUwordPair* gt_max_get_end(Gtmax *max);

GtUword gt_max_get_row_length(Gtmax *max);

GtUword gt_max_get_col_length(Gtmax *max);
#endif
