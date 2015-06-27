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

#include "core/unused_api.h"
#include "core/assert_api.h"
#include "core/error_api.h"
#include "core/ma_api.h"
#include "extended/max.h"

struct Gtmaxcoordvalue{
    GtWord value;
    GtUwordPair start;
    GtUwordPair end;
};

Gtmaxcoordvalue* gt_max_new(void)
{
  Gtmaxcoordvalue *max;
  max = gt_calloc((size_t) 1, sizeof (Gtmaxcoordvalue));
  max->value = 0;

  max->start.a=0;
  max->start.b=0;

  max->end.a=0;
  max->end.b=0;
  return max;
}

void gt_max_delete(Gtmaxcoordvalue *max)
{
  if (max != NULL)
    gt_free(max);
}

void gt_max_set_value(Gtmaxcoordvalue *max, GtWord value)
{
  gt_assert(max != NULL);
  max->value=value;
}

GtWord gt_max_get_value(Gtmaxcoordvalue *max)
{
  gt_assert(max != NULL);
  return(max->value);
}

void gt_max_set_start(Gtmaxcoordvalue *max, GtUwordPair start )
{
  gt_assert(max != NULL);
  max->start=start;
}

GtUwordPair gt_max_get_start(Gtmaxcoordvalue *max)
{
  gt_assert(max != NULL);
  return(max->start);
}

void gt_max_set_end_with_pair(Gtmaxcoordvalue *max, GtUwordPair end)
{
  gt_assert(max != NULL);
  max->end = end;
}

void gt_max_set_end(Gtmaxcoordvalue *max, GtUword a, GtUword b)
{
  gt_assert(max != NULL);
  max->end.a = a;
  max->end.b = b ;
}

GtUwordPair gt_max_get_end(Gtmaxcoordvalue *max)
{
  gt_assert(max != NULL);
  return(max->end);
}

GtUword gt_max_get_row_length(Gtmaxcoordvalue *max)
{
  gt_assert(max != NULL);

  GtUword end = (gt_max_get_end(max)).a;
  GtUword start = (gt_max_get_start(max)).a;

  gt_assert(end >= start);
  return end-start;
}

GtUword gt_max_get_col_length(Gtmaxcoordvalue*max)
{
  gt_assert(max != NULL);

  GtUword end = (gt_max_get_end(max)).b;
  GtUword start = (gt_max_get_start(max)).b;

  gt_assert(end >= start);
  return end-start;
}
