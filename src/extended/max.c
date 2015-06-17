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

struct Gtmax{
    GtWord value;
    GtUwordPair *start;
    GtUwordPair *end;
};

Gtmax* gt_max_new(void)
{
  Gtmax *max;
  max = gt_calloc((size_t) 1, sizeof (Gtmax));
  max->value = 0;

  max->start = gt_malloc(sizeof *(max->start));
  max->start->a=0;
  max->start->b=0;

  max->end = gt_malloc(sizeof *(max->end));
  max->end->a=0;
  max->end->b=0;
  return max;
}

void gt_max_delete(Gtmax *max)
{
  if (max != NULL)
    gt_free(max);
}

void gt_max_set_value(Gtmax *max, GtWord value)
{
  gt_assert(max!=NULL);
  max->value=value;
}

GtWord gt_max_get_value(Gtmax *max)
{
  gt_assert(max!=NULL);
  return(max->value);
}

void gt_max_set_start(Gtmax *max, GtUwordPair *start )
{
  gt_assert(max!=NULL);
  max->start=start;
}

GtUwordPair* gt_max_get_start(Gtmax *max)
{
  gt_assert(max!=NULL);
  return(max->start);
}

void gt_max_set_end_with_pair(Gtmax *max, GtUwordPair *end)
{
  gt_assert(max!=NULL && end!=NULL);
  max->end=end;
}

void gt_max_set_end(Gtmax *max, GtUword a, GtUword b)
{
  gt_assert(max!=NULL);
  max->end->a = a;
  max->end->b =b ;
}

GtUwordPair* gt_max_get_end(Gtmax *max)
{
  gt_assert(max!=NULL);
  return(max->end);
}

GtUword gt_max_get_row_length(Gtmax *max)
{
  GtUword length;
  gt_assert(max!=NULL);

  length = (gt_max_get_end(max))->a -(gt_max_get_start(max))->a;
  return(length);
}

GtUword gt_max_get_col_length(Gtmax *max)
{
  GtUword length;
  gt_assert(max!=NULL);

  length = (gt_max_get_end(max))->b -(gt_max_get_start(max))->b;
  return(length);
}
