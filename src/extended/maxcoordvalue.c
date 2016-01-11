/*
  Copyright (c) 2015 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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
#include "extended/maxcoordvalue.h"

struct GtMaxcoordvalue{
    GtWord value;
    GtUwordPair start;
    GtUwordPair end;
};

GtMaxcoordvalue* gt_maxcoordvalue_new(void)
{
  GtMaxcoordvalue *max;
  max = gt_calloc((size_t) 1, sizeof (GtMaxcoordvalue));
  max->value = 0;

  max->start.a=0;
  max->start.b=0;

  max->end.a=0;
  max->end.b=0;
  return max;
}

void gt_maxcoordvalue_delete(GtMaxcoordvalue *max)
{
  if (max != NULL)
    gt_free(max);
}

static void gt_maxcoordvalue_set_value(GtMaxcoordvalue *max, const GtWord value)
{
  gt_assert(max != NULL);
  max->value=value;
}

GtWord gt_maxcoordvalue_get_value(const GtMaxcoordvalue *max)
{
  gt_assert(max != NULL);
  return(max->value);
}

void gt_maxcoordvalue_set_start(GtMaxcoordvalue *max,
                                GtUword a, GtUword b)
{
  gt_assert(max != NULL);
  max->start.a = a;
  max->start.b = b ;
}

static void gt_maxcoordvalue_set_start_with_pair(GtMaxcoordvalue *max,
                                                 const GtUwordPair start )
{
  gt_assert(max != NULL);
  max->start=start;
}

GtUwordPair gt_maxcoordvalue_get_start(const GtMaxcoordvalue *max)
{
  gt_assert(max != NULL);
  return(max->start);
}

void gt_maxcoordvalue_set_end_with_pair(GtMaxcoordvalue *max,
                                        const GtUwordPair end)
{
  gt_assert(max != NULL);
  max->end = end;
}

static void gt_maxcoordvalue_set_end(GtMaxcoordvalue *max,
                                     GtUword a, GtUword b)
{
  gt_assert(max != NULL);
  max->end.a = a;
  max->end.b = b ;
}

GtUwordPair gt_maxcoordvalue_get_end(const GtMaxcoordvalue *max)
{
  gt_assert(max != NULL);
  return(max->end);
}

/*use this in linear space cases*/
void gt_maxcoordvalue_coord_update(GtMaxcoordvalue *max,
                                   GtWord value,
                                   GtUwordPair start,
                                   GtUword enda, GtUword endb)
{
  gt_assert(max != NULL);

  gt_maxcoordvalue_set_value(max, value);
  gt_maxcoordvalue_set_start_with_pair(max, start);
  gt_maxcoordvalue_set_end(max, enda, endb);
}

/*use this in square space cases*/
void gt_maxcoordvalue_coord_update_without_start (GtMaxcoordvalue *max,
                                                  GtWord value,
                                                  GtUword enda,
                                                  GtUword endb)
{
  gt_assert(max != NULL);

  gt_maxcoordvalue_set_value(max, value);
  gt_maxcoordvalue_set_end(max, enda, endb);
}

GtUword gt_maxcoordvalue_get_row_length(const GtMaxcoordvalue *max)
{
  gt_assert(max != NULL);

  GtUword end = (gt_maxcoordvalue_get_end(max)).a;
  GtUword start = (gt_maxcoordvalue_get_start(max)).a;

  gt_assert(end >= start);
  return end-start;
}

GtUword gt_maxcoordvalue_get_col_length(const GtMaxcoordvalue *max)
{
  gt_assert(max != NULL);

  GtUword end = (gt_maxcoordvalue_get_end(max)).b;
  GtUword start = (gt_maxcoordvalue_get_start(max)).b;

  gt_assert(end >= start);
  return end-start;
}

bool gt_maxcoordvalue_get_length_safe(const GtMaxcoordvalue *max)
{
  if (gt_maxcoordvalue_get_end(max).a == gt_maxcoordvalue_get_start(max).a &&
      gt_maxcoordvalue_get_end(max).b == gt_maxcoordvalue_get_start(max).b  )
    return false;
  return true;
}

void gt_maxcoordvalue_reset(GtMaxcoordvalue *max)
{
  gt_assert(max != NULL);

  max->value = 0;
  max->start.a=0;
  max->start.b=0;
  max->end.a=0;
  max->end.b=0;
}
