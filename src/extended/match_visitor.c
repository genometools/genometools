/*
  Copyright (c) 2011 Sascha Kastens <sascha.kastens@studium.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include "core/assert_api.h"
#include "core/class_alloc.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/match_visitor_rep.h"

struct GtMatchVisitorClass {
  size_t size;
  GtMatchVisitorFreeFunc free;
  GtMatchVisitorBlastFunc match_blast;
  GtMatchVisitorOpenFunc match_open;
};

const GtMatchVisitorClass*
gt_match_visitor_class_new(size_t size,
                           GtMatchVisitorFreeFunc free,
                           GtMatchVisitorBlastFunc match_blast,
                           GtMatchVisitorOpenFunc match_open)
{
  GtMatchVisitorClass *c_class;
  gt_assert(size);
  c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->free = free;
  c_class->match_blast = match_blast;
  c_class->match_open = match_open;
  return c_class;
}

GtMatchVisitor* gt_match_visitor_create(const GtMatchVisitorClass *mvc)
{
  GtMatchVisitor *mv;
  gt_assert(mvc && mvc->size);
  mv = gt_calloc(1, mvc->size);
  mv->c_class = mvc;
  return mv;
}

void* gt_match_visitor_cast(GT_UNUSED const GtMatchVisitorClass *mvc,
                            GtMatchVisitor *mv)
{
  gt_assert(mvc && mv && mv->c_class == mvc);
  return mv;
}

int gt_match_visitor_visit_match_blast(GtMatchVisitor *mv, GtMatchBlast *mb,
                                       GtError *err)
{
  gt_error_check(err);
  gt_assert(mv && mb && mv->c_class && mv->c_class->match_blast);
  return mv->c_class->match_blast(mv, mb, err);
}

int gt_match_visitor_visit_match_open(GtMatchVisitor *mv, GtMatchOpen *mo,
                                      GtError *err)
{
  gt_error_check(err);
  gt_assert(mv && mo && mv->c_class && mv->c_class->match_open);
  return mv->c_class->match_open(mv, mo, err);
}

void gt_match_visitor_delete(GtMatchVisitor *mv)
{
  if (!mv) return;
  gt_assert(mv->c_class);
  if (mv->c_class->free)
    mv->c_class->free(mv);
  gt_free(mv);
}
