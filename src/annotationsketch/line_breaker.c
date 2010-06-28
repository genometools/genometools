/*
  Copyright (c) 2008-2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2010 Center for Bioinformatics, University of Hamburg

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
#include "annotationsketch/line_breaker_rep.h"

struct GtLineBreakerMembers {
  unsigned int reference_count;
};

struct GtLineBreakerClass {
  size_t size;
  GtLineBreakerIsOccupiedFunc is_occupied;
  GtLineBreakerRegisterBlockFunc register_block;
  GtLineBreakerFreeFunc free;
};

const GtLineBreakerClass* gt_line_breaker_class_new(size_t size,
                                  GtLineBreakerIsOccupiedFunc is_occupied,
                                  GtLineBreakerRegisterBlockFunc register_block,
                                  GtLineBreakerFreeFunc free)
{
  GtLineBreakerClass *c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->is_occupied = is_occupied;
  c_class->register_block = register_block;
  c_class->free = free;
  return c_class;
}

GtLineBreaker* gt_line_breaker_create(const GtLineBreakerClass *lbc)
{
  GtLineBreaker *lb;
  gt_assert(lbc && lbc->size);
  lb = gt_calloc(1, lbc->size);
  lb->c_class = lbc;
  lb->pvt = gt_calloc(1, sizeof (GtLineBreakerMembers));
  return lb;
}

GtLineBreaker* gt_line_breaker_ref(GtLineBreaker *lb)
{
  gt_assert(lb);
  lb->pvt->reference_count++;
  return lb;
}

void gt_line_breaker_delete(GtLineBreaker *lb)
{
  if (!lb) return;
  if (lb->pvt->reference_count) {
    lb->pvt->reference_count--;
    return;
  }
  gt_assert(lb->c_class);
  if (lb->c_class->free)
    lb->c_class->free(lb);
  gt_free(lb->pvt);
  gt_free(lb);
}

int gt_line_breaker_line_is_occupied(GtLineBreaker *lb, bool *result,
                                     GtLine *line, GtBlock *block, GtError *err)
{
  gt_assert(lb && lb->c_class && result && line && block);
  return lb->c_class->is_occupied(lb, result, line, block, err);
}

int gt_line_breaker_register_block(GtLineBreaker *lb, GtLine *line,
                                    GtBlock *block, GtError *err)
{
  gt_assert(lb && lb->c_class && line && block);
  return lb->c_class->register_block(lb, line, block, err);
}

void* gt_line_breaker_cast(GT_UNUSED const GtLineBreakerClass *lbc,
                           GtLineBreaker *lb)
{
  gt_assert(lbc && lb && lb->c_class == lbc);
  return lb;
}
