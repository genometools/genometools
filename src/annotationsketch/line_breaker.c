/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
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

#include <assert.h>
#include "core/ma.h"
#include "core/unused_api.h"
#include "annotationsketch/line_breaker_rep.h"

GT_LineBreaker* gt_line_breaker_create(const GT_LineBreakerClass *lbc)
{
  GT_LineBreaker *lb;
  assert(lbc && lbc->size);
  lb = ma_calloc(1, lbc->size);
  lb->c_class = lbc;
  return lb;
}

GT_LineBreaker* gt_line_breaker_ref(GT_LineBreaker *lb)
{
  assert(lb);
  lb->reference_count++;
  return lb;
}

void gt_line_breaker_delete(GT_LineBreaker *lb)
{
  if (!lb) return;
  if (lb->reference_count) {
    lb->reference_count--;
    return;
  }
  assert(lb->c_class);
  if (lb->c_class->free)
    lb->c_class->free(lb);
  ma_free(lb);
}

bool gt_line_breaker_gt_line_is_occupied(GT_LineBreaker *lb, GT_Line *line, GT_Block *block)
{
  assert(lb && lb->c_class && line && block);
  return lb->c_class->is_occupied(lb, line, block);
}

void gt_line_breaker_register_block(GT_LineBreaker *lb, GT_Line *line, GT_Block *block)
{
  assert(lb && lb->c_class && line && block);
  lb->c_class->register_block(lb, line, block);
}

void* gt_line_breaker_cast(GT_UNUSED const GT_LineBreakerClass *lbc,
                        GT_LineBreaker *lb)
{
  assert(lbc && lb && lb->c_class == lbc);
  return lb;
}
