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
#include "libgtcore/ma.h"
#include "libgtcore/unused.h"
#include "libgtview/line_breaker_rep.h"

LineBreaker* line_breaker_create(const LineBreakerClass *lbc)
{
  LineBreaker *lb;
  assert(lbc && lbc->size);
  lb = ma_calloc(1, lbc->size);
  lb->c_class = lbc;
  return lb;
}

LineBreaker* line_breaker_ref(LineBreaker *lb)
{
  assert(lb);
  lb->reference_count++;
  return lb;
}

void line_breaker_delete(LineBreaker *lb)
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

bool line_breaker_line_is_occupied(LineBreaker *lb, Line *line, Block *block)
{
  assert(lb && lb->c_class && line && block);
  return lb->c_class->is_occupied(lb, line, block);
}

void line_breaker_register_block(LineBreaker *lb, Line *line, Block *block)
{
  assert(lb && lb->c_class && line && block);
  lb->c_class->register_block(lb, line, block);
}

void* line_breaker_cast(UNUSED const LineBreakerClass *lbc,
                        LineBreaker *lb)
{
  assert(lbc && lb && lb->c_class == lbc);
  return lb;
}
