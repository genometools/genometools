/*
  Copyright (c) 2010      Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>
  Copyright (c)      2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg

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

#include "core/class_alloc.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "match_iterator_api.h"
#include "match_iterator_rep.h"

const GtMatchIteratorClass* gt_match_iterator_class_new(size_t size,
                                           GtMatchIteratorNextFunc next,
                                           GtMatchIteratorFreeFunc free)
{
  GtMatchIteratorClass *c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->free = free;
  c_class->next = next;
  return c_class;
}

GtMatchIterator* gt_match_iterator_create(const GtMatchIteratorClass *mpc)
{
  GtMatchIterator *mp;
  gt_assert(mpc && mpc->size);
  mp = gt_calloc(1, mpc->size);
  mp->c_class = mpc;
  return mp;
}

GtMatchIteratorStatus gt_match_iterator_next(GtMatchIterator *mp,
                                             GtMatch **match, GtError *err)
{
  gt_assert(mp);
  if (mp->c_class->next)
    return mp->c_class->next(mp, match, err);
  return GT_MATCHER_STATUS_END;
}

void* gt_match_iterator_cast(GT_UNUSED const GtMatchIteratorClass *mc,
                             GtMatchIterator *m)
{
  gt_assert(mc && m && m->c_class == mc);
  return m;
}

void gt_match_iterator_delete(GtMatchIterator *mp)
{
  if (!mp) return;
  if (mp->c_class->free)
    mp->c_class->free(mp);
  gt_free(mp);
}
