/*
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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
#include "extended/cstr_iterator.h"
#include "extended/cstr_iterator_rep.h"

int gt_cstr_iterator_reset(GtCstrIterator *cstr_iterator,
                           GtError *err)
{
  gt_error_check(err);
  gt_assert(cstr_iterator && cstr_iterator->c_class);
  if (cstr_iterator->c_class->reset_func != NULL)
    return cstr_iterator->c_class->reset_func(cstr_iterator, err);
  return 0;
}

int gt_cstr_iterator_next(GtCstrIterator *cstr_iterator,
                          const char **string,
                          GtError *err)
{
  gt_error_check(err);
  gt_assert(cstr_iterator && cstr_iterator->c_class);
  if (cstr_iterator->c_class->next_func != NULL)
    return cstr_iterator->c_class->next_func(cstr_iterator,
                                             string,
                                             err);
  return 0;
}

void gt_cstr_iterator_delete(GtCstrIterator *cstr_iterator)
{
  if (cstr_iterator != NULL) {
    gt_assert(cstr_iterator->c_class != NULL);
    if (cstr_iterator->c_class->delete_func != NULL)
      cstr_iterator->c_class->delete_func(cstr_iterator);
    gt_free(cstr_iterator);
  }
}

GtCstrIteratorClass *
gt_cstr_iterator_class_new(size_t size,
                           GtCstrIteratorNextFunc next,
                           GtCstrIteratorResetFunc reset,
                           GtCstrIteratorDeleteFunc delete)
{
  GtCstrIteratorClass *cstr_iterator_c;
  gt_assert(size != 0);
  cstr_iterator_c = gt_class_alloc(sizeof (*cstr_iterator_c));
  cstr_iterator_c->size = size;
  cstr_iterator_c->next_func = next;
  cstr_iterator_c->reset_func = reset;
  cstr_iterator_c->delete_func = delete;
  return cstr_iterator_c;
}

GtCstrIterator *
gt_cstr_iterator_create(const GtCstrIteratorClass *cstr_iterator_c)
{
  GtCstrIterator *cstr_iterator;
  gt_assert(cstr_iterator_c && cstr_iterator_c->size);
  cstr_iterator = gt_calloc((size_t) 1, cstr_iterator_c->size);
  cstr_iterator->c_class = cstr_iterator_c;
  return cstr_iterator;
}

void *
gt_cstr_iterator_cast(GT_UNUSED const GtCstrIteratorClass *cstr_iterator_c,
                      GtCstrIterator *cstr_iterator)
{
  gt_assert(cstr_iterator_c && cstr_iterator);
  gt_assert(cstr_iterator->c_class == cstr_iterator_c);
  return cstr_iterator;
}
