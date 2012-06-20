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
#include "extended/string_iter.h"
#include "extended/string_iter_rep.h"

int gt_string_iter_reset(GtStringIter *string_iter,
                         GtError *err)
{
  gt_error_check(err);
  gt_assert(string_iter && string_iter->c_class);
  if (string_iter->c_class->reset_func != NULL)
    return string_iter->c_class->reset_func(string_iter, err);
  return 0;
}

int gt_string_iter_next(GtStringIter *string_iter,
                        const char **string,
                        GtError *err)
{
  gt_error_check(err);
  gt_assert(string_iter && string_iter->c_class);
  if (string_iter->c_class->next_func != NULL)
    return string_iter->c_class->next_func(string_iter,
                                           string,
                                           err);
  return 0;
}

void gt_string_iter_delete(GtStringIter *string_iter)
{
  if (string_iter != NULL) {
    gt_assert(string_iter->c_class != NULL);
    if (string_iter->c_class->delete_func != NULL)
      string_iter->c_class->delete_func(string_iter);
    gt_free(string_iter);
  }
}

GtStringIterClass *gt_string_iter_class_new(size_t size,
                                           GtStringIterNextFunc next,
                                           GtStringIterResetFunc reset,
                                           GtStringIterDeleteFunc delete)
{
  GtStringIterClass *str_iter_c;
  gt_assert(size != 0);
  str_iter_c = gt_class_alloc(sizeof (*str_iter_c));
  str_iter_c->size = size;
  str_iter_c->next_func = next;
  str_iter_c->reset_func = reset;
  str_iter_c->delete_func = delete;
  return str_iter_c;
}

GtStringIter *gt_string_iter_create(const GtStringIterClass *str_iter_c)
{
  GtStringIter *str_iter;
  gt_assert(str_iter_c && str_iter_c->size);
  str_iter = gt_calloc((size_t) 1, str_iter_c->size);
  str_iter->c_class = str_iter_c;
  return str_iter;
}

void *gt_string_iter_cast(GT_UNUSED const GtStringIterClass *str_iter_c,
                         GtStringIter *str_iter)
{
  gt_assert(str_iter_c && str_iter);
  gt_assert(str_iter->c_class == str_iter_c);
  return str_iter;
}
