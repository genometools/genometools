/*
  Copyright (c) 2008, 2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008       Center for Bioinformatics, University of Hamburg

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
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/type_checker_rep.h"

GtTypeChecker* gt_type_checker_create(const GtTypeCheckerClass
                                       *type_checker_class)
{
  GtTypeChecker *type_checker;
  gt_assert(type_checker_class && type_checker_class->size);
  type_checker = gt_calloc(1, type_checker_class->size);
  type_checker->c_class = type_checker_class;
  return type_checker;
}

GtTypeChecker* gt_type_checker_ref(GtTypeChecker *type_checker)
{
  gt_assert(type_checker);
  type_checker->reference_count++;
  return type_checker;
}

const char* gt_type_checker_description(GtTypeChecker *type_checker)
{
  gt_assert(type_checker && type_checker->c_class &&
            type_checker->c_class->description);
  return type_checker->c_class->description(type_checker);
}

bool gt_type_checker_is_valid(GtTypeChecker *type_checker, const char *type)
{
  gt_assert(type_checker && type_checker->c_class &&
            type_checker->c_class->is_valid);
  return type_checker->c_class->is_valid(type_checker, type);
}

bool gt_type_checker_is_partof(GtTypeChecker *type_checker,
                               const char *parent_type, const char *child_type)
{
  gt_assert(type_checker && type_checker->c_class &&
            type_checker->c_class->is_partof);
  return type_checker->c_class->is_partof(type_checker, parent_type,
                                          child_type);
}

void gt_type_checker_delete(GtTypeChecker *type_checker)
{
  if (!type_checker) return;
  if (type_checker->reference_count) {
    type_checker->reference_count--;
    return;
  }
  gt_assert(type_checker->c_class);
  if (type_checker->c_class->free)
    type_checker->c_class->free(type_checker);
  gt_free(type_checker);
}

void* gt_type_checker_cast(GT_UNUSED const GtTypeCheckerClass
                           *type_checker_class,
                           GtTypeChecker *type_checker)
{
  gt_assert(type_checker_class && type_checker &&
         type_checker->c_class == type_checker_class);
  return type_checker;
}
