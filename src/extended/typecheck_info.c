/*
  Copyright (c) 2013 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include "core/ma_api.h"
#include "extended/typecheck_info.h"
#include "extended/type_checker_builtin.h"
#include "extended/type_checker_obo.h"

struct GtTypecheckInfo {
  GtStr *typecheck;
  GtOption *typecheck_option;
  bool typecheck_built_in;
};

GtTypecheckInfo* gt_typecheck_info_new(void)
{
  GtTypecheckInfo *tci = gt_calloc(1, sizeof *tci);
  tci->typecheck = gt_str_new();
  return tci;
}

void gt_typecheck_info_delete(GtTypecheckInfo *tci)
{
  if (!tci) return;
  gt_option_delete(tci->typecheck_option);
  gt_str_delete(tci->typecheck);
  gt_free(tci);
}

void gt_typecheck_info_register_options(GtTypecheckInfo *tci,
                                        GtOptionParser *op)
{
  GtOption *typecheck_option, *built_in_option;
  gt_assert(tci && op);
  gt_assert(!tci->typecheck_option); /* can only called once */

  /* -typecheck */
  typecheck_option = gt_option_new_string("typecheck",
                                          "check GFF3 types against \"id\" and "
                                          "\"name\" tags in given OBO file and "
                                          "validate parent (part-of) "
                                          "relationships", tci->typecheck,
                                          NULL);
  gt_option_argument_is_optional(typecheck_option);
  gt_option_parser_add_option(op, typecheck_option);
  tci->typecheck_option = gt_option_ref(typecheck_option);

  /* -typecheck-built-in */
  built_in_option = gt_option_new_bool("typecheck-built-in",
                                       "use built-in type checker",
                                       &tci->typecheck_built_in, false);
  gt_option_is_development_option(built_in_option);
  gt_option_parser_add_option(op, built_in_option);
  gt_option_exclude(typecheck_option, built_in_option);
}

bool gt_typecheck_info_option_used(const GtTypecheckInfo *tci)
{
  gt_assert(tci && tci->typecheck_option);
  if (gt_option_is_set(tci->typecheck_option) || tci->typecheck_built_in)
    return true;
  return false;
}

GtTypeChecker* gt_typecheck_info_create_type_checker(const GtTypecheckInfo *tci,
                                                     GtError *err)
{
  GtTypeChecker *type_checker = NULL;
  gt_error_check(err);
  gt_assert(tci);
  if (tci->typecheck_built_in)
    type_checker = gt_type_checker_builtin_new();
  else {
    gt_assert(gt_option_is_set(tci->typecheck_option));
    type_checker = gt_type_checker_obo_new(gt_str_get(tci->typecheck), err);
  }
  return type_checker;
}
