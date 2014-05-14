/*
  Copyright (c) 2013 Gordon Gremme <gordon@gremme.org>

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

#include "core/cstr_api.h"
#include "core/error.h"
#include "core/fileutils.h"
#include "core/gtdatapath.h"
#include "core/ma_api.h"
#include "extended/typecheck_info.h"
#include "extended/type_checker_builtin_api.h"
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
  typecheck_option =
    gt_option_new_string("typecheck", "check GFF3 types against \"id\" and "
                         "\"name\" tags in given OBO file and validate parent "
                         "(part-of) relationships.\nIf no argument is given, "
                         "the sofa.obo file from the gtdata/obo_files "
                         "directory is used.\nIf an argument is given, it is "
                         "used as an OBO filename.\nIn the case that such a "
                         "file does not exist '.obo' is added to the argument "
                         "and loading the resulting filename from the "
                         "gtdata/obo_files directory is attempted.",
                         tci->typecheck, NULL);
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

static GtStr* get_obo_path(GtError *err)
{
  const char *progname;
  GtStr *obo_path, *prog;
  gt_error_check(err);
  progname = gt_error_get_progname(err);
  gt_assert(progname != NULL);
  prog = gt_str_new();
  gt_str_append_cstr_nt(prog, progname,
                        gt_cstr_length_up_to_char(progname, ' '));
  obo_path = gt_get_gtdata_path(gt_str_get(prog), err);
  if (obo_path)
    gt_str_append_cstr(obo_path, "/obo_files/");
  gt_str_delete(prog);
  return obo_path;
}

GtTypeChecker* gt_typecheck_info_create_type_checker(const GtTypecheckInfo *tci,
                                                     GtError *err)
{
  GtTypeChecker *type_checker = NULL;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(tci);
  if (tci->typecheck_built_in)
    type_checker = gt_type_checker_builtin_new();
  else {
    GtStr *obo_file;
    gt_assert(gt_option_is_set(tci->typecheck_option));
    if (!gt_str_length(tci->typecheck)) {
      /* a. */
      if (!(obo_file = get_obo_path(err)))
        had_err = -1;
      if (!had_err)
        gt_str_append_cstr(obo_file, "sofa.obo");
    }
    else if (gt_file_exists(gt_str_get(tci->typecheck))) {
      /* b. */
      obo_file = gt_str_new_cstr(gt_str_get(tci->typecheck));
    }
    else {
      /* c. */
      if (!(obo_file = get_obo_path(err)))
        had_err = -1;
      if (!had_err) {
        gt_str_append_str(obo_file, tci->typecheck);
        gt_str_append_cstr(obo_file, ".obo");
      }
    }

    if (!had_err)
      type_checker = gt_type_checker_obo_new(gt_str_get(obo_file), err);

    gt_str_delete(obo_file);
  }
  return type_checker;
}
