/*
  Copyright (c) 2013 Sascha Steinbiss <ss34@sanger.ac.uk>

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
#include "extended/xrfcheck_info.h"
#include "extended/xrf_checker_api.h"

struct GtXRFCheckInfo {
  GtStr *xrfcheck;
  GtOption *xrfcheck_option;
};

GtXRFCheckInfo* gt_xrfcheck_info_new(void)
{
  GtXRFCheckInfo *xci = gt_calloc(1, sizeof *xci);
  xci->xrfcheck = gt_str_new();
  return xci;
}

void gt_xrfcheck_info_delete(GtXRFCheckInfo *xci)
{
  if (!xci) return;
  gt_option_delete(xci->xrfcheck_option);
  gt_str_delete(xci->xrfcheck);
  gt_free(xci);
}

void gt_xrfcheck_info_register_options(GtXRFCheckInfo *xci,
                                        GtOptionParser *op)
{
  GtOption *xrfcheck_option;
  gt_assert(xci && op);
  gt_assert(!xci->xrfcheck_option); /* can only called once */

  /* -xrfcheck */
  xrfcheck_option =
    gt_option_new_string("xrfcheck", "check Dbxref and Ontology_term "
                         "attributes for correct syntax according to a "
                         "abbreviation definition file.\nIf no argument is "
                         "given, the GO.xrf_abbs file from the "
                         "gtdata/xrf_abbr directory is used.\nIf an argument "
                         "is given, it is used as an specific filename for an "
                         "abbreviation file.\nIn the case that such a "
                         "file does not exist, '.xrf_abbr' is added to the "
                         "argument and loading the resulting filename from the "
                         "gtdata/xrf_abbr directory is attempted.",
                         xci->xrfcheck, NULL);
  gt_option_argument_is_optional(xrfcheck_option);
  gt_option_parser_add_option(op, xrfcheck_option);
  xci->xrfcheck_option = gt_option_ref(xrfcheck_option);
}

bool gt_xrfcheck_info_option_used(const GtXRFCheckInfo *xci)
{
  gt_assert(xci && xci->xrfcheck_option);
  return gt_option_is_set(xci->xrfcheck_option);
}

static GtStr* get_xrf_path(GtError *err)
{
  const char *progname;
  GtStr *xrf_path, *prog;
  gt_error_check(err);
  progname = gt_error_get_progname(err);
  gt_assert(progname != NULL);
  prog = gt_str_new();
  gt_str_append_cstr_nt(prog, progname,
                        gt_cstr_length_up_to_char(progname, ' '));
  xrf_path = gt_get_gtdata_path(gt_str_get(prog), err);
  if (xrf_path)
    gt_str_append_cstr(xrf_path, "/xrf_abbr/");
  gt_str_delete(prog);
  return xrf_path;
}

GtXRFChecker* gt_xrfcheck_info_create_xrf_checker(const GtXRFCheckInfo *xci,
                                                  GtError *err)
{
  GtXRFChecker *xrf_checker = NULL;
  int had_err = 0;
  GtStr *xrf_file;
  gt_error_check(err);
  gt_assert(xci);
  gt_assert(gt_option_is_set(xci->xrfcheck_option));
  if (!gt_str_length(xci->xrfcheck)) {
    /* a. */
    if (!(xrf_file = get_xrf_path(err)))
      had_err = -1;
    if (!had_err)
      gt_str_append_cstr(xrf_file, "GO.xrf_abbr");
  }
  else if (gt_file_exists(gt_str_get(xci->xrfcheck))) {
    /* b. */
    xrf_file = gt_str_new_cstr(gt_str_get(xci->xrfcheck));
  }
  else {
    /* c. */
    if (!(xrf_file = get_xrf_path(err)))
      had_err = -1;
    if (!had_err) {
      gt_str_append_str(xrf_file, xci->xrfcheck);
      gt_str_append_cstr(xrf_file, ".xrf_abbr");
    }
  }

  if (!had_err)
    xrf_checker = gt_xrf_checker_new(gt_str_get(xrf_file), err);

  gt_str_delete(xrf_file);
  return xrf_checker;
}
