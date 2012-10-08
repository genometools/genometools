/*
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "core/logger.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "match/esa-seqread.h"
#include "match/rdj-errfind.h"
#include "match/rdj-version.h"
#include "tools/gt_readjoiner_correct.h"

#define GT_READJOINER_CORRECT_TAG ".corrected"

typedef struct {
  unsigned long k, c;
  GtStr *indexname;
  unsigned long debug_value;
  bool edit_twobitencoding;
} GtReadjoinerCorrectArguments;

static void* gt_readjoiner_correct_arguments_new(void)
{
  GtReadjoinerCorrectArguments *arguments = gt_calloc((size_t)1,
      sizeof *arguments);
  arguments->indexname = gt_str_new();
  return arguments;
}

static void gt_readjoiner_correct_arguments_delete(void *tool_arguments)
{
  GtReadjoinerCorrectArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->indexname);
  gt_free(arguments);
}

static GtOptionParser* gt_readjoiner_correct_option_parser_new(
    void *tool_arguments)
{
  GtReadjoinerCorrectArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  op = gt_option_parser_new("[options] -ii indexname",
      "Readjoiner k-mer based error correction.");

  option = gt_option_new_ulong_min("k","k-mer length",
                                  &arguments->k, 31UL, 1UL);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_ulong_min("c","minimal trusted count",
                                  &arguments->c, 3UL, 1UL);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("encseq","correct encoded sequence",
      &arguments->edit_twobitencoding, true);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_string("ii", "input index",
                                arguments->indexname, NULL);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_ulong("dbg","(debug value) "
      "see code for current meaning", &arguments->debug_value, GT_UNDEF_ULONG);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_version_func(op, gt_readjoiner_show_version);
  gt_option_parser_set_max_args(op, 0);

  return op;
}

static int gt_readjoiner_correct_runner(GT_UNUSED int argc,
    GT_UNUSED const char **argv, GT_UNUSED int parsed_args,
    void *tool_arguments, GT_UNUSED GtError *err)
{
  GtReadjoinerCorrectArguments *arguments = tool_arguments;
  Sequentialsuffixarrayreader *ssar;
  GtLogger *logger;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);
  logger = gt_logger_new(true, GT_LOGGER_DEFLT_PREFIX, stderr);

  ssar = gt_newSequentialsuffixarrayreaderfromfile(
      gt_str_get(arguments->indexname), SARR_LCPTAB | SARR_SUFTAB |
      SARR_ESQTAB | SARR_SSPTAB, SEQ_scan, logger, err);
  if (ssar == NULL)
    had_err = -1;
  else
  {
    GtStr *indexname_corrected = gt_str_clone(arguments->indexname);
    gt_str_append_cstr(arguments->indexname, GT_READJOINER_CORRECT_TAG);
    had_err = gt_errfind(ssar, gt_encseqSequentialsuffixarrayreader(ssar),
           arguments->k, arguments->c, arguments->debug_value,
           arguments->edit_twobitencoding, gt_str_get(indexname_corrected),
           err);
    gt_freeSequentialsuffixarrayreader(&ssar);
    gt_str_delete(indexname_corrected);
  }
  gt_logger_delete(logger);
  return had_err;
}

GtTool* gt_readjoiner_correct(void)
{
  return gt_tool_new(gt_readjoiner_correct_arguments_new,
                  gt_readjoiner_correct_arguments_delete,
                  gt_readjoiner_correct_option_parser_new,
                  NULL,
                  gt_readjoiner_correct_runner);
}
