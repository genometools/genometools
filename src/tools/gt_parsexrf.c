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

#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/xrf_abbr_parse_tree.h"
#include "extended/xrf_abbr_entry.h"
#include "tools/gt_parsexrf.h"

typedef struct {
  bool bool_option_parsexrf;
  GtStr  *str_option_parsexrf;
} GtParsexrfArguments;

static void* gt_parsexrf_arguments_new(void)
{
  GtParsexrfArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->str_option_parsexrf = gt_str_new();
  return arguments;
}

static void gt_parsexrf_arguments_delete(void *tool_arguments)
{
  GtParsexrfArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->str_option_parsexrf);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_parsexrf_option_parser_new(GT_UNUSED
                                                           void *tool_arguments)
{
  GtOptionParser *op;
  /* GtOption *option; */
  gt_assert(tool_arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [file]",
                            "Parse and validate XRF abbreviation files.");

  return op;
}

static int gt_parsexrf_runner(GT_UNUSED int argc, const char **argv,
                              int parsed_args, GT_UNUSED void *tool_arguments,
                              GtError *err)
{
  GtXRFAbbrParseTree *xpt;
  int had_err = 0;
  GtUword i;
  gt_error_check(err);

  xpt = gt_xrf_abbr_parse_tree_new(argv[parsed_args], err);
  if (!xpt)
    had_err = -1;

  if (!had_err) {
    for (i = 0; i < gt_xrf_abbr_parse_tree_num_of_entries(xpt); i++) {
      const GtXRFAbbrEntry *e = gt_xrf_abbr_parse_tree_get_entry(xpt, i);
      printf("%s\n", gt_xrf_abbr_entry_get_value(e, "abbreviation"));
    }
  }

  gt_xrf_abbr_parse_tree_delete(xpt);

  return had_err;
}

GtTool* gt_parsexrf(void)
{
  return gt_tool_new(gt_parsexrf_arguments_new,
                     gt_parsexrf_arguments_delete,
                     gt_parsexrf_option_parser_new,
                     NULL,
                     gt_parsexrf_runner);
}
