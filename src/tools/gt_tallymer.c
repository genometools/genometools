/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#include "core/cstr_array.h"
#include "core/error.h"
#include "core/option.h"
#include "core/str.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "extended/toolbox.h"
#include "tools/gt_tallymer.h"

static int gt_tallymer_mkindex(GT_UNUSED int argc, 
                               GT_UNUSED const char *argv[], 
                               GT_UNUSED GtError *err)
{
  gt_error_set(err,"tallymer mkindex not implemented yet");
  return -1;
}

static int gt_tallymer_occratio(GT_UNUSED int argc, 
                                GT_UNUSED const char *argv[], 
                                GtError *err)
{
  gt_error_set(err,"tallymer occratio not implemented yet");
  return -1;
}

static int gt_tallymer_search(GT_UNUSED int argc, 
                                GT_UNUSED const char *argv[], 
                                GtError *err)
{
  gt_error_set(err,"tallymer search not implemented yet");
  return -1;
}

static void *gt_tallymer_arguments_new(void)
{
  GtToolbox *tallymer_toolbox = gt_toolbox_new();
  gt_toolbox_add(tallymer_toolbox, "mkindex", gt_tallymer_mkindex);
  gt_toolbox_add(tallymer_toolbox, "occratio", gt_tallymer_occratio);
  gt_toolbox_add(tallymer_toolbox, "search", gt_tallymer_search);
  return tallymer_toolbox;
}

static void gt_tallymer_arguments_delete(void *tool_arguments)
{
  GtToolbox *index_toolbox = tool_arguments;
  if (!index_toolbox) return;
  gt_toolbox_delete(index_toolbox);
}

static GtOptionParser* gt_tallymer_option_parser_new(void *tool_arguments)
{
  GtToolbox *index_toolbox = tool_arguments;
  GtOptionParser *op;

  assert(index_toolbox != NULL);
  op = gt_option_parser_new(
                    "[option ...] [mkindex|occratio|search] [argument ...]",
                    "Call tallymer with specific tool and "
                    "pass argument(s) to it.");
  gt_option_parser_set_comment_func(op, gt_toolbox_show, index_toolbox);
  gt_option_parser_refer_to_manual(op);
  return op;
}

static int gt_tallymer_runner(int argc, const char **argv, int parsed_args,
                                 void *tool_arguments, GtError *err)
{
  GtToolbox *index_toolbox = tool_arguments;
  GtToolfunc toolfunc;
  GtTool *tool = NULL;
  char **nargv = NULL;
  int had_err = 0;

  gt_error_check(err);
  assert(index_toolbox != NULL);

  /* determine tool */
  if (!gt_toolbox_has_tool(index_toolbox, argv[parsed_args])) 
  {
    gt_error_set(err, "tallymer tool '%s' not found; option -help lists "
                   "possible tools", argv[parsed_args]);
    had_err = -1;
  }

  /* call sub-tool */
  if (!had_err) 
  {
    if (!(toolfunc = gt_toolbox_get(index_toolbox, argv[parsed_args]))) 
    {
      tool = gt_toolbox_get_tool(index_toolbox, argv[parsed_args]);
      assert(tool != NULL);
    }
    nargv = gt_cstr_array_prefix_first(argv + parsed_args,
                                       gt_error_get_progname(err));
    gt_error_set_progname(err, nargv[0]);
    if (toolfunc != NULL)
      had_err = toolfunc(argc - parsed_args, (const char**) nargv, err);
    else
      had_err = gt_tool_run(tool, argc - parsed_args, (const char**) nargv,
                            err);
  }
  gt_cstr_array_delete(nargv);
  return had_err;
}

GtTool* gt_tallymer(void)
{
  return gt_tool_new(gt_tallymer_arguments_new,
                     gt_tallymer_arguments_delete,
                     gt_tallymer_option_parser_new,
                     NULL,
                     gt_tallymer_runner);
}
