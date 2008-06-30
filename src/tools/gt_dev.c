/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/cstr_array.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtext/toolbox.h"
#include "tools/gt_dev.h"
#include "tools/gt_extracttarget.h"
#include "tools/gt_guessprot.h"
#include "tools/gt_magicmatch.h"
#include "tools/gt_maxpairs.h"
#include "tools/gt_mergeesa.h"
#include "tools/gt_paircmp.h"
#include "tools/gt_patternmatch.h"
#include "tools/gt_regioncov.h"
#include "tools/gt_seqiterator.h"
#include "tools/gt_sfxmap.h"
#include "tools/gt_skproto.h"
#include "tools/gt_trieins.h"

static void* gt_dev_arguments_new(void)
{
  Toolbox *dev_toolbox = toolbox_new();
  /* add development tools here with a function call like this:
     toolbox_add(dev_toolbox, "devtool", gt_devtool); */
  toolbox_add_tool(dev_toolbox, "extracttarget", gt_extracttarget());
  toolbox_add(dev_toolbox, "guessprot", gt_guessprot);
  toolbox_add_tool(dev_toolbox, "magicmatch", gt_magicmatch());
  toolbox_add(dev_toolbox, "maxpairs", gt_maxpairs);
  toolbox_add(dev_toolbox, "mergeesa", gt_mergeesa);
  toolbox_add(dev_toolbox, "paircmp", gt_paircmp);
  toolbox_add(dev_toolbox, "patternmatch", gt_patternmatch);
  toolbox_add(dev_toolbox, "regioncov", gt_regioncov);
  toolbox_add(dev_toolbox, "seqiterator", gt_seqiterator);
  toolbox_add(dev_toolbox, "sfxmap", gt_sfxmap);
  toolbox_add_tool(dev_toolbox, "skproto", gt_skproto());
  toolbox_add(dev_toolbox, "trieins", gt_trieins);
  return dev_toolbox;
}

static void gt_dev_arguments_delete(void *tool_arguments)
{
  Toolbox *dev_toolbox = tool_arguments;
  if (!dev_toolbox) return;
  toolbox_delete(dev_toolbox);
}

static OptionParser* gt_dev_option_parser_new(void *tool_arguments)
{
  Toolbox *dev_toolbox = tool_arguments;
  OptionParser *op;
  assert(dev_toolbox);
  op = option_parser_new("[option ...] dev_tool_name [argument ...]",
                         "Call development tool with name dev_tool_name and "
                         "pass argument(s) to it.");
  option_parser_set_comment_func(op, toolbox_show, dev_toolbox);
  option_parser_set_min_args(op, 1);
  return op;
}

static int gt_dev_runner(int argc, const char **argv, int parsed_args,
                         void *tool_arguments, Error *err)
{
  Toolbox *dev_toolbox = tool_arguments;
  Toolfunc toolfunc;
  Tool *tool = NULL;
  int had_err = 0;
  char **nargv = NULL;

  error_check(err);
  assert(dev_toolbox);

  /* get development tools */
  if (!toolbox_has_tool(dev_toolbox, argv[parsed_args])) {
    error_set(err, "development tool '%s' not found; option -help lists "
                   "possible tools", argv[parsed_args]);
    had_err = -1;
  }

  /* call development tool */
  if (!had_err) {
    if (!(toolfunc = toolbox_get(dev_toolbox, argv[parsed_args]))) {
      tool = toolbox_get_tool(dev_toolbox, argv[parsed_args]);
      assert(tool);
    }
    nargv = cstr_array_prefix_first(argv + parsed_args,
                                    error_get_progname(err));
    error_set_progname(err, nargv[0]);
    if (toolfunc)
      had_err = toolfunc(argc - parsed_args , (const char**) nargv, err);
    else
      had_err = tool_run(tool, argc - parsed_args , (const char**) nargv, err);
  }

  /* free */
  cstr_array_delete(nargv);

  return had_err;
}

Tool* gt_dev(void)
{
  return tool_new(gt_dev_arguments_new,
                  gt_dev_arguments_delete,
                  gt_dev_option_parser_new,
                  NULL,
                  gt_dev_runner);
}
