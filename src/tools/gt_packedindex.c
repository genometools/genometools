/*
  Copyright (c) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/error.h"
#include "libgtcore/option.h"
#include "libgtcore/str.h"
#include "libgtext/toolbox.h"
#include "libgtcore/versionfunc.h"
#include "libgtmatch/sfx-run.h"
#include "tools/gt_packedindex.h"
#include "tools/gt_packedindex_mkctxmap.h"
#include "tools/gt_packedindex_trsuftab.h"
#include "tools/gt_packedindex_chk_integrity.h"
#include "tools/gt_packedindex_chk_search.h"

/* rely on suffixerator for on the fly index construction */
static int gt_packedindex_make(int argc, const char *argv[], Error *err)
{
  return parseargsandcallsuffixerator(false, argc, argv, err);
}

static void* gt_packedindex_arguments_new(void)
{
  Toolbox *packedindex_toolbox = toolbox_new();
  toolbox_add(packedindex_toolbox, "mkindex", gt_packedindex_make);
  toolbox_add(packedindex_toolbox, "mkctxmap", gt_packedindex_mkctxmap);
  toolbox_add(packedindex_toolbox, "trsuftab", gt_packedindex_trsuftab);
  toolbox_add(packedindex_toolbox, "chkintegrity",
              gt_packedindex_chk_integrity );
  toolbox_add(packedindex_toolbox, "chksearch", gt_packedindex_chk_search);
  return packedindex_toolbox;
}

static void gt_packedindex_arguments_delete(void *tool_arguments)
{
  Toolbox *index_toolbox = tool_arguments;
  if (!index_toolbox) return;
  toolbox_delete(index_toolbox);
}

static OptionParser* gt_packedindex_option_parser_new(void *tool_arguments)
{
  Toolbox *index_toolbox = tool_arguments;
  OptionParser *op;
  assert(index_toolbox);
  op = option_parser_new("[option ...] index_tool [argument ...]",
                         "Call packed index tool with name index_tool and "
                         "pass argument(s) to it.");
  option_parser_set_comment_func(op, toolbox_show, index_toolbox);
  return op;
}

static int gt_packedindex_runner(int argc, const char **argv, int parsed_args,
                                 void *tool_arguments, Error *err)
{
  Toolbox *index_toolbox = tool_arguments;
  Toolfunc toolfunc;
  Tool *tool = NULL;
  char **nargv = NULL;
  int had_err = 0;

  error_check(err);
  assert(index_toolbox);

  /* determine tool */
  if (!toolbox_has_tool(index_toolbox, argv[parsed_args])) {
    error_set(err, "packedindex tool '%s' not found; option -help lists "
                   "possible tools", argv[parsed_args]);
    had_err = -1;
  }

  /* call sub-tool */
  if (!had_err) {
    if (!(toolfunc = toolbox_get(index_toolbox, argv[parsed_args]))) {
      tool = toolbox_get_tool(index_toolbox, argv[parsed_args]);
      assert(tool);
    }
    nargv = cstr_array_prefix_first(argv + parsed_args,
                                    error_get_progname(err));
    error_set_progname(err, nargv[0]);
    if (toolfunc)
      had_err = toolfunc(argc - parsed_args, (const char**) nargv, err);
    else
      had_err = tool_run(tool, argc - parsed_args, (const char**) nargv, err);
  }

  cstr_array_delete(nargv);
  return had_err;
}

Tool* gt_packedindex(void)
{
  return tool_new(gt_packedindex_arguments_new,
                  gt_packedindex_arguments_delete,
                  gt_packedindex_option_parser_new,
                  NULL,
                  gt_packedindex_runner);
}
