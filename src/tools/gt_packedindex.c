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
#include "libgtcore/env.h"
#include "libgtcore/option.h"
#include "libgtcore/str.h"
#include "libgtcore/cstr.h"
#include "libgtext/toolbox.h"
#include "libgtcore/versionfunc.h"
#include "libgtmatch/sfx-run.pr"
#include "gt_packedindex_chk_integrity.h"
#include "gt_packedindex_chk_search.h"

static void
register_packedindextools(Toolbox *packedindex_toolbox, Env *env);

static OPrval
parse_subtool_options(int *parsed_args, int argc, const char **argv,
                      Toolbox *index_toolbox, Env *env);

int
gt_packedindex(int argc, const char **argv, Env *env)
{
  Toolbox *index_toolbox;
  Tool indexTool;
  int parsed_args;
  bool had_err = false;
  char **nargv = NULL;
  env_error_check(env);

  index_toolbox = toolbox_new();
  register_packedindextools(index_toolbox, env);

  switch (parse_subtool_options(&parsed_args, argc, argv, index_toolbox, env))
  {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      toolbox_delete(index_toolbox);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      toolbox_delete(index_toolbox);
      return 0;
  }
  assert(parsed_args < argc);

  /* determine tool */
  if (!(indexTool = toolbox_get(index_toolbox, argv[1]))) {
    env_error_set(env, "packedindex tool '%s' not found; option -help lists "
                  "possible tools", argv[1]);
    had_err = true;
  }

  /* call sub-tool */
  if (!had_err) {
    nargv = cstr_array_prefix_first(argv+parsed_args, argv[0], env);
    env_error_set_progname(env, nargv[0]);
    had_err = indexTool(argc-parsed_args, (const char**) nargv, env);
  }

  cstr_array_delete(nargv, env);
  toolbox_delete(index_toolbox);
  return had_err?-1:0;
}

static OPrval
parse_subtool_options(int *parsed_args, int argc, const char **argv,
                      Toolbox *index_toolbox, Env *env)
{
  OptionParser *op;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] index_tool [argument ...]",
                         "Call packed index tool with name index_tool and "
                         "pass argument(s) to it.", env);
  option_parser_set_comment_func(op, toolbox_show, index_toolbox);
  oprval = option_parser_parse_min_args(op, parsed_args, argc, argv,
                                        versionfunc, 1, env);
  option_parser_delete(op, env);
  return oprval;
}

static int
gt_packedindex_make(int argc, const char *argv[], Env *env);

static void
register_packedindextools(Toolbox *packedindex_toolbox, Env *env)
{
  assert(packedindex_toolbox);
  toolbox_add(packedindex_toolbox, "mkindex", gt_packedindex_make);
  toolbox_add(packedindex_toolbox, "chkintegrity",
              gt_packedindex_chk_integrity );
  toolbox_add(packedindex_toolbox, "chksearch", gt_packedindex_chk_search);
}

/***************************************************************************
 * rely on suffixerator for on the fly index construction
 ***************************************************************************/

static int
gt_packedindex_make(int argc, const char *argv[], Env *env)
{
  return parseargsandcallsuffixerator(false, argc, argv, env);
}

