/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"
#include "gt_guessprot.h"
#include "gt_png.h"
#include "gt_regioncov.h"
#include "gt_sfxmap.h"

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            Toolbox *dev_toolbox, Env *env)
{
  OptionParser *op;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] dev_tool_name [argument ...]",
                         "Call development tool with name dev_tool_name and "
                         "pass argument(s) to it.", env);
  option_parser_set_comment_func(op, toolbox_show, dev_toolbox);
  oprval = option_parser_parse_min_args(op, parsed_args, argc, argv,
                                        versionfunc, 1, env);
  option_parser_delete(op, env);
  return oprval;
}

void register_devtools(Toolbox *dev_toolbox, Env *env)
{
  assert(dev_toolbox);
  /* add development tools here with a function call like this:
     toolbox_add(dev_toolbox, "devtool", gt_devtool, env); */
  toolbox_add(dev_toolbox, "guessprot", gt_guessprot, env);
  toolbox_add(dev_toolbox, "png", gt_png, env);
  toolbox_add(dev_toolbox, "regioncov", gt_regioncov, env);
  toolbox_add(dev_toolbox, "sfxmap", gt_sfxmap, env);
}

int gt_dev(int argc, const char **argv, Env *env)
{
  Toolbox *dev_toolbox;
  Tool devtool;
  int parsed_args, has_err = 0;
  char **nargv = NULL;
  env_error_check(env);

  /* option parsing */
  dev_toolbox = toolbox_new(env);
  register_devtools(dev_toolbox, env);
  switch (parse_options(&parsed_args, argc, argv, dev_toolbox, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      toolbox_delete(dev_toolbox, env);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      toolbox_delete(dev_toolbox, env);
      return 0;
  }
  assert(parsed_args < argc);

  /* get development tools */
  if (!(devtool = toolbox_get(dev_toolbox, argv[1]))) {
    env_error_set(env, "development tool '%s' not found; option -help lists "
                  "possible tools", argv[1]);
    has_err = -1;
  }

  /* call development tool */
  if (!has_err) {
    nargv = cstr_array_prefix_first(argv+parsed_args, argv[0], env);
    has_err = devtool(argc-parsed_args, (const char**) nargv, env);
  }

  /* free */
  cstr_array_delete(nargv, env);
  toolbox_delete(dev_toolbox, env);

  return has_err;
}
