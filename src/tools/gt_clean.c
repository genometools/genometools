/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/bioseq.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/xposix.h"

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            Env *env)
{
  OptionParser *op;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("", "Remove all files in the current directory which "
                         "are automatically created by gt.", env);
  oprval = option_parser_parse_max_args(op, parsed_args, argc, argv,
                                        versionfunc, 0, env);
  option_parser_delete(op, env);
  return oprval;
}

static void remove_pattern_in_current_dir(const char *pattern, Env *env)
{
  char **files_to_remove;
  Str *path;
  glob_t g;

  path = str_new_cstr("./*", env);
  str_append_cstr(path, pattern, env);
  xglob(str_get(path), GLOB_NOCHECK, NULL, &g);

  /* remove found files */
  if (g.gl_pathc) {
    files_to_remove = g.gl_pathv;
    if (strcmp(*files_to_remove, str_get(path))) {
      while (*files_to_remove) {
        xunlink(*files_to_remove);
        files_to_remove++;
      }
    }
  }

  /* free */
  globfree(&g);
  str_delete(path, env);
}

int gt_clean(int argc, const char **argv, Env *env)
{
  int parsed_args;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args == 1);

  /* remove GT_BIOSEQ_INDEX files */
  remove_pattern_in_current_dir(GT_BIOSEQ_INDEX, env);

  /* remove GT_BIOSEQ_RAW files */
  remove_pattern_in_current_dir(GT_BIOSEQ_RAW, env);

  return 0;
}
