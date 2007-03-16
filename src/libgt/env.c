/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgt/cstr.h>
#include <libgt/env.h>
#include <libgt/option.h>
#include <libgt/splitter.h>
#include <libgt/versionfunc.h>
#include <libgt/xansi.h>

struct Env {
  MA *ma;
  Error *error;
  Log *log;
  bool spacepeak;
};

static OPrval parse_env_options(int argc, const char **argv, Env *env)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;
  assert(env);
  op = option_parser_new("GT_ENV_OPTIONS='[option ...]' ...",
                         "Parse the options contained in the "
                         "environment variable GT_ENV_OPTIONS.", env);
  o = option_new_bool("spacepeak", "show space peak on stdout upon deletion",
                      &env->spacepeak, false, env);
  option_parser_add_option(op, o, env);
  oprval = option_parser_parse_max_args(op, NULL, argc, argv, versionfunc, 0,
                                        env);
  option_parser_delete(op, env);
  return oprval;
}

static void proc_gt_env_options(Env *env)
{
  int argc;
  char *env_options, **argv;
  Splitter *splitter;
  assert(env);
  /* construct argument vector from $GT_ENV_OPTIONS */
  env_options = getenv("GT_ENV_OPTIONS");
  if (!env_options)
    return;
  env_options = cstr_dup(env_options, env); /* make writeable copy */
  splitter = splitter_new(env);
  splitter_split(splitter, env_options, strlen(env_options), ' ', env);
  argc = splitter_size(splitter);
  argv = cstr_array_preprend((const char**) splitter_get_tokens(splitter),
                             "env", env);
  argc++;
  /* parse options contained in $GT_ENV_OPTIONS */
  switch (parse_env_options(argc, (const char**) argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      fprintf(stderr, "error parsing $GT_ENV_OPTIONS: %s\n",
              env_error_get(env));
      env_error_unset(env);
      break;
    case OPTIONPARSER_REQUESTS_EXIT: break;
  }
  env_ma_free(env_options, env);
  splitter_delete(splitter, env);
  cstr_array_delete(argv, env);
}

Env* env_new(void)
{
  Env *env = xcalloc(1, sizeof (Env));
  env->ma = ma_new();
  ma_init(env->ma, env);
  env->error = error_new(env->ma);
  proc_gt_env_options(env);
  return env;
}

MA* env_ma(const Env *env)
{
  assert(env && env->ma);
  return env->ma;
}

Error* env_error(const Env *env)
{
  assert(env);
  return env->error;
}

Log* env_log(const Env *env)
{
  assert(env);
  return env->log;
}

void env_set_log(Env *env, Log *log)
{
  assert(env);
  env->log = log;
}

int env_delete(Env *env)
{
  int rval;
  assert(env);
  log_delete(env->log, env->ma);
  error_delete(env->error, env->ma);
  env->error = NULL;
  if (env->spacepeak)
    ma_show_space_peak(env->ma, stdout);
  rval = ma_check_space_leak(env->ma, env);
  ma_clean(env->ma, env);
  ma_delete(env->ma);
  free(env);
  return rval;
}

void env_ma_free_func(void *ptr, Env *env)
{
  assert(env);
  if (!ptr) return;
  ma_free(ptr, env_ma(env));
}

void env_error_set(Env *env, const char *format, ...)
{
  va_list ap;
  assert(env && format);
  va_start(ap, format);
  error_vset(env_error(env), format, ap);
  va_end(ap);
}

void env_log_log(Env *env, const char *format, ...)
{
  va_list ap;
  assert(env && format);
  if (!env_log(env)) return;
  va_start(ap, format);
  log_vlog(env_log(env), format, ap);
  va_end(ap);
}
