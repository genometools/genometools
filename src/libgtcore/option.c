/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include "libgtcore/array.h"
#include "libgtcore/cstr.h"
#include "libgtcore/error.h"
#include "libgtcore/mailaddress.h"
#include "libgtcore/minmax.h"
#include "libgtcore/option.h"
#include "libgtcore/str.h"
#include "libgtcore/undef.h"
#include "libgtcore/xansi.h"

typedef enum {
  OPTION_BOOL,
  OPTION_CHOICE,
  OPTION_DOUBLE,
  OPTION_HELP,
  OPTION_HELPPLUS,
  OPTION_HELPDEV,
  OPTION_OUTPUTFILE,
  OPTION_INT,
  OPTION_UINT,
  OPTION_LONG,
  OPTION_ULONG,
  OPTION_STRING,
  OPTION_STRINGARRAY,
  OPTION_VERSION,
} OptionType;

typedef struct {
  OptionParserHookFunc hook;
  void *data;
} HookInfo;

struct OptionParser {
  char *progname,
       *synopsis,
       *one_liner;
  Array *options,
        *hooks;
  bool parser_called;
  ShowCommentFunc comment_func;
  void *comment_func_data;
  const char *mailaddress;
};

struct Option {
  OptionType option_type;
  Str *option_str,
      *description;
  void *value;
  union {
    bool b;
    double d;
    FILE *fp;
    int i;
    unsigned int ui;
    long l;
    unsigned long ul;
    const char *s;
  } default_value;
  const char** domain;
  union {
    double d;
    int i;
    unsigned int ui;
    unsigned long ul;
  } min_value;
  union {
    double d;
    int i;
    unsigned int ui;
    unsigned long ul;
  } max_value;
  bool is_set,
       is_mandatory,
       hide_default,
       min_value_set,
       max_value_set,
       is_extended_option,
       is_development_option,
       argument_is_optional;
  Array *implications, /* contains option arrays, from each array at least one
                          option needs to be set */
        *exclusions;
  const Option *mandatory_either_option;
};

static Option *option_new(const char *option_str, const char *description,
                          void *value, Env *env)
{
  Option *o = env_ma_calloc(env, 1, sizeof (Option));
  assert(option_str && strlen(option_str));
  assert("an option string should not start with '-', this is added "
         "automatically"  && option_str[0] != '-');
  o->option_str = str_new_cstr(option_str, env);
  o->description = str_new_cstr(description, env);
  o->value = value;
  return o;
}

static Option* option_new_help(bool has_extended_options, Env *env)
{
  Option *o;
  if (has_extended_options) {
    o = option_new("help", "display help for basic options and exit", NULL,
                   env);
  }
  else
    o = option_new("help", "display help and exit", NULL, env);
  o->option_type = OPTION_HELP;
  o->default_value.b = false;
  return o;
}

static Option* option_new_helpplus(Env *env)
{
  Option *o = option_new("help+", "display help for all options and exit", NULL,
                         env);
  o->option_type = OPTION_HELPPLUS;
  o->default_value.b = false;
  return o;
}

static Option* option_new_helpdev(Env *env)
{
  Option *o = option_new("helpdev", "display help for development options and "
                         "exit", NULL, env);
  o->option_type = OPTION_HELPDEV;
  o->default_value.b = false;
  o->is_development_option = true;
  return o;
}

static Option* option_new_version(ShowVersionFunc versionfunc, Env *env)
{
  Option *o = option_new("version", "display version information and exit",
                         versionfunc, env);
  o->option_type = OPTION_VERSION;
  return o;
}

OptionParser* option_parser_new(const char *synopsis, const char *one_liner,
                                Env *env)
{
  OptionParser *op = env_ma_malloc(env, sizeof (OptionParser));
  assert(synopsis && one_liner);
  assert("one_liner must have upper case letter at start and '.' at end" &&
         strlen(one_liner) && isupper((int) one_liner[0]));
  assert(one_liner[strlen(one_liner)-1] == '.');
  op->progname = NULL;
  op->synopsis = cstr_dup(synopsis, env);
  op->one_liner = cstr_dup(one_liner, env);
  op->options = array_new(sizeof (Option*), env);
  op->hooks = NULL;
  op->parser_called = false;
  op->comment_func = NULL;
  op->comment_func_data = NULL;
  op->mailaddress = NULL;
  return op;
}

void option_parser_add_option(OptionParser *op, Option *o, Env *env)
{
  assert(op && o);
  array_add(op->options, o, env);
}

void option_parser_set_comment_func(OptionParser *op,
                                    ShowCommentFunc comment_func, void *data)
{
  assert(op);
  op->comment_func = comment_func;
  op->comment_func_data = data;
}

void option_parser_register_hook(OptionParser *op, OptionParserHookFunc hook,
                                 void *data, Env *env)
{
  HookInfo hookinfo;
  env_error_check(env);
  assert(op && hook);
  if (!op->hooks)
    op->hooks = array_new(sizeof (HookInfo), env);
  hookinfo.hook = hook;
  hookinfo.data = data;
  array_add(op->hooks, hookinfo, env);
}

void option_parser_set_mailaddress(OptionParser *op, const char *address)
{
  assert(op && address);
  op->mailaddress = address;
}

static void show_description(unsigned long initial_space, const char *desc,
                             unsigned long len)
{
  const unsigned long width = TERMINAL_WIDTH - initial_space;
  const char *tmp_ptr, *desc_ptr = desc;
  unsigned long i;
  bool continue_while = false;

  /* got space to show option */
  assert(initial_space < TERMINAL_WIDTH);

  while (desc_ptr < desc + len) {
    /* break, if the rest of the description fits on one line */
    if (desc_ptr + width - 1 >= desc + len - 1)
      break;
    /* go backwards to find a point to break description */
    for (tmp_ptr = desc_ptr + width; tmp_ptr >= desc_ptr; tmp_ptr--) {
      if (*tmp_ptr == ' ' || *tmp_ptr == '\n')
        break;
    }
    /* break point found, show description up to that point */
    for (; desc_ptr < tmp_ptr; desc_ptr++) {
      xputchar(*desc_ptr);
      if (*desc_ptr == '\n') {
        /* show leading spaces */
        for  (i = 0; i < initial_space; i++)
          xputchar(' ');
        desc_ptr++;
        continue_while = true;
        break;
      }
    }
    if (continue_while) {
      continue_while = false;
      continue;
    }
    /* we are at the break point now */
    assert(*desc_ptr == ' ' || *desc_ptr == '\n');
    /* show newline for break point */
    desc_ptr++;
    xputchar('\n');
    /* show leading spaces */
    for  (i = 0; i < initial_space; i++)
      xputchar(' ');
  }
  /* show final line */
  while (desc_ptr < desc + len) {
    xputchar(*desc_ptr);
    if (*desc_ptr == '\n') {
      /* show leading spaces */
      for  (i = 0; i < initial_space; i++)
        xputchar(' ');
      desc_ptr++;
      continue;
    }
    desc_ptr++;
  }
  xputchar('\n');
}

static int show_help(OptionParser *op, OptionType optiontype, Env *env)
{
  unsigned long i, max_option_length = 0;
  Option *option;
  int had_err = 0;
  env_error_check(env);
  assert(optiontype == OPTION_HELP || optiontype == OPTION_HELPPLUS ||
         optiontype == OPTION_HELPDEV);

  /* determine maximum option length */
  for (i = 0; i < array_size(op->options); i++) {
    option = *(Option**) array_get(op->options, i);
    /* skip option if necessary */
    if ((optiontype == OPTION_HELP && option->is_extended_option) ||
        (optiontype == OPTION_HELPDEV && !option->is_development_option) ||
        (!(optiontype == OPTION_HELPDEV) && option->is_development_option)) {
      continue;
    }
    if (str_length(option->option_str) > max_option_length)
      max_option_length = str_length(option->option_str);
  }
  assert(max_option_length);

  printf("Usage: %s %s\n", op->progname, op->synopsis);
  printf("%s\n\n", op->one_liner);
  for (i = 0; i < array_size(op->options); i++) {
    option = *(Option**) array_get(op->options, i);
    /* skip option if necessary */
    if ((optiontype == OPTION_HELP && option->is_extended_option) ||
        (optiontype == OPTION_HELPDEV && !option->is_development_option) ||
        (!(optiontype == OPTION_HELPDEV) && option->is_development_option)) {
      continue;
    }
    printf("-%s%*s ", str_get(option->option_str),
           (int) (max_option_length - str_length(option->option_str)), "");
    show_description(max_option_length + 2, str_get(option->description),
                     str_length(option->description));

    /* show default value for some types of options */
    if (!option->hide_default) {
      if (option->option_type == OPTION_BOOL) {
        printf("%*s  default: %s\n", (int) max_option_length, "",
               option->default_value.b ? "yes" : "no");
      }
      else if (option->option_type == OPTION_CHOICE) {
        printf("%*s  default: ", (int) max_option_length, "");
        if (!option->default_value.s || !strlen(option->default_value.s))
          xputs("undefined");
        else
          xputs(option->default_value.s);
      }

      else if (option->option_type == OPTION_DOUBLE) {
        printf("%*s  default: ", (int) max_option_length, "");
        if (option->default_value.d == UNDEF_DOUBLE)
          xputs("undefined");
        else
          printf("%.2f\n", option->default_value.d);
      }
      else if (option->option_type == OPTION_INT) {
        printf("%*s  default: ", (int) max_option_length, "");
        if (option->default_value.i == UNDEF_INT)
          xputs("undefined");
        else
          printf("%d\n", option->default_value.i);
      }
      else if (option->option_type == OPTION_UINT) {
        printf("%*s  default: ", (int) max_option_length, "");
        if (option->default_value.ui == UNDEF_UINT)
          xputs("undefined");
        else
          printf("%u\n", option->default_value.ui);
      }
      else if (option->option_type == OPTION_LONG) {
        printf("%*s  default: ", (int) max_option_length, "");
        if (option->default_value.ul == UNDEF_LONG)
          xputs("undefined");
        else
          printf("%ld\n", option->default_value.l);
      }
      else if (option->option_type == OPTION_ULONG) {
        printf("%*s  default: ", (int) max_option_length, "");
        if (option->default_value.ul == UNDEF_ULONG)
          xputs("undefined");
        else
          printf("%lu\n", option->default_value.ul);
      }
      else if (option->option_type == OPTION_STRING) {
        printf("%*s  default: ", (int) max_option_length, "");
        if (!option->default_value.s || !strlen(option->default_value.s))
          xputs("undefined");
        else
          xputs(option->default_value.s);
      }
    }
  }
  if (op->comment_func)
    had_err = op->comment_func(op->progname, op->comment_func_data, env);
  if (!had_err) {
    printf("\nReport bugs to %s.\n",
           op->mailaddress ? op->mailaddress : MAILADDRESS);
  }
  return had_err;
}

static bool optional_arg(Option *o, int argnum, int argc, const char **argv)
{
  assert(o);
  if (o->argument_is_optional &&
      (argnum + 1 >= argc || argv[argnum + 1][0] == '-' ||
       !strcmp(argv[argnum + 1], "--"))) {
    return true;
  }
  return false;
}

static int check_missing_argument(int argnum, int argc, Str *option, Env *env)
{
  env_error_check(env);
  if (argnum + 1 >= argc) {
    env_error_set(env, "missing argument to option \"-%s\"", str_get(option));
    return -1;
  }
  return 0;
}

static int check_mandatory_options(OptionParser *op, Env *env)
{
  unsigned long i;
  Option *o;
  env_error_check(env);
  assert(op);
  for (i = 0; i < array_size(op->options); i++) {
    o = *(Option**) array_get(op->options, i);
    if (o->is_mandatory && !o->is_set) {
      env_error_set(env, "option \"-%s\" is mandatory", str_get(o->option_str));
      return -1;
    }
  }
  return 0;
}

static int check_option_implications(OptionParser *op, Env *env)
{
  unsigned long i, j, k, l;
  Array *implied_option_array;
  Option *o, *implied_option;
  unsigned int option_set;
  Str *error_str;
  env_error_check(env);

  for (i = 0; i < array_size(op->options); i++) {
    o = *(Option**) array_get(op->options, i);
    if (o->implications && o->is_set) {
      for (j = 0; j < array_size(o->implications); j++) {
        implied_option_array = *(Array**) array_get(o->implications, j);
        assert(array_size(implied_option_array));
        if (array_size(implied_option_array) == 1) {
          /* special case: option implies exactly one option */
          implied_option = *(Option**) array_get(implied_option_array, 0);
          if (!implied_option->is_set) {
            env_error_set(env, "option \"-%s\" requires option \"-%s\"",
                          str_get(o->option_str),
                          str_get(implied_option->option_str));
            return -1;
          }
        }
        else {
          /* ``either'' case: option implied at least one of the options given
             in array */
          option_set = 0;
          for (k = 0; k < array_size(implied_option_array); k++) {
            implied_option = *(Option**) array_get(implied_option_array, k);
            if (implied_option->is_set) {
              option_set = 1;
              break;
            }
          }
          if (!option_set) {
            error_str = str_new_cstr("option \"-", env);
            str_append_str(error_str, o->option_str, env);
            str_append_cstr(error_str, "\" requires option", env);
            for (l = 0; l < array_size(implied_option_array) - 1; l++) {
              str_append_cstr(error_str, " \"-", env);
              str_append_str(error_str, (*(Option**)
                             array_get(implied_option_array, l))->option_str,
                             env);
              str_append_cstr(error_str, "\"", env);
              if (array_size(implied_option_array) > 2)
                str_append_char(error_str, ',', env);
            }
            str_append_cstr(error_str, " or \"-", env);
            str_append_str(error_str, (*(Option**)
                           array_get(implied_option_array,
                                     array_size(implied_option_array) - 1))
                                     ->option_str, env);
            str_append_cstr(error_str, "\"", env);
            env_error_set(env, "%s", str_get(error_str));
            str_delete(error_str, env);
            return -1;
          }
        }
      }
    }
  }
  return 0;
}

static int check_option_exclusions(OptionParser *op, Env *env)
{
  unsigned long i, j;
  Option *o, *excluded_option;
  env_error_check(env);

  for (i = 0; i < array_size(op->options); i++) {
    o = *(Option**) array_get(op->options, i);
    if (o->exclusions && o->is_set) {
      for (j = 0; j < array_size(o->exclusions); j++) {
        excluded_option = *(Option**) array_get(o->exclusions, j);
        if (excluded_option->is_set) {
          env_error_set(env,
                        "option \"-%s\" and option \"-%s\" exclude each other",
                        str_get(o->option_str),
                        str_get(excluded_option->option_str));
          return -1;
        }
      }
    }
  }
  return 0;
}

static int check_mandatory_either_options(OptionParser *op, Env *env)
{
  unsigned long i;
  Option *o;
  env_error_check(env);

  for (i = 0; i < array_size(op->options); i++) {
    o = *(Option**) array_get(op->options, i);
    if (o->mandatory_either_option) {
      if (!o->is_set && !o->mandatory_either_option->is_set) {
        env_error_set(env,
                      "either option \"-%s\" or option \"-%s\" is mandatory",
                      str_get(o->option_str),
                      str_get(o->mandatory_either_option->option_str));
        return -1;
      }
    }
  }
  return 0;
}

static bool has_extended_option(Array *options)
{
  unsigned long i;
  Option *option;
  assert(options);
  for (i = 0; i < array_size(options); i++) {
    option = *(Option**) array_get(options, i);
    if (option->is_extended_option)
      return true;
  }
  return false;
}

static OPrval parse(OptionParser *op, int *parsed_args, int argc,
                    const char **argv,
                    ShowVersionFunc versionfunc,
                    unsigned int min_additional_arguments,
                    unsigned int max_additional_arguments, Env *env)
{
  int argnum, int_value;
  unsigned long i;
  double double_value;
  HookInfo *hookinfo;
  Option *option;
  bool has_extended_options, option_parsed;
  long long_value;
  int had_err = 0;
  Str *error_str;

  env_error_check(env);
  assert(op);
  assert(!op->parser_called); /* to avoid multiple adding of common options */

  op->progname = cstr_dup(argv[0], env);

  /* add common options */
  has_extended_options = has_extended_option(op->options);
  option = option_new_help(has_extended_options, env);
  option_parser_add_option(op, option, env);
  if (has_extended_options) {
    option = option_new_helpplus(env);
    option_parser_add_option(op, option, env);
  }
  option = option_new_helpdev(env);
  option_parser_add_option(op, option, env);
  option = option_new_version(versionfunc, env);
  option_parser_add_option(op, option, env);

  for (argnum = 1; argnum < argc; argnum++) {
    if (!(argv[argnum] && argv[argnum][0] == '-' && strlen(argv[argnum]) > 1) ||
        !strcmp(argv[argnum], "--")) {
      break;
    }

    /* look for matching option */
    option_parsed = false;
    for (i = 0; i < array_size(op->options); i++) {
      option = *(Option**) array_get(op->options, i);

      if (strcmp(argv[argnum]+1, str_get(option->option_str)) == 0) {
        /* make sure option has not been used before */
        if (option->is_set) {
          env_error_set(env, "option \"%s\" already set",
                        str_get(option->option_str));
          had_err = -1;
        }
        else
          option->is_set = true;
        if (!had_err) {
          switch (option->option_type) {
            case OPTION_BOOL:
              if (argv[argnum+1] && argv[argnum+1][0] != '-') {
                if (!strcmp(argv[argnum+1], "yes") ||
                    !strcmp(argv[argnum+1], "true")) {
                  argnum++;
                  *(bool*) option->value = true;
                  option_parsed = true;
                  break;
                }
                else if (!strcmp(argv[argnum+1], "no") ||
                         !strcmp(argv[argnum+1], "false")) {
                  argnum++;
                  *(bool*) option->value = false;
                  option_parsed = true;
                  break;
                }
              }
              *(bool*) option->value = true;
              option_parsed = true;
              break;
            case OPTION_CHOICE:
              if (optional_arg(option, argnum, argc, argv)) {
                option_parsed = true;
                break;
              }
              assert(option->domain[0]);
              had_err = check_missing_argument(argnum, argc, option->option_str,
                        env);
              if (!had_err) {
                argnum++;
                if (strcmp(argv[argnum], option->domain[0])) {
                  error_str = str_new_cstr(option->domain[0], env);
                  i = 1;
                  while (option->domain[i] != NULL) {
                    if (!strcmp(argv[argnum], option->domain[i])) {
                      str_set(option->value, option->domain[i], env);
                      break;
                    }
                    str_append_cstr(error_str, ", ", env);
                    str_append_cstr(error_str, option->domain[i], env);
                    i++;
                  }
                  if (option->domain[i] == NULL) {
                    env_error_set(env, "argument to option \"-%s\" must be one "
                                  "of: %s", str_get(option->option_str),
                                  str_get(error_str));
                    had_err = -1;
                  }
                  str_delete(error_str, env);
                }
                else {
                  str_set(option->value, option->domain[0], env);
                }
              }
              if (!had_err) {
                option_parsed = true;
              }
              break;
            case OPTION_DOUBLE:
              if (optional_arg(option, argnum, argc, argv)) {
                option_parsed = true;
                break;
              }
              had_err = check_missing_argument(argnum, argc, option->option_str,
                                               env);
              if (!had_err) {
                argnum++;
                if (sscanf(argv[argnum], "%lf", &double_value) != 1) {
                  env_error_set(env, "argument to option \"-%s\" must be "
                                     "floating-point number",
                                str_get(option->option_str));
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* minimum value check */
                if (option->min_value_set &&
                    double_value < option->min_value.d) {
                  env_error_set(env, "argument to option \"-%s\" must be a "
                                     "floating point value >= %f",
                                str_get(option->option_str),
                                option->min_value.d);
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* maximum value check */
                if (option->max_value_set &&
                    double_value > option->max_value.d) {
                  env_error_set(env, "argument to option \"-%s\" must be a "
                                     "floating point value <= %f",
                                str_get(option->option_str),
                                option->max_value.d);
                  had_err = -1;
                }
              }
              if (!had_err) {
                *(double*) option->value = double_value;
                option_parsed = true;
              }
              break;
            case OPTION_HELP:
              if (show_help(op, OPTION_HELP, env))
                return OPTIONPARSER_ERROR;
              return OPTIONPARSER_REQUESTS_EXIT;
            case OPTION_HELPPLUS:
              if (show_help(op, OPTION_HELPPLUS, env))
                return OPTIONPARSER_ERROR;
              return OPTIONPARSER_REQUESTS_EXIT;
            case OPTION_HELPDEV:
              if (show_help(op, OPTION_HELPDEV, env))
                return OPTIONPARSER_ERROR;
              return OPTIONPARSER_REQUESTS_EXIT;
            case OPTION_OUTPUTFILE:
              had_err = check_missing_argument(argnum, argc, option->option_str,
                                               env);
              if (!had_err) {
                argnum++;
                *(FILE**) option->value = env_fa_xfopen(env, argv[argnum], "w");
                option_parsed = true;
              }
              break;
            case OPTION_INT:
              assert(!option->argument_is_optional);
              had_err = check_missing_argument(argnum, argc, option->option_str,
                                               env);
              if (!had_err) {
                argnum++;
                if (sscanf(argv[argnum], "%d", &int_value) != 1) {
                  env_error_set(env, "argument to option \"-%s\" must be an "
                                "integer", str_get(option->option_str));
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* minimum value check */
                if (option->min_value_set && int_value < option->min_value.i) {
                  env_error_set(env, "argument to option \"-%s\" must be an "
                                "integer >= %d", str_get(option->option_str),
                                option->min_value.i);
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* maximum value check */
                if (option->max_value_set && int_value > option->max_value.i) {
                  env_error_set(env, "argument to option \"-%s\" must be an "
                                "integer <= %d", str_get(option->option_str),
                                option->max_value.i);
                  had_err = -1;
                }
              }
              if (!had_err) {
                *(int*) option->value = int_value;
                option_parsed = true;
              }
              break;
            case OPTION_UINT:
              if (optional_arg(option, argnum, argc, argv)) {
                option_parsed = true;
                break;
              }
              had_err = check_missing_argument(argnum, argc, option->option_str,
                                               env);
              if (!had_err) {
                argnum++;
                if (sscanf(argv[argnum], "%d", &int_value) != 1 ||
                    int_value < 0) {
                  env_error_set(env, "argument to option \"-%s\" must be a "
                                "non-negative integer",
                                str_get(option->option_str));
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* minimum value check */
                if (option->min_value_set && int_value < option->min_value.ui) {
                  env_error_set(env, "argument to option \"-%s\" must be an "
                                "integer >= %u", str_get(option->option_str),
                                option->min_value.ui);
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* maximum value check */
                if (option->max_value_set && int_value > option->max_value.ui) {
                  env_error_set(env, "argument to option \"-%s\" must be an "
                                "integer <= %u", str_get(option->option_str),
                                option->max_value.ui);
                  had_err = -1;
                }
              }
              if (!had_err) {
                *(unsigned int*) option->value = int_value;
                option_parsed = true;
              }
              break;
            case OPTION_LONG:
              assert(!option->argument_is_optional);
              had_err = check_missing_argument(argnum, argc, option->option_str,
                        env);
              if (!had_err) {
                argnum++;
                if (sscanf(argv[argnum], "%ld", &long_value) != 1) {
                  env_error_set(env, "argument to option \"-%s\" must be an "
                                "integer", str_get(option->option_str));
                  had_err = -1;
                }
              }
              if (!had_err) {
                *(long*) option->value = long_value;
                option_parsed = true;
              }
              break;
            case OPTION_ULONG:
              if (optional_arg(option, argnum, argc, argv)) {
                option_parsed = true;
                break;
              }
              had_err = check_missing_argument(argnum, argc, option->option_str,
                                               env);
              if (!had_err) {
                argnum++;
                if (sscanf(argv[argnum], "%ld", &long_value) != 1 ||
                    long_value < 0) {
                  env_error_set(env, "argument to option \"-%s\" must be a "
                                "non-negative integer",
                                str_get(option->option_str));
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* minimum value check */
                if (option->min_value_set &&
                    long_value < option->min_value.ul) {
                  env_error_set(env, "argument to option \"-%s\" must be an "
                                "integer >= %lu", str_get(option->option_str),
                                option->min_value.ul);
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* maximum value check */
                if (option->max_value_set &&
                    long_value > option->max_value.ul) {
                  env_error_set(env, "argument to option \"-%s\" must be an "
                                "integer <= %lu", str_get(option->option_str),
                                option->max_value.ul);
                  had_err = -1;
                }
              }
              if (!had_err) {
                *(unsigned long*) option->value = long_value;
                option_parsed = true;
              }
              break;
            case OPTION_STRING:
              if (optional_arg(option, argnum, argc, argv)) {
                option_parsed = true;
                break;
              }
              had_err = check_missing_argument(argnum, argc, option->option_str,
                        env);
              if (!had_err) {
                argnum++;
                str_set(option->value, argv[argnum], env);
                option_parsed = true;
              }
              break;
            case OPTION_STRINGARRAY:
              if (optional_arg(option, argnum, argc, argv)) {
                option_parsed = true;
                break;
              }
              had_err = check_missing_argument(argnum, argc, option->option_str,
                                               env);
              while (!had_err) {
                if (argnum + 1 < argc && argv[argnum+1][0] != '-') {
                  argnum++;
                  strarray_add_cstr(option->value, argv[argnum], env);
                  option_parsed = true;
                }
                else
                  break;
              }
              break;
            case OPTION_VERSION:
              ((ShowVersionFunc) option->value)(op->progname);
              return OPTIONPARSER_REQUESTS_EXIT;
            default: assert(0);
          }
        }
        if (had_err || option_parsed)
          break;
      }
    }

    if (had_err)
      break;
    if (option_parsed)
      continue;

    /* no matching option found -> error */
    assert(!had_err);
    env_error_set(env, "unknown option: %s (-help shows possible options)",
                  argv[argnum]);
    had_err = -1;
    break;
  }

  /* skip "--" if necessary */
  if (argnum < argc && !strcmp(argv[argnum], "--"))
    argnum++;

  /* check for minimum number of additional arguments, if necessary */
  if (!had_err && min_additional_arguments != UNDEF_UINT &&
      argc - argnum < min_additional_arguments) {
    env_error_set(env, "missing argument\nUsage: %s %s", op->progname,
                  op->synopsis);
    had_err = -1;
  }

  /* check for maximal number of additional arguments, if necessary */
  if (!had_err && max_additional_arguments != UNDEF_UINT &&
      argc - argnum > max_additional_arguments) {
    env_error_set(env, "superfluous argument \"%s\"\nUsage: %s %s",
                  argv[argnum + max_additional_arguments], op->progname,
                  op->synopsis);
    had_err = -1;
  }

  if (!had_err)
    had_err = check_mandatory_options(op, env);
  if (!had_err)
    had_err = check_option_implications(op, env);
  if (!had_err)
    had_err = check_option_exclusions(op, env);
  if (!had_err)
    had_err = check_mandatory_either_options(op, env);

  /* call hooks */
  for (i = 0; !had_err && i < array_size(op->hooks); i++) {
    hookinfo = array_get(op->hooks, i);
    had_err = hookinfo->hook(hookinfo->data, env);
  }

  op->parser_called = true;
  if (parsed_args)
    *parsed_args = argnum;

  if (had_err)
    return OPTIONPARSER_ERROR;
  return OPTIONPARSER_OK;
}

OPrval option_parser_parse(OptionParser *op, int *parsed_args, int argc,
                           const char **argv, ShowVersionFunc versionfunc,
                           Env *env)
{
  env_error_check(env);
  return parse(op, parsed_args, argc, argv, versionfunc, UNDEF_UINT, UNDEF_UINT,
               env);
}

OPrval option_parser_parse_min_args(OptionParser *op, int *parsed_args,
                                    int argc, const char **argv,
                                    ShowVersionFunc versionfunc,
                                    unsigned int min_additional_arguments,
                                    Env *env)
{
  env_error_check(env);
  return parse(op, parsed_args, argc, argv, versionfunc,
               min_additional_arguments, UNDEF_UINT, env);
}

OPrval option_parser_parse_max_args(OptionParser *op, int *parsed_args,
                                    int argc, const char **argv,
                                    ShowVersionFunc versionfunc,
                                    unsigned int max_additional_arguments,
                                    Env *env)
{
  env_error_check(env);
  return parse(op, parsed_args, argc, argv, versionfunc, UNDEF_UINT,
               max_additional_arguments, env);
}

OPrval option_parser_parse_min_max_args(OptionParser *op, int *parsed_args,
                                        int argc, const char **argv,
                                        ShowVersionFunc versionfunc,
                                        unsigned int min_additional_arguments,
                                        unsigned int max_additional_arguments,
                                        Env *env)
{
  env_error_check(env);
  return parse(op, parsed_args, argc, argv, versionfunc,
               min_additional_arguments, max_additional_arguments, env);
}

void option_parser_delete(OptionParser *op, Env *env)
{
  unsigned long i;
  if (!op) return;
  env_ma_free(op->progname, env);
  env_ma_free(op->synopsis, env);
  env_ma_free(op->one_liner, env);
  for (i = 0; i < array_size(op->options); i++)
    option_delete(*(Option**) array_get(op->options, i), env);
  array_delete(op->options, env);
  array_delete(op->hooks, env);
  env_ma_free(op, env);
}

Option* option_new_outputfile(FILE **outfp, Env *env)
{
  Option *o = option_new("o", "redirect output to specified file (will "
                         "overwrite existing file!)", outfp, env);
  o->option_type = OPTION_OUTPUTFILE;
  o->default_value.fp = stdout;
  *outfp = stdout;
  return o;
}

Option* option_new_verbose(bool *value, Env *env)
{
  return option_new_bool("v", "be verbose", value, false, env);
}

Option* option_new_debug(bool *value, Env *env)
{
  Option *o = option_new_bool("debug", "enable debugging output", value, false,
                              env);
  o->is_development_option = true;
  return o;
}

Option* option_new_bool(const char *option_str, const char *description,
                        bool *value, bool default_value, Env *env)
{
  Option *o = option_new(option_str, description, value, env);
  o->option_type = OPTION_BOOL;
  o->default_value.b = default_value;
  *value = default_value;
  return o;
}

Option* option_new_double(const char *option_str, const char *description,
                          double *value, double default_value, Env *env)
{
  Option *o = option_new(option_str, description, value, env);
  o->option_type = OPTION_DOUBLE;
  o->default_value.d = default_value;
  *value = default_value;
  return o;
}

Option *option_new_double_min(const char *option_str, const char *description,
                              double *value, double default_value,
                              double min_value, Env *env)
{
  Option *o = option_new_double(option_str, description, value, default_value,
                                env);
  o->min_value_set = true;
  o->min_value.d = min_value;
  return o;
}

Option *option_new_double_min_max(const char *option_str,
                                  const char *description, double *value,
                                  double default_value, double min_value,
                                  double max_value, Env *env)
{
  Option *o = option_new_double(option_str, description, value, default_value,
                                env);
  o->min_value_set = true;
  o->min_value.d = min_value;
  o->max_value_set = true;
  o->max_value.d = max_value;
  return o;
}

Option* option_new_probability(const char *option_str, const char *description,
                               double *value, double default_value, Env *env)
{
  return option_new_double_min_max(option_str, description, value,
                                   default_value, 0.0, 1.0, env);
}

Option* option_new_int(const char *option_str, const char *description,
                       int *value, int default_value, Env *env)
{
  Option *o = option_new(option_str, description, value, env);
  o->option_type = OPTION_INT;
  o->default_value.i = default_value;
  *value = default_value;
  return o;
}

Option* option_new_int_min(const char *option_str, const char *description,
                           int *value, int default_value, int min_value,
                           Env *env)
{
  Option *o = option_new_int(option_str, description, value, default_value,
                             env);
  o->min_value_set = true;
  o->min_value.i = min_value;
  return o;
}

Option* option_new_int_max(const char *option_str, const char *description,
                           int *value, int default_value, int max_value,
                           Env *env)
{
  Option *o = option_new_int(option_str, description, value, default_value,
                              env);
  o->max_value_set = true;
  o->max_value.i = max_value;
  return o;
}

Option* option_new_uint(const char *option_str, const char *description,
                        unsigned int *value, unsigned int default_value,
                        Env *env)
{
  Option *o = option_new(option_str, description, value, env);
  o->option_type = OPTION_UINT;
  o->default_value.ui = default_value;
  *value = default_value;
  return o;
}

Option* option_new_uint_min(const char *option_str, const char *description,
                            unsigned int *value, unsigned int default_value,
                            unsigned int min_value, Env *env)
{
  Option *o = option_new_uint(option_str, description, value, default_value,
                              env);
  o->min_value_set = true;
  o->min_value.ui = min_value;
  return o;
}

Option* option_new_uint_max(const char *option_str, const char *description,
                            unsigned int *value, unsigned int default_value,
                            unsigned int max_value, Env *env)
{
  Option *o = option_new_uint(option_str, description, value, default_value,
                              env);
  o->max_value_set = true;
  o->max_value.ui = max_value;
  return o;
}

Option *option_new_uint_min_max(const char *option_str, const char *description,
                                unsigned int *value, unsigned int default_value,
                                unsigned int min_value, unsigned int max_value,
                                Env *env)
{
  Option *o = option_new_uint(option_str, description, value, default_value,
                               env);
  o->min_value_set = true;
  o->min_value.i = min_value;
  o->max_value_set = true;
  o->max_value.i = max_value;
  return o;
}

Option* option_new_long(const char *option_str, const char *description,
                        long *value, long default_value, Env *env)
{
  Option *o = option_new(option_str, description, value, env);
  o->option_type = OPTION_LONG;
  o->default_value.l = default_value;
  *value = default_value;
  return o;
}

Option* option_new_ulong(const char *option_str, const char *description,
                         unsigned long *value, unsigned long default_value,
                         Env *env)
{
  Option *o = option_new(option_str, description, value, env);
  o->option_type = OPTION_ULONG;
  o->default_value.ul = default_value;
  *value = default_value;
  return o;
}

Option* option_new_ulong_min(const char *option_str, const char *description,
                             unsigned long *value, unsigned long default_value,
                             unsigned long min_value, Env *env)
{
  Option *o = option_new_ulong(option_str, description, value, default_value,
                               env);
  o->min_value_set = true;
  o->min_value.ul = min_value;
  return o;
}

Option *option_new_ulong_min_max(const char *option_str,
                                 const char *description, unsigned long *value,
                                 unsigned long default_value,
                                 unsigned long min_value,
                                 unsigned long max_value, Env *env)
{
  Option *o = option_new_ulong(option_str, description, value, default_value,
                               env);
  o->min_value_set = true;
  o->min_value.ul = min_value;
  o->max_value_set = true;
  o->max_value.ul = max_value;
  return o;
}

Option* option_new_string(const char *option_str, const char *description,
                          Str *value, const char *default_value, Env *env)
{
  Option *o = option_new(option_str, description, value, env);
  o->option_type = OPTION_STRING;
  o->default_value.s = default_value;
  str_set(value, default_value, env);
  return o;
}

Option* option_new_stringarray(const char *option_str,
                               const char *description, StrArray *value,
                               Env *env)
{
  Option *o = option_new(option_str, description, value, env);
  o->option_type = OPTION_STRINGARRAY;
  return o;
}

/* the following function would allow to handle files differently from strings
   later on (e.g., for CGI scripts), but for now the are implemented in the same
   way */
Option* option_new_filename(const char *option_str, const char *description,
                            Str*filename, Env *env)
{
  return option_new_string(option_str, description, filename, NULL, env);
}

/* the following function would allow to handle file arrays differently from
   string arrays later on (e.g., for CGI scripts) , but for now the are
   implemented in the same way */
Option* option_new_filenamearray(const char *option_str,
                                 const char *description, StrArray *filenames,
                                 Env *env)
{
  return option_new_stringarray(option_str, description, filenames, env);
}

Option* option_new_choice(const char *option_str, const char *description,
                          Str *value, const char* default_value,
                          const char** domain, Env *env)
{
  Option *o;
#ifndef NDEBUG
  unsigned long in_domain = 1;
  if (default_value) {
    while (domain[(in_domain - 1)] != NULL) {
      if (domain[(in_domain - 1)] == default_value) {
        in_domain = 0;
        break;
      }
      in_domain++;
    }
  }
  else
    in_domain = 0;
  assert(!in_domain);
#endif

  o = option_new_string(option_str, description, value, default_value, env);
  o->option_type = OPTION_CHOICE;
  o->domain = domain;

  return o;
}

void option_is_mandatory(Option *o)
{
  assert(o);
  o->is_mandatory = true;
}

void option_is_mandatory_either(Option *o, const Option *meo)
{
  assert(o && meo);
  assert(!o->mandatory_either_option);
  o->mandatory_either_option = meo;
}

void option_is_extended_option(Option *o)
{
  assert(o);
  o->is_extended_option = true;
}

void option_is_development_option(Option *o)
{
  assert(o);
  o->is_development_option = true;
}

void option_imply(Option *o, const Option *implied_option, Env *env)
{
  Array *option_array;
  assert(o && implied_option);
  if (!o->implications)
    o->implications = array_new(sizeof (Array*), env);
  option_array = array_new(sizeof (Option*), env);
  array_add(option_array, implied_option, env);
  array_add(o->implications, option_array, env);
}

void option_imply_either_2(Option *o, const Option *io1, const Option *io2,
                           Env *env)
{
  Array *option_array;
  assert(o && io1 && io2);
  if (!o->implications)
    o->implications = array_new(sizeof (Array*), env);
  option_array = array_new(sizeof (Option*), env);
  array_add(option_array, io1, env);
  array_add(option_array, io2, env);
  array_add(o->implications, option_array, env);
}

void option_exclude(Option *o_a, Option *o_b, Env *env)
{
  assert(o_a && o_b);
  if (!o_a->exclusions)
    o_a->exclusions = array_new(sizeof (Option*), env);
  if (!o_b->exclusions)
    o_b->exclusions = array_new(sizeof (Option*), env);
  array_add(o_a->exclusions, o_b, env);
  array_add(o_b->exclusions, o_a, env);
}

void option_hide_default(Option *o)
{
  assert(o);
  o->hide_default = true;
}

void option_argument_is_optional(Option *o)
{
  assert(o);
  o->argument_is_optional = true;
}

bool option_is_set(const Option *o)
{
  assert(o);
  return o->is_set;
}

void option_delete(Option *o, Env *env)
{
  unsigned long i;
  if (!o) return;
  str_delete(o->option_str, env);
  str_delete(o->description, env);
  for (i = 0; i < array_size(o->implications); i++)
    array_delete(*(Array**) array_get(o->implications, i), env);
  array_delete(o->implications, env);
  array_delete(o->exclusions, env);
  env_ma_free(o, env);
}
