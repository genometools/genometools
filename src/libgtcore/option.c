/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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
#include <limits.h>
#include <stdio.h>
#include "libgtcore/array.h"
#include "libgtcore/cstr.h"
#include "libgtcore/error.h"
#include "libgtcore/fa.h"
#include "libgtcore/ma.h"
#include "libgtcore/mailaddress.h"
#include "libgtcore/minmax.h"
#include "libgtcore/option.h"
#include "libgtcore/parseutils.h"
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
  OPTION_RANGE,
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
  unsigned int min_additional_arguments,
               max_additional_arguments;
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
    Range r;
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
  unsigned int reference_count;
};

static Option *option_new(const char *option_str, const char *description,
                          void *value)
{
  Option *o = ma_calloc(1, sizeof (Option));
  assert(option_str && strlen(option_str));
  assert("an option string should not start with '-', this is added "
         "automatically"  && option_str[0] != '-');
  o->option_str = str_new_cstr(option_str);
  o->description = str_new_cstr(description);
  o->value = value;
  return o;
}

static Option* option_new_help(bool has_extended_options)
{
  Option *o;
  if (has_extended_options)
    o = option_new("help", "display help for basic options and exit", NULL);
  else
    o = option_new("help", "display help and exit", NULL);
  o->option_type = OPTION_HELP;
  o->default_value.b = false;
  return o;
}

static Option* option_new_helpplus(void)
{
  Option *o = option_new("help+", "display help for all options and exit",
                         NULL);
  o->option_type = OPTION_HELPPLUS;
  o->default_value.b = false;
  return o;
}

static Option* option_new_helpdev(void)
{
  Option *o = option_new("helpdev", "display help for development options and "
                         "exit", NULL);
  o->option_type = OPTION_HELPDEV;
  o->default_value.b = false;
  o->is_development_option = true;
  return o;
}

static Option* option_new_version(ShowVersionFunc versionfunc)
{
  Option *o = option_new("version", "display version information and exit",
                         versionfunc);
  o->option_type = OPTION_VERSION;
  return o;
}

Option* option_ref(Option *o)
{
  assert(o);
  o->reference_count++;
  return o;
}

OptionParser* option_parser_new(const char *synopsis, const char *one_liner)
{
  OptionParser *op = ma_malloc(sizeof (OptionParser));
  assert(synopsis && one_liner);
  assert("one_liner must have upper case letter at start and '.' at end" &&
         strlen(one_liner) && isupper((int) one_liner[0]));
  assert(one_liner[strlen(one_liner)-1] == '.');
  op->progname = NULL;
  op->synopsis = cstr_dup(synopsis);
  op->one_liner = cstr_dup(one_liner);
  op->options = array_new(sizeof (Option*));
  op->hooks = NULL;
  op->parser_called = false;
  op->comment_func = NULL;
  op->comment_func_data = NULL;
  op->mailaddress = NULL;
  op->min_additional_arguments = UNDEF_UINT;
  op->max_additional_arguments = UNDEF_UINT;
  return op;
}

void option_parser_add_option(OptionParser *op, Option *o)
{
  assert(op && o);
  array_add(op->options, o);
}

void option_parser_set_comment_func(OptionParser *op,
                                    ShowCommentFunc comment_func, void *data)
{
  assert(op);
  op->comment_func = comment_func;
  op->comment_func_data = data;
}

void option_parser_register_hook(OptionParser *op, OptionParserHookFunc hook,
                                 void *data)
{
  HookInfo hookinfo;
  assert(op && hook);
  if (!op->hooks)
    op->hooks = array_new(sizeof (HookInfo));
  hookinfo.hook = hook;
  hookinfo.data = data;
  array_add(op->hooks, hookinfo);
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

static int show_help(OptionParser *op, OptionType optiontype, Error *err)
{
  unsigned long i, max_option_length = 0;
  Option *option;
  int had_err = 0;
  error_check(err);
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
      else if (option->option_type == OPTION_RANGE) {
        printf("%*s  default: ", (int) max_option_length, "");
        if (option->default_value.r.start == UNDEF_ULONG)
          xputs("undefined");
        else {
          printf("%lu %lu\n", option->default_value.r.start,
                 option->default_value.r.end);
        }
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
    had_err = op->comment_func(op->progname, op->comment_func_data, err);
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

static int check_missing_argument(int argnum, int argc, Str *option, Error *err)
{
  error_check(err);
  if (argnum + 1 >= argc) {
    error_set(err, "missing argument to option \"-%s\"", str_get(option));
    return -1;
  }
  return 0;
}

static int check_mandatory_options(OptionParser *op, Error *err)
{
  unsigned long i;
  Option *o;
  error_check(err);
  assert(op);
  for (i = 0; i < array_size(op->options); i++) {
    o = *(Option**) array_get(op->options, i);
    if (o->is_mandatory && !o->is_set) {
      error_set(err, "option \"-%s\" is mandatory", str_get(o->option_str));
      return -1;
    }
  }
  return 0;
}

static int check_option_implications(OptionParser *op, Error *err)
{
  unsigned long i, j, k, l;
  Array *implied_option_array;
  Option *o, *implied_option;
  unsigned int option_set;
  Str *error_str;
  error_check(err);

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
            error_set(err, "option \"-%s\" requires option \"-%s\"",
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
            error_str = str_new_cstr("option \"-");
            str_append_str(error_str, o->option_str);
            str_append_cstr(error_str, "\" requires option");
            for (l = 0; l < array_size(implied_option_array) - 1; l++) {
              str_append_cstr(error_str, " \"-");
              str_append_str(error_str, (*(Option**)
                             array_get(implied_option_array, l))->option_str);
              str_append_cstr(error_str, "\"");
              if (array_size(implied_option_array) > 2)
                str_append_char(error_str, ',');
            }
            str_append_cstr(error_str, " or \"-");
            str_append_str(error_str, (*(Option**)
                           array_get(implied_option_array,
                                     array_size(implied_option_array) - 1))
                                     ->option_str);
            str_append_cstr(error_str, "\"");
            error_set(err, "%s", str_get(error_str));
            str_delete(error_str);
            return -1;
          }
        }
      }
    }
  }
  return 0;
}

static int check_option_exclusions(OptionParser *op, Error *err)
{
  unsigned long i, j;
  Option *o, *excluded_option;
  error_check(err);

  for (i = 0; i < array_size(op->options); i++) {
    o = *(Option**) array_get(op->options, i);
    if (o->exclusions && o->is_set) {
      for (j = 0; j < array_size(o->exclusions); j++) {
        excluded_option = *(Option**) array_get(o->exclusions, j);
        if (excluded_option->is_set) {
          error_set(err, "option \"-%s\" and option \"-%s\" exclude each other",
                    str_get(o->option_str),
                    str_get(excluded_option->option_str));
          return -1;
        }
      }
    }
  }
  return 0;
}

static int check_mandatory_either_options(OptionParser *op, Error *err)
{
  unsigned long i;
  Option *o;
  error_check(err);

  for (i = 0; i < array_size(op->options); i++) {
    o = *(Option**) array_get(op->options, i);
    if (o->mandatory_either_option) {
      if (!o->is_set && !o->mandatory_either_option->is_set) {
        error_set(err, "either option \"-%s\" or option \"-%s\" is mandatory",
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

void option_parser_set_min_args(OptionParser *op,
                                unsigned int min_additional_arguments)
{
  assert(op);
  op->min_additional_arguments = min_additional_arguments;
}

void option_parser_set_max_args(OptionParser *op,
                                unsigned int max_additional_arguments)
{
  assert(op);
  op->max_additional_arguments = max_additional_arguments;
}

void option_parser_set_min_max_args(OptionParser *op,
                                    unsigned int min_additional_arguments,
                                    unsigned int max_additional_arguments)
{
  assert(op);
  op->min_additional_arguments = min_additional_arguments;
  op->max_additional_arguments = max_additional_arguments;
}

OPrval option_parser_parse(OptionParser *op, int *parsed_args, int argc,
                           const char **argv, ShowVersionFunc versionfunc,
                           Error *err)
{
  int argnum, int_value;
  unsigned int uint_value;
  unsigned long i;
  double double_value;
  HookInfo *hookinfo;
  Option *option;
  bool has_extended_options, option_parsed;
  long long_value;
  int minus_offset, had_err = 0;
  Str *error_str;

  error_check(err);
  assert(op);
  assert(!op->parser_called); /* to avoid multiple adding of common options */

  op->progname = cstr_dup(argv[0]);

  /* add common options */
  has_extended_options = has_extended_option(op->options);
  option = option_new_help(has_extended_options);
  option_parser_add_option(op, option);
  if (has_extended_options) {
    option = option_new_helpplus();
    option_parser_add_option(op, option);
  }
  option = option_new_helpdev();
  option_parser_add_option(op, option);
  option = option_new_version(versionfunc);
  option_parser_add_option(op, option);

  for (argnum = 1; argnum < argc; argnum++) {
    if (!(argv[argnum] && argv[argnum][0] == '-' && strlen(argv[argnum]) > 1) ||
        !strcmp(argv[argnum], "--")) {
      break;
    }

    /* look for matching option */
    option_parsed = false;
    for (i = 0; i < array_size(op->options); i++) {
      option = *(Option**) array_get(op->options, i);

      /* allow options to start with '--', too */
      minus_offset = argv[argnum][1] == '-' ? 1 : 0;
      if (!strcmp(argv[argnum]+1+minus_offset, str_get(option->option_str))) {
        /* make sure option has not been used before */
        if (option->is_set) {
          error_set(err, "option \"%s\" already set",
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
                                               err);
              if (!had_err) {
                argnum++;
                if (strcmp(argv[argnum], option->domain[0])) {
                  error_str = str_new_cstr(option->domain[0]);
                  i = 1;
                  while (option->domain[i] != NULL) {
                    if (!strcmp(argv[argnum], option->domain[i])) {
                      str_set(option->value, option->domain[i]);
                      break;
                    }
                    str_append_cstr(error_str, ", ");
                    str_append_cstr(error_str, option->domain[i]);
                    i++;
                  }
                  if (option->domain[i] == NULL) {
                    error_set(err, "argument to option \"-%s\" must be one "
                                   "of: %s", str_get(option->option_str),
                              str_get(error_str));
                    had_err = -1;
                  }
                  str_delete(error_str);
                }
                else {
                  str_set(option->value, option->domain[0]);
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
                                               err);
              if (!had_err) {
                argnum++;
                if (parse_double(&double_value, argv[argnum])) {
                  error_set(err, "argument to option \"-%s\" must be "
                                 "floating-point number",
                            str_get(option->option_str));
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* minimum value check */
                if (option->min_value_set &&
                    double_value < option->min_value.d) {
                  error_set(err, "argument to option \"-%s\" must be a "
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
                  error_set(err, "argument to option \"-%s\" must be a "
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
              if (show_help(op, OPTION_HELP, err))
                return OPTIONPARSER_ERROR;
              return OPTIONPARSER_REQUESTS_EXIT;
            case OPTION_HELPPLUS:
              if (show_help(op, OPTION_HELPPLUS, err))
                return OPTIONPARSER_ERROR;
              return OPTIONPARSER_REQUESTS_EXIT;
            case OPTION_HELPDEV:
              if (show_help(op, OPTION_HELPDEV, err))
                return OPTIONPARSER_ERROR;
              return OPTIONPARSER_REQUESTS_EXIT;
            case OPTION_OUTPUTFILE:
              had_err = check_missing_argument(argnum, argc, option->option_str,
                                               err);
              if (!had_err) {
                argnum++;
                *(FILE**) option->value = fa_xfopen(argv[argnum], "w");
                option_parsed = true;
              }
              break;
            case OPTION_INT:
              assert(!option->argument_is_optional);
              had_err = check_missing_argument(argnum, argc, option->option_str,
                                               err);
              if (!had_err) {
                argnum++;
                if (parse_int(&int_value, argv[argnum])) {
                  error_set(err, "argument to option \"-%s\" must be an "
                                 "integer", str_get(option->option_str));
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* minimum value check */
                if (option->min_value_set && int_value < option->min_value.i) {
                  error_set(err, "argument to option \"-%s\" must be an "
                                 "integer >= %d", str_get(option->option_str),
                            option->min_value.i);
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* maximum value check */
                if (option->max_value_set && int_value > option->max_value.i) {
                  error_set(err, "argument to option \"-%s\" must be an "
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
                                               err);
              if (!had_err) {
                argnum++;
                if (parse_uint(&uint_value, argv[argnum])) {
                  error_set(err, "argument to option \"-%s\" must be a "
                                 "non-negative integer <= %u",
                            str_get(option->option_str), UINT_MAX);
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* minimum value check */
                if (option->min_value_set
                    && uint_value < option->min_value.ui) {
                  error_set(err, "argument to option \"-%s\" must be an "
                                 "integer >= %u", str_get(option->option_str),
                            option->min_value.ui);
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* maximum value check */
                if (option->max_value_set
                    && uint_value > option->max_value.ui) {
                  error_set(err, "argument to option \"-%s\" must be an "
                                 "integer <= %u", str_get(option->option_str),
                            option->max_value.ui);
                  had_err = -1;
                }
              }
              if (!had_err) {
                *(unsigned int*) option->value = uint_value;
                option_parsed = true;
              }
              break;
            case OPTION_LONG:
              assert(!option->argument_is_optional);
              had_err = check_missing_argument(argnum, argc, option->option_str,
                                               err);
              if (!had_err) {
                argnum++;
                if (parse_long(&long_value, argv[argnum])) {
                  error_set(err, "argument to option \"-%s\" must be an "
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
                                               err);
              if (!had_err) {
                argnum++;
                if (parse_long(&long_value, argv[argnum]) || long_value < 0) {
                  error_set(err, "argument to option \"-%s\" must be a "
                                 "non-negative integer",
                            str_get(option->option_str));
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* minimum value check */
                if (option->min_value_set &&
                    long_value < option->min_value.ul) {
                  error_set(err, "argument to option \"-%s\" must be an "
                                 "integer >= %lu", str_get(option->option_str),
                            option->min_value.ul);
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* maximum value check */
                if (option->max_value_set &&
                    long_value > option->max_value.ul) {
                  error_set(err, "argument to option \"-%s\" must be an "
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
            case OPTION_RANGE:
              if (optional_arg(option, argnum, argc, argv)) {
                option_parsed = true;
                break;
              }
              /* parse first argument */
              had_err = check_missing_argument(argnum, argc, option->option_str,
                                               err);
              if (!had_err) {
                argnum++;
                if (parse_long(&long_value, argv[argnum]) || long_value < 0) {
                  error_set(err, "first argument to option \"-%s\" must be a "
                                 "non-negative integer",
                            str_get(option->option_str));
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* minimum value check */
                if (option->min_value_set &&
                    long_value < option->min_value.ul) {
                  error_set(err, "first argument to option \"-%s\" must be an "
                                 "integer >= %lu", str_get(option->option_str),
                            option->min_value.ul);
                  had_err = -1;
                }
                else
                  ((Range*) option->value)->start = long_value;
              }
              /* parse second argument */
              if (!had_err) {
                had_err = check_missing_argument(argnum, argc,
                                                 option->option_str, err);
              }
              if (!had_err) {
                argnum++;
                if (parse_long(&long_value, argv[argnum]) || long_value < 0) {
                  error_set(err, "second argument to option \"-%s\" must be a "
                                 "non-negative integer",
                            str_get(option->option_str));
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* maximum value check */
                if (option->max_value_set &&
                    long_value > option->max_value.ul) {
                  error_set(err, "second argument to option \"-%s\" must be an "
                                 "integer <= %lu", str_get(option->option_str),
                            option->max_value.ul);
                  had_err = -1;
                }
                else
                  ((Range*) option->value)->end = long_value;
              }
              /* check arguments */
              if (!had_err && (((Range*) option->value)->start >
                               ((Range*) option->value)->end)) {
                error_set(err, "first argument %lu to option \"-%s\" must be "
                               "<= than second argument %lu",
                          ((Range*) option->value)->start,
                          str_get(option->option_str),
                          ((Range*) option->value)->end);
                had_err = -1;
              }
              if (!had_err)
                option_parsed = true;
              break;
            case OPTION_STRING:
              if (optional_arg(option, argnum, argc, argv)) {
                option_parsed = true;
                break;
              }
              had_err = check_missing_argument(argnum, argc, option->option_str,
                                               err);
              if (!had_err) {
                argnum++;
                str_set(option->value, argv[argnum]);
                option_parsed = true;
              }
              break;
            case OPTION_STRINGARRAY:
              if (optional_arg(option, argnum, argc, argv)) {
                option_parsed = true;
                break;
              }
              had_err = check_missing_argument(argnum, argc, option->option_str,
                                               err);
              while (!had_err) {
                if (argnum + 1 < argc && argv[argnum+1][0] != '-') {
                  argnum++;
                  strarray_add_cstr(option->value, argv[argnum]);
                  option_parsed = true;
                }
                else {
                  if (!option_parsed) {
                    error_set(err, "missing argument to option \"-%s\"",
                              str_get(option->option_str));
                    had_err = -1;
                  }
                  break;
                }
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
    error_set(err, "unknown option: %s (-help shows possible options)",
              argv[argnum]);
    had_err = -1;
    break;
  }

  /* skip "--" if necessary */
  if (argnum < argc && !strcmp(argv[argnum], "--"))
    argnum++;

  /* check for minimum number of additional arguments, if necessary */
  if (!had_err && op->min_additional_arguments != UNDEF_UINT &&
      argc - argnum < op->min_additional_arguments) {
    error_set(err, "missing argument\nUsage: %s %s", op->progname,
              op->synopsis);
    had_err = -1;
  }

  /* check for maximal number of additional arguments, if necessary */
  if (!had_err && op->max_additional_arguments != UNDEF_UINT &&
      argc - argnum > op->max_additional_arguments) {
    error_set(err, "superfluous argument \"%s\"\nUsage: %s %s",
              argv[argnum + op->max_additional_arguments], op->progname,
              op->synopsis);
    had_err = -1;
  }

  if (!had_err)
    had_err = check_mandatory_options(op, err);
  if (!had_err)
    had_err = check_option_implications(op, err);
  if (!had_err)
    had_err = check_option_exclusions(op, err);
  if (!had_err)
    had_err = check_mandatory_either_options(op, err);

  /* call hooks */
  for (i = 0; !had_err && i < array_size(op->hooks); i++) {
    hookinfo = array_get(op->hooks, i);
    had_err = hookinfo->hook(hookinfo->data, err);
  }

  op->parser_called = true;
  if (parsed_args)
    *parsed_args = argnum;

  if (had_err)
    return OPTIONPARSER_ERROR;
  return OPTIONPARSER_OK;
}

void option_parser_delete(OptionParser *op)
{
  unsigned long i;
  if (!op) return;
  ma_free(op->progname);
  ma_free(op->synopsis);
  ma_free(op->one_liner);
  for (i = 0; i < array_size(op->options); i++)
    option_delete(*(Option**) array_get(op->options, i));
  array_delete(op->options);
  array_delete(op->hooks);
  ma_free(op);
}

Option* option_new_outputfile(FILE **outfp)
{
  Option *o = option_new("o", "redirect output to specified file (will "
                         "overwrite existing file!)", outfp);
  o->option_type = OPTION_OUTPUTFILE;
  o->default_value.fp = stdout;
  *outfp = stdout;
  return o;
}

Option* option_new_verbose(bool *value)
{
  return option_new_bool("v", "be verbose", value, false);
}

Option* option_new_debug(bool *value)
{
  Option *o = option_new_bool("debug", "enable debugging output", value, false);
  o->is_development_option = true;
  return o;
}

Option* option_new_bool(const char *option_str, const char *description,
                        bool *value, bool default_value)
{
  Option *o = option_new(option_str, description, value);
  o->option_type = OPTION_BOOL;
  o->default_value.b = default_value;
  *value = default_value;
  return o;
}

Option* option_new_double(const char *option_str, const char *description,
                          double *value, double default_value)
{
  Option *o = option_new(option_str, description, value);
  o->option_type = OPTION_DOUBLE;
  o->default_value.d = default_value;
  *value = default_value;
  return o;
}

Option *option_new_double_min(const char *option_str, const char *description,
                              double *value, double default_value,
                              double min_value)
{
  Option *o = option_new_double(option_str, description, value, default_value);
  o->min_value_set = true;
  o->min_value.d = min_value;
  return o;
}

Option *option_new_double_min_max(const char *option_str,
                                  const char *description, double *value,
                                  double default_value, double min_value,
                                  double max_value)
{
  Option *o = option_new_double(option_str, description, value, default_value);
  o->min_value_set = true;
  o->min_value.d = min_value;
  o->max_value_set = true;
  o->max_value.d = max_value;
  return o;
}

Option* option_new_probability(const char *option_str, const char *description,
                               double *value, double default_value)
{
  return option_new_double_min_max(option_str, description, value,
                                   default_value, 0.0, 1.0);
}

Option* option_new_int(const char *option_str, const char *description,
                       int *value, int default_value)
{
  Option *o = option_new(option_str, description, value);
  o->option_type = OPTION_INT;
  o->default_value.i = default_value;
  *value = default_value;
  return o;
}

Option* option_new_int_min(const char *option_str, const char *description,
                           int *value, int default_value, int min_value)
{
  Option *o = option_new_int(option_str, description, value, default_value);
  o->min_value_set = true;
  o->min_value.i = min_value;
  return o;
}

Option* option_new_int_max(const char *option_str, const char *description,
                           int *value, int default_value, int max_value)
{
  Option *o = option_new_int(option_str, description, value, default_value);
  o->max_value_set = true;
  o->max_value.i = max_value;
  return o;
}

Option* option_new_int_min_max(const char *option_str, const char *description,
                               int *value, int default_value,
                               int min_value, int max_value)
{
  Option *o = option_new_int(option_str, description, value, default_value);
  o->min_value_set = true;
  o->min_value.i = min_value;
  o->max_value_set = true;
  o->max_value.i = max_value;
  return o;
}

Option* option_new_uint(const char *option_str, const char *description,
                        unsigned int *value, unsigned int default_value)
{
  Option *o = option_new(option_str, description, value);
  o->option_type = OPTION_UINT;
  o->default_value.ui = default_value;
  *value = default_value;
  return o;
}

Option* option_new_uint_min(const char *option_str, const char *description,
                            unsigned int *value, unsigned int default_value,
                            unsigned int min_value)
{
  Option *o = option_new_uint(option_str, description, value, default_value);
  o->min_value_set = true;
  o->min_value.ui = min_value;
  return o;
}

Option* option_new_uint_max(const char *option_str, const char *description,
                            unsigned int *value, unsigned int default_value,
                            unsigned int max_value)
{
  Option *o = option_new_uint(option_str, description, value, default_value);
  o->max_value_set = true;
  o->max_value.ui = max_value;
  return o;
}

Option *option_new_uint_min_max(const char *option_str, const char *description,
                                unsigned int *value, unsigned int default_value,
                                unsigned int min_value, unsigned int max_value)
{
  Option *o = option_new_uint(option_str, description, value, default_value);
  o->min_value_set = true;
  o->min_value.i = min_value;
  o->max_value_set = true;
  o->max_value.i = max_value;
  return o;
}

Option* option_new_long(const char *option_str, const char *description,
                        long *value, long default_value)
{
  Option *o = option_new(option_str, description, value);
  o->option_type = OPTION_LONG;
  o->default_value.l = default_value;
  *value = default_value;
  return o;
}

Option* option_new_ulong(const char *option_str, const char *description,
                         unsigned long *value, unsigned long default_value)
{
  Option *o = option_new(option_str, description, value);
  o->option_type = OPTION_ULONG;
  o->default_value.ul = default_value;
  *value = default_value;
  return o;
}

Option* option_new_ulong_min(const char *option_str, const char *description,
                             unsigned long *value, unsigned long default_value,
                             unsigned long min_value)
{
  Option *o = option_new_ulong(option_str, description, value, default_value);
  o->min_value_set = true;
  o->min_value.ul = min_value;
  return o;
}

Option *option_new_ulong_min_max(const char *option_str,
                                 const char *description, unsigned long *value,
                                 unsigned long default_value,
                                 unsigned long min_value,
                                 unsigned long max_value)
{
  Option *o = option_new_ulong(option_str, description, value, default_value);
  o->min_value_set = true;
  o->min_value.ul = min_value;
  o->max_value_set = true;
  o->max_value.ul = max_value;
  return o;
}

Option* option_new_range(const char *option_str, const char *description,
                         Range *value, Range *default_value)
{
  Option *o = option_new(option_str, description, value);
  o->option_type = OPTION_RANGE;
  o->default_value.r.start = default_value ? default_value->start : UNDEF_ULONG;
  o->default_value.r.end   = default_value ? default_value->end   : UNDEF_ULONG;
  value->start = o->default_value.r.start;
  value->end   = o->default_value.r.end;
  return o;
}

Option* option_new_range_min_max(const char *option_str,
                                 const char *description, Range *value,
                                 Range *default_value,
                                 unsigned long min_value,
                                 unsigned long max_value)
{
   Option *o = option_new_range(option_str, description,
                                value, default_value);
   o->min_value_set = true;
   o->min_value.ul = min_value;
   o->max_value_set = true;
   o->max_value.ul = max_value;
   return o;
}

Option* option_new_string(const char *option_str, const char *description,
                          Str *value, const char *default_value)
{
  Option *o = option_new(option_str, description, value);
  o->option_type = OPTION_STRING;
  o->default_value.s = default_value;
  str_set(value, default_value);
  return o;
}

Option* option_new_stringarray(const char *option_str,
                               const char *description, StrArray *value)
{
  Option *o = option_new(option_str, description, value);
  o->option_type = OPTION_STRINGARRAY;
  return o;
}

/* the following function would allow to handle files differently from strings
   later on (e.g., for CGI scripts), but for now the are implemented in the same
   way */
Option* option_new_filename(const char *option_str, const char *description,
                            Str *filename)
{
  return option_new_string(option_str, description, filename, NULL);
}

/* the following function would allow to handle file arrays differently from
   string arrays later on (e.g., for CGI scripts) , but for now the are
   implemented in the same way */
Option* option_new_filenamearray(const char *option_str,
                                 const char *description, StrArray *filenames)
{
  return option_new_stringarray(option_str, description, filenames);
}

Option* option_new_choice(const char *option_str, const char *description,
                          Str *value, const char *default_value,
                          const char **domain)
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

  o = option_new_string(option_str, description, value, default_value);
  o->option_type = OPTION_CHOICE;
  o->domain = domain;

  return o;
}

const char* option_get_name(const Option *o)
{
  assert(o);
  return str_get(o->option_str);
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

void option_imply(Option *o, const Option *implied_option)
{
  Array *option_array;
  assert(o && implied_option);
  if (!o->implications)
    o->implications = array_new(sizeof (Array*));
  option_array = array_new(sizeof (Option*));
  array_add(option_array, implied_option);
  array_add(o->implications, option_array);
}

void option_imply_either_2(Option *o, const Option *io1, const Option *io2)
{
  Array *option_array;
  assert(o && io1 && io2);
  if (!o->implications)
    o->implications = array_new(sizeof (Array*));
  option_array = array_new(sizeof (Option*));
  array_add(option_array, io1);
  array_add(option_array, io2);
  array_add(o->implications, option_array);
}

void option_exclude(Option *o_a, Option *o_b)
{
  assert(o_a && o_b);
  if (!o_a->exclusions)
    o_a->exclusions = array_new(sizeof (Option*));
  if (!o_b->exclusions)
    o_b->exclusions = array_new(sizeof (Option*));
  array_add(o_a->exclusions, o_b);
  array_add(o_b->exclusions, o_a);
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

void option_delete(Option *o)
{
  unsigned long i;
  if (!o) return;
  if (o->reference_count) {
    o->reference_count--;
    return;
  }
  str_delete(o->option_str);
  str_delete(o->description);
  for (i = 0; i < array_size(o->implications); i++)
    array_delete(*(Array**) array_get(o->implications, i));
  array_delete(o->implications);
  array_delete(o->exclusions);
  ma_free(o);
}
