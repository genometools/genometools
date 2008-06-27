/*
  Copyright (c) CCYY YOUR NAME HERE <user@your.dom.ain>
  Copyright (c) CCYY Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/ma.h"
#include "libgtcore/unused.h"
#include "tools/gt_template.h"

typedef struct {
  bool bool_option_template;
  Str  *str_option_template;
} TemplateArguments;

static void* gt_template_arguments_new(void)
{
  TemplateArguments *arguments = ma_calloc(1, sizeof *arguments);
  arguments->str_option_template = str_new();
  return arguments;
}

static void gt_template_arguments_delete(void *tool_arguments)
{
  TemplateArguments *arguments = tool_arguments;
  if (!arguments) return;
  str_delete(arguments->str_option_template);
  ma_free(arguments);
}

static OptionParser* gt_template_option_parser_new(void *tool_arguments)
{
  TemplateArguments *arguments = tool_arguments;
  OptionParser *op;
  Option *option;
  assert(arguments);

  /* init */
  op = option_parser_new("[option ...] [file]", /* XXX */
                         "DESCRIBE YOUR TOOL IN ONE LINE HERE."); /* XXX */

  /* -bool */
  option = option_new_bool("bool", "bool option template",
                           &arguments->bool_option_template, false);
  option_parser_add_option(op, option);

  /* -str */
  option = option_new_string("str", "str option template",
                             arguments->str_option_template, NULL);
  option_parser_add_option(op, option);

  return op;
}

static int gt_template_arguments_check(UNUSED int rest_argc,
                                       void *tool_arguments, UNUSED Error *err)
{
  TemplateArguments *arguments = tool_arguments;
  int had_err = 0;
  error_check(err);
  assert(arguments);

  /* XXX: do some checking after the option have been parsed (usally this is not
     necessary and this function can be removed completely). */
  if (str_length(arguments->str_option_template))
    printf("%s\n", str_get(arguments->str_option_template));

  return had_err;
}

static int gt_template_runner(int argc, const char **argv, int parsed_args,
                              void *tool_arguments, UNUSED Error *err)
{
  TemplateArguments *arguments = tool_arguments;
  int had_err = 0;

  error_check(err);
  assert(arguments);

  /* XXX */
  if (arguments->bool_option_template)
    printf("argc=%d, parsed_args=%d\n", argc, parsed_args);
  printf("argv[0]=%s\n", argv[0]);

  return had_err;
}

Tool* gt_template(void)
{
  return tool_new(gt_template_arguments_new,
                  gt_template_arguments_delete,
                  gt_template_option_parser_new,
                  gt_template_arguments_check,
                  gt_template_runner);
}
