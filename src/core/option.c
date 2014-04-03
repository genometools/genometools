/*
  Copyright (c) 2006-2011 Gordon Gremme <gordon@gremme.org>
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

#include <ctype.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include "gt_config.h"
#include "core/cstr_api.h"
#include "core/error.h"
#include "core/fa.h"
#include "core/fileutils_api.h"
#include "core/hashmap_api.h"
#include "core/grep_api.h"
#include "core/ma.h"
#include "core/mail_address.h"
#include "core/option.h"
#include "core/parseutils.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "core/xansi_api.h"
#include "core/xposix.h"

typedef enum {
  OPTION_BOOL,
  OPTION_CHOICE,
  OPTION_DOUBLE,
  OPTION_FILENAME,
  OPTION_HELP,
  OPTION_HELPPLUS,
  OPTION_HELPDEV,
  OPTION_INT,
  OPTION_UINT,
  OPTION_LONG,
  OPTION_ULONG,
  OPTION_RANGE,
  OPTION_STRING,
  OPTION_STRING_ARRAY,
  OPTION_VERSION,
} GtOptionType;

typedef struct {
  GtOptionParserHookFunc hook;
  void *data;
} HookInfo;

struct GtOptionParser {
  char *progname,
       *synopsis,
       *one_liner;
  GtArray *options,
          *hooks;
  bool common_options_added,
       refer_to_manual;
  GtShowCommentFunc comment_func;
  void *comment_func_data;
  GtShowVersionFunc version_func;
  const char *mail_address;
  unsigned int min_additional_arguments,
               max_additional_arguments;
  GtHashmap *optionindex;
};

struct GtOption {
  GtOptionType option_type;
  GtStr *option_str,
      *description;
  void *value;
  union {
    bool b;
    double d;
    int i;
    unsigned int ui;
    GtWord l;
    GtUword ul;
    GtRange r;
    const char *s;
  } default_value;
  const char** domain;
  union {
    double d;
    int i;
    unsigned int ui;
    GtUword ul;
  } min_value;
  union {
    double d;
    int i;
    unsigned int ui;
    GtUword ul;
  } max_value;
  bool is_set,
       is_mandatory,
       hide_default,
       min_value_set,
       max_value_set,
       is_extended_option,
       is_development_option,
       argument_is_optional;
  GtArray *implications, /* contains option arrays, from each array at least
                            one option needs to be set */
          *exclusions,
          *mandatory_either_options;
  unsigned int reference_count;
};

static GtOption *gt_option_new(const char *option_str, const char *description,
                               void *value)
{
  GtOption *o = gt_calloc(1, sizeof (GtOption));
  gt_assert(option_str && strlen(option_str));
  gt_assert("an option string should not start with '-', this is added "
            "automatically"  && option_str[0] != '-');
  o->option_str = gt_str_new_cstr(option_str);
  o->description = gt_str_new_cstr(description);
  o->value = value;
  return o;
}

static GtOption* gt_option_new_help(bool has_extended_options)
{
  GtOption *o;
  if (has_extended_options)
    o = gt_option_new("help", "display help for basic options and exit", NULL);
  else
    o = gt_option_new("help", "display help and exit", NULL);
  o->option_type = OPTION_HELP;
  o->default_value.b = false;
  return o;
}

static GtOption* gt_option_new_helpplus(void)
{
  GtOption *o = gt_option_new("help+", "display help for all options and exit",
                              NULL);
  o->option_type = OPTION_HELPPLUS;
  o->default_value.b = false;
  return o;
}

static GtOption* gt_option_new_helpdev(void)
{
  GtOption *o = gt_option_new("helpdev", "display help for development options "
                              "and exit", NULL);
  o->option_type = OPTION_HELPDEV;
  o->default_value.b = false;
  o->is_development_option = true;
  return o;
}

static GtOption* gt_option_new_version(void )
{
  GtOption *o = gt_option_new("version", "display version information and exit",
                              NULL);
  o->option_type = OPTION_VERSION;
  return o;
}

GtOption* gt_option_ref(GtOption *o)
{
  gt_assert(o);
  o->reference_count++;
  return o;
}

static void gt_option_reset(GtOption *option)
{
  gt_assert(option);
  if (option->option_type == OPTION_BOOL) {
    *(bool*) option->value = option->default_value.b;
  }
  else if (option->option_type == OPTION_CHOICE
             || option->option_type == OPTION_STRING) {
    gt_str_set((GtStr*) option->value, option->default_value.s);
  }
  else if (option->option_type == OPTION_DOUBLE) {
    *(double*) option->value = option->default_value.d;
  }
  else if (option->option_type == OPTION_INT) {
    *(int*) option->value = option->default_value.i;
  }
  else if (option->option_type == OPTION_UINT) {
    *(unsigned int*) option->value = option->default_value.ui;
  }
  else if (option->option_type == OPTION_LONG) {
    *(GtWord*) option->value = option->default_value.l;
  }
  else if (option->option_type == OPTION_ULONG) {
    *(GtUword*) option->value = option->default_value.ul;
  }
  else if (option->option_type == OPTION_RANGE) {
   *(GtRange*) option->value = option->default_value.r;
  }
}

GtOptionParser* gt_option_parser_new(const char *synopsis,
                                     const char *one_liner)
{
  GtOptionParser *op;
  gt_assert(synopsis && one_liner);
  gt_assert("one_liner must have upper case letter at start and '.' at end" &&
            strlen(one_liner) && isupper((int) one_liner[0]));
  gt_assert(one_liner[strlen(one_liner)-1] == '.');
  op = gt_calloc(1, sizeof *op);
  op->synopsis = gt_cstr_dup(synopsis);
  op->one_liner = gt_cstr_dup(one_liner);
  op->options = gt_array_new(sizeof (GtOption*));
  op->min_additional_arguments = GT_UNDEF_UINT;
  op->max_additional_arguments = GT_UNDEF_UINT;
  op->optionindex = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  return op;
}

void gt_option_parser_add_option(GtOptionParser *op, GtOption *o)
{
  gt_assert(op && o);
  gt_array_add(op->options, o);
  gt_hashmap_add(op->optionindex, gt_str_get(o->option_str), o);
}

GtOption* gt_option_parser_get_option(GtOptionParser *op,
                                      const char *option_str)
{
  gt_assert(op && option_str);
  return gt_hashmap_get(op->optionindex, option_str);
}

void gt_option_parser_refer_to_manual(GtOptionParser *op)
{
  gt_assert(op);
  op->refer_to_manual = true;
}

void gt_option_parser_set_comment_func(GtOptionParser *op,
                                       GtShowCommentFunc comment_func,
                                       void *data)
{
  gt_assert(op);
  op->comment_func = comment_func;
  op->comment_func_data = data;
}

void gt_option_parser_set_version_func(GtOptionParser *op,
                                       GtShowVersionFunc version_func)
{
  gt_assert(op && version_func);
  op->version_func = version_func;
}

void gt_option_parser_register_hook(GtOptionParser *op,
                                    GtOptionParserHookFunc hook,
                                    void *data)
{
  HookInfo hookinfo;
  gt_assert(op && hook);
  if (!op->hooks)
    op->hooks = gt_array_new(sizeof (HookInfo));
  hookinfo.hook = hook;
  hookinfo.data = data;
  gt_array_add(op->hooks, hookinfo);
}

static int reset_option(GT_UNUSED void *key, void *value, GT_UNUSED void *data,
                        GT_UNUSED GtError *err)
{
  GtOption *o = (GtOption*) value;
  gt_option_reset(o);
  return 0;
}

void gt_option_parser_reset(GtOptionParser *op)
{
  GT_UNUSED int rval;
  gt_assert(op);
  rval = gt_hashmap_foreach(op->optionindex, reset_option, NULL, NULL);
  gt_assert(!rval); /* reset_option() is sane */
}

void gt_option_parser_set_mail_address(GtOptionParser *op, const char *address)
{
  gt_assert(op && address);
  op->mail_address = address;
}

static void show_description(GtUword initial_space, const char *desc,
                             GtUword len)
{
  const GtUword width = GT_OPTION_PARSER_TERMINAL_WIDTH - initial_space;
  const char *tmp_ptr, *desc_ptr = desc;
  GtUword i;
  bool continue_while = false;

  /* got space to show option */
  gt_assert(initial_space < GT_OPTION_PARSER_TERMINAL_WIDTH);

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
      gt_xputchar(*desc_ptr);
      if (*desc_ptr == '\n') {
        /* show leading spaces */
        for  (i = 0; i < initial_space; i++)
          gt_xputchar(' ');
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
    gt_assert(*desc_ptr == ' ' || *desc_ptr == '\n');
    /* show newline for break point */
    desc_ptr++;
    gt_xputchar('\n');
    /* show leading spaces */
    for  (i = 0; i < initial_space; i++)
      gt_xputchar(' ');
  }
  /* show final line */
  while (desc_ptr < desc + len) {
    gt_xputchar(*desc_ptr);
    if (*desc_ptr == '\n') {
      /* show leading spaces */
      for  (i = 0; i < initial_space; i++)
        gt_xputchar(' ');
      desc_ptr++;
      continue;
    }
    desc_ptr++;
  }
  gt_xputchar('\n');
}

static int show_help(GtOptionParser *op, GtOptionType optiontype, GtError *err)
{
  GtUword i, max_option_length = 0;
  GtOption *option;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(optiontype == OPTION_HELP || optiontype == OPTION_HELPPLUS ||
            optiontype == OPTION_HELPDEV);

  /* determine maximum option length */
  for (i = 0; i < gt_array_size(op->options); i++) {
    option = *(GtOption**) gt_array_get(op->options, i);
    /* skip option if necessary */
    if ((optiontype == OPTION_HELP && option->is_extended_option) ||
        (optiontype == OPTION_HELPDEV && !option->is_development_option) ||
        (!(optiontype == OPTION_HELPDEV) && option->is_development_option)) {
      continue;
    }
    if (gt_str_length(option->option_str) > max_option_length)
      max_option_length = gt_str_length(option->option_str);
  }
  gt_assert(max_option_length);

  printf("Usage: %s %s\n", op->progname, op->synopsis);
  printf("%s\n\n", op->one_liner);
  for (i = 0; i < gt_array_size(op->options); i++) {
    option = *(GtOption**) gt_array_get(op->options, i);
    /* skip option if necessary */
    if ((optiontype == OPTION_HELP && option->is_extended_option) ||
        (optiontype == OPTION_HELPDEV && !option->is_development_option) ||
        (!(optiontype == OPTION_HELPDEV) && option->is_development_option)) {
      continue;
    }
    printf("-%s%*s ", gt_str_get(option->option_str),
           (int) (max_option_length - gt_str_length(option->option_str)), "");
    show_description(max_option_length + 2, gt_str_get(option->description),
                     gt_str_length(option->description));

    /* show default value for some types of options */
    if (!option->hide_default) {
      if (option->option_type == OPTION_BOOL) {
        printf("%*s  default: %s\n", (int) max_option_length, "",
               option->default_value.b ? "yes" : "no");
      }
      else if (option->option_type == OPTION_CHOICE) {
        printf("%*s  default: ", (int) max_option_length, "");
        if (!option->default_value.s || !strlen(option->default_value.s))
          gt_xputs("undefined");
        else
          gt_xputs(option->default_value.s);
      }

      else if (option->option_type == OPTION_DOUBLE) {
        printf("%*s  default: ", (int) max_option_length, "");
        if (option->default_value.d == GT_UNDEF_DOUBLE)
          gt_xputs("undefined");
        else
          printf("%.2f\n", option->default_value.d);
      }
      else if (option->option_type == OPTION_INT) {
        printf("%*s  default: ", (int) max_option_length, "");
        if (option->default_value.i == GT_UNDEF_INT)
          gt_xputs("undefined");
        else
          printf("%d\n", option->default_value.i);
      }
      else if (option->option_type == OPTION_UINT) {
        printf("%*s  default: ", (int) max_option_length, "");
        if (option->default_value.ui == GT_UNDEF_UINT)
          gt_xputs("undefined");
        else
          printf("%u\n", option->default_value.ui);
      }
      else if (option->option_type == OPTION_LONG) {
        printf("%*s  default: ", (int) max_option_length, "");
        if (option->default_value.ul == GT_UNDEF_WORD)
          gt_xputs("undefined");
        else
          printf(""GT_WD"\n", option->default_value.l);
      }
      else if (option->option_type == OPTION_ULONG) {
        printf("%*s  default: ", (int) max_option_length, "");
        if (option->default_value.ul == GT_UNDEF_UWORD)
          gt_xputs("undefined");
        else
          printf(""GT_WU"\n", option->default_value.ul);
      }
      else if (option->option_type == OPTION_RANGE) {
        printf("%*s  default: ", (int) max_option_length, "");
        if (option->default_value.r.start == GT_UNDEF_UWORD)
          gt_xputs("undefined");
        else {
          printf(""GT_WU" "GT_WU"\n", option->default_value.r.start,
                 option->default_value.r.end);
        }
      }
      else if (option->option_type == OPTION_FILENAME ||
               option->option_type == OPTION_STRING) {
        printf("%*s  default: ", (int) max_option_length, "");
        if (!option->default_value.s || !strlen(option->default_value.s))
          gt_xputs("undefined");
        else
          gt_xputs(option->default_value.s);
      }
    }
  }
  if (op->comment_func)
    had_err = op->comment_func(op->progname, op->comment_func_data, err);
  if (!had_err) {
    if (op->refer_to_manual) {
      printf("\nFor detailed information, please refer to the manual of%s.",
             op->progname + gt_cstr_length_up_to_char(op->progname, ' '));
    }
    printf("\nReport bugs to %s.\n",
           op->mail_address ? op->mail_address : GT_MAIL_ADDRESS);
  }
  return had_err;
}

static void print_asciidoc_header(const char *hdr, GtStr *outstr) {
  GtUword i;
  gt_str_append_cstr(outstr, hdr);
  gt_str_append_char(outstr, '\n');
  for (i = 0; i < strlen(hdr); i++)
    gt_str_append_char(outstr, '-');
  gt_str_append_cstr(outstr, "\n\n");
}

static void print_toolname(const char *toolname, GtStr *outstr, bool upper) {
  GtUword i;
  for (i = 0; i < strlen(toolname); i++) {
    char c = toolname[i];
    if (c == ' ')
      gt_str_append_char(outstr, '-');
    else {
      if (upper)
        gt_str_append_char(outstr, toupper(c));
      else
        gt_str_append_char(outstr, c);
    }
  }
}

static bool has_extended_option(GtArray *options)
{
  GtUword i;
  GtOption *option;
  gt_assert(options);
  for (i = 0; i < gt_array_size(options); i++) {
    option = *(GtOption**) gt_array_get(options, i);
    if (option->is_extended_option)
      return true;
  }
  return false;
}

static void add_common_options(GtOptionParser *op)
{
  bool has_extended_options;
  GtOption *option;
  gt_assert(op);
  has_extended_options = has_extended_option(op->options);
  option = gt_option_new_help(has_extended_options);
  gt_option_parser_add_option(op, option);
  if (has_extended_options) {
    option = gt_option_new_helpplus();
    gt_option_parser_add_option(op, option);
  }
  option = gt_option_new_helpdev();
  gt_option_parser_add_option(op, option);
  option = gt_option_new_version();
  gt_option_parser_add_option(op, option);
}

int gt_option_parser_manpage(GtOptionParser *op, const char *toolname,
                             GtStr *outstr, GtError *err)
{
  GtUword i;
  GtOption *option;
  GtStr *default_string;
  int had_err = 0;
  gt_assert(op);
  gt_error_check(err);
  gt_assert(strlen(toolname) > 0);

  /* add common options, if necessary */
  if (!op->common_options_added) {
    op->common_options_added = true;
    add_common_options(op);
  }

  /* title and section */
  print_toolname(toolname, outstr, true);
  gt_str_append_cstr(outstr, "(1)\n");
  for (i = 0; i < strlen(toolname) + 3; i++)
    gt_str_append_char(outstr, '=');
  gt_str_append_char(outstr, '\n');
  gt_str_append_cstr(outstr, ":man source:   GenomeTools\n");
  gt_str_append_cstr(outstr, ":man version:  ");
  gt_str_append_cstr(outstr, GT_VERSION);
  gt_str_append_char(outstr, '\n');
  gt_str_append_cstr(outstr, ":man manual:   GenomeTools Manual\n");
  gt_str_append_char(outstr, '\n');

  /* name */
  print_asciidoc_header("NAME", outstr);
  print_toolname(toolname, outstr, false);
  gt_str_append_cstr(outstr, " - ");
  gt_str_append_cstr(outstr, op->one_liner);
  gt_str_append_cstr(outstr, "\n\n");

  /* synopsis */
  print_asciidoc_header("SYNOPSIS", outstr);
  gt_str_append_cstr(outstr, "*");
  gt_str_append_cstr(outstr, toolname);
  gt_str_append_cstr(outstr, "* ");
  gt_str_append_cstr(outstr, op->synopsis);
  gt_str_append_cstr(outstr, "\n\n");

  /* description */
  if (gt_array_size(op->options)) {
    print_asciidoc_header("DESCRIPTION", outstr);
    default_string = gt_str_new();
    for (i = 0; i < gt_array_size(op->options); i++) {
      option = *(GtOption**) gt_array_get(op->options, i);
      /* skip dev options */
      if (option->is_development_option) continue;
      gt_str_append_cstr(outstr, "*-");
      gt_str_append_cstr(outstr, gt_str_get(option->option_str));
      gt_str_append_cstr(outstr, "* ");
      if (option->option_type == OPTION_BOOL) {
        gt_str_append_cstr(outstr, "['yes|no']");
        gt_str_append_cstr(default_string,
                              option->default_value.b ? "yes" : "no");
      }
      else if (option->option_type == OPTION_CHOICE) {
        gt_str_append_cstr(outstr, "['...']");
        if (!option->default_value.s || !strlen(option->default_value.s))
          gt_str_append_cstr(default_string, "undefined");
        else
          gt_str_append_cstr(default_string, option->default_value.s);
      }
      else if (option->option_type == OPTION_DOUBLE) {
        gt_str_append_cstr(outstr, "['value']");
        if (option->default_value.d == GT_UNDEF_DOUBLE)
          gt_str_append_cstr(default_string, "undefined");
        else
          gt_str_append_double(default_string, option->default_value.d, 6);
      }
      else if (option->option_type == OPTION_FILENAME) {
        gt_str_append_cstr(outstr, "['filename']");
        if (!option->default_value.s || !strlen(option->default_value.s))
          gt_str_append_cstr(default_string, "undefined");
        else
          gt_str_append_cstr(default_string, option->default_value.s);
      }
      else if (option->option_type == OPTION_INT) {
        gt_str_append_cstr(outstr, "['value']");
        if (option->default_value.i == GT_UNDEF_INT)
          gt_str_append_cstr(default_string, "undefined");
        else
          gt_str_append_int(default_string, option->default_value.i);
      }
      else if (option->option_type == OPTION_UINT) {
        gt_str_append_cstr(outstr, "['value']");
        if (option->default_value.ui == GT_UNDEF_UINT)
          gt_str_append_cstr(default_string, "undefined");
        else
          gt_str_append_uint(default_string, option->default_value.ui);
      }
      else if (option->option_type == OPTION_LONG) {
        gt_str_append_cstr(outstr, "['value']");
      }
      else if (option->option_type == OPTION_ULONG) {
        gt_str_append_cstr(outstr, "['value']");
        if (option->default_value.ul == GT_UNDEF_UWORD)
          gt_str_append_cstr(default_string, "undefined");
        else
          gt_str_append_ulong(default_string, option->default_value.ul);
      }
      else if (option->option_type == OPTION_RANGE) {
        gt_str_append_cstr(outstr, "['start' 'end']");
        if (option->default_value.r.start == GT_UNDEF_UWORD)
          gt_str_append_cstr(default_string, "undefined");
        else {
          gt_str_append_char(default_string, '[');
          gt_str_append_ulong(default_string, option->default_value.r.start);
          gt_str_append_cstr(default_string, "..");
          gt_str_append_ulong(default_string, option->default_value.r.end);
          gt_str_append_char(default_string, ']');
        }
      }
      else if (option->option_type == OPTION_STRING) {
        gt_str_append_cstr(outstr, "['string']");
        if (!option->default_value.s || !strlen(option->default_value.s))
          gt_str_append_cstr(default_string, "undefined");
        else
          gt_str_append_cstr(default_string, option->default_value.s);
      }
      gt_str_append_cstr(outstr, "::\n");
      gt_str_append_cstr(outstr, gt_str_get(option->description));
      if (gt_str_length(default_string) > 0) {
        gt_str_append_cstr(outstr, " (default: ");
        gt_str_append_cstr(outstr, gt_str_get(default_string));
        gt_str_append_cstr(outstr, ")");
        }
      gt_str_append_cstr(outstr, "\n\n");
      gt_str_reset(default_string);
    }
    gt_str_delete(default_string);
  }

  if (op->comment_func) {
    int old_stdout, out_pipe[2];
#ifndef _WIN32
    int rval;
    GtWord flags;
#endif
    char c,
         prognamebuf[BUFSIZ];

    /* Lua docu scripts print to stdout by themselves, so temporarily
       redirect stdout to a pipe. */
#ifndef _WIN32
    fflush(stdout);
    old_stdout = dup(STDOUT_FILENO);
    if ((rval = pipe(out_pipe)) == -1) {
      perror("pipe");
      exit(EXIT_FAILURE);  /* XXX */
    }

    flags = fcntl(out_pipe[0], F_GETFL);
    flags |= O_NONBLOCK;
    (void) fcntl(out_pipe[0], F_SETFL, flags);

    dup2(out_pipe[1], STDOUT_FILENO);
    close(out_pipe[1]);
#else
    /* XXX */
    fprintf(stderr, "pipe() stuff not implemented\n");
    exit(EXIT_FAILURE);
#endif

    if (!strcmp(toolname, "gt"))
      (void) snprintf(prognamebuf, BUFSIZ, "%s", gt_error_get_progname(err));
    else {
      (void) snprintf(prognamebuf, BUFSIZ, "%s %s",
                      gt_error_get_progname(err),
                      toolname + gt_cstr_length_up_to_char(toolname, ' '));
    }
    had_err = op->comment_func(prognamebuf, op->comment_func_data, err);
    fflush(stdout);

    while (read(out_pipe[0], &c, sizeof (char)) > 0)
      gt_str_append_char(outstr, c);

    dup2(old_stdout, STDOUT_FILENO);
    close(old_stdout);
  }

  gt_str_append_cstr(outstr, "\n");
  if (!had_err) {
    if (op->refer_to_manual) {
      print_asciidoc_header("ADDITIONAL INFORMATION", outstr);
      gt_str_append_cstr(outstr, "For detailed information, please refer to "
                                 "the manual of");
      gt_str_append_cstr(outstr, toolname +
                                 gt_cstr_length_up_to_char(toolname, ' '));
      gt_str_append_cstr(outstr, ".\n\n");
    }
    print_asciidoc_header("REPORTING BUGS", outstr);
    gt_str_append_cstr(outstr, "Report bugs to ");
    gt_str_append_cstr(outstr,
                       op->mail_address ? op->mail_address : GT_MAIL_ADDRESS);
    gt_str_append_cstr(outstr, ".\n");
  }
  return had_err;
}

static bool optional_arg(GtOption *o, int argnum, int argc, const char **argv)
{
  gt_assert(o);
  if (o->argument_is_optional &&
      (argnum + 1 >= argc || argv[argnum + 1][0] == '-' ||
       !strcmp(argv[argnum + 1], "--"))) {
    return true;
  }
  return false;
}

static int check_missing_argument(int argnum, int argc, GtStr *option,
                                  GtError *err)
{
  gt_error_check(err);
  if (argnum + 1 >= argc) {
    gt_error_set(err, "missing argument to option \"-%s\"", gt_str_get(option));
    return -1;
  }
  return 0;
}

static int check_missing_argument_and_minus_sign(int argnum, int argc,
                                                 GtStr *option,
                                                 const char **argv,
                                                 GtError *err)
{
  gt_error_check(err);
  if (argnum + 1 >= argc || (argv[argnum+1][0] == '-'
                               && argv[argnum+1][1] != '\0')) {
    gt_error_set(err, "missing argument to option \"-%s\"", gt_str_get(option));
    return -1;
  }
  return 0;
}

static int check_mandatory_options(GtOptionParser *op, GtError *err)
{
  GtUword i;
  GtOption *o;
  gt_error_check(err);
  gt_assert(op);
  for (i = 0; i < gt_array_size(op->options); i++) {
    o = *(GtOption**) gt_array_get(op->options, i);
    if (o->is_mandatory && !o->is_set) {
      gt_error_set(err, "option \"-%s\" is mandatory",
                   gt_str_get(o->option_str));
      return -1;
    }
  }
  return 0;
}

static int check_option_implications(GtOptionParser *op, GtError *err)
{
  GtUword i, j, k, l;
  GtArray *implied_option_array;
  GtOption *o, *implied_option;
  unsigned int option_set;
  GtStr *gt_error_str;
  gt_error_check(err);

  for (i = 0; i < gt_array_size(op->options); i++) {
    o = *(GtOption**) gt_array_get(op->options, i);
    if (o->implications && o->is_set) {
      for (j = 0; j < gt_array_size(o->implications); j++) {
        implied_option_array = *(GtArray**) gt_array_get(o->implications, j);
        gt_assert(gt_array_size(implied_option_array));
        if (gt_array_size(implied_option_array) == 1) {
          /* special case: option implies exactly one option */
          implied_option = *(GtOption**) gt_array_get(implied_option_array, 0);
          if (!implied_option->is_set) {
            gt_error_set(err, "option \"-%s\" requires option \"-%s\"",
                      gt_str_get(o->option_str),
                      gt_str_get(implied_option->option_str));
            return -1;
          }
        }
        else {
          /* ``either'' case: option implied at least one of the options given
             in array */
          option_set = 0;
          for (k = 0; k < gt_array_size(implied_option_array); k++) {
            implied_option = *(GtOption**) gt_array_get(implied_option_array,
                                                        k);
            if (implied_option->is_set) {
              option_set = 1;
              break;
            }
          }
          if (!option_set) {
            gt_error_str = gt_str_new_cstr("option \"-");
            gt_str_append_str(gt_error_str, o->option_str);
            gt_str_append_cstr(gt_error_str, "\" requires option");
            for (l = 0; l < gt_array_size(implied_option_array) - 1; l++) {
              gt_str_append_cstr(gt_error_str, " \"-");
              gt_str_append_str(gt_error_str,
                                (*(GtOption**)
                                 gt_array_get(implied_option_array, l))
                                ->option_str);
              gt_str_append_cstr(gt_error_str, "\"");
              if (gt_array_size(implied_option_array) > 2)
                gt_str_append_char(gt_error_str, ',');
            }
            gt_str_append_cstr(gt_error_str, " or \"-");
            gt_str_append_str(gt_error_str,
                              (*(GtOption**)
                               gt_array_get(implied_option_array,
                                            gt_array_size(implied_option_array)
                                            - 1))
                              ->option_str);
            gt_str_append_cstr(gt_error_str, "\"");
            gt_error_set(err, "%s", gt_str_get(gt_error_str));
            gt_str_delete(gt_error_str);
            return -1;
          }
        }
      }
    }
  }
  return 0;
}

static int check_option_exclusions(GtOptionParser *op, GtError *err)
{
  GtUword i, j;
  GtOption *o, *excluded_option;
  gt_error_check(err);

  for (i = 0; i < gt_array_size(op->options); i++) {
    o = *(GtOption**) gt_array_get(op->options, i);
    if (o->exclusions && o->is_set) {
      for (j = 0; j < gt_array_size(o->exclusions); j++) {
        excluded_option = *(GtOption**) gt_array_get(o->exclusions, j);
        if (excluded_option->is_set) {
          gt_error_set(err, "option \"-%s\" and option \"-%s\" exclude each "
                       "other", gt_str_get(o->option_str),
                       gt_str_get(excluded_option->option_str));
          return -1;
        }
      }
    }
  }
  return 0;
}

static int check_mandatory_either_options(GtOptionParser *op, GtError *err)
{
  GtUword i;
  GtOption *o, *meo_a, *meo_b, *meo_c;
  gt_error_check(err);

  for (i = 0; i < gt_array_size(op->options); i++) {
    o = *(GtOption**) gt_array_get(op->options, i);
    if (o->mandatory_either_options) {
      if (gt_array_size(o->mandatory_either_options) == 1) {
        meo_a = *(GtOption**) gt_array_get_first(o->mandatory_either_options);
        if (!o->is_set && !meo_a->is_set) {
          gt_error_set(err, "either option \"-%s\" or option \"-%s\" is "
                       "mandatory", gt_str_get(o->option_str),
                       gt_str_get(meo_a->option_str));
          return -1;
        }
      }
      else if (gt_array_size(o->mandatory_either_options) == 2) {
        meo_a = *(GtOption**) gt_array_get(o->mandatory_either_options, 0);
        meo_b = *(GtOption**) gt_array_get(o->mandatory_either_options, 1);
        if (!o->is_set && !meo_a->is_set && !meo_b->is_set) {
          gt_error_set(err, "either option \"-%s\", option \"-%s\" or option "
                       "\"-%s\" is mandatory", gt_str_get(o->option_str),
                       gt_str_get(meo_a->option_str),
                       gt_str_get(meo_b->option_str));
          return -1;
        }
      }
      else {
        /* XXX: this code needs to be generalized fo more than 3 mandatory
           options (if the corresponding methods are added) */
        if (gt_array_size(o->mandatory_either_options) == 3) {
          meo_a = *(GtOption**) gt_array_get(o->mandatory_either_options, 0);
          meo_b = *(GtOption**) gt_array_get(o->mandatory_either_options, 1);
          meo_c = *(GtOption**) gt_array_get(o->mandatory_either_options, 2);
          if (!o->is_set && !meo_a->is_set && !meo_b->is_set
                && !meo_c->is_set) {
            gt_error_set(err, "either option \"-%s\", option \"-%s\", option "
                              "\"-%s\" or option \"-%s\" is mandatory",
                              gt_str_get(o->option_str),
                              gt_str_get(meo_a->option_str),
                              gt_str_get(meo_b->option_str),
                              gt_str_get(meo_c->option_str));
            return -1;
          }
        }
      }
    }
  }
  return 0;
}

void gt_option_parser_set_min_args(GtOptionParser *op,
                                   unsigned int min_additional_arguments)
{
  gt_assert(op);
  op->min_additional_arguments = min_additional_arguments;
}

void gt_option_parser_set_max_args(GtOptionParser *op,
                                   unsigned int max_additional_arguments)
{
  gt_assert(op);
  op->max_additional_arguments = max_additional_arguments;
}

void gt_option_parser_set_min_max_args(GtOptionParser *op,
                                       unsigned int min_additional_arguments,
                                       unsigned int max_additional_arguments)
{
  gt_assert(op);
  op->min_additional_arguments = min_additional_arguments;
  op->max_additional_arguments = max_additional_arguments;
}

static void gt_option_parser_warn_on_nonascii(const char *str)
{
  GtUword j = 1;
  while (*str != '\0') {
    if (!isprint(*str)) {
      gt_warning("option '%s' contains non-ASCII character at "
                 "offset "GT_WU" -- possible copy&paste from "
                 "non-ASCII context?",
                 str, j);
      break;
    }
    str++;
    j++;
  }
}

GtOPrval gt_option_parser_parse(GtOptionParser *op, int *parsed_args, int argc,
                                const char **argv,
                                GtShowVersionFunc version_func, GtError *err)
{
  int argnum, int_value;
  unsigned int uint_value;
  GtUword i;
  double double_value;
  HookInfo *hookinfo;
  GtOption *option;
  bool option_parsed;
  GtWord long_value;
  GtUword ulong_value;
  int minus_offset, had_err = 0;
  GtStr *gt_error_str;

  gt_error_check(err);
  gt_assert(op);

  /* add common options, if necessary */
  if (!op->common_options_added) {
    op->common_options_added = true;
    add_common_options(op);
  }

  /* set progname */
  if (op->progname)
    gt_free(op->progname);
  op->progname = gt_cstr_dup(argv[0]);

  for (argnum = 1; argnum < argc; argnum++) {
    gt_option_parser_warn_on_nonascii(argv[argnum]);
  }

  for (argnum = 1; argnum < argc; argnum++) {
    if (!(argv[argnum] && argv[argnum][0] == '-' && strlen(argv[argnum]) > 1) ||
        !strcmp(argv[argnum], "--")) {
      break;
    }

    /* look for matching option */
    option_parsed = false;
    for (i = 0; i < gt_array_size(op->options); i++) {
      option = *(GtOption**) gt_array_get(op->options, i);

      /* allow options to start with '--', too */
      minus_offset = argv[argnum][1] == '-' ? 1 : 0;
      if (!strcmp(argv[argnum] + 1 + minus_offset,
                  gt_str_get(option->option_str))) {
        /* make sure option has not been used before */
        if (option->is_set) {
          gt_error_set(err, "option \"%s\" already set",
                       gt_str_get(option->option_str));
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
              gt_assert(option->domain[0]);
              had_err = check_missing_argument_and_minus_sign(argnum, argc,
                                                              option
                                                              ->option_str,
                                                              argv, err);
              if (!had_err) {
                argnum++;
                if (strcmp(argv[argnum], option->domain[0])) {
                  gt_error_str = gt_str_new_cstr(option->domain[0]);
                  i = 1;
                  while (option->domain[i] != NULL) {
                    if (!strcmp(argv[argnum], option->domain[i])) {
                      gt_str_set(option->value, option->domain[i]);
                      break;
                    }
                    gt_str_append_cstr(gt_error_str, ", ");
                    gt_str_append_cstr(gt_error_str, option->domain[i]);
                    i++;
                  }
                  if (option->domain[i] == NULL) {
                    gt_error_set(err, "argument to option \"-%s\" must be one "
                                      "of: %s", gt_str_get(option->option_str),
                                 gt_str_get(gt_error_str));
                    had_err = -1;
                  }
                  gt_str_delete(gt_error_str);
                }
                else {
                  gt_str_set(option->value, option->domain[0]);
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
                if (gt_parse_double(&double_value, argv[argnum])) {
                  gt_error_set(err, "argument to option \"-%s\" must be "
                                    "floating-point number",
                            gt_str_get(option->option_str));
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* minimum value check */
                if (option->min_value_set &&
                    double_value < option->min_value.d) {
                  gt_error_set(err, "argument to option \"-%s\" must be a "
                                    "floating point value >= %f",
                               gt_str_get(option->option_str),
                               option->min_value.d);
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* maximum value check */
                if (option->max_value_set &&
                    double_value > option->max_value.d) {
                  gt_error_set(err, "argument to option \"-%s\" must be a "
                                    "floating point value <= %f",
                               gt_str_get(option->option_str),
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
                return GT_OPTION_PARSER_ERROR;
              return GT_OPTION_PARSER_REQUESTS_EXIT;
            case OPTION_HELPPLUS:
              if (show_help(op, OPTION_HELPPLUS, err))
                return GT_OPTION_PARSER_ERROR;
              return GT_OPTION_PARSER_REQUESTS_EXIT;
            case OPTION_HELPDEV:
              if (show_help(op, OPTION_HELPDEV, err))
                return GT_OPTION_PARSER_ERROR;
              return GT_OPTION_PARSER_REQUESTS_EXIT;
            case OPTION_INT:
              gt_assert(!option->argument_is_optional);
              had_err = check_missing_argument(argnum, argc, option->option_str,
                                               err);
              if (!had_err) {
                argnum++;
                if (gt_parse_int(&int_value, argv[argnum])) {
                  gt_error_set(err, "argument to option \"-%s\" must be an "
                                    "integer", gt_str_get(option->option_str));
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* minimum value check */
                if (option->min_value_set && int_value < option->min_value.i) {
                  gt_error_set(err, "argument to option \"-%s\" must be an "
                               "integer >= %d", gt_str_get(option->option_str),
                               option->min_value.i);
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* maximum value check */
                if (option->max_value_set && int_value > option->max_value.i) {
                  gt_error_set(err, "argument to option \"-%s\" must be an "
                               "integer <= %d", gt_str_get(option->option_str),
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
              had_err = check_missing_argument_and_minus_sign(argnum, argc,
                                                             option->option_str,
                                                             argv, err);
              if (!had_err) {
                argnum++;
                if (gt_parse_uint(&uint_value, argv[argnum])) {
                  gt_error_set(err, "argument to option \"-%s\" must be a "
                                    "non-negative integer <= %u",
                               gt_str_get(option->option_str), UINT_MAX);
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* minimum value check */
                if (option->min_value_set
                    && uint_value < option->min_value.ui) {
                  gt_error_set(err, "argument to option \"-%s\" must be an "
                               "integer >= %u", gt_str_get(option->option_str),
                               option->min_value.ui);
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* maximum value check */
                if (option->max_value_set
                    && uint_value > option->max_value.ui) {
                  gt_error_set(err, "argument to option \"-%s\" must be an "
                               "integer <= %u", gt_str_get(option->option_str),
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
              gt_assert(!option->argument_is_optional);
              had_err = check_missing_argument(argnum, argc, option->option_str,
                                               err);
              if (!had_err) {
                argnum++;
                if (gt_parse_long(&long_value, argv[argnum])) {
                  gt_error_set(err, "argument to option \"-%s\" must be an "
                                    "integer", gt_str_get(option->option_str));
                  had_err = -1;
                }
              }
              if (!had_err) {
                *(GtWord*) option->value = long_value;
                option_parsed = true;
              }
              break;
            case OPTION_ULONG:
              if (optional_arg(option, argnum, argc, argv)) {
                option_parsed = true;
                break;
              }
              had_err = check_missing_argument_and_minus_sign(argnum, argc,
                                                              option
                                                              ->option_str,
                                                              argv, err);
              if (!had_err) {
                argnum++;
                if (gt_parse_ulong(&ulong_value, argv[argnum])) {
                  gt_error_set(err, "argument to option \"-%s\" is out of "
                               "range", gt_str_get(option->option_str));
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* minimum value check */
                if (option->min_value_set &&
                    ulong_value < option->min_value.ul) {
                  gt_error_set(err, "argument to option \"-%s\" must be an "
                               "integer >= "GT_WU"",
                               gt_str_get(option->option_str),
                               option->min_value.ul);
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* maximum value check */
                if (option->max_value_set &&
                    ulong_value > option->max_value.ul) {
                  gt_error_set(err, "argument to option \"-%s\" must be an "
                               "integer <= "GT_WU"",
                               gt_str_get(option->option_str),
                               option->max_value.ul);
                  had_err = -1;
                }
              }
              if (!had_err) {
                *(GtUword*) option->value = ulong_value;
                option_parsed = true;
              }
              break;
            case OPTION_RANGE:
              if (optional_arg(option, argnum, argc, argv)) {
                option_parsed = true;
                break;
              }
              /* parse first argument */
              had_err = check_missing_argument_and_minus_sign(argnum, argc,
                                                              option
                                                              ->option_str,
                                                              argv, err);
              if (!had_err) {
                argnum++;
                if (gt_parse_long(&long_value, argv[argnum]) ||
                    long_value < 0) {
                  gt_error_set(err, "first argument to option \"-%s\" must be "
                               "a non-negative integer",
                               gt_str_get(option->option_str));
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* minimum value check */
                if (option->min_value_set &&
                    long_value < option->min_value.ul) {
                  gt_error_set(err, "first argument to option \"-%s\" must be "
                               "an integer >= "GT_WU"",
                               gt_str_get(option->option_str),
                               option->min_value.ul);
                  had_err = -1;
                }
                else
                  ((GtRange*) option->value)->start = long_value;
              }
              /* parse second argument */
              if (!had_err) {
                had_err = check_missing_argument_and_minus_sign(argnum, argc,
                                                                option
                                                                ->option_str,
                                                                argv, err);
              }
              if (!had_err) {
                argnum++;
                if (gt_parse_long(&long_value, argv[argnum]) ||
                    long_value < 0) {
                  gt_error_set(err, "second argument to option \"-%s\" must be "
                               "a non-negative integer",
                               gt_str_get(option->option_str));
                  had_err = -1;
                }
              }
              if (!had_err) {
                /* maximum value check */
                if (option->max_value_set &&
                    long_value > option->max_value.ul) {
                  gt_error_set(err, "second argument to option \"-%s\" must be "
                               "an integer <= "GT_WU"",
                               gt_str_get(option->option_str),
                               option->max_value.ul);
                  had_err = -1;
                }
                else
                  ((GtRange*) option->value)->end = long_value;
              }
              /* check arguments */
              if (!had_err && (((GtRange*) option->value)->start >
                               ((GtRange*) option->value)->end)) {
                gt_error_set(err,
                             "first argument "GT_WU" to option \"-%s\" must "
                             "be <= than second argument "GT_WU"",
                             ((GtRange*) option->value)->start,
                             gt_str_get(option->option_str),
                             ((GtRange*) option->value)->end);
                had_err = -1;
              }
              if (!had_err)
                option_parsed = true;
              break;
            case OPTION_FILENAME:
            case OPTION_STRING:
              if (optional_arg(option, argnum, argc, argv)) {
                option_parsed = true;
                break;
              }
              had_err = check_missing_argument_and_minus_sign(argnum, argc,
                                                              option
                                                              ->option_str,
                                                              argv, err);
              if (!had_err) {
                argnum++;
                gt_str_set(option->value, argv[argnum]);
                option_parsed = true;
              }
              break;
            case OPTION_STRING_ARRAY:
              if (optional_arg(option, argnum, argc, argv)) {
                option_parsed = true;
                break;
              }
              had_err = check_missing_argument_and_minus_sign(argnum, argc,
                                                             option->option_str,
                                                             argv, err);
              while (!had_err) {
                if (argnum + 1 < argc && argv[argnum+1][0] != '-') {
                  argnum++;
                  gt_str_array_add_cstr(option->value, argv[argnum]);
                  option_parsed = true;
                }
                else {
                  if (!option_parsed) {
                    gt_error_set(err, "missing argument to option \"-%s\"",
                                 gt_str_get(option->option_str));
                    had_err = -1;
                  }
                  break;
                }
              }
              break;
            case OPTION_VERSION:
              if (op->version_func)
                op->version_func(op->progname);
              else
                version_func(op->progname);
              return GT_OPTION_PARSER_REQUESTS_EXIT;
            default: gt_assert(0);
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
    gt_assert(!had_err);
    gt_error_set(err, "unknown option: %s (-help shows possible options)",
                 argv[argnum]);
    had_err = -1;
    break;
  }

  /* skip "--" if necessary */
  if (argnum < argc && !strcmp(argv[argnum], "--"))
    argnum++;

  /* check for minimum number of additional arguments, if necessary */
  if (!had_err && op->min_additional_arguments != GT_UNDEF_UINT &&
      argc - argnum < op->min_additional_arguments) {
    gt_error_set(err, "missing argument\nUsage: %s %s", op->progname,
                 op->synopsis);
    had_err = -1;
  }

  /* check for maximal number of additional arguments, if necessary */
  if (!had_err && op->max_additional_arguments != GT_UNDEF_UINT &&
      argc - argnum > op->max_additional_arguments) {
    gt_error_set(err, "superfluous argument \"%s\"\nUsage: %s %s",
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
  for (i = 0; !had_err && i < gt_array_size(op->hooks); i++) {
    hookinfo = gt_array_get(op->hooks, i);
    had_err = hookinfo->hook(hookinfo->data, err);
  }

  if (parsed_args)
    *parsed_args = argnum;

  if (had_err)
    return GT_OPTION_PARSER_ERROR;
  return GT_OPTION_PARSER_OK;
}

void gt_option_parser_delete(GtOptionParser *op)
{
  GtUword i;
  if (!op) return;
  gt_free(op->progname);
  gt_free(op->synopsis);
  gt_free(op->one_liner);
  for (i = 0; i < gt_array_size(op->options); i++)
    gt_option_delete(*(GtOption**) gt_array_get(op->options, i));
  gt_array_delete(op->options);
  gt_array_delete(op->hooks);
  gt_hashmap_delete(op->optionindex);
  gt_free(op);
}

GtOption* gt_option_new_debug(bool *value)
{
  GtOption *o = gt_option_new_bool("debug", "enable debugging output", value,
                                   false);
  o->is_development_option = true;
  return o;
}

GtOption* gt_option_new_verbose(bool *value)
{
  return gt_option_new_bool("v", "be verbose", value, false);
}

GtOption* gt_option_new_width(GtUword *value)
{
  return gt_option_new_uword("width", "set output width for FASTA sequence "
                             "printing\n(0 disables formatting)", value, 0);
}

GtOption* gt_option_new_bool(const char *option_str, const char *description,
                             bool *value, bool default_value)
{
  GtOption *o = gt_option_new(option_str, description, value);
  o->option_type = OPTION_BOOL;
  o->default_value.b = default_value;
  *value = default_value;
  return o;
}

GtOption* gt_option_new_double(const char *option_str, const char *description,
                               double *value, double default_value)
{
  GtOption *o = gt_option_new(option_str, description, value);
  o->option_type = OPTION_DOUBLE;
  o->default_value.d = default_value;
  *value = default_value;
  return o;
}

GtOption *gt_option_new_double_min(const char *option_str,
                                   const char *description,
                                   double *value, double default_value,
                                   double min_value)
{
  GtOption *o = gt_option_new_double(option_str, description, value,
                                     default_value);
  o->min_value_set = true;
  o->min_value.d = min_value;
  return o;
}

GtOption *gt_option_new_double_min_max(const char *option_str,
                                       const char *description, double *value,
                                       double default_value, double min_value,
                                       double max_value)
{
  GtOption *o = gt_option_new_double(option_str, description, value,
                                     default_value);
  o->min_value_set = true;
  o->min_value.d = min_value;
  o->max_value_set = true;
  o->max_value.d = max_value;
  return o;
}

GtOption* gt_option_new_probability(const char *option_str,
                                    const char *description,
                                    double *value, double default_value)
{
  return gt_option_new_double_min_max(option_str, description, value,
                                   default_value, 0.0, 1.0);
}

GtOption* gt_option_new_int(const char *option_str, const char *description,
                            int *value, int default_value)
{
  GtOption *o = gt_option_new(option_str, description, value);
  o->option_type = OPTION_INT;
  o->default_value.i = default_value;
  *value = default_value;
  return o;
}

GtOption* gt_option_new_int_min(const char *option_str, const char *description,
                                int *value, int default_value, int min_value)
{
  GtOption *o = gt_option_new_int(option_str, description, value,
                                  default_value);
  o->min_value_set = true;
  o->min_value.i = min_value;
  return o;
}

GtOption* gt_option_new_int_max(const char *option_str, const char *description,
                                int *value, int default_value, int max_value)
{
  GtOption *o = gt_option_new_int(option_str, description, value,
                                  default_value);
  o->max_value_set = true;
  o->max_value.i = max_value;
  return o;
}

GtOption* gt_option_new_int_min_max(const char *option_str,
                                    const char *description,
                                    int *value, int default_value,
                                    int min_value, int max_value)
{
  GtOption *o = gt_option_new_int(option_str, description, value,
                                  default_value);
  o->min_value_set = true;
  o->min_value.i = min_value;
  o->max_value_set = true;
  o->max_value.i = max_value;
  return o;
}

GtOption* gt_option_new_uint(const char *option_str, const char *description,
                             unsigned int *value, unsigned int default_value)
{
  GtOption *o = gt_option_new(option_str, description, value);
  o->option_type = OPTION_UINT;
  o->default_value.ui = default_value;
  *value = default_value;
  return o;
}

GtOption* gt_option_new_uint_min(const char *option_str,
                                 const char *description,
                                 unsigned int *value,
                                 unsigned int default_value,
                                 unsigned int min_value)
{
  GtOption *o = gt_option_new_uint(option_str, description, value,
                                   default_value);
  o->min_value_set = true;
  o->min_value.ui = min_value;
  return o;
}

GtOption* gt_option_new_uint_max(const char *option_str,
                                 const char *description,
                                 unsigned int *value,
                                 unsigned int default_value,
                                 unsigned int max_value)
{
  GtOption *o = gt_option_new_uint(option_str, description, value,
                                   default_value);
  o->max_value_set = true;
  o->max_value.ui = max_value;
  return o;
}

GtOption *gt_option_new_uint_min_max(const char *option_str,
                                     const char *description,
                                     unsigned int *value,
                                     unsigned int default_value,
                                     unsigned int min_value,
                                     unsigned int max_value)
{
  GtOption *o = gt_option_new_uint(option_str, description, value,
                                   default_value);
  o->min_value_set = true;
  o->min_value.i = min_value;
  o->max_value_set = true;
  o->max_value.i = max_value;
  return o;
}

GtOption* gt_option_new_word(const char *option_str, const char *description,
                             GtWord *value, GtWord default_value)
{
  GtOption *o = gt_option_new(option_str, description, value);
  o->option_type = OPTION_LONG;
  o->default_value.l = default_value;
  *value = default_value;
  return o;
}

GtOption* gt_option_new_uword(const char *option_str, const char *description,
                              GtUword *value, GtUword default_value)
{
  GtOption *o = gt_option_new(option_str, description, value);
  o->option_type = OPTION_ULONG;
  o->default_value.ul = default_value;
  *value = default_value;
  return o;
}

GtOption* gt_option_new_uword_min(const char *option_str,
                                  const char *description,
                                  GtUword *value,
                                  GtUword default_value,
                                  GtUword min_value)
{
  GtOption *o = gt_option_new_uword(option_str, description, value,
                                    default_value);
  o->min_value_set = true;
  o->min_value.ul = min_value;
  return o;
}

GtOption *gt_option_new_uword_min_max(const char *option_str,
                                      const char *description,
                                      GtUword *value,
                                      GtUword default_value,
                                      GtUword min_value,
                                      GtUword max_value)
{
  GtOption *o = gt_option_new_uword(option_str, description, value,
                                    default_value);
  o->min_value_set = true;
  o->min_value.ul = min_value;
  o->max_value_set = true;
  o->max_value.ul = max_value;
  return o;
}

GtOption *gt_option_new_long(const char *option_str,
                             const char *description,
                             GtWord *value,
                             GtWord default_value)
{
  return gt_option_new_word(option_str, description, value,
                            default_value);
}

GtOption *gt_option_new_ulong(const char *option_str,
                              const char *description,
                              GtUword *value,
                              GtUword default_value)
{
  return gt_option_new_uword(option_str, description, value,
                             default_value);
}

GtOption *gt_option_new_ulong_min(const char *option_str,
                                  const char *description,
                                  GtUword *value,
                                  GtUword default_value,
                                  GtUword min_value)
{
  return gt_option_new_uword_min(option_str, description, value,
                                 default_value, min_value);
}

GtOption *gt_option_new_ulong_min_max(const char *option_str,
                                      const char *description,
                                      GtUword *value,
                                      GtUword default_value,
                                      GtUword min_value,
                                      GtUword max_value)
{
  return gt_option_new_uword_min_max(option_str, description, value,
                                     default_value, min_value,
                                     max_value);
}

GtOption* gt_option_new_range(const char *option_str, const char *description,
                              GtRange *value, GtRange *default_value)
{
  GtOption *o = gt_option_new(option_str, description, value);
  o->option_type = OPTION_RANGE;
  o->default_value.r.start = default_value ? default_value->start
                                           : GT_UNDEF_UWORD;
  o->default_value.r.end   = default_value ? default_value->end
                                           : GT_UNDEF_UWORD;
  value->start = o->default_value.r.start;
  value->end   = o->default_value.r.end;
  return o;
}

GtOption* gt_option_new_range_min_max(const char *option_str,
                                      const char *description, GtRange *value,
                                      GtRange *default_value,
                                      GtUword min_value,
                                      GtUword max_value)
{
   GtOption *o = gt_option_new_range(option_str, description, value,
                                     default_value);
   o->min_value_set = true;
   o->min_value.ul = min_value;
   o->max_value_set = true;
   o->max_value.ul = max_value;
   return o;
}

GtOption* gt_option_new_string(const char *option_str, const char *description,
                               GtStr *value, const char *default_value)
{
  GtOption *o = gt_option_new(option_str, description, value);
  o->option_type = OPTION_STRING;
  o->default_value.s = default_value;
  gt_str_set(value, default_value);
  return o;
}

GtOption* gt_option_new_string_array(const char *option_str,
                                     const char *description, GtStrArray *value)
{
  GtOption *o = gt_option_new(option_str, description, value);
  o->option_type = OPTION_STRING_ARRAY;
  return o;
}

/* the following function would allow to handle files differently from strings
   later on (e.g., for CGI scripts), but for now the are implemented in the same
   way */
GtOption* gt_option_new_filename(const char *option_str,
                                 const char *description,
                                 GtStr *filename)
{
  GtOption *o = gt_option_new(option_str, description, filename);
  o->option_type = OPTION_FILENAME;
  o->default_value.s = NULL;
  gt_str_reset(filename);
  return o;
}

/* the following function would allow to handle file arrays differently from
   string arrays later on (e.g., for CGI scripts) , but for now the are
   implemented in the same way */
GtOption* gt_option_new_filename_array(const char *option_str,
                                       const char *description,
                                       GtStrArray *filenames)
{
  return gt_option_new_string_array(option_str, description, filenames);
}

GtOption* gt_option_new_choice(const char *option_str, const char *description,
                               GtStr *value, const char *default_value,
                               const char **domain)
{
  GtOption *o;
#ifndef NDEBUG
  GtUword in_domain = 1;
  if (default_value) {
    while (domain[in_domain - 1] != NULL) {
      if (domain[in_domain - 1] == default_value) {
        in_domain = 0;
        break;
      }
      in_domain++;
    }
  }
  else
    in_domain = 0;
  gt_assert(!in_domain);
#endif

  o = gt_option_new_string(option_str, description, value, default_value);
  o->option_type = OPTION_CHOICE;
  o->domain = domain;

  return o;
}

const char* gt_option_get_name(const GtOption *o)
{
  gt_assert(o);
  return gt_str_get(o->option_str);
}

void gt_option_is_mandatory(GtOption *o)
{
  gt_assert(o);
  o->is_mandatory = true;
}

void gt_option_is_mandatory_either(GtOption *o, const GtOption *meo)
{
  gt_assert(o && meo);
  gt_assert(!o->mandatory_either_options);
  o->mandatory_either_options = gt_array_new(sizeof (GtOption*));
  gt_array_add(o->mandatory_either_options, meo);
}

void gt_option_is_mandatory_either_3(GtOption *o, const GtOption *meo_a,
                                     const GtOption *meo_b)
{
  gt_assert(o && meo_a && meo_b);
  gt_assert(!o->mandatory_either_options);
  o->mandatory_either_options = gt_array_new(sizeof (GtOption*));
  gt_array_add(o->mandatory_either_options, meo_a);
  gt_array_add(o->mandatory_either_options, meo_b);
}

void gt_option_is_mandatory_either_4(GtOption *o, const GtOption *meo_a,
                                     const GtOption *meo_b,
                                     const GtOption *meo_c)
{
  gt_assert(o && meo_a && meo_b);
  gt_assert(!o->mandatory_either_options);
  o->mandatory_either_options = gt_array_new(sizeof (GtOption*));
  gt_array_add(o->mandatory_either_options, meo_a);
  gt_array_add(o->mandatory_either_options, meo_b);
  gt_array_add(o->mandatory_either_options, meo_c);
}

void gt_option_is_extended_option(GtOption *o)
{
  gt_assert(o);
  o->is_extended_option = true;
}

void gt_option_is_development_option(GtOption *o)
{
  gt_assert(o);
  o->is_development_option = true;
}

void gt_option_imply(GtOption *o, const GtOption *implied_option)
{
  GtArray *option_array;
  gt_assert(o && implied_option);
  if (!o->implications)
    o->implications = gt_array_new(sizeof (GtArray*));
  option_array = gt_array_new(sizeof (GtOption*));
  gt_array_add(option_array, implied_option);
  gt_array_add(o->implications, option_array);
}

void gt_option_imply_either_2(GtOption *o, const GtOption *io1,
                              const GtOption *io2)
{
  GtArray *option_array;
  gt_assert(o && io1 && io2);
  if (!o->implications)
    o->implications = gt_array_new(sizeof (GtArray*));
  option_array = gt_array_new(sizeof (GtOption*));
  gt_array_add(option_array, io1);
  gt_array_add(option_array, io2);
  gt_array_add(o->implications, option_array);
}

void gt_option_imply_either_3(GtOption *o, const GtOption *io1,
                              const GtOption *io2, const GtOption *io3)
{
  GtArray *option_array;
  gt_assert(o && io1 && io2);
  if (!o->implications)
    o->implications = gt_array_new(sizeof (GtArray*));
  option_array = gt_array_new(sizeof (GtOption*));
  gt_array_add(option_array, io1);
  gt_array_add(option_array, io2);
  gt_array_add(option_array, io3);
  gt_array_add(o->implications, option_array);
}

void gt_option_exclude(GtOption *o_a, GtOption *o_b)
{
  gt_assert(o_a && o_b);
  if (!o_a->exclusions)
    o_a->exclusions = gt_array_new(sizeof (GtOption*));
  if (!o_b->exclusions)
    o_b->exclusions = gt_array_new(sizeof (GtOption*));
  gt_array_add(o_a->exclusions, o_b);
  gt_array_add(o_b->exclusions, o_a);
}

void gt_option_hide_default(GtOption *o)
{
  gt_assert(o);
  o->hide_default = true;
}

void gt_option_argument_is_optional(GtOption *o)
{
  gt_assert(o);
  o->argument_is_optional = true;
}

bool gt_option_is_set(const GtOption *o)
{
  gt_assert(o);
  return o->is_set;
}

void gt_option_delete(GtOption *o)
{
  GtUword i;
  if (!o) return;
  if (o->reference_count) {
    o->reference_count--;
    return;
  }
  gt_str_delete(o->option_str);
  gt_str_delete(o->description);
  for (i = 0; i < gt_array_size(o->implications); i++)
    gt_array_delete(*(GtArray**) gt_array_get(o->implications, i));
  gt_array_delete(o->implications);
  gt_array_delete(o->exclusions);
  gt_array_delete(o->mandatory_either_options);
  gt_free(o);
}

int gt_option_parse_spacespec(GtUword *maximumspace,
                              const char *optname,
                              const GtStr *memlimit,
                              GtError *err)
{
  int had_err = 0;
  bool match = false;

  had_err = gt_grep(&match, "^[0-9]+(MB|GB)$", gt_str_get(memlimit), err);
  if (had_err || !match)
  {
    gt_error_set(err,"option -%s must have one positive "
                     "integer argument followed by one of "
                     "the keywords MB and GB",optname);
    had_err = -1;
  }
  if (!had_err)
  {
    int readint;
    char buffer[2+1];

    (void) sscanf(gt_str_get(memlimit), "%d%s", &readint, buffer);
    *maximumspace = (GtUword) readint;
    if (strcmp(buffer, "GB") == 0)
    {
      if (sizeof (GtUword) == (size_t) 4 && *maximumspace > 3UL)
      {
        gt_error_set(err,"option -%s: for 32bit binaries one cannot "
                         "specify more than 3 GB as maximum space",optname);
        had_err = -1;
      }
      if (had_err != 1)
      {
        *maximumspace <<= 30;
      }
    } else
    {
      if (strcmp(buffer, "MB") == 0)
      {
        if (sizeof (GtUword) == (size_t) 4 && *maximumspace > 4095UL)
        {
          gt_error_set(err,"option -%s: for 32bit binaries one cannot "
                           "specify more than 4095 MB as maximum space",
                           optname);
          had_err = -1;
        }
        if (!had_err)
        {
          *maximumspace <<= 20;
        }
      }
    }
  }
  return had_err;
}

const char *gt_option_parser_synopsis(const GtOptionParser *option_parser)
{
  gt_assert(option_parser);
  return option_parser->synopsis;
}

const char *gt_option_parser_one_liner(const GtOptionParser *option_parser)
{
  gt_assert(option_parser);
  return option_parser->one_liner;
}
