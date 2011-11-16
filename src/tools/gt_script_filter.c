/*
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/script_filter.h"
#include "tools/gt_script_filter.h"

typedef struct {
  bool oneline,
       showinfo,
       scriptname,
       validate;
} GtScriptFilterArguments;

static void* gt_script_filter_arguments_new(void)
{
  GtScriptFilterArguments *arguments = gt_calloc(1, sizeof *arguments);
  return arguments;
}

static void gt_script_filter_arguments_delete(void *tool_arguments)
{
  GtScriptFilterArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_free(arguments);
}

static GtOptionParser* gt_script_filter_option_parser_new(void *tool_arguments)
{
  GtScriptFilterArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [script file(s)]",
                            "Get info about and validate Lua script filters.");

  /* -showinfo */
  option = gt_option_new_bool("showinfo", "show information about filter",
                              &arguments->showinfo, true);
  gt_option_parser_add_option(op, option);

  /* -validate */
  option = gt_option_new_bool("validate", "validate filter function",
                              &arguments->validate, true);
  gt_option_parser_add_option(op, option);

  /* -oneline */
  option = gt_option_new_bool("oneline", "show compact information on one line",
                             &arguments->oneline, false);
  gt_option_parser_add_option(op, option);

  /* -scriptname */
  option = gt_option_new_bool("scriptname", "show script name",
                             &arguments->scriptname, true);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_min_args(op, 1);

  return op;
}

static int gt_script_filter_runner(int argc, const char **argv, int parsed_args,
                                   void *tool_arguments, GT_UNUSED GtError *err)
{
  GtScriptFilterArguments *arguments = tool_arguments;
  int had_err = 0;
  unsigned long i;

  gt_error_check(err);
  gt_assert(arguments);

  for (i = parsed_args; !had_err && i < argc; i++) {
    GtScriptFilter *sf = NULL;
    sf = gt_script_filter_new(argv[i], err);
    if (!sf) {
      had_err = -1;
      break;
    } else {
      if (arguments->showinfo) {
        if (arguments->oneline) {
          const char *name,
                     *version,
                     *author;

          if (!(name = gt_script_filter_get_name(sf, err)))
            had_err = -1;

          if (!had_err) {
            if (!(version = gt_script_filter_get_version(sf, err)))
              had_err = -1;
          }

          if (!had_err) {
            if (!(author = gt_script_filter_get_author(sf, err)))
              had_err = -1;
          }

          if (!had_err) {
            printf("%s v%s (by %s)\n", name, version, author);
          }
        } else {
          const char *out;

          if (arguments->scriptname)
            printf("script name:\t%s\n", argv[i]);

          out = gt_script_filter_get_name(sf, err);
          if (out)
            printf("filter name:\t%s\n", out);
          else
            had_err = -1;

          if (!had_err) {
            out = gt_script_filter_get_version(sf, err);
            if (out)
              printf("version:\t%s\n", out);
            else
              had_err = -1;
          }

          if (!had_err) {
            out = gt_script_filter_get_author(sf, err);
            if (out)
              printf("author:\t\t%s\n", out);
            else
              had_err = -1;
          }

          if (!had_err) {
            out = gt_script_filter_get_email(sf, err);
            if (out)
              printf("email:\t\t%s\n", out);
            else
              had_err = -1;
          }

          if (!had_err) {
            out = gt_script_filter_get_description(sf, err);
            if (out)
              printf("description:\t%s\n", out);
            else
              had_err = -1;
          }

          if (i != argc-1)
            printf("\n");
        }
      }
      if (arguments->validate) {
        GT_UNUSED bool select;
        GtStr *seqid = gt_str_new_cstr("foo");
        GtFeatureNode *fn = (GtFeatureNode*) gt_feature_node_new(seqid, "gene",
                                                             23, 42,
                                                             GT_STRAND_FORWARD);
        had_err = gt_script_filter_run(sf, fn, &select, err);
        gt_genome_node_delete((GtGenomeNode*) fn);
        gt_str_delete(seqid);
      }
      gt_script_filter_delete(sf);
    }
  }

  return had_err;
}

GtTool* gt_script_filter(void)
{
  return gt_tool_new(gt_script_filter_arguments_new,
                  gt_script_filter_arguments_delete,
                  gt_script_filter_option_parser_new,
                  NULL,
                  gt_script_filter_runner);
}
