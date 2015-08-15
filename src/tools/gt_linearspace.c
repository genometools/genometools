/*
  Copyright (c) 2015 Annika <annika.seidel@studium.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/fa.h"
#include "core/ma.h"
#include "core/str.h"
#include "core/str_api.h"
#include "core/str_array.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "extended/linearalign_affinegapcost.h"
#include "extended/linearalign.h"
#include "tools/gt_linearspace.h"

typedef struct {
  GtStr *outputfile;
  GtStrArray *strings,
             *linearcosts,
             *affinecosts;
  bool global,
       local,
       show;
} GtLinearspaceArguments;

static void* gt_linearspace_arguments_new(void)
{
  GtLinearspaceArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->outputfile = gt_str_new();
  arguments->strings = gt_str_array_new();
  arguments->linearcosts = gt_str_array_new();
  arguments->affinecosts = gt_str_array_new();
  return arguments;
}

static void gt_linearspace_arguments_delete(void *tool_arguments)
{
  GtLinearspaceArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->outputfile);
    gt_str_array_delete(arguments->strings);
    gt_str_array_delete(arguments->linearcosts);
    gt_str_array_delete(arguments->affinecosts);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_linearspace_option_parser_new(void *tool_arguments)
{
  GtLinearspaceArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *optionstrings, *optionglobal, *optionlocal,
  *optionlinearcosts, *optionaffinecosts,*optionshow,*optionoutputfile;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("options",
                            "Apply function to compute alignment.");
  gt_option_parser_set_mail_address(op,
                          "<annika.seidel@studium.uni-hamburg.de>");

  /* -bool */
  optionglobal = gt_option_new_bool("global", "global alignment",
                              &arguments->global, false);
  gt_option_parser_add_option(op, optionglobal);

  optionlocal = gt_option_new_bool("local", "local alignment",
                              &arguments->local, false);
  gt_option_parser_add_option(op, optionlocal);

  optionshow = gt_option_new_bool("show", "show alignment",
                              &arguments->show, false);
  gt_option_parser_add_option(op, optionshow);

  /* -str */
  optionstrings = gt_option_new_string_array("ss", "use two strings",
                                             arguments->strings);
  gt_option_parser_add_option(op, optionstrings);

  optionlinearcosts = gt_option_new_string_array("l", "lineargapcosts, "
                                                 "use three values",
                                                arguments->linearcosts);
  gt_option_parser_add_option(op, optionlinearcosts);

  optionaffinecosts = gt_option_new_string_array("a", "affinegapcosts, "
                                                 "use four values",
                                                 arguments->affinecosts);
  gt_option_parser_add_option(op, optionaffinecosts);

  optionoutputfile = gt_option_new_string("o", "use outputfile",
                                          arguments->outputfile, "stdout");
  gt_option_parser_add_option(op, optionoutputfile);

  /* dependencies*/
  gt_option_is_mandatory(optionstrings);
  gt_option_exclude(optionlocal, optionglobal);
  gt_option_exclude(optionlinearcosts, optionaffinecosts);
  gt_option_imply(optionoutputfile, optionshow);
  gt_option_imply_either_2(optionstrings, optionglobal, optionlocal);
  gt_option_imply_either_2(optionlocal, optionlinearcosts, optionaffinecosts);
  gt_option_imply_either_2(optionglobal, optionlinearcosts, optionaffinecosts);

  return op;
}

static int gt_linearspace_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtLinearspaceArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if ((gt_str_array_size(arguments->strings) > 0) &&
     (gt_str_array_size(arguments->strings) != 2UL))
  {
    gt_error_set(err, "option -ss requires two string arguments");
    had_err = 1;
  }
  if ((gt_str_array_size(arguments->linearcosts) > 0) &&
     (gt_str_array_size(arguments->linearcosts) != 3UL))
  {
    gt_error_set(err, "option -l requires "
                      "match, mismatch, gap costs/scores");
    had_err = 1;
  }
  if ((gt_str_array_size(arguments->affinecosts) > 0) &&
     (gt_str_array_size(arguments->affinecosts) != 4UL))
  {
    gt_error_set(err, "option -a requires match, mismatch, "
                      "gap_opening, gap_extending costs/scores");
    had_err = 1;
  }

  return had_err;
}

static void show(const GtLinearspaceArguments *arguments,
                 const GtAlignment *align, FILE *fp)
{
  if (gt_str_array_size(arguments->strings) > 0)
  {
    if (fp != NULL)
    {
      fprintf(fp,"# two strings \"%s\" \"%s\"\n",
              gt_str_array_get(arguments->strings,0),
              gt_str_array_get(arguments->strings,1UL));
      gt_alignment_show(align, fp, 80);
    }
  }
}

static GtWord* select_costs(const GtStrArray *arr,
                            GT_UNUSED GtError *err)
{
  GtWord *evalues, check;
  GtUword size, i;

  size = gt_str_array_size(arr);
  evalues = malloc(sizeof(*evalues)*size);
  for (i = 0; i < size; i++)
  {
    check = sscanf(gt_str_array_get(arr,i),GT_WD, &evalues[i]);
    if (check != 1)
    {
      gt_error_set(err, "find invalid cost or score");
      gt_free(evalues);
      return NULL;
    }
  }

  return evalues;
}

static int gt_linearspace_runner(GT_UNUSED int argc,
                                 GT_UNUSED const char **argv,
                                 GT_UNUSED int parsed_args,
                                 void *tool_arguments,
                                 GtError *err)
{
  GtLinearspaceArguments *arguments = tool_arguments;
  int had_err = 0;
  GtWord *linearcosts, *affinecosts;
  GtAlignment *align = NULL;
  FILE *fp;

  gt_error_check(err);
  gt_assert(arguments);

  /* linear gap costs */
  if (gt_str_array_size(arguments->linearcosts) > 0)
  {
    gt_assert(gt_str_array_size(arguments->linearcosts) == 3UL);
    linearcosts = select_costs(arguments->linearcosts, err);
    if (linearcosts == NULL)
      return 1;

    if (arguments->global)
    {
      align = gt_computelinearspace(
              (const GtUchar *) gt_str_array_get(arguments->strings,0),0,
              (GtUword) strlen(gt_str_array_get(arguments->strings,0)),
              (const GtUchar *) gt_str_array_get(arguments->strings,1UL),0,
              (GtUword) strlen(gt_str_array_get(arguments->strings,1UL)),
              linearcosts[0],linearcosts[1],linearcosts[2]);
    }
    else if (arguments->local)
    {
      align = gt_computelinearspace_local(
              (const GtUchar *) gt_str_array_get(arguments->strings,0),0,
              (GtUword) strlen(gt_str_array_get(arguments->strings,0)),
              (const GtUchar *) gt_str_array_get(arguments->strings,1UL),0,
              (GtUword) strlen(gt_str_array_get(arguments->strings,1UL)),
              linearcosts[0],linearcosts[1],linearcosts[2]);
    }
    gt_free(linearcosts);
  }/* affine gap costs */
  else if (gt_str_array_size(arguments->affinecosts) > 0)
  {
    gt_assert(gt_str_array_size(arguments->affinecosts) == 4UL);
    affinecosts = select_costs(arguments->affinecosts, err);
    if (affinecosts == NULL)
      return 1;

    if (arguments->global)
    {
      align = gt_computeaffinelinearspace(
              (const GtUchar *) gt_str_array_get(arguments->strings,0),0,
              (GtUword) strlen(gt_str_array_get(arguments->strings,0)),
              (const GtUchar *) gt_str_array_get(arguments->strings,1UL),0,
              (GtUword) strlen(gt_str_array_get(arguments->strings,1UL)),
              affinecosts[0],affinecosts[1],affinecosts[2], affinecosts[3]);
    }
    else if (arguments->local)
    {
      align = gt_computeaffinelinearspace_local(
              (const GtUchar *) gt_str_array_get(arguments->strings,0),0,
              (GtUword) strlen(gt_str_array_get(arguments->strings,0)),
              (const GtUchar *) gt_str_array_get(arguments->strings,1UL),0,
              (GtUword) strlen(gt_str_array_get(arguments->strings,1UL)),
              affinecosts[0],affinecosts[1],affinecosts[2], affinecosts[3]);
    }
     gt_free(affinecosts);
  }

  /* show */
  if (arguments->show)
  {
    if (!strcmp(gt_str_get(arguments->outputfile),"stdout"))
      show(arguments,align,stdout);
    else{
      fp = gt_fa_fopen_func(gt_str_get(arguments->outputfile),
                            "a", __FILE__,__LINE__,err);
      gt_error_check(err);
      show(arguments,align,fp);
      gt_fa_fclose(fp);
    }
  }
  gt_alignment_delete(align);

  return had_err;
}

GtTool* gt_linearspace(void)
{
  return gt_tool_new(gt_linearspace_arguments_new,
                     gt_linearspace_arguments_delete,
                     gt_linearspace_option_parser_new,
                     gt_linearspace_arguments_check,
                     gt_linearspace_runner);
}
