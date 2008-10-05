/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "core/cstr_array.h"
#include "core/error.h"
#include "core/ma_api.h"
#include "core/option.h"
#include "core/str.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "extended/toolbox.h"
#include "tools/gt_tallymer.h"
#include "match/tyr-mkindex.h"
#include "match/verbose-def.h"

typedef enum
{
  Autoprefixlength,
  Undeterminedprefixlength,
  Determinedprefixlength
} Prefixlengthflag;

typedef struct
{
  unsigned long value;
  Prefixlengthflag flag;
} Prefixlengthvalue;

typedef struct
{
  unsigned long searchlength,
                userdefinedminocc,
                userdefinedmaxocc,
                userdefinedprefixlength;
  Prefixlengthvalue prefixlength;
  GtOption *refoptionpl;
  GtStr *str_storeindex,
        *str_inputindex;
  bool storecounts,
       performtest,
       verbose;
} Tallymer_mkindex_options;

static void *gt_tallymer_mkindex_arguments_new(void)
{
  Tallymer_mkindex_options *arguments
    = gt_malloc(sizeof (Tallymer_mkindex_options));
  arguments->str_storeindex = gt_str_new();
  arguments->str_inputindex = gt_str_new();
  return arguments;
}

static void gt_tallymer_mkindex_arguments_delete(void *tool_arguments)
{
  Tallymer_mkindex_options *arguments = tool_arguments;

  if (!arguments)
  {
    return;
  }
  gt_str_delete(arguments->str_storeindex);
  gt_str_delete(arguments->str_inputindex);
  gt_option_delete(arguments->refoptionpl);
  gt_free(arguments);
}

static GtOptionParser
            *gt_tallymer_mkindex_option_parser_new(void *tool_arguments)
{
  GtOptionParser *op;
  GtOption *option,
           *optionminocc,
           *optionmaxocc,
           *optionpl,
           *optionstoreindex,
           *optionstorecounts;
  Tallymer_mkindex_options *arguments = tool_arguments;

  op = gt_option_parser_new("[options] enhanced-suffix-array",
                            "Count and index k-mers in the given enhanced "
                            "suffix array for a fixed value of k.");
  gt_option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");

  option = gt_option_new_ulong("mersize",
                               "Specify the mer size.",
                               &arguments->searchlength,
                               20UL);
  gt_option_parser_add_option(op, option);

  optionminocc
    = gt_option_new_ulong("minocc",
                          "Specify the minimum occurrence number for "
                          "the mers to output/index",
                          &arguments->userdefinedminocc,
                          0);
  gt_option_parser_add_option(op, optionminocc);

  optionmaxocc
    = gt_option_new_ulong("maxocc",
                          "Specify the maximum occurrence number for "
                          "the mers to output/index",
                          &arguments->userdefinedmaxocc,
                          0);
  gt_option_parser_add_option(op, optionmaxocc);

  optionpl = gt_option_new_ulong_min("pl",
                 "specify prefix length for bucket boundary construction\n"
                 "recommendation: use without argument;\n"
                 "then a reasonable prefix length is automatically determined",
                 &arguments->userdefinedprefixlength,
                 0,
                 1UL);
  gt_option_argument_is_optional(optionpl);
  gt_option_parser_add_option(op, optionpl);
  arguments->refoptionpl = gt_option_ref(optionpl);

  optionstoreindex = gt_option_new_string("indexname",
                                          "store the mers specified by options "
                                          "-maxocc and -minocc in an index",
                                          arguments->str_storeindex, NULL);
  gt_option_parser_add_option(op, optionstoreindex);

  optionstorecounts = gt_option_new_bool("counts", "store counts of the mers",
                                         &arguments->storecounts,false);
  gt_option_parser_add_option(op, optionstorecounts);

  option = gt_option_new_bool("test", "perform tests to verify program "
                                      "correctness", &arguments->performtest,
                                      false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  gt_option_imply(optionpl, optionstoreindex);
  gt_option_imply(optionstorecounts, optionstoreindex);
  gt_option_imply_either_2(optionstoreindex,optionminocc,optionmaxocc);
  return op;
}

static int gt_tallymer_mkindex_arguments_check(int rest_argc,
                                               void *tool_arguments,
                                               GtError *err)
{
  Tallymer_mkindex_options *arguments = tool_arguments;

  if (gt_option_is_set(arguments->refoptionpl))
  {
    if (arguments->userdefinedprefixlength == 0)
    {
      arguments->prefixlength.flag = Autoprefixlength;
      arguments->prefixlength.value = 0;
    } else
    {
      arguments->prefixlength.flag = Determinedprefixlength;
      arguments->prefixlength.value = arguments->userdefinedprefixlength;
    }
  } else
  {
    arguments->prefixlength.flag = Undeterminedprefixlength;
    arguments->prefixlength.value = 0;
  }
  if (rest_argc != 1)
  {
    gt_error_set(err,"missing name of enhanced suffix array index");
    return -1;
  }
  return 0;
}

static int gt_tallymer_mkindex_runner(int argc,
                                      const char **argv,
                                      int parsed_args,
                                      void *tool_arguments,
                                      GtError *err)
{
  Tallymer_mkindex_options *arguments = tool_arguments;
  Verboseinfo *verboseinfo;
  bool haserr = false;

  assert(parsed_args + 1 == argc);
  gt_str_set(arguments->str_inputindex,argv[parsed_args]);
  verboseinfo = newverboseinfo(arguments->verbose);
  if (arguments->verbose)
  {
    printf("# mersize=%lu\n",arguments->searchlength);
    if (arguments->userdefinedminocc > 0)
    {
      printf("# minocc=%lu\n",arguments->userdefinedminocc);
    } else
    {
      printf("# minocc=undefined\n");
    }
    if (arguments->userdefinedmaxocc > 0)
    {
      printf("# maxocc=%lu\n",arguments->userdefinedmaxocc);
    } else
    {
      printf("# maxocc=undefined\n");
    }
    printf("# prefixlength=");
    if (arguments->prefixlength.flag == Autoprefixlength)
    {
      printf("automatic");
    } else
    {
      if (arguments->prefixlength.flag == Determinedprefixlength)
      {
        printf("%lu",arguments->prefixlength.value);
      } else
      {
        printf("undefined");
      }
    }
    printf("\n");
    if (gt_str_length(arguments->str_storeindex) > 0)
    {
      printf("# storeindex=%s\n",gt_str_get(arguments->str_storeindex));
    }
    printf("# inputindex=%s\n",gt_str_get(arguments->str_inputindex));
  }
  if (merstatistics(arguments->str_inputindex,
                    arguments->searchlength,
                    arguments->userdefinedminocc,
                    arguments->userdefinedmaxocc,
                    arguments->str_storeindex,
                    true,
                    arguments->storecounts,
                    arguments->performtest,
                    verboseinfo,
                    err) != 0)
  {
    haserr = true;
  }
  freeverboseinfo(&verboseinfo);
  return haserr ? - 1 : 0;
}

static GtTool* gt_tallymer_mkindex(void)
{
  return gt_tool_new(gt_tallymer_mkindex_arguments_new,
                     gt_tallymer_mkindex_arguments_delete,
                     gt_tallymer_mkindex_option_parser_new,
                     gt_tallymer_mkindex_arguments_check,
                     gt_tallymer_mkindex_runner);
}

static int gt_tallymer_occratio(GT_UNUSED int argc,
                                GT_UNUSED const char *argv[],
                                GtError *err)
{
  gt_error_set(err,"tallymer occratio not implemented yet");
  return -1;
}

static int gt_tallymer_search(GT_UNUSED int argc,
                              GT_UNUSED const char *argv[],
                              GtError *err)
{
  gt_error_set(err,"tallymer search not implemented yet");
  return -1;
}

static void *gt_tallymer_arguments_new(void)
{
  GtToolbox *tallymer_toolbox = gt_toolbox_new();
  gt_toolbox_add_tool(tallymer_toolbox, "mkindex", gt_tallymer_mkindex());
  gt_toolbox_add(tallymer_toolbox, "occratio", gt_tallymer_occratio);
  gt_toolbox_add(tallymer_toolbox, "search", gt_tallymer_search);
  return tallymer_toolbox;
}

static void gt_tallymer_arguments_delete(void *tool_arguments)
{
  GtToolbox *index_toolbox = tool_arguments;
  if (!index_toolbox) return;
  gt_toolbox_delete(index_toolbox);
}

static GtOptionParser* gt_tallymer_option_parser_new(void *tool_arguments)
{
  GtToolbox *index_toolbox = tool_arguments;
  GtOptionParser *op;

  assert(index_toolbox != NULL);
  op = gt_option_parser_new(
                    "[option ...] [mkindex|occratio|search] [argument ...]",
                    "Call tallymer with specific tool and "
                    "pass argument(s) to it.");
  gt_option_parser_set_comment_func(op, gt_toolbox_show, index_toolbox);
  gt_option_parser_refer_to_manual(op);
  return op;
}

static int gt_tallymer_runner(int argc, const char **argv, int parsed_args,
                                 void *tool_arguments, GtError *err)
{
  GtToolbox *index_toolbox = tool_arguments;
  GtToolfunc toolfunc;
  GtTool *tool = NULL;
  char **nargv = NULL;
  int had_err = 0;

  gt_error_check(err);
  assert(index_toolbox != NULL);

  /* determine tool */
  if (!gt_toolbox_has_tool(index_toolbox, argv[parsed_args]))
  {
    gt_error_set(err, "tallymer tool '%s' not found; option -help lists "
                   "possible tools", argv[parsed_args]);
    had_err = -1;
  }

  /* call sub-tool */
  if (!had_err)
  {
    if (!(toolfunc = gt_toolbox_get(index_toolbox, argv[parsed_args])))
    {
      tool = gt_toolbox_get_tool(index_toolbox, argv[parsed_args]);
      assert(tool != NULL);
    }
    nargv = gt_cstr_array_prefix_first(argv + parsed_args,
                                       gt_error_get_progname(err));
    gt_error_set_progname(err, nargv[0]);
    if (toolfunc != NULL)
      had_err = toolfunc(argc - parsed_args, (const char**) nargv, err);
    else
      had_err = gt_tool_run(tool, argc - parsed_args, (const char**) nargv,
                            err);
  }
  gt_cstr_array_delete(nargv);
  return had_err;
}

GtTool* gt_tallymer(void)
{
  return gt_tool_new(gt_tallymer_arguments_new,
                     gt_tallymer_arguments_delete,
                     gt_tallymer_option_parser_new,
                     NULL,
                     gt_tallymer_runner);
}
