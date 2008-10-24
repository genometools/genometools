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
#include "match/tyr-search.h"
#include "match/tyr-mersplit.h"
#include "match/verbose-def.h"
#include "match/optionargmode.h"
#include "match/defined-types.h"

typedef enum
{
  Autoprefixlength,
  Undeterminedprefixlength,
  Determinedprefixlength
} Prefixlengthflag;

typedef struct
{
  unsigned int value;
  Prefixlengthflag flag;
} Prefixlengthvalue;

typedef struct
{
  unsigned long mersize,
                userdefinedminocc,
                userdefinedmaxocc;
  unsigned int userdefinedprefixlength;
  Prefixlengthvalue prefixlength;
  GtOption *refoptionpl;
  GtStr *str_storeindex,
        *str_inputindex;
  bool storecounts,
       performtest,
       verbose;
} Tyr_mkindex_options;

static void *gt_tyr_mkindex_arguments_new(void)
{
  Tyr_mkindex_options *arguments
    = gt_malloc(sizeof (Tyr_mkindex_options));
  arguments->str_storeindex = gt_str_new();
  arguments->str_inputindex = gt_str_new();
  return arguments;
}

static void gt_tyr_mkindex_arguments_delete(void *tool_arguments)
{
  Tyr_mkindex_options *arguments = tool_arguments;

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
            *gt_tyr_mkindex_option_parser_new(void *tool_arguments)
{
  GtOptionParser *op;
  GtOption *option,
           *optionminocc,
           *optionmaxocc,
           *optionpl,
           *optionstoreindex,
           *optionstorecounts;
  Tyr_mkindex_options *arguments = tool_arguments;

  op = gt_option_parser_new("[options] enhanced-suffix-array",
                            "Count and index k-mers in the given enhanced "
                            "suffix array for a fixed value of k.");
  gt_option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");

  option = gt_option_new_ulong("mersize",
                               "Specify the mer size.",
                               &arguments->mersize,
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

  optionpl = gt_option_new_uint_min("pl",
                 "specify prefix length for bucket boundary construction\n"
                 "recommendation: use without argument;\n"
                 "then a reasonable prefix length is automatically determined",
                 &arguments->userdefinedprefixlength,
                 0,
                 1U);
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

static int gt_tyr_mkindex_arguments_check(int rest_argc,
                                               void *tool_arguments,
                                               GtError *err)
{
  Tyr_mkindex_options *arguments = tool_arguments;

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

static int gt_tyr_mkindex_runner(GT_UNUSED int argc,
                                 const char **argv,
                                 int parsed_args,
                                 void *tool_arguments,
                                 GtError *err)
{
  Tyr_mkindex_options *arguments = tool_arguments;
  Verboseinfo *verboseinfo;
  bool haserr = false;

  gt_assert(parsed_args + 1 == argc);
  gt_str_set(arguments->str_inputindex,argv[parsed_args]);
  verboseinfo = newverboseinfo(arguments->verbose);
  if (arguments->verbose)
  {
    printf("# mersize=%lu\n",arguments->mersize);
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
        printf("%u",arguments->prefixlength.value);
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
                    arguments->mersize,
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
  if (!haserr &&
      gt_str_length(arguments->str_inputindex) > 0 &&
      arguments->prefixlength.flag != Undeterminedprefixlength)
  {
    Definedunsignedint callprefixlength;

    if (arguments->prefixlength.flag == Determinedprefixlength)
    {
      callprefixlength.defined = true;
      callprefixlength.valueunsignedint = arguments->prefixlength.value;
    } else
    {
      callprefixlength.defined = false;
    }
    if (constructmerbuckets(arguments->str_inputindex,&callprefixlength) < 0)
    {
      haserr = true;
    }
  }
  freeverboseinfo(&verboseinfo);
  return haserr ? - 1 : 0;
}

static GtTool* gt_tyr_mkindex(void)
{
  return gt_tool_new(gt_tyr_mkindex_arguments_new,
                     gt_tyr_mkindex_arguments_delete,
                     gt_tyr_mkindex_option_parser_new,
                     gt_tyr_mkindex_arguments_check,
                     gt_tyr_mkindex_runner);
}

static int gt_tyr_occratio(GT_UNUSED int argc,
                                GT_UNUSED const char *argv[],
                                GtError *err)
{
  gt_error_set(err,"tyr occratio not implemented yet");
  return -1;
}

typedef struct
{
  GtStr *str_indexname;
  GtStrArray *queryfilenames;
  GtStr *strandspec;
  GtStrArray *showmodespec;
  unsigned int strand,
               showmode;
  bool verbose,
       performtest;
} Tyr_search_options;

static void *gt_tyr_search_arguments_new(void)
{
  Tyr_search_options *arguments
    = gt_malloc(sizeof (Tyr_search_options));
  arguments->str_indexname = gt_str_new();
  arguments->strandspec = gt_str_new();
  arguments->queryfilenames = gt_str_array_new();
  arguments->showmodespec = gt_str_array_new();
  arguments->showmode = 0;
  arguments->strand = 0;
  return arguments;
}

static void gt_tyr_search_arguments_delete(void *tool_arguments)
{
  Tyr_search_options *arguments = tool_arguments;

  if (!arguments)
  {
    return;
  }
  gt_str_delete(arguments->str_indexname);
  gt_str_delete(arguments->strandspec);
  gt_str_array_delete(arguments->queryfilenames);
  gt_str_array_delete(arguments->showmodespec);
  gt_free(arguments);
}

static GtOptionParser
          *gt_tyr_search_option_parser_new(void *tool_arguments)
{
  GtOptionParser *op;
  GtOption *option;
  Tyr_search_options *arguments = tool_arguments;

  op = gt_option_parser_new("[options] tyr-indexname queryfile0 "
                            "[queryfile1..]",
                            "Search a set of k-mers in an index constructed "
                            "by \"gt tyr mkindex\".");
  gt_option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");

  option = gt_option_new_string("strand",
                                "specify the strand to be searched: "
                                "use f (for forward strand) or "
                                "p (for reverse complemented strand) or "
                                "fp (for both); default is f",
                                arguments->strandspec,
                                "f");
  gt_option_parser_add_option(op, option);

  option = gt_option_new_stringarray("output",
                                     "specify output flags "
                                     "(qseqnum, qpos, counts, sequence)",
                                     arguments->showmodespec);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("test", "perform tests to verify program "
                                      "correctness", &arguments->performtest,
                                      false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_tyr_search_arguments_check(int rest_argc,
                                              void *tool_arguments,
                                              GtError *err)
{
  Optionargmodedesc showmodedesctable[] =
  {
    {"qseqnum",SHOWQSEQNUM},
    {"qpos",SHOWQPOS},
    {"counts",SHOWCOUNTS},
    {"sequence",SHOWSEQUENCE}
  };

  Optionargmodedesc stranddesctable[] =
  {
    {"f",STRAND_FORWARD},
    {"p",STRAND_REVERSE},
    {"fp",STRAND_FORWARD | STRAND_REVERSE}
  };

  unsigned long idx;
  Tyr_search_options *arguments = tool_arguments;

  if (rest_argc < 1)
  {
    gt_error_set(err,"missing tyr-indexnames and queryfilenames");
    return -1;
  }
  if (rest_argc < 2)
  {
    gt_error_set(err,"missing queryfilenames");
    return -1;
  }
  for (idx=0; idx<gt_str_array_size(arguments->showmodespec); idx++)
  {
    if (optionaddbitmask(showmodedesctable,
                         sizeof (showmodedesctable)/
                         sizeof (showmodedesctable[0]),
                         &arguments->showmode,
                         "-output",
                         gt_str_array_get(arguments->showmodespec,idx),
                         err) != 0)
    {
      return -1;
    }
  }
  if (optionaddbitmask(stranddesctable,
                       sizeof (stranddesctable)/
                       sizeof (stranddesctable[0]),
                       &arguments->strand,
                       "-output",
                       gt_str_get(arguments->strandspec),err) != 0)
  {
    return -1;
  }
  return 0;
}

static int gt_tyr_search_runner(int argc,
                                     const char **argv,
                                     int parsed_args,
                                     void *tool_arguments,
                                     GtError *err)
{
  int idx;
  Tyr_search_options *arguments = tool_arguments;

  gt_assert(parsed_args + 2 <= argc);
  gt_str_set(arguments->str_indexname,argv[parsed_args]);
  for (idx=parsed_args+1; idx<argc; idx++)
  {
    gt_str_array_add_cstr(arguments->queryfilenames,argv[idx]);
  }
  if (tyrsearch(arguments->str_indexname,
                     arguments->queryfilenames,
                     arguments->showmode,
                     arguments->strand,
                     arguments->verbose,
                     arguments->performtest,
                     err) != 0)
  {
    return -1;
  }
  return 0;
}

static GtTool *gt_tyr_search(void)
{
  return gt_tool_new(gt_tyr_search_arguments_new,
                     gt_tyr_search_arguments_delete,
                     gt_tyr_search_option_parser_new,
                     gt_tyr_search_arguments_check,
                     gt_tyr_search_runner);
}

static void *gt_tyr_arguments_new(void)
{
  GtToolbox *tyr_toolbox = gt_toolbox_new();
  gt_toolbox_add_tool(tyr_toolbox, "mkindex", gt_tyr_mkindex());
  gt_toolbox_add(tyr_toolbox, "occratio", gt_tyr_occratio);
  gt_toolbox_add_tool(tyr_toolbox, "search", gt_tyr_search());
  return tyr_toolbox;
}

static void gt_tyr_arguments_delete(void *tool_arguments)
{
  GtToolbox *index_toolbox = tool_arguments;
  if (!index_toolbox) return;
  gt_toolbox_delete(index_toolbox);
}

static GtOptionParser* gt_tyr_option_parser_new(void *tool_arguments)
{
  GtToolbox *index_toolbox = tool_arguments;
  GtOptionParser *op;

  gt_assert(index_toolbox != NULL);
  op = gt_option_parser_new(
                    "[option ...] [mkindex|occratio|search] [argument ...]",
                    "Call tyr with specific tool and "
                    "pass argument(s) to it.");
  gt_option_parser_set_comment_func(op, gt_toolbox_show, index_toolbox);
  gt_option_parser_refer_to_manual(op);
  return op;
}

static int gt_tyr_runner(int argc, const char **argv, int parsed_args,
                                 void *tool_arguments, GtError *err)
{
  GtToolbox *index_toolbox = tool_arguments;
  GtToolfunc toolfunc;
  GtTool *tool = NULL;
  char **nargv = NULL;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(index_toolbox != NULL);

  /* determine tool */
  if (!gt_toolbox_has_tool(index_toolbox, argv[parsed_args]))
  {
    gt_error_set(err, "tyr tool '%s' not found; option -help lists "
                   "possible tools", argv[parsed_args]);
    had_err = -1;
  }

  /* call sub-tool */
  if (!had_err)
  {
    if (!(toolfunc = gt_toolbox_get(index_toolbox, argv[parsed_args])))
    {
      tool = gt_toolbox_get_tool(index_toolbox, argv[parsed_args]);
      gt_assert(tool != NULL);
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
  return gt_tool_new(gt_tyr_arguments_new,
                     gt_tyr_arguments_delete,
                     gt_tyr_option_parser_new,
                     NULL,
                     gt_tyr_runner);
}
