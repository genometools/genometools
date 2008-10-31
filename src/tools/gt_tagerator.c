/*
  Copyright (c) 2006-2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#include "core/option.h"
#include "core/ma.h"
#include "core/str_array.h"
#include "core/unused_api.h"
#include "core/tool.h"
#include "match/tagerator.h"
#include "tools/gt_tagerator.h"

/*
  Remark Gordon from July 2008: How to have access to an option in the
  GtTool argument parser framework:

  typedef struct
  {
    GtOption *foo;
  } Arguments;

  foo = gt_option_parser_new();
  optionparseraddoption(op,foo);
  arguments->foo = gt_option_ref(foo);

  ...

  arguments_delete()
  {
    gt_option_delete(arguments_foo);
  }
*/

static void *gt_tagerator_arguments_new(void)
{
  return gt_malloc(sizeof (TageratorOptions));
}

static void gt_tagerator_arguments_delete(void *tool_arguments)
{
  TageratorOptions *arguments = tool_arguments;

  if (!arguments)
  {
    return;
  }
  gt_str_delete(arguments->esaindexname);
  gt_str_delete(arguments->pckindexname);
  gt_str_array_delete(arguments->tagfiles);
  gt_free(arguments);
}

static GtOptionParser* gt_tagerator_option_parser_new(void *tool_arguments)
{
  TageratorOptions *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionrw, *optiononline, *optioncmp, *optionesaindex,
           *optionpckindex, *optionmaxdepth;

  gt_assert(arguments != NULL);
  arguments->esaindexname = gt_str_new();
  arguments->pckindexname = gt_str_new();
  arguments->tagfiles = gt_str_array_new();
  op = gt_option_parser_new("[options] -t tagfile [-esa|-pck] indexname",
                         "Map short sequence tags in given index.");
  gt_option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  option = gt_option_new_filenamearray("q",
                                    "Specify files containing the short "
                                    "sequence tags",
                                    arguments->tagfiles);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  option = gt_option_new_long("e",
                           "Specify the allowed number of differences "
                           "(replacements/insertions/deletions)",
                           &arguments->maxdistance,
                           -1L);
  gt_option_parser_add_option(op, option);

  optionesaindex = gt_option_new_string("esa",
                                     "Specify index (enhanced suffix array)",
                                     arguments->esaindexname, NULL);
  gt_option_parser_add_option(op, optionesaindex);

  optionpckindex = gt_option_new_string("pck",
                                     "Specify index (packed index)",
                                     arguments->pckindexname, NULL);
  gt_option_parser_add_option(op, optionpckindex);
  gt_option_exclude(optionesaindex,optionpckindex);
  gt_option_is_mandatory_either(optionesaindex,optionpckindex);

  optionmaxdepth = gt_option_new_int("maxdepth",
                                  "Use the data in the .pbt file only up to "
                                  "this depth (only relevant with option -pck)",
                                  &arguments->userdefinedmaxdepth,
                                  -1);
  gt_option_parser_add_option(op, optionmaxdepth);
  gt_option_is_development_option(optionmaxdepth);

  optiononline = gt_option_new_bool("online","Perform online searches",
                            &arguments->online, false);
  gt_option_parser_add_option(op, optiononline);
  gt_option_is_development_option(optiononline);

  optioncmp = gt_option_new_bool("cmp","compare results of offline and online "
                              "searches",
                            &arguments->docompare, false);
  gt_option_parser_add_option(op, optioncmp);
  gt_option_exclude(optiononline,optioncmp);
  gt_option_is_development_option(optioncmp);

  optionrw = gt_option_new_bool("rw","Replace wildcard in tag by random char",
                             &arguments->replacewildcard, false);
  gt_option_parser_add_option(op, optionrw);
  gt_option_is_development_option(optionrw);

  option = gt_option_new_bool("nod","Do not compute direct matches",
                           &arguments->nofwdmatch, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("nop","Do not compute palindromic matches "
                           "(i.e. no reverse complemented matches.)",
                             &arguments->norcmatch, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_ulong_min("maxocc",
                                "specify max number of match-occurrences",
                                &arguments->maxintervalwidth,0,1UL);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("skpp",
                           "Skip prefix of pattern (only in pdiff mode)",
                           &arguments->skpp, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("withwildcards","output matches containing "
                              "wildcard characters (e.g. N); only relevant for "
                              "approximate matching",
                              &arguments->nowildcards, true);
  gt_option_parser_add_option(op, option);
  return op;
}

static int gt_tagerator_arguments_check(GT_UNUSED int rest_argc,
                                        void *tool_arguments,
                                        GtError *err)
{
  TageratorOptions *arguments = tool_arguments;

  if (arguments->maxdistance < 0)
  {
    if (arguments->online)
    {
      gt_error_set(err,"option -online requires option -e");
      return -1;
    }
    if (!arguments->nowildcards)
    {
      arguments->nowildcards = true;
    }
    if (arguments->maxintervalwidth == 0)
    {
      gt_error_set(err,
                   "if option -e is not used then option -maxocc is required");
      return -1;
    }
  } else
  {
    if (arguments->skpp &&
        (arguments->maxdistance == 0 || arguments->maxintervalwidth == 0))
    {
      gt_error_set(err,"option -skpp only works in pdiff mode");
      return -1;
    }
  }
  return 0;
}

static int gt_tagerator_runner(GT_UNUSED int argc,
                               GT_UNUSED const char **argv,
                               GT_UNUSED int parsed_args,
                               void *tool_arguments, GtError *err)
{
  TageratorOptions *arguments = tool_arguments;
  bool haserr = false;
  unsigned long idx;

  gt_error_check(err);
  gt_assert(arguments != NULL);

  gt_assert(parsed_args == argc);
  if (arguments->maxdistance == -1L)
  {
    printf("# computing matching statistics\n");
  } else
  {
    if (arguments->maxintervalwidth == 0)
    {
      printf("# computing complete matches");
    } else
    {
      printf("# computing prefix matches");
    }
    if (arguments->maxdistance == 0)
    {
      printf(" without differences (exact matches)");
    } else
    {
      printf(" with up to %ld differences",arguments->maxdistance);
    }
    if (arguments->maxintervalwidth > 0)
    {
      printf(" and at most %lu occurrences in the subject sequences",
             arguments->maxintervalwidth);
    }
    printf("\n");
  }
  if (gt_str_length(arguments->esaindexname) > 0)
  {
    printf("# indexname(esa)=%s\n",gt_str_get(arguments->esaindexname));
  } else
  {
    gt_assert(gt_str_length(arguments->pckindexname) > 0);
    printf("# indexname(pck)=%s\n",gt_str_get(arguments->pckindexname));
  }
  for (idx=0; idx<gt_str_array_size(arguments->tagfiles); idx++)
  {
    printf("# queryfile=%s\n",gt_str_array_get(arguments->tagfiles,idx));
  }
  if (runtagerator(arguments,err) != 0)
  {
    haserr = true;
  }
  return haserr ? -1 : 0;
}

GtTool* gt_tagerator(void)
{
  return gt_tool_new(gt_tagerator_arguments_new,
                     gt_tagerator_arguments_delete,
                     gt_tagerator_option_parser_new,
                     gt_tagerator_arguments_check,
                     gt_tagerator_runner);
}
