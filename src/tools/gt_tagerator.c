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

#include "libgtcore/option.h"
#include "libgtcore/ma.h"
#include "libgtcore/strarray.h"
#include "libgtcore/unused.h"
#include "libgtcore/tool.h"
#include "libgtmatch/tagerator.h"
#include "tools/gt_tagerator.h"

/*
  Remark Gordon from July 2008: How to have access to an option in the
  Tool argument parser framework:

  typedef struct
  {
    Option *foo;
  } Arguments;

  foo = option_parser_new();
  optionparseraddoption(op,foo);
  arguments->foo = option_ref(foo);

  ...

  arguments_delete()
  {
    option_delete(arguments_foo);
  }
*/

static void *gt_tagerator_arguments_new(void)
{
  return ma_malloc(sizeof (TageratorOptions));
}

static void gt_tagerator_arguments_delete(void *tool_arguments)
{
  TageratorOptions *arguments = tool_arguments;

  if (!arguments)
  {
    return;
  }
  str_delete(arguments->esaindexname);
  str_delete(arguments->pckindexname);
  strarray_delete(arguments->tagfiles);
  ma_free(arguments);
}

static OptionParser* gt_tagerator_option_parser_new(void *tool_arguments)
{
  TageratorOptions *arguments = tool_arguments;
  OptionParser *op;
  Option *option, *optionrw, *optiononline, *optioncmp, *optionesaindex,
         *optionpckindex, *optionmaxdepth;

  assert(arguments != NULL);
  arguments->esaindexname = str_new();
  arguments->pckindexname = str_new();
  arguments->tagfiles = strarray_new();
  op = option_parser_new("[options] -t tagfile [-esa|-pck] indexname",
                         "Map short sequence tags in given index.");
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  option = option_new_filenamearray("q",
                                    "Specify files containing the short "
                                    "sequence tags",
                                    arguments->tagfiles);
  option_parser_add_option(op, option);
  option_is_mandatory(option);

  option = option_new_long("e",
                           "Specify the allowed number of differences "
                           "(replacements/insertions/deletions)",
                           &arguments->maxdistance,
                           -1L);
  option_parser_add_option(op, option);

  optionesaindex = option_new_string("esa",
                                     "Specify index (enhanced suffix array)",
                                     arguments->esaindexname, NULL);
  option_parser_add_option(op, optionesaindex);

  optionpckindex = option_new_string("pck",
                                     "Specify index (packed index)",
                                     arguments->pckindexname, NULL);
  option_parser_add_option(op, optionpckindex);
  option_exclude(optionesaindex,optionpckindex);
  option_is_mandatory_either(optionesaindex,optionpckindex);

  optionmaxdepth = option_new_int("maxdepth",
                                  "Use the data in the .pbt file only up to "
                                  "this depth (only relevant with option -pck)",
                                  &arguments->userdefinedmaxdepth,
                                  -1);
  option_parser_add_option(op, optionmaxdepth);
  option_is_development_option(optionmaxdepth);

  optiononline = option_new_bool("online","Perform online searches",
                            &arguments->online, false);
  option_parser_add_option(op, optiononline);
  option_is_development_option(optiononline);

  optioncmp = option_new_bool("cmp","compare results of offline and online "
                              "searches",
                            &arguments->docompare, false);
  option_parser_add_option(op, optioncmp);
  option_exclude(optiononline,optioncmp);
  option_is_development_option(optioncmp);

  optionrw = option_new_bool("rw","Replace wildcard in tag by random char",
                             &arguments->replacewildcard, false);
  option_parser_add_option(op, optionrw);
  option_is_development_option(optionrw);

  option = option_new_bool("nod","Do not compute direct matches",
                           &arguments->nofwdmatch, false);
  option_parser_add_option(op, option);

  option = option_new_bool("nop","Do not compute palindromic matches "
                           "(i.e. no reverse complemented matches.)",
                             &arguments->norcmatch, false);
  option_parser_add_option(op, option);

  option = option_new_ulong_min("maxocc",
                                "specify max number of match-occurrencs",
                                &arguments->maxintervalwidth,0,1UL);
  option_parser_add_option(op, option);

  option = option_new_bool("nowildcards","do not output matches containing "
                           "wildcard characters (e.g. N); only relevant for "
                           "approximate matching",
                           &arguments->nowildcards, false);
  option_parser_add_option(op, option);
  return op;
}

static int gt_tagerator_runner(UNUSED int argc,
                               UNUSED const char **argv,
                               UNUSED int parsed_args,
                               void *tool_arguments, Error *err)
{
  TageratorOptions *arguments = tool_arguments;
  bool haserr = false;
  unsigned long idx;

  error_check(err);
  assert(arguments != NULL);

  assert(parsed_args == argc);
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
  if (str_length(arguments->esaindexname) > 0)
  {
    printf("# indexname(esa)=%s\n",str_get(arguments->esaindexname));
  } else
  {
    assert(str_length(arguments->pckindexname) > 0);
    printf("# indexname(pck)=%s\n",str_get(arguments->pckindexname));
  }
  for (idx=0; idx<strarray_size(arguments->tagfiles); idx++)
  {
    printf("# queryfile=%s\n",strarray_get(arguments->tagfiles,idx));
  }
  if (runtagerator(arguments,err) != 0)
  {
    haserr = true;
  }
  return haserr ? -1 : 0;
}

static int gt_tagerator_arguments_check(UNUSED int rest_argc,
                                        void *tool_arguments,
                                        Error *err)
{
  TageratorOptions *arguments = tool_arguments;

  if (arguments->maxdistance < 0)
  {
    if (arguments->online)
    {
      error_set(err,"option -online requires option -e");
      return -1;
    }
    if (!arguments->nowildcards)
    {
      arguments->nowildcards = true;
    }
    if (arguments->maxintervalwidth == 0)
    {
      error_set(err,"if option -e is not used then option -maxocc is required");
      return -1;
    }
  }
  return 0;
}

Tool* gt_tagerator(void)
{
  return tool_new(gt_tagerator_arguments_new,
                  gt_tagerator_arguments_delete,
                  gt_tagerator_option_parser_new,
                  gt_tagerator_arguments_check,
                  gt_tagerator_runner);
}
