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
#include "tools/gt_tagerator.h"
#include "libgtmatch/tagerator.h"

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
         *optionpckindex;

  assert(arguments != NULL);
  arguments->esaindexname = str_new();
  arguments->pckindexname = str_new();
  arguments->tagfiles = strarray_new();
  op = option_parser_new("[options] -t tagfile [-esa|-pck] indexname",
                         "Map short sequence tags in given index.");
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  option = option_new_filenamearray("t","Specify files containing the tags",
                                    arguments->tagfiles);
  option_parser_add_option(op, option);
  option_is_mandatory(option);
  option = option_new_ulong("k",
                            "Specify the allowed number of difference",
                            &arguments->maxdistance,
                            0);
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

  optiononline = option_new_bool("online","Perform online searches",
                            &arguments->online, false);
  option_parser_add_option(op, optiononline);

  optioncmp = option_new_bool("cmp","compare results of offline and online "
                                 "searches",
                            &arguments->docompare, false);
  option_parser_add_option(op, optioncmp);
  option_exclude(optiononline,optioncmp);

  optionrw = option_new_bool("rw","Replace wildcard in tag by random char",
                             &arguments->replacewildcard, false);
  option_parser_add_option(op, optionrw);
  option_is_development_option(optionrw);

  option = option_new_bool("d","Compute direct matches (default)",
                             &arguments->fwdmatch, true);
  option_parser_add_option(op, option);

  option = option_new_bool("p","Compute palindromic "
                           "(i.e. reverse complemented matches)",
                             &arguments->rcmatch, false);
  option_parser_add_option(op, option);

  option = option_new_bool("nospecials","do not output matches containing "
                           "wildcard characters (e.g. N)",
                           &arguments->nospecials, false);
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
  for (idx=0; idx<strarray_size(arguments->tagfiles); idx++)
  {
    printf("# tagfile=%s\n",strarray_get(arguments->tagfiles,idx));
  }
  printf("# maxdifference=%lu\n",arguments->maxdistance);
  if (str_length(arguments->esaindexname) > 0)
  {
    printf("# indexname(esa)=%s\n",str_get(arguments->esaindexname));
  } else
  {
    assert(str_length(arguments->pckindexname) > 0);
    printf("# indexname(pck)=%s\n",str_get(arguments->pckindexname));
  }
  if (runtagerator(arguments,err) != 0)
  {
    haserr = true;
  }
  return haserr ? -1 : 0;
}

Tool* gt_tagerator(void)
{
  return tool_new(gt_tagerator_arguments_new,
                  gt_tagerator_arguments_delete,
                  gt_tagerator_option_parser_new,
                  NULL,
                  gt_tagerator_runner);
}
