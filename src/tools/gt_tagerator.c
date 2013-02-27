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

#include "core/option_api.h"
#include "core/ma.h"
#include "core/str_array.h"
#include "core/unused_api.h"
#include "core/tool_api.h"
#include "match/tagerator.h"
#include "match/optionargmode.h"
#include "tools/gt_tagerator.h"

static const Optionargmodedesc outputmodedesctable[] =
{
  {"tagnum","ordinal number of tag",TAGOUT_TAGNUM},
  {"tagseq","tag sequence",TAGOUT_TAGSEQ},
  {"dblength","length of match in database",TAGOUT_DBLENGTH},
  {"dbstartpos","start position of match in database",TAGOUT_DBSTARTPOS},
  {"abspos","absolute value of dbstartpos",TAGOUT_DBABSPOS},
  {"dbsequence","sequence of match",TAGOUT_DBSEQUENCE},
  {"strand","strand",TAGOUT_STRAND},
  {"edist","edit distance",TAGOUT_EDIST},
  {"tagstartpos","start position of match in tag (only for -maxocc)",
                 TAGOUT_TAGSTARTPOS},
  {"taglength","length of match in tag (only for -maxocc)",TAGOUT_TAGLENGTH},
  {"tagsuffixseq","suffix tag involved in match (only for -maxocc)",
                  TAGOUT_TAGSUFFIXSEQ}
};

static void *gt_tagerator_arguments_new(void)
{
  TageratorOptions *arguments;

  arguments = gt_malloc(sizeof (TageratorOptions));
  arguments->indexname = gt_str_new();
  arguments->tagfiles = gt_str_array_new();
  arguments->outputspec = gt_str_array_new();
  arguments->outputmode = 0;
  arguments->numberofmodedescentries = sizeof (outputmodedesctable)/
                                       sizeof (outputmodedesctable[0]);
  arguments->outputhelp
    = gt_getargmodekeywords(outputmodedesctable,
                            arguments->numberofmodedescentries,
                            "output");
  arguments->modedesc = outputmodedesctable;
  return arguments;
}

static void gt_tagerator_arguments_delete(void *tool_arguments)
{
  TageratorOptions *arguments = tool_arguments;

  if (!arguments)
  {
    return;
  }
  gt_str_delete(arguments->indexname);
  gt_str_delete(arguments->outputhelp);
  gt_str_array_delete(arguments->tagfiles);
  gt_str_array_delete(arguments->outputspec);
  gt_option_delete(arguments->refoptionesaindex);
  gt_option_delete(arguments->refoptionpckindex);
  gt_free(arguments);
}

static GtOptionParser* gt_tagerator_option_parser_new(void *tool_arguments)
{
  TageratorOptions *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionrw, *optiononline, *optioncmp, *optionesaindex,
           *optionpckindex, *optionmaxdepth, *optionbest;

  gt_assert(arguments != NULL);
  op = gt_option_parser_new("[options] -q tagfile [-esa|-pck] indexname",
                         "Map short sequence tags in given index.");
  gt_option_parser_set_mail_address(op,"<kurtz@zbh.uni-hamburg.de>");
  option = gt_option_new_filename_array("q",
                                    "Specify files containing the short "
                                    "sequence tags",
                                    arguments->tagfiles);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  option = gt_option_new_long("e",
                           "Specify the allowed number of differences "
                           "(replacements/insertions/deletions)",
                           &arguments->userdefinedmaxdistance,
                           -1L);
  gt_option_parser_add_option(op, option);

  optionesaindex = gt_option_new_string("esa",
                                     "Specify index (enhanced suffix array)",
                                     arguments->indexname, NULL);
  gt_option_parser_add_option(op, optionesaindex);
  arguments->refoptionesaindex = gt_option_ref(optionesaindex);

  optionpckindex = gt_option_new_string("pck",
                                     "Specify index (packed index)",
                                     arguments->indexname, NULL);
  gt_option_parser_add_option(op, optionpckindex);
  arguments->refoptionpckindex = gt_option_ref(optionpckindex);
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
                                    &arguments->doonline, false);
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

  optionbest = gt_option_new_bool("best","Compute only best matches, i.e. only "
                                  "for smallest edit distance with matches",
                                  &arguments->best, false);
  gt_option_parser_add_option(op, optionbest);
  gt_option_exclude(optiononline,optionbest);
  gt_option_exclude(optioncmp,optionbest);

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

  option = gt_option_new_string_array("output",
                                      gt_str_get(arguments->outputhelp),
                                      arguments->outputspec);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_tagerator_arguments_check(GT_UNUSED int rest_argc,
                                        void *tool_arguments,
                                        GtError *err)
{
  TageratorOptions *arguments = tool_arguments;
  unsigned long idx;

  if (!arguments->nowildcards && arguments->userdefinedmaxdistance <= 0)
  {
    arguments->nowildcards = true;
  }
  if (gt_option_is_set(arguments->refoptionesaindex))
  {
    arguments->withesa = true;
  } else
  {
    gt_assert(gt_option_is_set(arguments->refoptionpckindex));
    arguments->withesa = false;
  }
  if (arguments->userdefinedmaxdistance < 0)
  {
    if (arguments->doonline)
    {
      gt_error_set(err,"option -online requires option -e");
      return -1;
    }
    if (arguments->maxintervalwidth == 0)
    {
      gt_error_set(err,
                   "if option -e is not used then option -maxocc is required");
      return -1;
    }
    if (arguments->best)
    {
      gt_error_set(err,"option -best requires option -e");
      return -1;
    }
  } else
  {
    if (arguments->skpp &&
        (arguments->userdefinedmaxdistance == 0 ||
         arguments->maxintervalwidth == 0))
    {
      gt_error_set(err,"option -skpp only works in pdiff mode");
      return -1;
    }
  }
  for (idx=0; idx<gt_str_array_size(arguments->outputspec); idx++)
  {
    if (gt_optionargaddbitmask(outputmodedesctable,
                            sizeof (outputmodedesctable)/
                            sizeof (outputmodedesctable[0]),
                            &arguments->outputmode,
                            "-output",
                            gt_str_array_get(arguments->outputspec,idx),
                            err) != 0)
    {
      return -1;
    }
  }
  if (arguments->outputmode == 0)
  {
    arguments->outputmode = TAGOUT_TAGNUM | TAGOUT_TAGSEQ |
                            TAGOUT_DBLENGTH | TAGOUT_DBSTARTPOS |
                            TAGOUT_STRAND;
    if (arguments->maxintervalwidth > 0)
    {
      arguments->outputmode |= TAGOUT_TAGLENGTH;
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
  if (arguments->userdefinedmaxdistance == -1L)
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
    if (arguments->userdefinedmaxdistance == 0)
    {
      printf(" without differences (exact matches)");
    } else
    {
      printf(" with up to %ld differences",arguments->userdefinedmaxdistance);
    }
    if (arguments->maxintervalwidth > 0)
    {
      printf(" and at most %lu occurrences in the subject sequences",
             arguments->maxintervalwidth);
    }
    printf("\n");
  }
  printf("# indexname(%s)=%s\n",arguments->withesa ? "esa" : "pck",
                                gt_str_get(arguments->indexname));
  for (idx=0; idx<gt_str_array_size(arguments->tagfiles); idx++)
  {
    printf("# queryfile=%s\n",gt_str_array_get(arguments->tagfiles,idx));
  }
  if (gt_runtagerator(arguments,err) != 0)
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
