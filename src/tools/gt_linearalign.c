#include <string.h>
#include <inttypes.h>
#include "core/error.h"
#include "core/ma_api.h"
#include "core/versionfunc.h"
#include "core/option_api.h"
#include "core/str.h"
#include "core/str_array.h"
#include "core/types_api.h"
#include "extended/linearedist.h"
#include "extended/linearspaceAlignment.h"
#include "extended/localAlignment.h"
#include "tools/gt_linearalign.h"

typedef struct
{
  GtStrArray *strings,
             *costs,
             *scores;
  bool showedist,
       global,
       local;
} Alignopt;


static void showsimpleoptions(const Alignopt *opt)
{
  if (gt_str_array_size(opt->strings) > 0)
  {
    if (!opt->showedist)
      printf("# two strings \"%s\" \"%s\"\n", gt_str_array_get(opt->strings,0),
             gt_str_array_get(opt->strings,1UL));
    //return;
  }
}

static void freesimpleoption(Alignopt *opt)
{
  gt_str_array_delete(opt->strings);
  gt_str_array_delete(opt->costs);
  gt_str_array_delete(opt->scores);
}

static GtOPrval parse_options(int *parsed_args,
                              Alignopt *aop,
                              int argc, const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOption *optionstrings,
           *optionshowedist,
           *optioncosts,
           *optionscores;
  GtOPrval oprval;

  gt_error_check(err);
  aop->strings = gt_str_array_new();
  aop->costs = gt_str_array_new();
  aop->scores = gt_str_array_new();
  //op->files = gt_str_array_new();
  aop->showedist = false;
  aop->local = false;
  aop->global = true;
  
  op = gt_option_parser_new("options", "Apply function to compute alignment.");
  gt_option_parser_set_mail_address(op, "<annika.seidel@studium.uni-hamburg.de>");

  optionstrings = gt_option_new_string_array("ss", "use two strings",
                                             aop->strings);
  gt_option_parser_add_option(op, optionstrings);

  optionshowedist = gt_option_new_bool("e", "output unit edit distance",
                      &aop->showedist, false);
  gt_option_parser_add_option(op, optionshowedist);

  //optionprint = gt_option_new_bool("p", "print alignments",
    //                  &aop->print, false);
                       
 // gt_option_parser_add_option(op, optionprint);
  
  optioncosts = gt_option_new_string_array("g", "use three costs", aop->costs);
  
  gt_option_parser_add_option(op, optioncosts);

  optionscores = gt_option_new_string_array("l", "use three scores", aop->scores);
  
  gt_option_parser_add_option(op, optionscores);

  //gt_option_exclude(optionstrings, optionfiles);
  //gt_option_exclude(optionlocal, optionglobal);
  gt_option_imply(optionshowedist, optionstrings);
  //gt_option_imply(optionprint, optionstrings);
  gt_option_imply(optioncosts, optionstrings);
  gt_option_imply(optionscores, optionstrings);
  
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  if (oprval == GT_OPTION_PARSER_OK)
  {
    if (gt_option_is_set(optionstrings))
    {
      if (gt_str_array_size(aop->strings) != 2UL)
      {
        gt_error_set(err, "option -ss requires two string arguments");
        oprval = GT_OPTION_PARSER_ERROR;
      }
      
      if(gt_option_is_set(optionscores))
      {
          aop->local=true;
          if(gt_str_array_size(aop->scores) != 3UL)
          {
            gt_error_set(err, "option -l requires matchscore, mismatchscore, gapscore");
            oprval = GT_OPTION_PARSER_ERROR;
          }  
      }
      if(gt_option_is_set(optioncosts))
      {
          if(gt_str_array_size(aop->costs) != 3UL)
          {
            gt_error_set(err, "option -g requires matchcost, mismatchcost, gapcosts");
            oprval = GT_OPTION_PARSER_ERROR;
          }  
      }
    }else
    {
      gt_error_set(err,"use exactly one of the options -ss");//-ff
      oprval = GT_OPTION_PARSER_ERROR;
    }
  }
  gt_option_parser_delete(op);
  if (oprval == GT_OPTION_PARSER_OK && *parsed_args != argc)
  {
    gt_error_set(err, "superfluous program parameters");
    oprval = GT_OPTION_PARSER_ERROR;
  }
  return oprval;
}

int gt_linearalign(int argc, const char **argv, GtError *err)
{
  int parsed_args;
  Alignopt aop;
  GtOPrval oprval;

  gt_error_check(err);

  oprval = parse_options(&parsed_args, &aop, argc, argv, err);
  if (oprval == GT_OPTION_PARSER_OK)
  {
    gt_assert(parsed_args == argc);
    showsimpleoptions(&aop);
    if (aop.showedist)//global alignment
    {
      GtUword edist;
      edist = gt_calc_linearedist(
        (const GtUchar *) gt_str_array_get(aop.strings,0),
        (GtUword) strlen(gt_str_array_get(aop.strings,0)),
        (const GtUchar *) gt_str_array_get(aop.strings,1UL),
        (GtUword) strlen(gt_str_array_get(aop.strings,1UL)));
      printf(GT_WU "\n", edist);
    }
   /* else if (cmppairwise.print)
    {
      gt_computelinearspace(
        (const GtUchar *) gt_str_array_get(cmppairwise.strings,0),
        (GtUword) strlen(gt_str_array_get(cmppairwise.strings,0)),
        (const GtUchar *) gt_str_array_get(cmppairwise.strings,1UL),
        (GtUword) strlen(gt_str_array_get(cmppairwise.strings,1UL)));
    }*/
    else if (aop.local)
    {
      GtWord matchscore, mismatchscore, gapscore, check;
      check = sscanf(gt_str_array_get(aop.scores,0),GT_WD, &matchscore);
      gt_assert(check == 1);
      check = sscanf(gt_str_array_get(aop.scores,1UL),GT_WD, &mismatchscore);
      gt_assert(check == 1);
      check = sscanf(gt_str_array_get(aop.scores,2UL),GT_WD, &gapscore);
      gt_assert(check == 1);
      
      gt_computelinearspace_local(
      (const GtUchar *) gt_str_array_get(aop.strings,0),
      (GtUword) strlen(gt_str_array_get(aop.strings,0)),
      (const GtUchar *) gt_str_array_get(aop.strings,1UL),
      (GtUword) strlen(gt_str_array_get(aop.strings,1UL)),
      matchscore,mismatchscore,gapscore);  
    }else
    {
      GtWord matchcost = 0, mismatchcost = 1, gapcost = 1, check;
      if (gt_str_array_size(aop.costs))
      {
        check = sscanf(gt_str_array_get(aop.costs,0),GT_WD, &matchcost);
        gt_assert(check == 1);
        check = sscanf(gt_str_array_get(aop.costs,1UL),GT_WD, &mismatchcost);
        gt_assert(check == 1);
        check = sscanf(gt_str_array_get(aop.costs,2UL),GT_WD, &gapcost);
        gt_assert(check == 1);
      }

      gt_computelinearspace_with_output(
      (const GtUchar *) gt_str_array_get(aop.strings,0),
      (GtUword) strlen(gt_str_array_get(aop.strings,0)),
      (const GtUchar *) gt_str_array_get(aop.strings,1UL),
      (GtUword) strlen(gt_str_array_get(aop.strings,1UL)),
      matchcost,mismatchcost,gapcost);
    }
  }
  freesimpleoption(&aop);
  if (oprval == GT_OPTION_PARSER_REQUESTS_EXIT)
  {
    return 0;
  }
  if (oprval == GT_OPTION_PARSER_ERROR)
  {
    return -1;
  }
  return 0;
}
