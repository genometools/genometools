#include <string.h>
#include <inttypes.h>
#include "core/error.h"
#include "core/fa.h"
#include "core/ma_api.h"
#include "core/versionfunc.h"
#include "core/option_api.h"
#include "core/str.h"
#include "core/str_api.h"
#include "core/str_array.h"
#include "core/types_api.h"
#include "extended/affinealign_linear.h"
#include "extended/affinealign_linear_local.h"
#include "extended/linearedist.h"
#include "extended/linearspace.h"
#include "extended/linearspace_local.h"
#include "tools/gt_linearalign.h"

typedef struct
{
  GtStrArray *strings,
             *costs,
             *scores;
  GtStr *outfilename;
  bool showedist,
       showevalue,
       global,
       local,
       affine,
       outfile;
} Alignopt;

static void showsimpleoptions(const Alignopt *opt)
{
  if (gt_str_array_size(opt->strings) > 0)
  {
    if (!opt->showedist)
      printf("# two strings \"%s\" \"%s\"\n", gt_str_array_get(opt->strings,0),
             gt_str_array_get(opt->strings,1UL));
  }
}

static void freesimpleoption(Alignopt *opt)
{
  gt_str_array_delete(opt->strings);
  gt_str_array_delete(opt->costs);
  gt_str_array_delete(opt->scores);
  gt_str_delete(opt->outfilename);
  
}

static GtOPrval parse_options(int *parsed_args,
                              Alignopt *aop,
                              int argc, const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOption *optionstrings,
           *optionshowedist,
           *optionshowevalue,
           *optioncosts,
           *optionscores,
           *optionaffine,
           *optionoutfilename;
  GtOPrval oprval;

  gt_error_check(err);
  aop->strings = gt_str_array_new();
  aop->costs = gt_str_array_new();
  aop->scores = gt_str_array_new();
  aop->outfilename = gt_str_new();
  aop->showedist = false;
  aop->showevalue = false;
  aop->local = false;
  aop->global = true;
  aop->outfile = false;
  aop->affine = false;

  op = gt_option_parser_new("options", "Apply function to compute alignment.");
  gt_option_parser_set_mail_address(op, "<>");

  optionstrings = gt_option_new_string_array("ss", "use two strings",
                                             aop->strings);
  gt_option_parser_add_option(op, optionstrings);

  optionshowedist = gt_option_new_bool("e", "output unit edit distance",
                      &aop->showedist, false);
  gt_option_parser_add_option(op, optionshowedist);

  optionoutfilename = gt_option_new_string("o", "use outputfile",
                      aop->outfilename, "stdout");
  gt_option_parser_add_option(op, optionoutfilename);

  optioncosts = gt_option_new_string_array("g", "global alignment, "
                                           "use three costs", aop->costs);

  gt_option_parser_add_option(op, optioncosts);

  optionscores = gt_option_new_string_array("l", "local alignment, "
                                            "use three scores",
                                            aop->scores);

  gt_option_parser_add_option(op, optionscores);

  optionaffine = gt_option_new_bool("a", "use g or l with affine gapcosts",
                      &aop->affine, false);
  gt_option_parser_add_option(op, optionaffine);
  
  optionshowevalue = gt_option_new_bool("v", "output costs",
                      &aop->showevalue, false);
  gt_option_parser_add_option(op, optionshowevalue);

  gt_option_exclude(optioncosts, optionscores);
  gt_option_exclude(optionshowedist, optionaffine);
  gt_option_exclude(optionshowedist, optionshowevalue);
  gt_option_exclude(optionshowedist, optioncosts);
  gt_option_exclude(optionshowedist, optionscores);
  gt_option_imply(optionshowedist, optionstrings);
  gt_option_imply(optioncosts, optionstrings);
  gt_option_imply(optionscores, optionstrings);
  gt_option_imply(optionshowevalue, optionstrings);
  //gt_option_imply_either_2(optionshowevalue, optioncosts, optionscores);
  gt_option_imply_either_2(optionaffine, optioncosts, optionscores);

  oprval = gt_option_parser_parse(op, parsed_args, argc, argv,
                                  gt_versionfunc, err);
  if (oprval == GT_OPTION_PARSER_OK)
  {
    if (gt_option_is_set(optionstrings))
    {
      if (gt_str_array_size(aop->strings) != 2UL)
      {
        gt_error_set(err, "option -ss requires two string arguments");
        oprval = GT_OPTION_PARSER_ERROR;
      }

      if (gt_option_is_set(optionscores))
      {
          aop->local=true;
          if (gt_str_array_size(aop->scores) != 3UL)
          {
            gt_error_set(err, "option -l requires matchscore, "
                              "mismatchscore, gapscore");
            oprval = GT_OPTION_PARSER_ERROR;
          }
      }
      if (gt_option_is_set(optioncosts))
      {
          if (gt_str_array_size(aop->costs) != 3UL)
          {
            if (gt_option_is_set(optionaffine))
              gt_error_set(err, "option -g related to -a requires "
               "mismatchcost, gap_opening_cost, gap_extension_cost");
            else
              gt_error_set(err, "option -g requires matchcost, "
                                "mismatchcost, gapcosts");
            oprval = GT_OPTION_PARSER_ERROR;
          }
      }
      if (gt_option_is_set(optionoutfilename))
        aop->outfile = true;
      if (gt_option_is_set(optionaffine))
        aop->affine = true;
    }
    else
    {
      gt_error_set(err,"use exactly one of the options -ss");/*-ff*/
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

static void select_evalues(GtStrArray *evalues,
                           GtWord *val1,
                           GtWord *val2,
                           GtWord *val3)
{
  GtWord check;

  gt_assert(gt_str_array_size(evalues) == 3UL);

  check = sscanf(gt_str_array_get(evalues,0),GT_WD, val1);
        gt_assert(check == 1);
        check = sscanf(gt_str_array_get(evalues,1UL),GT_WD, val2);
        gt_assert(check == 1);
        check = sscanf(gt_str_array_get(evalues,2UL),GT_WD, val3);
        gt_assert(check == 1);
}

int gt_linearalign(int argc, const char **argv, GtError *err)
{
  int parsed_args;
  Alignopt aop;
  GtOPrval oprval;
  FILE *fp;
  GtWord matchcost = 0, mismatchcost = 1, gapcost = 1,
  matchscore,mismatchscore,gapscore;/*TODO:evalue*/

  gt_error_check(err);

  oprval = parse_options(&parsed_args, &aop, argc, argv, err);
  if (oprval == GT_OPTION_PARSER_OK)
  {
    gt_assert(parsed_args == argc);
   if (!aop.outfile)
      showsimpleoptions(&aop);

    fp = (aop.outfile? gt_fa_fopen_func(
                       gt_str_get(aop.outfilename), "a", __FILE__,__LINE__,err)
                       :stdout);
    gt_error_check(err);

    if (aop.showedist)/*global alignment*/
    {
      GtUword edist;
      edist = gt_calc_linearedist(
        (const GtUchar *) gt_str_array_get(aop.strings,0),
        (GtUword) strlen(gt_str_array_get(aop.strings,0)),
        (const GtUchar *) gt_str_array_get(aop.strings,1UL),
        (GtUword) strlen(gt_str_array_get(aop.strings,1UL)));
      fprintf(fp, "edist: "GT_WU "\n", edist);
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
      select_evalues(aop.scores, &matchscore, &mismatchscore, &gapscore);
      if (aop.affine)
      {
        gt_computeaffinelinearspace_local(//aop.showevalue,
      (const GtUchar *) gt_str_array_get(aop.strings,0),0,
      (GtUword) strlen(gt_str_array_get(aop.strings,0)),
      (const GtUchar *) gt_str_array_get(aop.strings,1UL),0,
      (GtUword) strlen(gt_str_array_get(aop.strings,1UL)),
      6,matchscore,mismatchscore,gapscore,fp);
        /*fprintf(stderr,"-l -a is not implemented\n");
        exit(GT_EXIT_PROGRAMMING_ERROR);*/
        
      }
      else{
      

      gt_computelinearspace_local(//aop.showevalue,
      (const GtUchar *) gt_str_array_get(aop.strings,0),0,
      (GtUword) strlen(gt_str_array_get(aop.strings,0)),
      (const GtUchar *) gt_str_array_get(aop.strings,1UL),0,
      (GtUword) strlen(gt_str_array_get(aop.strings,1UL)),
      matchscore,mismatchscore,gapscore,fp);}
    }else /* global */
    {
      if (gt_str_array_size(aop.costs))
      {
        select_evalues(aop.costs, &matchcost, &mismatchcost, &gapcost);
        if (matchcost < 0 || mismatchcost < 0 || gapcost < 0)
        {
          gt_error_set(err, "invalid cost value to option -g");
          oprval = GT_OPTION_PARSER_ERROR;
         return -1;
        }
      }
      if (aop.affine)
      {
        /* matchcost=replacement_cost,
         * mismatchcost=gap_opening_cost,
         * gapcost=gap_extension_cost*/
        gt_computeaffinelinearspace(//aop.showevalue,
        (const GtUchar *) gt_str_array_get(aop.strings,0),0,
        (GtUword) strlen(gt_str_array_get(aop.strings,0)),
        (const GtUchar *) gt_str_array_get(aop.strings,1UL),0,
        (GtUword) strlen(gt_str_array_get(aop.strings,1UL)),
         matchcost,mismatchcost, gapcost,fp);/*TODO:variablebennenung???*/
      }else
      {
        gt_computelinearspace2(//aop.showevalue,
        (const GtUchar *) gt_str_array_get(aop.strings,0),0,
        (GtUword) strlen(gt_str_array_get(aop.strings,0)),
        (const GtUchar *) gt_str_array_get(aop.strings,1UL),0,
        (GtUword) strlen(gt_str_array_get(aop.strings,1UL)),
        matchcost,mismatchcost,gapcost,fp);
      }
    }
    if (aop.outfile)
      gt_fa_fclose(fp);
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
