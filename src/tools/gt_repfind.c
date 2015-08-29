/*
  Copyright (c) 2007-2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2015 Center for Bioinformatics, University of Hamburg

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

#include <inttypes.h>
#include <stdarg.h>
#include "core/error_api.h"
#include "core/format64.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/ma_api.h"
#include "core/option_api.h"
#include "core/str_api.h"
#include "core/tool_api.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "core/minmax.h"
#include "core/encseq.h"
#include "core/showtime.h"
#include "core/timer_api.h"
#include "match/esa-maxpairs.h"
#include "match/esa-mmsearch.h"
#include "match/querymatch.h"
#include "match/test-maxpairs.h"
#include "match/seed-extend.h"
#include "tools/gt_repfind.h"

typedef struct
{
  unsigned int userdefinedleastlength, /* minalilen */
               seedlength; /* like kmerlen, but here MEMs are used */
  GtXdropscore xdropbelowscore;
  GtUword samples,
          extendxdrop,
          maxfreq,
          history,        /* number of bits used for history of alignments */
          perc_mat_history, /* percent of matches in history */
          minidentity, /* We prefer to specify the minidentity. The use of the
           notion of error percentage may be misleading, as in Myers paper it
           refers to the percentage of errors in a read. If (for compatibility
           reasons, the option -err is used, then the minidentity contains the
           error rate, a value in the range from 1 to 30. */
          maxalignedlendifference, /* maxfrontdist */
          extendgreedy, /* determines which of the tables in
                           seed-extend-params.h is used */
          alignmentwidth; /* 0 for no alignment display and otherwidth number
                             of columns of alignment per line displayed. */
  bool scanfile, beverbose, forward, reverse, searchspm,
       check_extend_symmetry, silent, trimstat;
  GtStr *indexname, *cam_string; /* parse this using
                                    gt_greedy_extend_char_access*/
  GtStrArray *queryfiles;
  GtOption *refforwardoption, *refseedlengthoption,
           *refuserdefinedleastlengthoption,
           *refextendxdropoption,
           *refextendgreedyoption,
           *refalignmentwidthoption;
} Maxpairsoptions;

static int gt_simpleexactselfmatchoutput(void *info,
                                         const GtGenericEncseq *genericencseq,
                                         GtUword len,
                                         GtUword pos1,
                                         GtUword pos2,
                                         GtError *err)
{
  GtUword queryseqnum, seqstartpos, seqlength;
  const GtEncseq *encseq;

  gt_assert(pos1 < pos2);
  gt_assert(genericencseq != NULL && genericencseq->hasencseq);
  encseq = genericencseq->seqptr.encseq;
  queryseqnum = gt_encseq_seqnum(encseq,pos2);
  seqstartpos = gt_encseq_seqstartpos(encseq,queryseqnum);
  seqlength = gt_encseq_seqlength(encseq,queryseqnum);
  gt_assert(pos2 >= seqstartpos);
  return gt_querymatch_fill_and_output(
                     len,
                     pos1,
                     GT_READMODE_FORWARD,
                     false,
                     0,
                     0,
                     true,
                     (uint64_t) queryseqnum,
                     len,
                     pos2 - seqstartpos,
                     (GtQuerymatchoutoptions *) info,
                     encseq,
                     NULL,
                     seqlength,
                     false,
                     err);
}

static int gt_simplesuffixprefixmatchoutput(GT_UNUSED void *info,
                                            const GtGenericEncseq
                                              *genericencseq,
                                            GtUword matchlen,
                                            GtUword pos1,
                                            GtUword pos2,
                                            GT_UNUSED GtError *err)
{
  GtUword seqnum1, relpos1, seqnum2, relpos2, seqstartpos;
  const GtEncseq *encseq;

  gt_assert(pos1 < pos2);
  gt_assert(genericencseq != NULL && genericencseq->hasencseq);
  encseq = genericencseq->seqptr.encseq;
  seqnum1 = gt_encseq_seqnum(encseq,pos1);
  seqstartpos = gt_encseq_seqstartpos(encseq, seqnum1);
  gt_assert(seqstartpos <= pos1);
  relpos1 = pos1 - seqstartpos;
  seqnum2 = gt_encseq_seqnum(encseq,pos2);
  seqstartpos = gt_encseq_seqstartpos(encseq, seqnum2);
  gt_assert(seqstartpos <= pos2);
  relpos2 = pos2 - seqstartpos;
  if (relpos1 == 0)
  {
    GtUword seqlen2 = gt_encseq_seqlength(encseq,seqnum2);

    if (relpos2 + matchlen == seqlen2)
    {
      printf(GT_WU " " GT_WU " " GT_WU "\n",seqnum2,seqnum1,matchlen);
    }
  } else
  {
    if (relpos2 == 0)
    {
      GtUword seqlen1 = gt_encseq_seqlength(encseq,seqnum1);

      if (relpos1 + matchlen == seqlen1)
      {
        printf(GT_WU " " GT_WU " " GT_WU "\n",seqnum1,seqnum2,matchlen);
      }
    }
  }
  return 0;
}

static void *gt_repfind_arguments_new(void)
{
  Maxpairsoptions *arguments;

  arguments = gt_malloc(sizeof (*arguments));
  arguments->indexname = gt_str_new();
  arguments->cam_string = gt_str_new();
  arguments->queryfiles = gt_str_array_new();
  return arguments;
}

static void gt_repfind_arguments_delete(void *tool_arguments)
{
  Maxpairsoptions *arguments = tool_arguments;

  if (!arguments)
  {
    return;
  }
  gt_str_delete(arguments->indexname);
  gt_str_delete(arguments->cam_string);
  gt_str_array_delete(arguments->queryfiles);
  gt_option_delete(arguments->refforwardoption);
  gt_option_delete(arguments->refseedlengthoption);
  gt_option_delete(arguments->refuserdefinedleastlengthoption);
  gt_option_delete(arguments->refextendxdropoption);
  gt_option_delete(arguments->refextendgreedyoption);
  gt_option_delete(arguments->refalignmentwidthoption);
  gt_free(arguments);
}

static GtOptionParser *gt_repfind_option_parser_new(void *tool_arguments)
{
  GtOptionParser *op;
  GtOption *option, *reverseoption, *queryoption, *extendxdropoption,
           *extendgreedyoption, *scanoption, *sampleoption, *forwardoption,
           *spmoption, *seedlengthoption, *minidentityoption,
           *maxalilendiffoption, *leastlength_option, *char_access_mode_option,
           *check_extend_symmetry_option, *xdropbelowoption, *historyoption,
           *percmathistoryoption, *errorpercentageoption, *optiontrimstat,
           *optionwithalignment;
  Maxpairsoptions *arguments = tool_arguments;

  op = gt_option_parser_new("[options] -ii indexname",
                                  "Compute maximal repeats (and more).");
  gt_option_parser_set_mail_address(op,"<kurtz@zbh.uni-hamburg.de>");

  leastlength_option
    = gt_option_new_uint("l","Specify minimum length of repeats",
                         &arguments->userdefinedleastlength,
                         0U);
  gt_option_parser_add_option(op, leastlength_option);
  arguments->refuserdefinedleastlengthoption
    = gt_option_ref(leastlength_option);

  forwardoption = gt_option_new_bool("f","Compute maximal forward repeats",
                                     &arguments->forward,
                                     true);
  gt_option_parser_add_option(op, forwardoption);
  arguments->refforwardoption = gt_option_ref(forwardoption);

  reverseoption = gt_option_new_bool("r","Compute maximal reverse matches",
                                     &arguments->reverse,
                                     false);
  gt_option_parser_add_option(op, reverseoption);

  seedlengthoption = gt_option_new_uint_min("seedlength",
                                             "Specify minimum length of seed",
                                             &arguments->seedlength,
                                             0,
                                             1UL);
  gt_option_parser_add_option(op, seedlengthoption);
  arguments->refseedlengthoption = gt_option_ref(seedlengthoption);

  option = gt_option_new_uword_min("maxfreq",
                                   "Specify maximal frequency of maximal exact "
                                   "matches in reference sequence",
                                   &arguments->maxfreq,
                                   0,
                                   2UL);
  gt_option_parser_add_option(op, option);

  extendxdropoption
    = gt_option_new_uword_min_max("extendxdrop",
                         "Extend seed to both sides using xdrop algorithm,"
                         "optional parameter specifies sensitivity",
                         &arguments->extendxdrop,
                         93,90,100);
  gt_option_argument_is_optional(extendxdropoption);
  gt_option_parser_add_option(op, extendxdropoption);
  arguments->refextendxdropoption = gt_option_ref(extendxdropoption);

  xdropbelowoption = gt_option_new_word("xdropbelow",
                                        "Specify xdrop cutoff score "
                                        "(argument 0 means undefined). If "
                                        "undefined an optimal value is "
                                        "determined automatically depending "
                                        "on the error "
                                        "rate",
                                        &arguments->xdropbelowscore,
                                        0);
  gt_option_parser_add_option(op, xdropbelowoption);

  extendgreedyoption
    = gt_option_new_uword_min_max("extendgreedy",
                                  "Extend seed to both sides using "
                                  "greedy algorithm with trimming of waves, "
                                  "optional parameter specifies sensitivity",
                                  &arguments->extendgreedy, 93, 90, 100);
  gt_option_argument_is_optional(extendgreedyoption);
  gt_option_parser_add_option(op, extendgreedyoption);
  arguments->refextendgreedyoption = gt_option_ref(extendgreedyoption);

  errorpercentageoption
    = gt_option_new_uword_min_max("err","Specify error percentage of matches "
                         "as integer in the range from 1 to 30 "
                         "(for xdrop and greedy extension) [deprecated option, "
                         "kept for backwards compatibility]",
                         &arguments->minidentity,
                         10,
                         1,
                         100 - GT_EXTEND_MIN_IDENTITY_PERCENTAGE);
  gt_option_is_development_option(errorpercentageoption);
  gt_option_parser_add_option(op, errorpercentageoption);

  minidentityoption
    = gt_option_new_uword_min_max("minidentity",
                          "Specify minimum identity of matches\n"
                          "as integer in the range from 70 to 99 "
                          "(for xdrop and greedy extension)",
                          &arguments->minidentity,
                          80,
                          GT_EXTEND_MIN_IDENTITY_PERCENTAGE,99);
  gt_option_parser_add_option(op, minidentityoption);

  maxalilendiffoption
    = gt_option_new_uword("maxalilendiff","Specify maximum difference of "
                          "alignment length (trimming for greedy extension), "
                          "if option is not used or parameter 0 is specified, "
                          "then good value is automatically chosen",
                          &arguments->maxalignedlendifference,
                          0);
  gt_option_parser_add_option(op, maxalilendiffoption);
  gt_option_is_development_option(maxalilendiffoption);

  historyoption
    = gt_option_new_uword_min_max("history",
                                  "Specify size of history in range [1..64] "
                                  "(trimming for greedy extension)",
                                  &arguments->history,
                                  60,
                                  0,
                                  64);
  gt_option_parser_add_option(op, historyoption);
  gt_option_is_development_option(historyoption);

  percmathistoryoption
    = gt_option_new_uword_min_max("percmathistory",
                                  "percentage of matches required in history",
                                  &arguments->perc_mat_history,
                                  0,
                                  1,
                                  100);
  gt_option_parser_add_option(op, percmathistoryoption);
  gt_option_is_development_option(percmathistoryoption);

  optionwithalignment
    = gt_option_new_uword_min("a",
                              "show alignments/sequences for exact matches "
                              "(optional argument is number of columns per "
                              "line)",
                              &arguments->alignmentwidth,
                              70,
                              20);
  gt_option_argument_is_optional(optionwithalignment);
  gt_option_parser_add_option(op, optionwithalignment);
  arguments->refalignmentwidthoption = gt_option_ref(optionwithalignment);

  char_access_mode_option = gt_option_new_string("cam",
                                                 gt_cam_extendgreedy_comment(),
                                                 arguments->cam_string,"");
  gt_option_parser_add_option(op, char_access_mode_option);

  option = gt_option_new_string("ii",
                                "Specify input index",
                                arguments->indexname, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  char_access_mode_option = gt_option_new_string("cam",
                                                 gt_cam_extendgreedy_comment(),
                                                 arguments->cam_string,"");
  gt_option_parser_add_option(op, char_access_mode_option);
  gt_option_is_development_option(char_access_mode_option);

  option = gt_option_new_bool("silent","do not report matches",
                               &arguments->silent, false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  optiontrimstat = gt_option_new_bool("trimstat","show trimming statistics",
                                      &arguments->trimstat, false);
  gt_option_parser_add_option(op, optiontrimstat);
  gt_option_is_development_option(optiontrimstat);

  /* the following option are options special to repfind */

  queryoption = gt_option_new_filename_array("q",
                                             "Specify query files",
                                             arguments->queryfiles);
  gt_option_is_development_option(queryoption);
  gt_option_parser_add_option(op, queryoption);

  sampleoption = gt_option_new_uword_min("samples","Specify number of samples",
                                         &arguments->samples,
                                         0,
                                         1UL);
  gt_option_is_development_option(sampleoption);
  gt_option_parser_add_option(op, sampleoption);

  spmoption = gt_option_new_bool("spm","Search for suffix prefix matches",
                                       &arguments->searchspm,
                                       false);
  gt_option_is_development_option(spmoption);
  gt_option_parser_add_option(op, spmoption);

  check_extend_symmetry_option
    = gt_option_new_bool("check_extend_symmetry",
                         "check that left/right greedy extension is symmetric "
                         "for sequences mirror around seed",
                         &arguments->check_extend_symmetry,
                         false);
  gt_option_parser_add_option(op, check_extend_symmetry_option);
  gt_option_is_development_option(check_extend_symmetry_option);

  scanoption = gt_option_new_bool("scan","scan index rather than map "
                                         "it to main memory",
                                  &arguments->scanfile,
                                  false);
  gt_option_parser_add_option(op, scanoption);

  option = gt_option_new_verbose(&arguments->beverbose);
  gt_option_parser_add_option(op, option);

  gt_option_exclude(queryoption,sampleoption);
  gt_option_exclude(queryoption,scanoption);
  gt_option_exclude(queryoption,reverseoption);
  gt_option_exclude(queryoption,spmoption);
  gt_option_exclude(sampleoption,spmoption);
  gt_option_exclude(reverseoption,spmoption);
  gt_option_exclude(extendgreedyoption,extendxdropoption);
  gt_option_exclude(errorpercentageoption,minidentityoption);
  gt_option_imply(xdropbelowoption,extendxdropoption);
  gt_option_imply(char_access_mode_option,extendgreedyoption);
  gt_option_imply(historyoption,extendgreedyoption);
  gt_option_imply(maxalilendiffoption,extendgreedyoption);
  gt_option_imply(percmathistoryoption,extendgreedyoption);
  gt_option_imply(optiontrimstat,extendgreedyoption);
  gt_option_imply_either_2(seedlengthoption,extendxdropoption,
                           extendgreedyoption);
  gt_option_imply_either_2(minidentityoption,extendxdropoption,
                           extendgreedyoption);
  gt_option_imply_either_2(errorpercentageoption,extendxdropoption,
                           extendgreedyoption);
  return op;
}

static int gt_repfind_arguments_check(GT_UNUSED int rest_argc,
                                      void *tool_arguments,
                                      GT_UNUSED GtError *err)
{
  Maxpairsoptions *arguments = tool_arguments;

  if (!gt_option_is_set(arguments->refforwardoption) && arguments->reverse)
  {
    arguments->forward = false;
  }
  if (!gt_option_is_set(arguments->refuserdefinedleastlengthoption))
  {
    if (!gt_option_is_set(arguments->refseedlengthoption))
    {
      arguments->seedlength = arguments->userdefinedleastlength = 20U;
    } else
    {
      arguments->userdefinedleastlength = arguments->seedlength;
    }
  } else
  {
    if (!gt_option_is_set(arguments->refseedlengthoption))
    {
      arguments->seedlength = arguments->userdefinedleastlength;
    } else
    {
      if (arguments->seedlength > arguments->userdefinedleastlength)
      {
        arguments->seedlength = arguments->userdefinedleastlength;
      }
    }
  }
  return 0;
}

static int gt_generic_simplexdropselfmatchoutput(
                                           void *info,
                                           const GtGenericEncseq *genericencseq,
                                           GtUword len,
                                           GtUword pos1,
                                           GtUword pos2,
                                           GtError *err)
{
  gt_assert(genericencseq != NULL && genericencseq->hasencseq);
  return gt_simplexdropselfmatchoutput(info,
                                       genericencseq->seqptr.encseq,
                                       len,
                                       pos1,
                                       pos2,
                                       err);
}

static int gt_generic_simplegreedyselfmatchoutput(
                                           void *processinfo,
                                           const GtGenericEncseq *genericencseq,
                                           GtUword len,
                                           GtUword pos1,
                                           GtUword pos2,
                                           GtError *err)
{
  gt_assert(genericencseq != NULL && genericencseq->hasencseq);
  return gt_simplegreedyselfmatchoutput(processinfo,
                                        genericencseq->seqptr.encseq,
                                        len,
                                        pos1,
                                        pos2,
                                        err);
}

static GtUword minidentity2errorpercentage(GtUword minidentity)
{
  if (minidentity >= 1 &&
      minidentity <= 100 - GT_EXTEND_MIN_IDENTITY_PERCENTAGE)
  {
    return minidentity;
  } else
  {
    gt_assert(minidentity >= GT_EXTEND_MIN_IDENTITY_PERCENTAGE);
    return 100 - minidentity;
  }
}

static int gt_repfind_runner(int argc,
                             GT_UNUSED const char **argv,
                             int parsed_args,
                             void *tool_arguments, GtError *err)
{
  bool haserr = false;
  Maxpairsoptions *arguments = (Maxpairsoptions *) tool_arguments;
  GtLogger *logger = NULL;
  GtXdropmatchinfo *xdropmatchinfo = NULL;
  GtGreedyextendmatchinfo *greedyextendmatchinfo = NULL;
  GtTimer *repfindtimer = NULL;
  GtExtendCharAccess extend_char_access = GT_EXTEND_CHAR_ACCESS_ANY;

  gt_error_check(err);
  logger = gt_logger_new(arguments->beverbose, GT_LOGGER_DEFLT_PREFIX, stdout);
  if (gt_showtime_enabled())
  {
    repfindtimer = gt_timer_new();
    gt_timer_start(repfindtimer);
  }
  if (parsed_args < argc)
  {
    gt_error_set(err,"superfluous arguments: \"%s\"",argv[argc-1]);
    haserr = true;
  }
  if (!haserr && !gt_option_is_set(arguments->refalignmentwidthoption))
  {
    arguments->alignmentwidth = 0;
  }
  if (!haserr && gt_option_is_set(arguments->refextendxdropoption))
  {
    xdropmatchinfo
      = gt_xdrop_matchinfo_new(arguments->userdefinedleastlength,
                               minidentity2errorpercentage(
                                            arguments->minidentity),
                               arguments->xdropbelowscore,
                               arguments->extendxdrop,
                               gt_str_array_size(arguments->queryfiles) == 0
                                             ? true : false);
    gt_assert(xdropmatchinfo != NULL);
    if (arguments->beverbose)
    {
      gt_xdrop_matchinfo_verbose_set(xdropmatchinfo);
    }
    if (arguments->silent)
    {
      gt_xdrop_matchinfo_silent_set(xdropmatchinfo);
    }
  }
  if (!haserr && gt_option_is_set(arguments->refextendgreedyoption))
  {
    extend_char_access
      = gt_greedy_extend_char_access(gt_str_get(arguments->cam_string),err);

    if ((int) extend_char_access == -1)
    {
      haserr = true;
    }
    if (!haserr)
    {
      greedyextendmatchinfo
        = gt_greedy_extend_matchinfo_new(minidentity2errorpercentage(
                                             arguments->minidentity),
                                         arguments->maxalignedlendifference,
                                         arguments->history,
                                         arguments->perc_mat_history,
                                         arguments->userdefinedleastlength,
                                         extend_char_access,
                                         arguments->extendgreedy);
      if (arguments->beverbose)
      {
        gt_greedy_extend_matchinfo_verbose_set(greedyextendmatchinfo);
      }
      if (arguments->check_extend_symmetry)
      {
        gt_greedy_extend_matchinfo_check_extend_symmetry_set(
                                                  greedyextendmatchinfo);
      }
      if (arguments->silent)
      {
        gt_greedy_extend_matchinfo_silent_set(greedyextendmatchinfo);
      }
      if (arguments->trimstat)
      {
        gt_greedy_extend_matchinfo_trimstat_set(greedyextendmatchinfo);
      }
    }
  }
  if (!haserr)
  {
    GtProcessinfo_and_outoptions processinfo_and_outoptions;
    GtUword sensitivity = gt_option_is_set(arguments->refextendgreedyoption)
                            ? arguments->extendgreedy
                            : 100;

    processinfo_and_outoptions.processinfo = NULL;
    if (arguments->alignmentwidth > 0)
    {
      processinfo_and_outoptions.querymatchoutoptions
        = gt_querymatchoutoptions_new(arguments->alignmentwidth,
                                      minidentity2errorpercentage(
                                               arguments->minidentity),
                                      arguments->maxalignedlendifference,
                                      arguments->history,
                                      arguments->perc_mat_history,
                                      extend_char_access,
                                      sensitivity);
    } else
    {
      processinfo_and_outoptions.querymatchoutoptions = NULL;
    }
    if (gt_str_array_size(arguments->queryfiles) == 0)
    {
      if (arguments->samples == 0)
      {
        if (arguments->forward)
        {
          GtProcessmaxpairs processmaxpairs;
          void *processmaxpairsdata;

          if (arguments->searchspm)
          {
            processmaxpairs = gt_simplesuffixprefixmatchoutput;
            processmaxpairsdata = NULL;
          } else
          {
            if (gt_option_is_set(arguments->refextendxdropoption))
            {
              processmaxpairs = gt_generic_simplexdropselfmatchoutput;
              processinfo_and_outoptions.processinfo = (void *) xdropmatchinfo;
              processmaxpairsdata = (void *) &processinfo_and_outoptions;
            } else
            {
              if (gt_option_is_set(arguments->refextendgreedyoption))
              {
                processmaxpairs = gt_generic_simplegreedyselfmatchoutput;
                processinfo_and_outoptions.processinfo
                  = (void *) greedyextendmatchinfo;
                processmaxpairsdata = (void *) &processinfo_and_outoptions;
              } else
              {
                processmaxpairs = gt_simpleexactselfmatchoutput;
                processmaxpairsdata
                  = (void *) processinfo_and_outoptions.querymatchoutoptions;
              }
            }
          }
          if (gt_callenummaxpairs(gt_str_get(arguments->indexname),
                                  arguments->seedlength,
                                  arguments->maxfreq,
                                  arguments->scanfile,
                                  processmaxpairs,
                                  processmaxpairsdata,
                                  logger,
                                  err) != 0)
          {
            haserr = true;
          }
        }
        if (!haserr && arguments->reverse)
        {
          if (gt_callenumselfmatches(gt_str_get(arguments->indexname),
                                     GT_READMODE_REVERSE,
                                     arguments->seedlength,
                                     gt_querymatch_output,
                                     processinfo_and_outoptions.
                                         querymatchoutoptions,
                                     logger,
                                     err) != 0)
          {
            haserr = true;
          }
        }
      } else
      {
        if (gt_testmaxpairs(gt_str_get(arguments->indexname),
                            arguments->samples,
                            arguments->seedlength,
                            (GtUword)
                            (100 * arguments->seedlength),
                            logger,
                            err) != 0)
        {
          haserr = true;
        }
      }
    } else
    {
      GtProcessquerymatch processquerymatch = NULL;
      void *processquerymatch_data = NULL;

      if (gt_option_is_set(arguments->refextendxdropoption))
      {
        processquerymatch = gt_processxdropquerymatches;
        processinfo_and_outoptions.processinfo = xdropmatchinfo;
        processquerymatch_data = (void *) &processinfo_and_outoptions;
      } else
      {
        if (gt_option_is_set(arguments->refextendgreedyoption))
        {
          gt_assert(false);
        } else
        {
          processquerymatch = gt_querymatch_output;
          processquerymatch_data
            = processinfo_and_outoptions.querymatchoutoptions;
        }
      }
      if (gt_callenumquerymatches(gt_str_get(arguments->indexname),
                                  arguments->queryfiles,
                                  false,
                                  true,
                                  false,
                                  arguments->seedlength,
                                  NULL,
                                  processquerymatch,
                                  processquerymatch_data,
                                  logger,
                                  err) != 0)
      {
        haserr = true;
      }
    }
    gt_querymatchoutoptions_delete(processinfo_and_outoptions.
                                   querymatchoutoptions);
  }
  gt_xdrop_matchinfo_delete(xdropmatchinfo);
  gt_greedy_extend_matchinfo_delete(greedyextendmatchinfo);
  gt_logger_delete(logger);
  if (repfindtimer != NULL)
  {
    printf("# TIME repfind");
    if (gt_option_is_set(arguments->refextendgreedyoption))
    {
      printf("-greedy");
    } else
    {
      if (gt_option_is_set(arguments->refextendxdropoption))
      {
        printf("-xdrop");
      }
    }
    printf("-%u-%u",arguments->seedlength,arguments->userdefinedleastlength);
    if (gt_option_is_set(arguments->refextendgreedyoption) ||
        gt_option_is_set(arguments->refextendxdropoption))
    {
      printf("-" GT_WU,
             100 - minidentity2errorpercentage(arguments->minidentity));
    }
    if (gt_option_is_set(arguments->refextendgreedyoption))
    {
      GtUword maxalignedlendifference, perc_mat_history;
      gt_optimal_maxalilendiff_perc_mat_history(
                &maxalignedlendifference,
                &perc_mat_history,
                arguments->maxalignedlendifference,
                arguments->perc_mat_history,
                minidentity2errorpercentage(arguments->minidentity),
                arguments->extendgreedy);
      printf("-" GT_WU,maxalignedlendifference);
      printf("-" GT_WU,perc_mat_history);
    } else
    {
      if (gt_option_is_set(arguments->refextendxdropoption))
      {
        if (arguments->xdropbelowscore == 0)
        {
          printf("-" GT_WD,
                 gt_optimalxdropbelowscore(
                       minidentity2errorpercentage(arguments->minidentity),
                       arguments->extendxdrop));
        } else
        {
          printf("-" GT_WD,arguments->xdropbelowscore);
        }
      }
    }
    gt_timer_show_formatted(repfindtimer," overall " GT_WD ".%02ld\n",stdout);
    gt_timer_delete(repfindtimer);
  }
  return haserr ? -1 : 0;
}

GtTool* gt_repfind(void)
{
  return gt_tool_new(gt_repfind_arguments_new,
                     gt_repfind_arguments_delete,
                     gt_repfind_option_parser_new,
                     gt_repfind_arguments_check,
                     gt_repfind_runner);
}
