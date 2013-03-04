/*
  Copyright (c) 2007-2013 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2013 Center for Bioinformatics, University of Hamburg

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
#include "core/error_api.h"
#include "core/str_api.h"
#include "core/logger.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "core/option_api.h"
#include "core/tool_api.h"
#include "core/versionfunc.h"
#include "core/format64.h"
#include "match/esa-seqread.h"
#include "match/esa-mmsearch.h"
#include "match/esa-maxpairs.h"
#include "match/test-maxpairs.pr"
#include "match/querymatch.h"
#include "match/greedyedist.h"
#include "match/xdrop.h"
#include "tools/gt_maxpairs.h"

typedef struct
{
  unsigned int userdefinedleastlength;
  unsigned long samples;
  bool scanfile, beverbose, forward, reverse, searchspm, extendseed;
  GtStr *indexname;
  GtStrArray *queryfiles;
  GtOption *refforwardoption;
} Maxpairsoptions;

static int gt_simpleexactselfmatchoutput(void *info,
                                         const GtEncseq *encseq,
                                         unsigned long len,
                                         unsigned long pos1,
                                         unsigned long pos2,
                                         GT_UNUSED GtError *err)
{
  unsigned long queryseqnum, seqstartpos, seqlength;
  GtQuerymatch *querymatch = (GtQuerymatch *) info;

  if (pos1 > pos2)
  {
    unsigned long tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }
  queryseqnum = gt_encseq_seqnum(encseq,pos2);
  seqstartpos = gt_encseq_seqstartpos(encseq, queryseqnum);
  seqlength = gt_encseq_seqlength(encseq, queryseqnum);
  gt_assert(pos2 >= seqstartpos);
  gt_querymatch_fill(querymatch,
                     len,
                     pos1,
                     GT_READMODE_FORWARD,
                     false,
                     0,
                     0,
                     true,
                     (uint64_t) queryseqnum,
                     len,
                     pos2 - seqstartpos);
  return gt_querymatch_output(info, encseq, querymatch, NULL, seqlength, err);
}

typedef struct
{
  GtQuerymatch *querymatchspaceptr;
  GtXdropArbitraryscores arbitscores;
  GtXdropresources *res;
  GtFrontResource *frontresource;
  GtXdropbest best_left;
  GtXdropbest best_right;
  GtXdropscore belowscore;
  GtSeqabstract *useq, *vseq;
  const GtUchar *query_sequence;
  unsigned long query_totallength;
} GtXdropmatchinfo;

static int gt_simplexdropselfmatchoutput(void *info,
                                         const GtEncseq *encseq,
                                         unsigned long len,
                                         unsigned long pos1,
                                         unsigned long pos2,
                                         GtError *err)
{
  GtXdropmatchinfo *xdropmatchinfo = (GtXdropmatchinfo *) info;
  GtXdropscore score;
  unsigned long dbseqnum, dbseqstartpos, dbseqlength, dbstart, dblen,
                querystart, queryseqnum, querylen, queryseqlength,
                queryseqstartpos;
  const unsigned long dbtotallength = gt_encseq_total_length(encseq);

  if (pos1 > pos2)
  {
    unsigned long tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }
  dbseqnum = gt_encseq_seqnum(encseq,pos1),
  dbseqstartpos = gt_encseq_seqstartpos(encseq,dbseqnum),
  dbseqlength = gt_encseq_seqlength(encseq,dbseqnum);

  if (pos2 < dbseqstartpos + dbseqlength)
  {
    queryseqnum = dbseqnum;
    queryseqstartpos = dbseqstartpos;
    queryseqlength = dbseqlength;
  } else
  {
    queryseqnum = gt_encseq_seqnum(encseq,pos2);
    gt_assert(dbseqnum < queryseqnum);
    queryseqstartpos = gt_encseq_seqstartpos(encseq,queryseqnum);
    queryseqlength = gt_encseq_seqlength(encseq,queryseqnum);
  }
  if (pos1 > 0 && pos2 > 0)
  {
    gt_assert(pos1 >= dbseqstartpos && pos2 >= queryseqstartpos);
    gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,encseq,
                                 pos1 - dbseqstartpos,0);
    gt_seqabstract_reinit_encseq(xdropmatchinfo->vseq,encseq,
                                 pos2 - queryseqstartpos,0);
    gt_evalxdroparbitscoresextend(false,
                                  &xdropmatchinfo->best_left,
                                  xdropmatchinfo->res,
                                  xdropmatchinfo->useq,
                                  xdropmatchinfo->vseq,
                                  pos1,
                                  pos2,
                                  xdropmatchinfo->belowscore);
  } else
  {
    xdropmatchinfo->best_left.ivalue = 0;
    xdropmatchinfo->best_left.jvalue = 0;
    xdropmatchinfo->best_left.score = 0;
  }
  if (pos1 + len < dbtotallength && pos2 + len < dbtotallength)
  {
    const unsigned long seqend1 = dbseqstartpos + dbseqlength;
    const unsigned long seqend2 = queryseqstartpos + queryseqlength;

    gt_assert(seqend1 >= pos1 + len && seqend2 >= pos2 + len);
    gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,
                                 encseq,seqend1 - (pos1 + len),0);
    gt_seqabstract_reinit_encseq(xdropmatchinfo->vseq,
                                 encseq,seqend2 - (pos2 + len),0);
    gt_evalxdroparbitscoresextend(true,
                                  &xdropmatchinfo->best_right,
                                  xdropmatchinfo->res,
                                  xdropmatchinfo->useq,
                                  xdropmatchinfo->vseq,
                                  pos1 + len,
                                  pos2 + len,
                                  xdropmatchinfo->belowscore);
  } else
  {
    xdropmatchinfo->best_right.ivalue = 0;
    xdropmatchinfo->best_right.jvalue = 0;
    xdropmatchinfo->best_right.score = 0;
  }
  gt_assert(pos1 >= (unsigned long) xdropmatchinfo->best_left.ivalue &&
            pos2 >= (unsigned long) xdropmatchinfo->best_left.jvalue);
  querystart = pos2 - xdropmatchinfo->best_left.jvalue;
  gt_assert(querystart >= queryseqstartpos);
  dblen = len + xdropmatchinfo->best_left.ivalue
              + xdropmatchinfo->best_right.ivalue;
  dbstart = pos1 - xdropmatchinfo->best_left.ivalue;
  querylen = len + xdropmatchinfo->best_left.jvalue
                 + xdropmatchinfo->best_right.jvalue,
  score = (GtXdropscore) len * xdropmatchinfo->arbitscores.mat +
          xdropmatchinfo->best_left.score +
          xdropmatchinfo->best_right.score;
  gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,
                               encseq,
                               dblen,
                               dbstart);
  gt_seqabstract_reinit_encseq(xdropmatchinfo->vseq,
                               encseq,
                               querylen,
                               querystart);
  gt_querymatch_fill(xdropmatchinfo->querymatchspaceptr,
                     dblen,
                     dbstart,
                     GT_READMODE_FORWARD,
                     false,
                     score,
                     greedyunitedist(xdropmatchinfo->frontresource,
                                     xdropmatchinfo->useq,xdropmatchinfo->vseq),
                     true,
                     (uint64_t) queryseqnum,
                     querylen,
                     querystart - queryseqstartpos);
  return gt_querymatch_output(info, encseq, xdropmatchinfo->querymatchspaceptr,
                              NULL, gt_encseq_seqlength(encseq, queryseqnum),
                              err);
}

static int gt_processxdropquerymatches(void *info,
                                       const GtEncseq *encseq,
                                       const GtQuerymatch *querymatch,
                                       const GtUchar *query,
                                       unsigned long query_totallength,
                                       GtError *err)
{
  GtXdropmatchinfo *xdropmatchinfo = (GtXdropmatchinfo *) info;
  GtXdropscore score;
  unsigned long querystart, dblen, dbstart, querylen;
  unsigned long pos1 = gt_querymatch_dbstart(querymatch);
  unsigned long pos2 = gt_querymatch_querystart(querymatch);
  unsigned long len = gt_querymatch_querylen(querymatch);
  const unsigned long dbtotallength = gt_encseq_total_length(encseq);
  uint64_t queryseqnum;
  unsigned long dbseqnum, dbseqstartpos, dbseqlength;

  dbseqnum = gt_encseq_seqnum(encseq,pos1);
  dbseqstartpos = gt_encseq_seqstartpos(encseq,dbseqnum);
  dbseqlength = gt_encseq_seqlength(encseq,dbseqnum);
  if (pos1 > 0 && pos2 > 0)
  {
    gt_assert(dbseqstartpos < pos1);
    gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,encseq,
                                 pos1 - dbseqstartpos,0);
    gt_seqabstract_reinit_ptr(xdropmatchinfo->vseq,query,pos2,0);
    gt_evalxdroparbitscoresextend(false,
                                  &xdropmatchinfo->best_left,
                                  xdropmatchinfo->res,
                                  xdropmatchinfo->useq,
                                  xdropmatchinfo->vseq,
                                  pos1,
                                  pos2,
                                  xdropmatchinfo->belowscore);
  } else
  {
    xdropmatchinfo->best_left.ivalue = 0;
    xdropmatchinfo->best_left.jvalue = 0;
    xdropmatchinfo->best_left.score = 0;
  }
  if (pos1 + len < dbtotallength && pos2 + len < query_totallength)
  {
    gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,
                                 encseq,dbseqstartpos + dbseqlength -
                                        (pos1 + len),0);
    gt_seqabstract_reinit_ptr(xdropmatchinfo->vseq,
                              query,query_totallength - (pos2 + len),0);
    gt_evalxdroparbitscoresextend(true,
                                  &xdropmatchinfo->best_right,
                                  xdropmatchinfo->res,
                                  xdropmatchinfo->useq,
                                  xdropmatchinfo->vseq,
                                  pos1 + len,
                                  pos2 + len,
                                  xdropmatchinfo->belowscore);
  } else
  {
    xdropmatchinfo->best_right.ivalue = 0;
    xdropmatchinfo->best_right.jvalue = 0;
    xdropmatchinfo->best_right.score = 0;
  }
  gt_assert(pos1 >= (unsigned long) xdropmatchinfo->best_left.ivalue &&
            pos2 >= (unsigned long) xdropmatchinfo->best_left.jvalue);
  querystart = pos2 - xdropmatchinfo->best_left.jvalue;
  queryseqnum = gt_querymatch_queryseqnum(querymatch);
  dblen = len + xdropmatchinfo->best_left.ivalue
              + xdropmatchinfo->best_right.ivalue;
  dbstart = pos1 - xdropmatchinfo->best_left.ivalue;
  querylen = len + xdropmatchinfo->best_left.jvalue
                 + xdropmatchinfo->best_right.jvalue,
  score = (GtXdropscore) len * xdropmatchinfo->arbitscores.mat +
          xdropmatchinfo->best_left.score +
          xdropmatchinfo->best_right.score;
  gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,
                               encseq,
                               dblen,
                               dbstart);
  gt_seqabstract_reinit_ptr(xdropmatchinfo->vseq,
                            query,
                            querylen,
                            querystart);
  gt_querymatch_fill(xdropmatchinfo->querymatchspaceptr,
                     dblen,
                     dbstart,
                     GT_READMODE_FORWARD,
                     false,
                     score,
                     greedyunitedist(xdropmatchinfo->frontresource,
                                     xdropmatchinfo->useq,xdropmatchinfo->vseq),
                     false,
                     queryseqnum,
                     querylen,
                     querystart);
  return gt_querymatch_output(info, encseq, xdropmatchinfo->querymatchspaceptr,
                              query, query_totallength,
                              err);
}

static int gt_simplesuffixprefixmatchoutput(GT_UNUSED void *info,
                                            const GtEncseq *encseq,
                                            unsigned long matchlen,
                                            unsigned long pos1,
                                            unsigned long pos2,
                                            GT_UNUSED GtError *err)
{
  unsigned long seqnum1, relpos1, seqnum2, relpos2, seqstartpos;

  if (pos1 > pos2)
  {
    unsigned long tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }
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
    unsigned long seqlen2 = gt_encseq_seqlength(encseq,seqnum2);

    if (relpos2 + matchlen == seqlen2)
    {
      printf("%lu %lu %lu\n",seqnum2,seqnum1,matchlen);
    }
  } else
  {
    if (relpos2 == 0)
    {
      unsigned long seqlen1 = gt_encseq_seqlength(encseq,seqnum1);

      if (relpos1 + matchlen == seqlen1)
      {
        printf("%lu %lu %lu\n",seqnum1,seqnum2,matchlen);
      }
    }
  }
  return 0;
}

static int callenummaxpairs(const char *indexname,
                            unsigned int userdefinedleastlength,
                            bool scanfile,
                            Processmaxpairs processmaxpairs,
                            void *processmaxpairsinfo,
                            GtLogger *logger,
                            GtError *err)
{
  bool haserr = false;
  Sequentialsuffixarrayreader *ssar;

  gt_error_check(err);
  ssar = gt_newSequentialsuffixarrayreaderfromfile(indexname,
                                                   SARR_LCPTAB |
                                                   SARR_SUFTAB |
                                                   SARR_ESQTAB |
                                                   SARR_SSPTAB,
                                                   scanfile ? SEQ_scan
                                                            : SEQ_mappedboth,
                                                   logger,
                                                   err);
  if (ssar == NULL)
  {
    haserr = true;
  }
  if (!haserr &&
      gt_enumeratemaxpairs(ssar,
                           gt_encseqSequentialsuffixarrayreader(ssar),
                           gt_readmodeSequentialsuffixarrayreader(ssar),
                           userdefinedleastlength,
                           processmaxpairs,
                           processmaxpairsinfo,
                           err) != 0)
  {
    haserr = true;
  }
  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}

static void *gt_repfind_arguments_new(void)
{
  Maxpairsoptions *arguments;

  arguments = gt_malloc(sizeof (*arguments));
  arguments->indexname = gt_str_new();
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
  gt_str_array_delete(arguments->queryfiles);
  gt_option_delete(arguments->refforwardoption);
  gt_free(arguments);
}

static GtOptionParser *gt_repfind_option_parser_new(void *tool_arguments)
{
  GtOptionParser *op;
  GtOption *option, *reverseoption, *queryoption, *extendseedoption,
           *scanoption, *sampleoption, *forwardoption, *spmoption;
  Maxpairsoptions *arguments = tool_arguments;

  op = gt_option_parser_new("[options] -ii indexname",
                            "Compute maximal repeats.");
  gt_option_parser_set_mail_address(op,"<kurtz@zbh.uni-hamburg.de>");

  option = gt_option_new_uint_min("l","Specify minimum length of repeats",
                                  &arguments->userdefinedleastlength,
                                  20U,
                                  1U);
  gt_option_parser_add_option(op, option);

  forwardoption = gt_option_new_bool("f","Compute maximal forward repeats",
                                     &arguments->forward,
                                     true);
  gt_option_parser_add_option(op, forwardoption);
  arguments->refforwardoption = gt_option_ref(forwardoption);

  reverseoption = gt_option_new_bool("r","Compute maximal reverse matches",
                                     &arguments->reverse,
                                     false);
  gt_option_parser_add_option(op, reverseoption);

  sampleoption = gt_option_new_ulong_min("samples","Specify number of samples",
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

  extendseedoption = gt_option_new_bool("extend","Extend seed to both sides",
                                        &arguments->extendseed,
                                    false);
  gt_option_is_development_option(extendseedoption);
  gt_option_parser_add_option(op, extendseedoption);

  scanoption = gt_option_new_bool("scan","scan index rather than mapping "
                                         "it to main memory",
                                  &arguments->scanfile,
                                  false);
  gt_option_parser_add_option(op, scanoption);

  option = gt_option_new_string("ii",
                                "Specify input index",
                                arguments->indexname, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  queryoption = gt_option_new_filename_array("q",
                                             "Specify query files",
                                             arguments->queryfiles);
  gt_option_is_development_option(queryoption);
  gt_option_parser_add_option(op, queryoption);

  option = gt_option_new_bool("v",
                              "be verbose ",
                              &arguments->beverbose,
                              false);
  gt_option_parser_add_option(op, option);

  gt_option_exclude(queryoption,sampleoption);
  gt_option_exclude(queryoption,scanoption);
  gt_option_exclude(queryoption,reverseoption);
  gt_option_exclude(queryoption,spmoption);
  gt_option_exclude(reverseoption,spmoption);
  gt_option_exclude(queryoption,spmoption);
  gt_option_exclude(sampleoption,spmoption);
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
  return 0;
}

static int gt_repfind_runner(GT_UNUSED int argc,
                             GT_UNUSED const char **argv,
                             GT_UNUSED int parsed_args,
                             void *tool_arguments, GtError *err)
{
  bool haserr = false;
  Maxpairsoptions *arguments = tool_arguments;
  GtLogger *logger = NULL;
  GtQuerymatch *querymatchspaceptr = gt_querymatch_new();
  GtXdropmatchinfo xdropmatchinfo;

  gt_error_check(err);
  xdropmatchinfo.querymatchspaceptr = querymatchspaceptr;
  xdropmatchinfo.useq = gt_seqabstract_new_empty();
  xdropmatchinfo.vseq = gt_seqabstract_new_empty();
  xdropmatchinfo.arbitscores.mat = 2;
  xdropmatchinfo.arbitscores.mis = -2;
  xdropmatchinfo.arbitscores.ins = -3;
  xdropmatchinfo.arbitscores.del = -3;
  xdropmatchinfo.frontresource = gt_frontresource_new(100UL);
  xdropmatchinfo.res = gt_xdrop_resources_new(&xdropmatchinfo.arbitscores);
  xdropmatchinfo.belowscore = 5L;
  logger = gt_logger_new(arguments->beverbose, GT_LOGGER_DEFLT_PREFIX, stdout);
  if (parsed_args < argc)
  {
    gt_error_set(err,"superfluous arguments: \"%s\"",argv[argc-1]);
    haserr = true;
  }
  if (!haserr)
  {
    if (gt_str_array_size(arguments->queryfiles) == 0)
    {
      if (arguments->samples == 0)
      {
        if (arguments->forward)
        {
          if (callenummaxpairs(gt_str_get(arguments->indexname),
                               arguments->userdefinedleastlength,
                               arguments->scanfile,
                               arguments->searchspm
                                 ? gt_simplesuffixprefixmatchoutput
                                 : (arguments->extendseed
                                     ? gt_simplexdropselfmatchoutput
                                     : gt_simpleexactselfmatchoutput),
                               arguments->extendseed
                                 ? (void *) &xdropmatchinfo
                                 : (void *) querymatchspaceptr,
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
                                     arguments->userdefinedleastlength,
                                     /*arguments->extendseed
                                       ? gt_processxdropquerymatches
                                       :*/ gt_querymatch_output,
                                     /*arguments->extendseed
                                       ? (void *) &xdropmatchinfo
                                       :*/ NULL,
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
                            arguments->userdefinedleastlength,
                            (unsigned long)
                            (100 * arguments->userdefinedleastlength),
                            logger,
                            err) != 0)
        {
          haserr = true;
        }
      }
    } else
    {
      if (gt_callenumquerymatches(gt_str_get(arguments->indexname),
                                  arguments->queryfiles,
                                  false,
                                  true,
                                  false,
                                  arguments->userdefinedleastlength,
                                  NULL,
                                  arguments->extendseed
                                    ? gt_processxdropquerymatches
                                    : gt_querymatch_output,
                                  arguments->extendseed
                                    ? (void *) &xdropmatchinfo
                                    : NULL,
                                  logger,
                                  err) != 0)
      {
        haserr = true;
      }
    }
  }
  gt_querymatch_delete(querymatchspaceptr);
  gt_seqabstract_delete(xdropmatchinfo.useq);
  gt_seqabstract_delete(xdropmatchinfo.vseq);
  gt_xdrop_resources_delete(xdropmatchinfo.res);
  gt_frontresource_delete(xdropmatchinfo.frontresource);
  gt_logger_delete(logger);
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
