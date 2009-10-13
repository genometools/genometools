/*
  Copyright (c) 2007-2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2009 Center for Bioinformatics, University of Hamburg

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
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "core/option.h"
#include "core/tool.h"
#include "core/versionfunc.h"
#include "match/esa-seqread.h"
#include "match/esa-mmsearch.h"
#include "match/format64.h"
#include "match/verbose-def.h"
#include "match/esa-maxpairs.h"
#include "match/test-maxpairs.pr"
#include "tools/gt_maxpairs.h"

typedef struct
{
  unsigned int userdefinedleastlength;
  unsigned long samples;
  bool scanfile, beverbose, forward, reverse;
  GtStr *indexname;
  GtStrArray *queryfiles;
  GtOption *refforwardoption;
} Maxpairsoptions;

static int simpleexactmatchoutput(GT_UNUSED void *info,
                                  const Encodedsequence *encseq,
                                  Seqpos len,
                                  Seqpos absdbstart,
                                  Readmode readmode,
                                  uint64_t queryunitnum,
                                  Seqpos querystart,
                                  Seqpos querytotallength,
                                  GT_UNUSED GtError *err)
{
  const char *outflag = "FRCP";
  Seqinfo seqinfo;
  unsigned long dbseqnum;

  gt_assert(encseq != NULL);
  dbseqnum = getencseqfrompos2seqnum(encseq,absdbstart);
  getencseqSeqinfo(&seqinfo,encseq,dbseqnum);
  gt_assert((int) readmode < 4);
  if (readmode == Reversemode)
  {
    gt_assert(querystart + len <= querytotallength);
    /*
    printf("totallength = %lu, start=%lu\n",(unsigned long) querytotallength,
                                           (unsigned long) querystart);
    */
    querystart = querytotallength - querystart - len;
  }
  printf(FormatSeqpos " %lu " FormatSeqpos " %c " FormatSeqpos
         " " Formatuint64_t " " FormatSeqpos "\n",
         PRINTSeqposcast(len),
         dbseqnum,
         PRINTSeqposcast(absdbstart - seqinfo.seqstartpos),
         outflag[readmode],
         PRINTSeqposcast(len),
         PRINTuint64_tcast(queryunitnum),
         PRINTSeqposcast(querystart));
  return 0;
}

static int simpleexactselfmatchoutput(void *info,
                                      const Encodedsequence *encseq,
                                      Seqpos len,
                                      Seqpos pos1,
                                      Seqpos pos2,
                                      GT_UNUSED GtError *err)
{
  Seqinfo seqinfo;
  unsigned long queryseqnum;

  if (pos1 > pos2)
  {
    Seqpos tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }
  queryseqnum = getencseqfrompos2seqnum(encseq,pos2);
  getencseqSeqinfo(&seqinfo,encseq,queryseqnum);
  gt_assert(pos2 >= seqinfo.seqstartpos);
  return simpleexactmatchoutput(info,
                                encseq,
                                len,
                                pos1,
                                Forwardmode,
                                (uint64_t) queryseqnum,
                                pos2 - seqinfo.seqstartpos,
                                seqinfo.seqlength,
                                err);
}

static int callenummaxpairs(const GtStr *indexname,
                            unsigned int userdefinedleastlength,
                            bool scanfile,
                            Processmaxpairs processmaxpairs,
                            void *processmaxpairsinfo,
                            Verboseinfo *verboseinfo,
                            GtError *err)
{
  bool haserr = false;
  Sequentialsuffixarrayreader *ssar;

  gt_error_check(err);
  ssar = newSequentialsuffixarrayreaderfromfile(indexname,
                                                SARR_LCPTAB |
                                                SARR_SUFTAB |
                                                SARR_ESQTAB |
                                                SARR_SSPTAB,
                                                scanfile
                                                  ? SEQ_scan : SEQ_mappedboth,
                                                err);
  if (ssar == NULL)
  {
    haserr = true;
  }
  if (!haserr &&
      enumeratemaxpairs(ssar,
                        encseqSequentialsuffixarrayreader(ssar),
                        readmodeSequentialsuffixarrayreader(ssar),
                        userdefinedleastlength,
                        processmaxpairs,
                        processmaxpairsinfo,
                        verboseinfo,
                        err) != 0)
  {
    haserr = true;
  }
  if (ssar != NULL)
  {
    freeSequentialsuffixarrayreader(&ssar);
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
  GtOption *option, *reverseoption, *queryoption,
           *scanoption, *sampleoption, *forwardoption;
  Maxpairsoptions *arguments = tool_arguments;

  op = gt_option_parser_new("[options] -ii indexname",
                            "Compute maximal repeats.");
  gt_option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");

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

  queryoption = gt_option_new_filenamearray("q",
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
  Verboseinfo *verboseinfo;

  gt_error_check(err);
  verboseinfo = newverboseinfo(arguments->beverbose);
  if (parsed_args < argc)
  {
    gt_error_set(err,"superfluous arguments: \"%s\"\n",argv[argc-1]);
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
          if (callenummaxpairs(arguments->indexname,
                               arguments->userdefinedleastlength,
                               arguments->scanfile,
                               simpleexactselfmatchoutput,
                               NULL,
                               verboseinfo,
                               err) != 0)
          {
            haserr = true;
          }
        }
        if (!haserr && arguments->reverse)
        {
          if (callenumselfmatches(arguments->indexname,
                                  Reversemode,
                                  arguments->userdefinedleastlength,
                                  simpleexactmatchoutput,
                                  NULL,
                                  verboseinfo,
                                  err) != 0)
          {
            haserr = true;
          }
        }
      } else
      {
        if (testmaxpairs(arguments->indexname,
                         arguments->samples,
                         arguments->userdefinedleastlength,
                         (Seqpos) (100 * arguments->userdefinedleastlength),
                         verboseinfo,
                         err) != 0)
        {
          haserr = true;
        }
      }
    } else
    {
      if (callenumquerymatches(arguments->indexname,
                               arguments->queryfiles,
                               false,
                               arguments->userdefinedleastlength,
                               simpleexactmatchoutput,
                               NULL,
                               verboseinfo,
                               err) != 0)
      {
        haserr = true;
      }
    }
  }
  freeverboseinfo(&verboseinfo);
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
