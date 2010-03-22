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
#include "core/logger.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "core/option.h"
#include "core/tool.h"
#include "core/versionfunc.h"
#include "match/esa-seqread.h"
#include "match/esa-mmsearch.h"
#include "core/format64.h"
#include "match/esa-maxpairs.h"
#include "match/test-maxpairs.pr"
#include "match/querymatch.h"
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

static int simpleexactselfmatchoutput(void *info,
                                      const GtEncodedsequence *encseq,
                                      unsigned long len,
                                      unsigned long pos1,
                                      unsigned long pos2,
                                      GT_UNUSED GtError *err)
{
  GtSeqinfo seqinfo;
  unsigned long queryseqnum;
  Querymatch *querymatch = (Querymatch *) info;

  if (pos1 > pos2)
  {
    unsigned long tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }
  queryseqnum = getencseqfrompos2seqnum(encseq,pos2);
  gt_encodedsequence_seqinfo(encseq,&seqinfo,queryseqnum);
  gt_assert(pos2 >= seqinfo.seqstartpos);
  querymatch_fill(querymatch,
                  len,
                  pos1,
                  GT_READMODE_FORWARD,
                  true,
                  (uint64_t) queryseqnum,
                  pos2 - seqinfo.seqstartpos,
                  seqinfo.seqlength);
  return querymatch_output(info, encseq, querymatch, err);
}

static int callenummaxpairs(const GtStr *indexname,
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
                        logger,
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
  GtLogger *logger = NULL;
  Querymatch *querymatchspaceptr = querymatch_new();

  gt_error_check(err);
  logger = gt_logger_new(arguments->beverbose, GT_LOGGER_DEFLT_PREFIX, stdout);
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
                               querymatchspaceptr,
                               logger,
                               err) != 0)
          {
            haserr = true;
          }
        }
        if (!haserr && arguments->reverse)
        {
          if (callenumselfmatches(arguments->indexname,
                                  GT_READMODE_REVERSE,
                                  arguments->userdefinedleastlength,
                                  querymatch_output,
                                  NULL,
                                  logger,
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
                      (unsigned long) (100 * arguments->userdefinedleastlength),
                      logger,
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
                               querymatch_output,
                               NULL,
                               logger,
                               err) != 0)
      {
        haserr = true;
      }
    }
  }
  querymatch_delete(querymatchspaceptr);
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
