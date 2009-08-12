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
#include "match/esa-mmsearch-def.h"
#include "match/format64.h"
#include "match/verbose-def.h"

#include "match/esa-maxpairs.h"
#include "tools/gt_maxpairs.h"
#include "match/test-maxpairs.pr"

typedef struct
{
  unsigned int userdefinedleastlength;
  unsigned long samples;
  bool scanfile;
  GtStr *indexname;
  GtStrArray *queryfiles;
} Maxpairsoptions;

static int simpleexactselfmatchoutput(GT_UNUSED void *info,
                                      Seqpos len,
                                      Seqpos pos1,
                                      Seqpos pos2,
                                      GT_UNUSED GtError *err)
{
  if (pos1 > pos2)
  {
    Seqpos tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }
  printf(FormatSeqpos " " FormatSeqpos " " FormatSeqpos "\n",
            PRINTSeqposcast(len),
            PRINTSeqposcast(pos1),
            PRINTSeqposcast(pos2));
  return 0;
}

static int simpleexactquerymatchoutput(GT_UNUSED void *info,
                                       unsigned long len,
                                       Seqpos dbstart,
                                       uint64_t unitnum,
                                       unsigned long querystart,
                                       GT_UNUSED GtError *err)
{
  printf("%lu " FormatSeqpos " " Formatuint64_t " %lu\n",
           len,
           PRINTSeqposcast(dbstart),
           PRINTuint64_tcast(unitnum),
           querystart);
  return 0;
}

static int callenummaxpairs(const GtStr *indexname,
                            unsigned int userdefinedleastlength,
                            bool scanfile,
                            int(*processmaxpairs)(void *,Seqpos,Seqpos,
                                                  Seqpos,GtError *),
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
                                                SARR_ESQTAB,
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
  gt_free(arguments);
}

static GtOptionParser *gt_repfind_option_parser_new(void *tool_arguments)
{
  GtOptionParser *op;
  GtOption *option, *queryoption, *scanoption, *sampleoption;
  Maxpairsoptions *arguments = tool_arguments;

  op = gt_option_parser_new("[options] -ii indexname",
                            "Compute maximal repeats.");
  gt_option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");

  option = gt_option_new_uint_min("l","Specify minimum length",
                                  &arguments->userdefinedleastlength,
                                  20U,
                                  1U);
  gt_option_parser_add_option(op, option);

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

  gt_option_exclude(queryoption,sampleoption);
  gt_option_exclude(queryoption,scanoption);
  return op;
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
  gt_assert(parsed_args == argc);
  verboseinfo = newverboseinfo(false);
  if (gt_str_array_size(arguments->queryfiles) == 0)
  {
    if (arguments->samples == 0)
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
                             simpleexactquerymatchoutput,
                             NULL,
                             verboseinfo,
                             err) != 0)
    {
      haserr = true;
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
                     NULL,
                     gt_repfind_runner);
}
