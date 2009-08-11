/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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
#include "core/error.h"
#include "core/str.h"
#include "core/option.h"
#include "core/unused_api.h"
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

static OPrval parse_options(Maxpairsoptions *maxpairsoptions,
                            int *parsed_args,
                            int argc,
                            const char **argv,
                            GtError *err)
{
  GtOptionParser *op;
  GtOption *option, *queryoption, *scanoption, *sampleoption;
  OPrval oprval;

  gt_error_check(err);
  op = gt_option_parser_new("[options] -ii indexname",
                         "Perform Substring matches with or without query.");
  gt_option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");

  option = gt_option_new_uint_min("l","Specify minimum length",
                               &maxpairsoptions->userdefinedleastlength,
                               (unsigned int) 20,
                               (unsigned int) 1);
  gt_option_parser_add_option(op, option);

  sampleoption = gt_option_new_ulong_min("samples","Specify number of samples",
                                 &maxpairsoptions->samples,
                                 (unsigned long) 0,
                                 (unsigned long) 1);
  gt_option_parser_add_option(op, sampleoption);

  scanoption = gt_option_new_bool("scan","scan index",
                               &maxpairsoptions->scanfile,
                               false);
  gt_option_parser_add_option(op, scanoption);

  option = gt_option_new_string("ii",
                             "Specify input index",
                             maxpairsoptions->indexname, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  queryoption = gt_option_new_filenamearray("q",
                             "Specify query files",
                             maxpairsoptions->queryfiles);
  gt_option_parser_add_option(op, queryoption);

  oprval = gt_option_parser_parse(op, parsed_args, argc, argv,
                               gt_versionfunc, err);
  if (gt_option_is_set(queryoption))
  {
    if (gt_option_is_set(sampleoption))
    {
      gt_error_set(err, "option -samples cannot be combined with option -q");
      oprval = OPTIONPARSER_ERROR;
    } else
    {
      if (gt_option_is_set(scanoption))
      {
        gt_error_set(err, "option -scan cannot be combined with option -q");
        oprval = OPTIONPARSER_ERROR;
      }
    }
  }
  gt_option_parser_delete(op);
  return oprval;
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

int gt_maxpairs(int argc, const char **argv, GtError *err)
{
  bool haserr = false;
  int parsed_args;
  Maxpairsoptions maxpairsoptions;
  OPrval oprval;

  gt_error_check(err);

  maxpairsoptions.indexname = gt_str_new();
  maxpairsoptions.queryfiles = gt_str_array_new();
  oprval = parse_options(&maxpairsoptions,&parsed_args, argc, argv,err);
  if (oprval == OPTIONPARSER_OK)
  {
    Verboseinfo *verboseinfo = newverboseinfo(false);
    gt_assert(parsed_args == argc);
    if (gt_str_array_size(maxpairsoptions.queryfiles) == 0)
    {
      if (maxpairsoptions.samples == 0)
      {
        if (callenummaxpairs(maxpairsoptions.indexname,
                             maxpairsoptions.userdefinedleastlength,
                             maxpairsoptions.scanfile,
                             simpleexactselfmatchoutput,
                             NULL,
                             verboseinfo,
                             err) != 0)
        {
          haserr = true;
        }
      } else
      {
        if (testmaxpairs(maxpairsoptions.indexname,
                         maxpairsoptions.samples,
                         maxpairsoptions.userdefinedleastlength,
                         (Seqpos) (100 *
                                   maxpairsoptions.userdefinedleastlength),
                         verboseinfo,
                         err) != 0)
        {
          haserr = true;
        }
      }
    } else
    {
      if (callenumquerymatches(maxpairsoptions.indexname,
                               maxpairsoptions.queryfiles,
                               false,
                               maxpairsoptions.userdefinedleastlength,
                               simpleexactquerymatchoutput,
                               NULL,
                               verboseinfo,
                               err) != 0)
      {
        haserr = true;
      }
    }
    freeverboseinfo(&verboseinfo);
  }
  gt_str_delete(maxpairsoptions.indexname);
  gt_str_array_delete(maxpairsoptions.queryfiles);

  if (oprval == OPTIONPARSER_REQUESTS_EXIT)
  {
    return 0;
  }
  if (oprval == OPTIONPARSER_ERROR)
  {
    return -1;
  }
  return haserr ? -1 : 0;
}
