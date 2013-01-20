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

#include "core/basename_api.h"
#include "core/divmodmul.h"
#include "core/option_api.h"
#include "core/versionfunc.h"
#include "core/logger.h"
#include "core/ma_api.h"
#include "fmindex.h"
#include "fmi-save.pr"
#include "fmi-keyval.pr"
#include "fmi-sufbwtstream.pr"

typedef struct
{
  bool noindexpos;
  GtStrArray *indexnametab;
  GtStr *leveldesc,
      *outfmindex;
} Mkfmcallinfo;

typedef struct
{
  const char *name;
  unsigned int log2bsize,
               log2markdist;
} Indexleveldesc;

static Indexleveldesc indexlevel[] =
{
  {"tiny",  7U, 6U},
  {"small", 7U, 4U},
  {"medium",5U, 3U},
  {"big",   4U, 2U}
};

static GtOPrval parsemkfmindex(Mkfmcallinfo *mkfmcallinfo,
                               int argc,
                               const char **argv,
                               GtError *err)
{
  GtOptionParser *op;
  GtOption *option, *optionfmout;
  GtOPrval oprval;
  int parsed_args;

  gt_error_check(err);
  mkfmcallinfo->indexnametab = gt_str_array_new();
  mkfmcallinfo->outfmindex = gt_str_new();
  mkfmcallinfo->leveldesc = gt_str_new();
  op = gt_option_parser_new("[option ...] -ii indexfile [...]",
                         "Compute FM-index.");
  gt_option_parser_set_mail_address(op, "<kurtz@zbh.uni-hamburg.de>");
  optionfmout = gt_option_new_string("fmout",
                             "specify name of FM-index to be generated\n"
                             "(mandatory if more than one input index "
                             "is specified)",
                             mkfmcallinfo->outfmindex, NULL);
  gt_option_parser_add_option(op, optionfmout);

  option = gt_option_new_filename_array("ii", "specify indices to be used",
                                        mkfmcallinfo->indexnametab);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_string("size",
                             "specify size (tiny, small, medium, big)",
                             mkfmcallinfo->leveldesc, "medium");
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("noindexpos",
                           "store no index positions (hence the positions of\n"
                           "matches in the index cannot be retrieved)",
                           &mkfmcallinfo->noindexpos,false);
  gt_option_parser_add_option(op, option);

  oprval = gt_option_parser_parse(op, &parsed_args, argc, argv, gt_versionfunc,
                                  err);
  if (oprval == GT_OPTION_PARSER_OK)
  {
    if (!gt_option_is_set(optionfmout))
    {
      if (gt_str_array_size(mkfmcallinfo->indexnametab) > 1UL)
      {
        gt_error_set(err,"if more than one index is given, then "
                          "option -fmout is mandatory");
        oprval = GT_OPTION_PARSER_ERROR;
      } else
      {
        char *basenameptr;

        basenameptr = gt_basename(gt_str_array_get(mkfmcallinfo->indexnametab,
                                  0));
        gt_str_set(mkfmcallinfo->outfmindex,basenameptr);
        gt_free(basenameptr);
      }
    }
  }
  gt_option_parser_delete(op);
  if (oprval == GT_OPTION_PARSER_OK && parsed_args != argc)
  {
    gt_error_set(err,"superfluous program parameters");
    oprval = GT_OPTION_PARSER_ERROR;
  }
  return oprval;
}

static void freemkfmcallinfo(Mkfmcallinfo *mkfmcallinfo)
{
  gt_str_array_delete(mkfmcallinfo->indexnametab);
  gt_str_delete(mkfmcallinfo->outfmindex);
  gt_str_delete(mkfmcallinfo->leveldesc);
}

static int levedescl2levelnum(const char *name,
                              unsigned int *log2bsize,
                              unsigned int *log2markdist)
{
  size_t i;

  gt_assert(name != NULL);
  for (i=0; i<sizeof (indexlevel)/sizeof (indexlevel[0]); i++)
  {
    if (strcmp(name,indexlevel[i].name) == 0)
    {
      *log2bsize = indexlevel[i].log2bsize;
      *log2markdist = indexlevel[i].log2markdist;
      return 0;
    }
  }
  return -1;
}

static void freeconstructedfmindex(Fmindex *fm)
{
  GT_FREEARRAY (&fm->specpos, GtPairBwtidx);
  gt_free (fm->bfreq);
  gt_free (fm->superbfreq);
  gt_free (fm->tfreq);
  gt_free (fm->markpostable);
  if (fm->suffixlength > 0)
  {
    gt_free(fm->boundarray);
  }
}

static int mkfmindexoptions(Mkfmcallinfo *mkfmcallinfo,
                            int argc,const char **argv,GtError *err)
{
  GtOPrval rval;
  int retval = 0;

  gt_error_check(err);
  rval = parsemkfmindex(mkfmcallinfo,argc,argv,err);
  if (rval == GT_OPTION_PARSER_ERROR)
  {
    retval = -1;
  } else
  {
    if (rval == GT_OPTION_PARSER_REQUESTS_EXIT)
    {
      retval = 2;
    }
  }
  return retval;
}

static int runmkfmindex(Mkfmcallinfo *mkfmcallinfo,GtLogger *logger,
                        GtError *err)
{
  Fmindex fm;
  unsigned int log2bsize,
               log2markdist;
  bool haserr = false;
  GtSpecialcharinfo specialcharinfo;

  gt_error_check(err);
  GT_INITARRAY(&fm.specpos, GtPairBwtidx);
  fm.bfreq = NULL;
  fm.superbfreq = NULL;
  fm.tfreq = NULL;
  fm.markpostable = NULL;
  fm.boundarray = NULL;
  fm.suffixlength = 0;

  if (levedescl2levelnum(gt_str_get(mkfmcallinfo->leveldesc),
                        &log2bsize,
                        &log2markdist) != 0)
  {
    gt_error_set(err,"undefined level \"%s\"",
                      gt_str_get(mkfmcallinfo->leveldesc));
    haserr = true;
  }
  if (!haserr && gt_sufbwt2fmindex(&fm,
                                   &specialcharinfo,
                                   log2bsize,
                                   log2markdist,
                                   gt_str_get(mkfmcallinfo->outfmindex),
                                   mkfmcallinfo->indexnametab,
                                   mkfmcallinfo->noindexpos ? false : true,
                                   logger,
                                   err) != 0)
  {
    haserr = true;
  }
  if (!haserr && gt_saveFmindex(gt_str_get(mkfmcallinfo->outfmindex),
                                &fm,
                                &specialcharinfo,
                                mkfmcallinfo->noindexpos ? false : true,
                                err) < 0)
  {
    haserr = true;
  }
  freeconstructedfmindex(&fm);
  return haserr ? -1 : 0;
}

int gt_parseargsandcallmkfmindex(int argc,const char **argv,GtError *err)
{
  Mkfmcallinfo mkfmcallinfo;
  int retval;
  bool haserr = false;

  retval = mkfmindexoptions(&mkfmcallinfo,argc,argv,err);
  if (retval == 0)
  {
    GtLogger *logger = gt_logger_new(false, GT_LOGGER_DEFLT_PREFIX, stdout);
    if (runmkfmindex(&mkfmcallinfo,logger,err) < 0)
    {
      haserr = true;
    }
    gt_logger_delete(logger);
    logger = NULL;
  } else
  {
    if (retval < 0)
    {
      haserr = true;
    }
  }
  freemkfmcallinfo(&mkfmcallinfo);
  return haserr ? -1 : 0;
}
