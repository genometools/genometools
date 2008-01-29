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

#include "libgtcore/getbasename.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "divmodmul.h"
#include "fmindex.h"
#include "verbose-def.h"
#include "spacedef.h"
#include "stamp.h"

#include "fmi-save.pr"
#include "fmi-keyval.pr"
#include "fmi-sufbwtstream.pr"

typedef struct
{
  bool noindexpos;
  StrArray *indexnametab;
  Str *leveldesc,
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

static OPrval parsemkfmindex(Mkfmcallinfo *mkfmcallinfo,
                             int argc,
                             const char **argv,
                             Error *err)
{
  OptionParser *op;
  Option *option, *optionfmout;
  OPrval oprval;
  int parsed_args;

  error_check(err);
  mkfmcallinfo->indexnametab = strarray_new();
  mkfmcallinfo->outfmindex = str_new();
  mkfmcallinfo->leveldesc = str_new();
  op = option_parser_new("[option ...] -ii indexfile [...]",
                         "Compute FMindex.");
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  optionfmout = option_new_string("fmout",
                             "specify name of FMindex to be generated\n"
                             "(mandatory if more than one input index "
                             "is specified)",
                             mkfmcallinfo->outfmindex, NULL);
  option_parser_add_option(op, optionfmout);

  option = option_new_filenamearray("ii", "specify indices to be used",
                                    mkfmcallinfo->indexnametab);
  option_is_mandatory(option);
  option_parser_add_option(op, option);

  option = option_new_string("size",
                             "specify size (tiny, small, medium, big)",
                             mkfmcallinfo->leveldesc, "medium");
  option_parser_add_option(op, option);

  option = option_new_bool("noindexpos",
                           "store no index positions (hence the positions of\n"
                           "matches in the index cannot be retrieved)",
                           &mkfmcallinfo->noindexpos,false);
  option_parser_add_option(op, option);

  oprval = option_parser_parse(op, &parsed_args, argc, argv, versionfunc,err);
  if (oprval == OPTIONPARSER_OK)
  {
    if (!option_is_set(optionfmout))
    {
      if (strarray_size(mkfmcallinfo->indexnametab) > 1UL)
      {
        error_set(err,"if more than one index is given, then "
                          "option -fmout is mandatory");
        oprval = OPTIONPARSER_ERROR;
      } else
      {
        char *basenameptr;

        basenameptr = getbasename(strarray_get(mkfmcallinfo->indexnametab,0));
        str_set(mkfmcallinfo->outfmindex,basenameptr);
        ma_free(basenameptr);
      }
    }
  }
  option_parser_delete(op);
  if (oprval == OPTIONPARSER_OK && parsed_args != argc)
  {
    error_set(err,"superfluous program parameters");
    oprval = OPTIONPARSER_ERROR;
  }
  return oprval;
}

static void freemkfmcallinfo(Mkfmcallinfo *mkfmcallinfo)
{
  strarray_delete(mkfmcallinfo->indexnametab);
  str_delete(mkfmcallinfo->outfmindex);
  str_delete(mkfmcallinfo->leveldesc);
}

static int levedescl2levelnum(const char *name,
                              unsigned int *log2bsize,
                              unsigned int *log2markdist)
{
  size_t i;

  assert(name != NULL);
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
  FREEARRAY (&fm->specpos, PairBwtidx);
  FREESPACE (fm->bfreq);
  FREESPACE (fm->superbfreq);
  FREESPACE (fm->tfreq);
  FREESPACE (fm->markpostable);
  if (fm->suffixlength > 0)
  {
    FREESPACE(fm->boundarray);
  }
}

static int mkfmindexoptions(Mkfmcallinfo *mkfmcallinfo,
                            int argc,const char **argv,Error *err)
{
  OPrval rval;
  int retval = 0;

  error_check(err);
  rval = parsemkfmindex(mkfmcallinfo,argc,argv,err);
  if (rval == OPTIONPARSER_ERROR)
  {
    retval = -1;
  } else
  {
    if (rval == OPTIONPARSER_REQUESTS_EXIT)
    {
      retval = 2;
    }
  }
  return retval;
}

static int runmkfmindex(Mkfmcallinfo *mkfmcallinfo,Verboseinfo *verboseinfo,
                        Error *err)
{
  Fmindex fm;
  unsigned int log2bsize,
               log2markdist;
  bool haserr = false;

  error_check(err);
  INITARRAY(&fm.specpos, PairBwtidx);
  fm.bfreq = NULL;
  fm.superbfreq = NULL;
  fm.tfreq = NULL;
  fm.markpostable = NULL;
  fm.boundarray = NULL;
  fm.suffixlength = 0;

  if (levedescl2levelnum(str_get(mkfmcallinfo->leveldesc),
                        &log2bsize,
                        &log2markdist) != 0)
  {
    error_set(err,"undefined level \"%s\"",
                      str_get(mkfmcallinfo->leveldesc));
    haserr = true;
  }
  if (!haserr && sufbwt2fmindex(&fm,
                                log2bsize,
                                log2markdist,
                                mkfmcallinfo->outfmindex,
                                mkfmcallinfo->indexnametab,
                                mkfmcallinfo->noindexpos ? false : true,
                                verboseinfo,
                                err) != 0)
  {
    haserr = true;
  }
  if (!haserr && saveFmindex(mkfmcallinfo->outfmindex,
                            &fm,
                            mkfmcallinfo->noindexpos ? false : true,
                            err) < 0)
  {
    haserr = true;
  }
  freeconstructedfmindex(&fm);
  return haserr ? -1 : 0;
}

int parseargsandcallmkfmindex(int argc,const char **argv,Error *err)
{
  Mkfmcallinfo mkfmcallinfo;
  int retval;
  bool haserr = false;

  retval = mkfmindexoptions(&mkfmcallinfo,argc,argv,err);
  if (retval == 0)
  {
    Verboseinfo *verboseinfo = newverboseinfo(false);
    if (runmkfmindex(&mkfmcallinfo,verboseinfo,err) < 0)
    {
      haserr = true;
    }
    freeverboseinfo(&verboseinfo);
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
