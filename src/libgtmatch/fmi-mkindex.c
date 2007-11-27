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
  {"tiny",  (unsigned int) 7, (unsigned int) 6},
  {"small", (unsigned int) 7, (unsigned int) 4},
  {"medium",(unsigned int) 5, (unsigned int) 3},
  {"big",   (unsigned int) 4, (unsigned int) 2}
};

static OPrval parsemkfmindex(Mkfmcallinfo *mkfmcallinfo,
                             int argc,
                             const char **argv,
                             Env *env)
{
  OptionParser *op;
  Option *option, *optionfmout;
  OPrval oprval;
  int parsed_args;

  env_error_check(env);
  mkfmcallinfo->indexnametab = strarray_new(env);
  mkfmcallinfo->outfmindex = str_new();
  mkfmcallinfo->leveldesc = str_new();
  op = option_parser_new("[option ...] -ii indexfile [...]",
                         "Compute FMindex.", env);
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  optionfmout = option_new_string("fmout",
                             "specify name of FMindex to be generated\n"
                             "(mandatory if more than one input index "
                             "is specified)",
                             mkfmcallinfo->outfmindex, NULL, env);
  option_parser_add_option(op, optionfmout, env);

  option = option_new_filenamearray("ii", "specify indices to be used",
                                    mkfmcallinfo->indexnametab,env);
  option_is_mandatory(option);
  option_parser_add_option(op, option, env);

  option = option_new_string("size",
                             "specify size (tiny, small, medium, big)",
                             mkfmcallinfo->leveldesc, "medium", env);
  option_parser_add_option(op, option, env);

  option = option_new_bool("noindexpos",
                           "store no index positions (hence the positions of\n"
                           "matches in the index cannot be retrieved)",
                           &mkfmcallinfo->noindexpos,false,env);
  option_parser_add_option(op, option, env);

  oprval = option_parser_parse(op, &parsed_args, argc, argv, versionfunc, env);
  if (oprval == OPTIONPARSER_OK)
  {
    if (!option_is_set(optionfmout))
    {
      if (strarray_size(mkfmcallinfo->indexnametab) > (unsigned long) 1)
      {
        env_error_set(env,"if more than one index is given, then "
                          "option -fmout is mandatory");
        oprval = OPTIONPARSER_ERROR;
      } else
      {
        char *basenameptr;

        basenameptr = getbasename(strarray_get(mkfmcallinfo->indexnametab,0),
                                  env);
        str_set(mkfmcallinfo->outfmindex,basenameptr,env);
        ma_free(basenameptr);
      }
    }
  }
  option_parser_delete(op, env);
  if (oprval == OPTIONPARSER_OK && parsed_args != argc)
  {
    env_error_set(env,"superfluous program parameters");
    oprval = OPTIONPARSER_ERROR;
  }
  return oprval;
}

static void freemkfmcallinfo(Mkfmcallinfo *mkfmcallinfo,Env *env)
{
  strarray_delete(mkfmcallinfo->indexnametab,env);
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

static void freeconstructedfmindex(Fmindex *fm,Env *env)
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
                            int argc,const char **argv,Env *env)
{
  OPrval rval;
  int retval = 0;

  env_error_check(env);
  rval = parsemkfmindex(mkfmcallinfo,argc,argv,env);
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
                        Env *env)
{
  Fmindex fm;
  unsigned int log2bsize,
               log2markdist;
  bool haserr = false;

  env_error_check(env);
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
    env_error_set(env,"undefined level \"%s\"",
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
                                env) != 0)
  {
    haserr = true;
  }
  if (!haserr && saveFmindex(mkfmcallinfo->outfmindex,
                            &fm,
                            mkfmcallinfo->noindexpos ? false : true,
                            env) < 0)
  {
    haserr = true;
  }
  freeconstructedfmindex(&fm,env);
  return haserr ? -1 : 0;
}

int parseargsandcallmkfmindex(int argc,const char **argv,Env *env)
{
  Mkfmcallinfo mkfmcallinfo;
  int retval;
  bool haserr = false;

  retval = mkfmindexoptions(&mkfmcallinfo,argc,argv,env);
  if (retval == 0)
  {
    Verboseinfo *verboseinfo = newverboseinfo(false,env);
    if (runmkfmindex(&mkfmcallinfo,verboseinfo,env) < 0)
    {
      haserr = true;
    }
    freeverboseinfo(&verboseinfo,env);
  } else
  {
    if (retval < 0)
    {
      haserr = true;
    }
  }
  freemkfmcallinfo(&mkfmcallinfo,env);
  return haserr ? -1 : 0;
}

