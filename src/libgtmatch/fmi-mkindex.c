/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "divmodmul.h"
#include "fmindex.h"
#include "stamp.h"

#include "fmi-save.pr"
#include "fmi-keyval.pr"
#include "fmi-sufbwtstream.pr"
#include "getbasename.pr"

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
  mkfmcallinfo->outfmindex = str_new(env);
  mkfmcallinfo->leveldesc = str_new(env);
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
        env_ma_free(basenameptr,env);
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
  str_delete(mkfmcallinfo->outfmindex,env);
  str_delete(mkfmcallinfo->leveldesc,env);
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

static int runmkfmindex(Mkfmcallinfo *mkfmcallinfo,Env *env)
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
    if (runmkfmindex(&mkfmcallinfo,env) < 0)
    {
      haserr = true;
    }
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

