/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <inttypes.h>
#include "libgtcore/option.h"
#include "fmindex.h"
#include "stamp.h"

#include "fmi-map.pr"

#define SHOWSEQUENCE   ((uint32_t) 1)
#define SHOWSUBJECTPOS (((uint32_t) 1) << 1)
#define SHOWQUERYPOS   (((uint32_t) 1) << 2)

typedef struct
{
  void *subjectposinfo;
  uint32_t showmode;
  Definedunsignedlong minlength,
                      maxlength;
} Rangespecinfo;

typedef struct
{
  Definedunsignedlong minlength,
                      maxlength;
  uint32_t showmode;
  Str *fmindexname,
      *queryfilename;
} Uniquesubcallinfo;

typedef struct
{
  Seqpos subjectpos,
         querysubstringstart;
} Uniqueposinfo;

 DECLAREARRAYSTRUCT(Uniqueposinfo);

typedef struct
{
  char *name;
  uint32_t bitmask;
} Optionargmodedesc;

static void versionfunc(const char *progname)
{
  printf("%s version 0.1\n",progname);
}

static int optionaddbitmask(Optionargmodedesc *modedesc,
                            size_t numberofentries,
                            uint32_t *mode,
                            const char *optname,
                            const char *optionargument,
                            Env *env)
{
  size_t modecount;
  
  for(modecount=0; modecount < numberofentries; modecount++)
  {
    if(strcmp(optionargument,modedesc[modecount].name) == 0)
    {
      if(*mode & modedesc[modecount].bitmask)
      {
        env_error_set(env,"argument \"%s\" to option %s already specified",
                      modedesc[modecount].name,optname);
        return -1;
      }
      *mode |= modedesc[modecount].bitmask;
      return 0;
    }
  }
  env_error_set(env,"illegal argument \"%s\" to option %s",
                optionargument,optname);
  return -2;
}

static OPrval parseuniquesub(Uniquesubcallinfo *uniquesubcallinfo,
                             int argc, 
                             const char **argv, 
                             Env *env)
{
  OptionParser *op;
  Option *optionmin, *optionmax, *optionoutput, *optionfmindex, *optionquery;
  OPrval oprval;
  StrArray *flagsoutputoption;
  int parsed_args;
  Optionargmodedesc uniquesubmodedesctable[]
    = {
      {"sequence",SHOWSEQUENCE},
      {"subjectpos",SHOWSUBJECTPOS},
      {"querypos",SHOWQUERYPOS}
  };

  env_error_check(env);
  uniquesubcallinfo->minlength.defined = false;
  uniquesubcallinfo->maxlength.defined = false;
  uniquesubcallinfo->showmode = 0;
  uniquesubcallinfo->fmindexname = str_new(env);
  uniquesubcallinfo->queryfilename = str_new(env);
  flagsoutputoption = strarray_new(env);

  op = option_parser_new("options",
                         "Compute minumum unique prefixlengths.", env);
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  optionmin = option_new_ulong_min("min",
                                   "only output shortest unique prefixlength "
                                   "if >= given minimum length",
                                   &uniquesubcallinfo->minlength.
                                          valueunsignedlong,
                                   0,(unsigned long) 1,env);
  option_parser_add_option(op, optionmin, env);

  optionmax = option_new_ulong_min("max",
                                   "only output shortest unique prefixlength "
                                   "if <= given maximum length",
                                   &uniquesubcallinfo->maxlength.
                                          valueunsignedlong,
                                   0,(unsigned long) 1,env);
  option_parser_add_option(op, optionmax, env);

  optionoutput = option_new_filenamearray("output",
                          "set output flags (sequence, subjectpos, querypos)",
                          flagsoutputoption,env);
  option_parser_add_option(op, optionoutput, env);

  optionfmindex = option_new_string("fmi",
                                    "specify fmindex (mandatory)",
                                    uniquesubcallinfo->fmindexname,NULL,env);
  option_is_mandatory(optionfmindex);
  option_parser_add_option(op, optionfmindex, env);

  optionquery = option_new_string("query",
                                  "specify queryfiles (mandatory)",
                                  uniquesubcallinfo->queryfilename,
                                  NULL,env);
  option_is_mandatory(optionquery);
  option_parser_add_option(op, optionquery, env);

  oprval = option_parser_parse(op, &parsed_args, argc, argv, versionfunc, env);
  if (oprval == OPTIONPARSER_OK)
  {
    if (option_is_set(optionmin))
    {
       uniquesubcallinfo->minlength.defined = true;
    }
    if (option_is_set(optionmax))
    {
       uniquesubcallinfo->maxlength.defined = true;
    }
    if(uniquesubcallinfo->minlength.defined &&
       uniquesubcallinfo->maxlength.defined)
    {
      if(uniquesubcallinfo->maxlength.valueunsignedlong < 
         uniquesubcallinfo->minlength.valueunsignedlong)
      {
        env_error_set(env,"minvalue must be smaller or equal than maxvalue");
        oprval = OPTIONPARSER_ERROR;
      }
    }
    if (oprval != OPTIONPARSER_ERROR && option_is_set(optionoutput))
    {
      if(strarray_size(flagsoutputoption) == 0)
      {
        env_error_set(env,"missing arguments to option -output");
        oprval = OPTIONPARSER_ERROR;
      } else
      {
        unsigned long i;

        for(i=0; i<strarray_size(flagsoutputoption); i++)
        {
          if(optionaddbitmask(uniquesubmodedesctable,
                              sizeof(uniquesubmodedesctable)/
                              sizeof(uniquesubmodedesctable[0]),
                              &uniquesubcallinfo->showmode,
                              "-output",
                              strarray_get(flagsoutputoption,i),
                              env) != 0)
          {
            oprval = OPTIONPARSER_ERROR;
            break;
          }
        }
      }
    }
  }
  strarray_delete(flagsoutputoption,env);
  option_parser_delete(op, env);
  if(oprval == OPTIONPARSER_OK && parsed_args != argc)
  {
    env_error_set(env,"superfluous program parameters");
    oprval = OPTIONPARSER_ERROR;
  }
  return oprval;
}

int findminuniquesubstrings(int argc,const char **argv,Env *env)
{
  Uniquesubcallinfo uniquesubcallinfo;
  OPrval oprval;
  int retval = 0;
  bool haserr = false;

  oprval = parseuniquesub(&uniquesubcallinfo,
                          argc, 
                          argv, 
                          env);
  if (oprval == OPTIONPARSER_ERROR)
  {
    retval = -1;
  } else
  {
    if (oprval == OPTIONPARSER_REQUESTS_EXIT)
    {
      retval = 2;
    }
  }
  if (retval == 0)
  {
    Fmindex fmindex;

    if (mapfmindex (&fmindex, uniquesubcallinfo.fmindexname,env) != 0)
    {
      haserr = true;
    } else
    {
      freefmindex(&fmindex,env);
    }
  }
  str_delete(uniquesubcallinfo.fmindexname,env);
  str_delete(uniquesubcallinfo.queryfilename,env);
  return haserr ? -1 : 0;
}
