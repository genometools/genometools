/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <inttypes.h>
#include "libgtcore/option.h"
#include "fmindex.h"

#include "alphabet.pr"
#include "fmi-map.pr"
#include "fmi-fwduni.pr"

#define SHOWSEQUENCE   ((uint32_t) 1)
#define SHOWQUERYPOS   (((uint32_t) 1) << 1)

typedef struct
{
  bool defined;
  unsigned long valueunsignedlong;
} Definedunsignedlong;

typedef struct
{
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
  const char *name;
  uint32_t bitmask;
} Optionargmodedesc;

typedef int (*Preprocessuniquelength)(unsigned long,void *,const char *);
typedef int (*Processuniquelength)(const Alphabet *,
                                   const Uchar *,
                                   unsigned long,
                                   unsigned long,
                                   void *);
typedef int (*Postprocessuniquelength)(const Alphabet *,
                                       const char *,
                                       unsigned long,
                                       const Uchar *,
                                       unsigned long,
                                       void *);

typedef struct
{
  Fmindex *fmindex;
  Preprocessuniquelength preprocessuniquelength;
  Processuniquelength processuniquelength;
  Postprocessuniquelength postprocessuniquelength;
  void *processinfo;
} Substringinfo;

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
                          "set output flags (sequence, querypos)",
                          flagsoutputoption,env);
  option_parser_add_option(op, optionoutput, env);

  optionfmindex = option_new_string("fmi",
                                    "specify fmindex (mandatory)",
                                    uniquesubcallinfo->fmindexname,NULL,env);
  option_is_mandatory(optionfmindex);
  option_parser_add_option(op, optionfmindex, env);

  optionquery = option_new_string("query",
                                  "specify queryfile (mandatory)",
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
    if(!option_is_set(optionmin) && !option_is_set(optionmax))
    {
      env_error_set(env,"one of the options -min or -max must be set");
      oprval = OPTIONPARSER_ERROR;
    }
    if(oprval != OPTIONPARSER_ERROR)
    {
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

int uniqueposinsinglesequence(void *info, 
                                     unsigned long unitnum,
                                     const Uchar *start,
                                     const char *desc,
                                     unsigned long len)
{
  Substringinfo *substringinfo = (Substringinfo *) info;
  const Uchar *query;
  unsigned long uniquelength, remaining;

  if(substringinfo->preprocessuniquelength != NULL &&
     substringinfo->preprocessuniquelength(unitnum,
                                           substringinfo->processinfo,
                                           desc) != 0)
  {
    return -1;
  }
  for(query = start, remaining = len; remaining > 0; query++, remaining--)
  {
    uniquelength = skfmuniqueforward (substringinfo->fmindex,query,start+len);
    if(uniquelength > 0)
    {
      if(substringinfo->processuniquelength(substringinfo->fmindex->alphabet,
                                            start,
                                            uniquelength,
                                            (unsigned long) (query-start),
                                            substringinfo->processinfo) != 0)
      {
        return -2;
      }
    }
  }
  if(substringinfo->postprocessuniquelength != NULL &&
     substringinfo->postprocessuniquelength(substringinfo->fmindex->alphabet,
                                            desc,
                                            unitnum,
                                            start,
                                            len,
                                            substringinfo->processinfo) != 0)
  {
    return -3;
  }
  return 0;
}

static int showunitnum(unsigned long unitnum,
                       /*@unused@*/ void *info,
                       const char *desc)
{
  printf("unitnum=%lu",unitnum);
  if(desc != NULL)
  {
    printf(" (%s)",desc);
  }
  printf("\n");
  return 0;
}

static int showifinlengthrange(const Alphabet *alphabet,
                               const Uchar *start,
                               unsigned long uniquelength,
                               unsigned long querystart,
                               void *info)
{
  Rangespecinfo *rangespecinfo = (Rangespecinfo *) info;

  if((!rangespecinfo->minlength.defined ||
      uniquelength >= rangespecinfo->minlength.valueunsignedlong) &&
     (!rangespecinfo->maxlength.defined ||
      uniquelength <= rangespecinfo->maxlength.valueunsignedlong))
  {
    if(rangespecinfo->showmode & SHOWQUERYPOS)
    {
      printf("%lu ",querystart);
    }
    printf("%lu",uniquelength);
    if(rangespecinfo->showmode & SHOWSEQUENCE)
    {
      (void) putchar(' ');
      showsymbolstring(alphabet,
                       start + querystart,
                       uniquelength);
    }
    (void) putchar('\n');
  }
  return 0;
}

static int findsubquerymatch(Fmindex *fmindex,
                             /*@unused@*/ const Str *queryfilename,
                             Definedunsignedlong minlength,
                             Definedunsignedlong maxlength,
                             uint32_t showmode,
                             /*@unused@*/ Env *env)
{
  Substringinfo substringinfo;
  Rangespecinfo rangespecinfo;

  substringinfo.fmindex = fmindex;
  rangespecinfo.minlength = minlength;
  rangespecinfo.maxlength = maxlength;
  rangespecinfo.showmode = showmode;
  substringinfo.preprocessuniquelength = showunitnum;
  substringinfo.processuniquelength = showifinlengthrange;
  substringinfo.postprocessuniquelength = NULL;
  substringinfo.processinfo = &rangespecinfo;
  /*
  if(overallsequences(False,
                      querymultiseq,
                      &substringinfo,
                      uniqueposinsinglesequence) != 0)
  {
    return (Sint) -2;
  }
  */
  return 0;
}

/*
int overallquerysequence(
        const StrArray *filenametab,
        const Uchar *symbolmap,
        Env *env)
{
  Fastabufferstate fbs;
  Uchar charcode;
  int retval;
  unsigned long idx;

  *numofsequences = 0;
  specialcharinfo->specialcharacters = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofspecialsuffix = 0;

  initformatbufferstate(&fbs,
                        filenametab,
                        symbolmap,
                        plainformat,
                        filelengthtab,
                        env);
  specialrangelengths = initdistribution(env);
  for (pos = 0; ; pos++)
  {
    retval = readnextUchar(&charcode,&fbs,env);
    if (retval < 0)
    {
      return -1;
    }
    if (retval == 0)
    {
      if (lastspeciallength > 0)
      {
        idx = CALLCASTFUNC(Seqpos,unsigned_long,lastspeciallength);
        adddistribution(specialrangelengths,idx,env);
      }
      break;
    }
    if (ISSPECIAL(charcode))
    {
      if (specialprefix)
      {
        specialcharinfo->lengthofspecialprefix++;
      }
      specialcharinfo->specialcharacters++;
      if (lastspeciallength == 0)
      {
        lastspeciallength = (Seqpos) 1;
      } else
      {
        lastspeciallength++;
      }
      if (charcode == (Uchar) SEPARATOR)
      {
        (*numofsequences)++;
      }
    } else
    {
      if (specialprefix)
      {
        specialprefix = false;
      }
      if (lastspeciallength > 0)
      {
        idx = CALLCASTFUNC(Seqpos,unsigned_long,lastspeciallength);
        adddistribution(specialrangelengths,idx,env);
        lastspeciallength = 0;
      }
    }
  }
  specialcharinfo->specialranges = 0;
  (void) foreachdistributionvalue(specialrangelengths,updatesumranges,
                                  &specialcharinfo->specialranges,env);
  freedistribution(&specialrangelengths,env);
  specialcharinfo->lengthofspecialsuffix = lastspeciallength;
  (*numofsequences)++;
  *totallength = pos;
  return 0;
}
*/

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
      if(findsubquerymatch(&fmindex,
                           uniquesubcallinfo.queryfilename,
                           uniquesubcallinfo.minlength,
                           uniquesubcallinfo.maxlength,
                           uniquesubcallinfo.showmode,
                           env) != 0)
      {
        haserr = true;
      }
      freefmindex(&fmindex,env);
    }
  }
  str_delete(uniquesubcallinfo.fmindexname,env);
  str_delete(uniquesubcallinfo.queryfilename,env);
  return haserr ? -1 : 0;
}
