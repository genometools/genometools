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
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "fmindex.h"
#include "format64.h"
#include "seqdesc.h"
#include "gqueue-def.h"
#include "iterseq.h"

#include "fmi-map.pr"
#include "fmi-fwduni.pr"
#include "genericqueue.pr"

#define SHOWSEQUENCE   ((unsigned int) 1)
#define SHOWQUERYPOS   (SHOWSEQUENCE << 1)

typedef struct
{
  bool defined;
  unsigned long valueunsignedlong;
} Definedunsignedlong;

typedef struct
{
  unsigned int showmode;
  Definedunsignedlong minlength,
                      maxlength;
} Rangespecinfo;

typedef struct
{
  Definedunsignedlong minlength,
                      maxlength;
  unsigned int showmode;
  Str *fmindexname;
  StrArray *queryfilenames;
} Uniquesubcallinfo;

typedef struct
{
  const char *name;
  unsigned int bitmask;
} Optionargmodedesc;

typedef int (*Preprocessuniquelength)(uint64_t,
                                      const char *,
                                      void *,
                                      Env *env);
typedef int (*Processuniquelength)(const Alphabet *,
                                   const Uchar *,
                                   unsigned long,
                                   unsigned long,
                                   void *,
                                   Env *env);
typedef int (*Postprocessuniquelength)(const Alphabet *,
                                       uint64_t,
                                       const char *,
                                       const Uchar *,
                                       unsigned long,
                                       void *,
                                       Env *env);

typedef struct
{
  Fmindex *fmindex;
  Preprocessuniquelength preprocessuniquelength;
  Processuniquelength processuniquelength;
  Postprocessuniquelength postprocessuniquelength;
  void *processinfo;
} Substringinfo;

static int optionaddbitmask(Optionargmodedesc *modedesc,
                            size_t numberofentries,
                            unsigned int *mode,
                            const char *optname,
                            const char *optionargument,
                            Env *env)
{
  size_t modecount;

  env_error_check(env);
  for (modecount=0; modecount < numberofentries; modecount++)
  {
    if (strcmp(optionargument,modedesc[modecount].name) == 0)
    {
      if (*mode & modedesc[modecount].bitmask)
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
  uniquesubcallinfo->queryfilenames = strarray_new(env);
  flagsoutputoption = strarray_new(env);

  op = option_parser_new("[option ...] -fm fmindex -quer queryfile [...]",
                         "Compute length of minumum unique prefixes.", env);
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  optionmin = option_new_ulong_min("min",
                                   "only output length "
                                   "if >= given minimum length",
                                   &uniquesubcallinfo->minlength.
                                          valueunsignedlong,
                                   0,(unsigned long) 1,env);
  option_parser_add_option(op, optionmin, env);

  optionmax = option_new_ulong_min("max",
                                   "only output length "
                                   "if <= given maximum length",
                                   &uniquesubcallinfo->maxlength.
                                          valueunsignedlong,
                                   0,(unsigned long) 1,env);
  option_parser_add_option(op, optionmax, env);

  optionoutput = option_new_stringarray("output",
                          "set output flags (sequence, querypos)",
                          flagsoutputoption,env);
  option_parser_add_option(op, optionoutput, env);

  optionfmindex = option_new_string("fmi", "specify fmindex",
                                    uniquesubcallinfo->fmindexname,NULL,env);
  option_is_mandatory(optionfmindex);
  option_parser_add_option(op, optionfmindex, env);

  optionquery = option_new_filenamearray("query", "specify queryfiles",
                                         uniquesubcallinfo->queryfilenames,env);
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
    if (!option_is_set(optionmin) && !option_is_set(optionmax))
    {
      env_error_set(env,"one of the options -min or -max must be set");
      oprval = OPTIONPARSER_ERROR;
    }
    if (oprval != OPTIONPARSER_ERROR)
    {
      if (uniquesubcallinfo->minlength.defined &&
         uniquesubcallinfo->maxlength.defined)
      {
        if (uniquesubcallinfo->maxlength.valueunsignedlong <
           uniquesubcallinfo->minlength.valueunsignedlong)
        {
          env_error_set(env,"minvalue must be smaller or equal than maxvalue");
          oprval = OPTIONPARSER_ERROR;
        }
      }
    }
    if (oprval != OPTIONPARSER_ERROR && option_is_set(optionoutput))
    {
      if (strarray_size(flagsoutputoption) == 0)
      {
        env_error_set(env,"missing arguments to option -output");
        oprval = OPTIONPARSER_ERROR;
      } else
      {
        unsigned long i;

        for (i=0; i<strarray_size(flagsoutputoption); i++)
        {
          if (optionaddbitmask(uniquesubmodedesctable,
                              sizeof (uniquesubmodedesctable)/
                              sizeof (uniquesubmodedesctable[0]),
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
  if (oprval == OPTIONPARSER_OK && parsed_args != argc)
  {
    env_error_set(env,"superfluous program parameters");
    oprval = OPTIONPARSER_ERROR;
  }
  return oprval;
}

static int uniqueposinsinglesequence(void *info,
                                     uint64_t unitnum,
                                     const Uchar *start,
                                     unsigned long seqlen,
                                     const char *desc,
                                     Env *env)
{
  Substringinfo *substringinfo = (Substringinfo *) info;
  const Uchar *query;
  unsigned long uniquelength, remaining;

  env_error_check(env);
  if (substringinfo->preprocessuniquelength != NULL &&
     substringinfo->preprocessuniquelength(unitnum,
                                           desc,
                                           substringinfo->processinfo,
                                           env) != 0)
  {
    return -1;
  }
  for (query = start, remaining = seqlen; remaining > 0; query++, remaining--)
  {
    uniquelength = skfmuniqueforward (substringinfo->fmindex,
                                      query,start+seqlen);
    if (uniquelength > 0)
    {
      if (substringinfo->processuniquelength(substringinfo->fmindex->alphabet,
                                            start,
                                            uniquelength,
                                            (unsigned long) (query-start),
                                            substringinfo->processinfo,
                                            env) != 0)
      {
        return -2;
      }
    }
  }
  if (substringinfo->postprocessuniquelength != NULL &&
     substringinfo->postprocessuniquelength(substringinfo->fmindex->alphabet,
                                            unitnum,
                                            desc,
                                            start,
                                            seqlen,
                                            substringinfo->processinfo,
                                            env) != 0)
  {
    return -3;
  }
  return 0;
}

static int showunitnum(uint64_t unitnum,
                       const char *desc,
                       /*@unused@*/ void *info,
                       /*@unused@*/ Env *env)
{
  printf("unit " Formatuint64_t, PRINTuint64_tcast(unitnum));
  if (desc != NULL && desc[0] != '\0')
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
                               void *info,
                                /*@unused@*/ Env *env)
{
  Rangespecinfo *rangespecinfo = (Rangespecinfo *) info;

  if ((!rangespecinfo->minlength.defined ||
      uniquelength >= rangespecinfo->minlength.valueunsignedlong) &&
     (!rangespecinfo->maxlength.defined ||
      uniquelength <= rangespecinfo->maxlength.valueunsignedlong))
  {
    if (rangespecinfo->showmode & SHOWQUERYPOS)
    {
      printf("%lu ",querystart);
    }
    printf("%lu",uniquelength);
    if (rangespecinfo->showmode & SHOWSEQUENCE)
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
                             const StrArray *queryfilenames,
                             Definedunsignedlong minlength,
                             Definedunsignedlong maxlength,
                             unsigned int showmode,
                             Env *env)
{
  Substringinfo substringinfo;
  Rangespecinfo rangespecinfo;
  ArrayUchar sequencebuffer;
  Sequencedescription sequencedescription;
  bool haserr = false;

  env_error_check(env);
  substringinfo.fmindex = fmindex;
  rangespecinfo.minlength = minlength;
  rangespecinfo.maxlength = maxlength;
  rangespecinfo.showmode = showmode;
  substringinfo.preprocessuniquelength = showunitnum;
  substringinfo.processuniquelength = showifinlengthrange;
  substringinfo.postprocessuniquelength = NULL;
  substringinfo.processinfo = &rangespecinfo;
  INITARRAY(&sequencebuffer,Uchar);
  INITARRAY(&sequencedescription.headerbuffer,char);
  sequencedescription.descptr = emptyqueuegeneric(env);
  if (overallquerysequences(uniqueposinsinglesequence,
                           &substringinfo,
                           &sequencebuffer,
                           queryfilenames,
                           &sequencedescription,
                           getsymbolmapAlphabet(fmindex->alphabet),
                           env) != 0)
  {
    haserr = true;
  }
  wrapqueuegeneric(true,&sequencedescription.descptr,env);
  FREEARRAY(&sequencebuffer,Uchar);
  FREEARRAY(&sequencedescription.headerbuffer,char);
  return haserr ? -1 : 0;
}

int findminuniquesubstrings(int argc,const char **argv,Env *env)
{
  Uniquesubcallinfo uniquesubcallinfo;
  Fmindex fmindex;
  Verboseinfo *verboseinfo;
  int had_err = 0;

  env_error_check(env);

  switch (parseuniquesub(&uniquesubcallinfo, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      str_delete(uniquesubcallinfo.fmindexname,env);
      strarray_delete(uniquesubcallinfo.queryfilenames,env);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      str_delete(uniquesubcallinfo.fmindexname,env);
      strarray_delete(uniquesubcallinfo.queryfilenames,env);
      return 0;
  }
  verboseinfo = newverboseinfo(false,env);
  if (mapfmindex (&fmindex, uniquesubcallinfo.fmindexname,
                  verboseinfo,env) != 0)
  {
    had_err = -1;
  } else
  {
    if (findsubquerymatch(&fmindex,
                          uniquesubcallinfo.queryfilenames,
                          uniquesubcallinfo.minlength,
                          uniquesubcallinfo.maxlength,
                          uniquesubcallinfo.showmode,
                          env) != 0)
    {
      had_err = -1;
    }
    freefmindex(&fmindex,env);
  }
  freeverboseinfo(&verboseinfo,env);
  str_delete(uniquesubcallinfo.fmindexname,env);
  strarray_delete(uniquesubcallinfo.queryfilenames,env);
  return had_err;
}
