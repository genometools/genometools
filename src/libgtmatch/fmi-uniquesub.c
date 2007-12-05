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
#include <string.h>
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/seqiterator.h"
#include "fmindex.h"
#include "spacedef.h"
#include "format64.h"

#include "fmi-map.pr"
#include "fmi-fwduni.pr"

#define SHOWSEQUENCE   1U
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
                                      Error *err);
typedef int (*Processuniquelength)(const Alphabet *,
                                   const Uchar *,
                                   unsigned long,
                                   unsigned long,
                                   void *,
                                   Error *err);
typedef int (*Postprocessuniquelength)(const Alphabet *,
                                       uint64_t,
                                       const char *,
                                       const Uchar *,
                                       unsigned long,
                                       void *,
                                       Error *err);

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
                            Error *err)
{
  size_t modecount;

  error_check(err);
  for (modecount=0; modecount < numberofentries; modecount++)
  {
    if (strcmp(optionargument,modedesc[modecount].name) == 0)
    {
      if (*mode & modedesc[modecount].bitmask)
      {
        error_set(err,"argument \"%s\" to option %s already specified",
                      modedesc[modecount].name,optname);
        return -1;
      }
      *mode |= modedesc[modecount].bitmask;
      return 0;
    }
  }
  error_set(err,"illegal argument \"%s\" to option %s",
                optionargument,optname);
  return -2;
}

static OPrval parseuniquesub(Uniquesubcallinfo *uniquesubcallinfo,
                             int argc,
                             const char **argv,
                             Error *err)
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

  error_check(err);
  uniquesubcallinfo->minlength.defined = false;
  uniquesubcallinfo->maxlength.defined = false;
  uniquesubcallinfo->showmode = 0;
  uniquesubcallinfo->fmindexname = str_new();
  uniquesubcallinfo->queryfilenames = strarray_new();
  flagsoutputoption = strarray_new();

  op = option_parser_new("[option ...] -fm fmindex -query queryfile [...]",
                         "Compute length of minumum unique prefixes.");
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  optionmin = option_new_ulong_min("min",
                                   "only output length "
                                   "if >= given minimum length",
                                   &uniquesubcallinfo->minlength.
                                          valueunsignedlong,
                                   0,(unsigned long) 1);
  option_parser_add_option(op, optionmin);

  optionmax = option_new_ulong_min("max",
                                   "only output length "
                                   "if <= given maximum length",
                                   &uniquesubcallinfo->maxlength.
                                          valueunsignedlong,
                                   0,(unsigned long) 1);
  option_parser_add_option(op, optionmax);

  optionoutput = option_new_stringarray("output",
                          "set output flags (sequence, querypos)",
                          flagsoutputoption);
  option_parser_add_option(op, optionoutput);

  optionfmindex = option_new_string("fmi", "specify fmindex",
                                    uniquesubcallinfo->fmindexname,NULL);
  option_is_mandatory(optionfmindex);
  option_parser_add_option(op, optionfmindex);

  optionquery = option_new_filenamearray("query", "specify queryfiles",
                                         uniquesubcallinfo->queryfilenames);
  option_is_mandatory(optionquery);
  option_parser_add_option(op, optionquery);

  oprval = option_parser_parse(op, &parsed_args, argc, argv, versionfunc,
                               err);
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
      error_set(err,"one of the options -min or -max must be set");
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
          error_set(err,"minvalue must be smaller or equal than maxvalue");
          oprval = OPTIONPARSER_ERROR;
        }
      }
    }
    if (oprval != OPTIONPARSER_ERROR && option_is_set(optionoutput))
    {
      if (strarray_size(flagsoutputoption) == 0)
      {
        error_set(err,"missing arguments to option -output");
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
                              err) != 0)
          {
            oprval = OPTIONPARSER_ERROR;
            break;
          }
        }
      }
    }
  }
  strarray_delete(flagsoutputoption);
  option_parser_delete(op);
  if (oprval == OPTIONPARSER_OK && parsed_args != argc)
  {
    error_set(err,"superfluous program parameters");
    oprval = OPTIONPARSER_ERROR;
  }
  return oprval;
}

static int uniqueposinsinglesequence(Substringinfo *substringinfo,
                                     uint64_t unitnum,
                                     const Uchar *query,
                                     unsigned long querylen,
                                     const char *desc,
                                     Error *err)
{
  const Uchar *qptr;
  unsigned long uniquelength, remaining;

  error_check(err);
  if (substringinfo->preprocessuniquelength != NULL &&
      substringinfo->preprocessuniquelength(unitnum,
                                            desc,
                                            substringinfo->processinfo,
                                            err) != 0)
  {
    return -1;
  }
  for (qptr = query, remaining = querylen; remaining > 0; qptr++, remaining--)
  {
    uniquelength = skfmuniqueforward (substringinfo->fmindex,
                                      qptr,query+querylen);
    if (uniquelength > 0)
    {
      if (substringinfo->processuniquelength(substringinfo->fmindex->alphabet,
                                             query,
                                             uniquelength,
                                             (unsigned long) (qptr-query),
                                             substringinfo->processinfo,
                                             err) != 0)
      {
        return -2;
      }
    }
  }
  if (substringinfo->postprocessuniquelength != NULL &&
      substringinfo->postprocessuniquelength(substringinfo->fmindex->alphabet,
                                             unitnum,
                                             desc,
                                             query,
                                             querylen,
                                             substringinfo->processinfo,
                                             err) != 0)
  {
    return -3;
  }
  return 0;
}

static int showunitnum(uint64_t unitnum,
                       const char *desc,
                       /*@unused@*/ void *info,
                       /*@unused@*/ Error *err)
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
                                /*@unused@*/ Error *err)
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
                             Error *err)
{
  Substringinfo substringinfo;
  Rangespecinfo rangespecinfo;
  bool haserr = false;
  SeqIterator *seqit;
  const Uchar *query;
  unsigned long querylen;
  char *desc = NULL;
  int retval;
  uint64_t unitnum;

  error_check(err);
  substringinfo.fmindex = fmindex;
  rangespecinfo.minlength = minlength;
  rangespecinfo.maxlength = maxlength;
  rangespecinfo.showmode = showmode;
  substringinfo.preprocessuniquelength = showunitnum;
  substringinfo.processuniquelength = showifinlengthrange;
  substringinfo.postprocessuniquelength = NULL;
  substringinfo.processinfo = &rangespecinfo;
  seqit = seqiterator_new(queryfilenames,
                          getsymbolmapAlphabet(fmindex->alphabet),
                          true);
  for (unitnum = 0; /* Nothing */; unitnum++)
  {
    retval = seqiterator_next(seqit,
                              &query,
                              &querylen,
                              &desc,
                              err);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      break;
    }
    if (uniqueposinsinglesequence(&substringinfo,
                                  unitnum,
                                  query,
                                  querylen,
                                  desc,
                                  err) != 0)
    {
      haserr = true;
      break;
    }
    FREESPACE(desc);
  }
  seqiterator_delete(seqit);
  return haserr ? -1 : 0;
}

int findminuniquesubstrings(int argc,const char **argv,Error *err)
{
  Uniquesubcallinfo uniquesubcallinfo;
  Fmindex fmindex;
  Verboseinfo *verboseinfo;
  int had_err = 0;

  error_check(err);

  switch (parseuniquesub(&uniquesubcallinfo, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      str_delete(uniquesubcallinfo.fmindexname);
      strarray_delete(uniquesubcallinfo.queryfilenames);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      str_delete(uniquesubcallinfo.fmindexname);
      strarray_delete(uniquesubcallinfo.queryfilenames);
      return 0;
  }
  verboseinfo = newverboseinfo(false);
  if (mapfmindex (&fmindex, uniquesubcallinfo.fmindexname,
                  verboseinfo,err) != 0)
  {
    had_err = -1;
  } else
  {
    if (findsubquerymatch(&fmindex,
                          uniquesubcallinfo.queryfilenames,
                          uniquesubcallinfo.minlength,
                          uniquesubcallinfo.maxlength,
                          uniquesubcallinfo.showmode,
                          err) != 0)
    {
      had_err = -1;
    }
    freefmindex(&fmindex);
  }
  freeverboseinfo(&verboseinfo);
  str_delete(uniquesubcallinfo.fmindexname);
  strarray_delete(uniquesubcallinfo.queryfilenames);
  return had_err;
}
