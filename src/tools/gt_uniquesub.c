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

#include "libgtcore/error.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtmatch/fmindex.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/defined-types.h"
#include "libgtmatch/optionargmode.h"
#include "libgtmatch/uniquesub.h"

#include "tools/gt_uniquesub.h"
#include "libgtmatch/stamp.h"
#include "libgtmatch/fmi-fwduni.pr"
#include "libgtmatch/fmi-map.pr"
#include "libgtmatch/esa-map.pr"
#include "libgtmatch/esa-minunique.pr"

#define SHOWSEQUENCE   1U
#define SHOWQUERYPOS   (SHOWSEQUENCE << 1)

typedef enum
{
  Fmindextype,
  Saindextype
} Indextype;

typedef struct
{
  Definedunsignedlong minlength,
                      maxlength;
  unsigned int showmode;
  Str *indexname;
  StrArray *queryfilenames;
  Indextype indextype;
} Uniquesubcallinfo;

static OPrval parseuniquesub(Uniquesubcallinfo *uniquesubcallinfo,
                             int argc,
                             const char **argv,
                             Error *err)
{
  OptionParser *op;
  Option *optionmin, *optionmax, *optionoutput, *optionfmindex, 
         *optionesaindex, *optionquery;
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
  uniquesubcallinfo->indexname = str_new();
  uniquesubcallinfo->queryfilenames = strarray_new();
  flagsoutputoption = strarray_new();

  op = option_parser_new("[option ...] -query queryfile [...]",
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
                                    uniquesubcallinfo->indexname,NULL);
  option_parser_add_option(op, optionfmindex);
  
  optionesaindex = option_new_string("esa", "specify suffix array",
                                     uniquesubcallinfo->indexname,NULL);
  option_parser_add_option(op, optionesaindex);

  option_exclude(optionfmindex,optionesaindex);

  optionquery = option_new_filenamearray("query", "specify queryfiles",
                                         uniquesubcallinfo->queryfilenames);
  option_is_mandatory(optionquery);
  option_parser_add_option(op, optionquery);

  oprval = option_parser_parse(op, &parsed_args, argc, argv, versionfunc,
                               err);
  if (oprval == OPTIONPARSER_OK)
  {
    if (option_is_set(optionfmindex))
    {
      assert(!option_is_set(optionesaindex));
      uniquesubcallinfo->indextype = Fmindextype;
    } else
    {
      if (option_is_set(optionesaindex))
      {
        uniquesubcallinfo->indextype = Saindextype;
      } else
      {
        assert("unexpected case" == NULL);
      }
    }
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

static int findminuniquesubstrings(int argc,const char **argv,Error *err)
{
  Uniquesubcallinfo uniquesubcallinfo;
  Fmindex fmindex;
  Suffixarray suffixarray;
  Verboseinfo *verboseinfo;
  bool haserr = false;

  error_check(err);
  switch (parseuniquesub(&uniquesubcallinfo, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      str_delete(uniquesubcallinfo.indexname);
      strarray_delete(uniquesubcallinfo.queryfilenames);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      str_delete(uniquesubcallinfo.indexname);
      strarray_delete(uniquesubcallinfo.queryfilenames);
      return 0;
  }
  verboseinfo = newverboseinfo(false);
  if (uniquesubcallinfo.indextype == Fmindextype)
  {
    if (mapfmindex (&fmindex,uniquesubcallinfo.indexname,verboseinfo,err) != 0)
    {
      haserr = true;
    }
  } else
  {
    Seqpos totallength;

    if (mapsuffixarray(&suffixarray,
                       &totallength,
                       SARR_ESQTAB | SARR_SUFTAB,
                       uniquesubcallinfo.indexname,verboseinfo,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (uniquesubcallinfo.indextype == Fmindextype)
    {
      if (findsubqueryuniqueforward((const void *) &fmindex,
                                    skfmuniqueforward,
                                    fmindex.alphabet,
                                    uniquesubcallinfo.queryfilenames,
                                    uniquesubcallinfo.minlength,
                                    uniquesubcallinfo.maxlength,
                                    (uniquesubcallinfo.showmode & SHOWSEQUENCE)
                                       ? true : false,
                                    (uniquesubcallinfo.showmode & SHOWQUERYPOS)
                                       ? true : false,
                                    err) != 0)
      {
        haserr = true;
      }
    } else
    {
      if (findsubqueryuniqueforward((const void *) &suffixarray,
                                    suffixarrayuniqueforward,
                                    suffixarray.alpha,
                                    uniquesubcallinfo.queryfilenames,
                                    uniquesubcallinfo.minlength,
                                    uniquesubcallinfo.maxlength,
                                    (uniquesubcallinfo.showmode & SHOWSEQUENCE)
                                       ? true : false,
                                    (uniquesubcallinfo.showmode & SHOWQUERYPOS)
                                       ? true : false,
                                    err) != 0)
      {
        haserr = true;
      }
    }
  }
  if (uniquesubcallinfo.indextype == Fmindextype)
  {
    freefmindex(&fmindex);
  } else
  { 
    freesuffixarray(&suffixarray);
  }
  freeverboseinfo(&verboseinfo);
  str_delete(uniquesubcallinfo.indexname);
  strarray_delete(uniquesubcallinfo.queryfilenames);
  return haserr ? -1 : 0;
}

int gt_uniquesub(int argc, const char **argv, Error *err)
{
  error_check(err);
  return findminuniquesubstrings(argc, argv, err);
}
