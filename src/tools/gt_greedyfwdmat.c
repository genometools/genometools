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
#include "libgtcore/unused.h"
#include "libgtcore/versionfunc.h"
#include "libgtmatch/fmindex.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/defined-types.h"
#include "libgtmatch/optionargmode.h"
#include "libgtmatch/greedyfwdmat.h"

#include "libgtmatch/stamp.h"
#include "libgtmatch/fmi-fwduni.pr"
#include "libgtmatch/fmi-map.pr"
#include "libgtmatch/esa-map.pr"
#include "libgtmatch/esa-minunique.pr"
#include "libgtmatch/eis-voiditf.h"
#include "tools/gt_uniquesub.h"
#include "tools/gt_matchingstatistics.h"

#define SHOWSEQUENCE   1U
#define SHOWQUERYPOS   (SHOWSEQUENCE << 1)
#define SHOWSUBJECTPOS (SHOWSEQUENCE << 2)

typedef enum
{
  Fmindextype,
  Esaindextype,
  Packedindextype
} Indextype;

typedef struct
{
  Definedunsignedlong minlength,
                      maxlength;
  unsigned int showmode;
  bool verifywitnesspos;
  Str *indexname;
  StrArray *queryfilenames;
  Indextype indextype;
} Gfmsubcallinfo;

static OPrval parsegfmsub(bool doms,
                          Gfmsubcallinfo *gfmsubcallinfo,
                          int argc,
                          const char **argv,
                          Error *err)
{
  OptionParser *op;
  Option *optionmin, *optionmax, *optionoutput, *optionfmindex,
         *optionesaindex, *optionpckindex, *optionquery, *optionverify;
  OPrval oprval;
  StrArray *flagsoutputoption;
  int parsed_args;
  Optionargmodedesc msgfmsubmodedesctable[] =
  {
    {"sequence",SHOWSEQUENCE},
    {"querypos",SHOWQUERYPOS},
    {"subjectpos",SHOWSUBJECTPOS}
  };
  Optionargmodedesc gfmsubmodedesctable[] =
  {
    {"sequence",SHOWSEQUENCE},
    {"querypos",SHOWQUERYPOS}
  };

  error_check(err);
  gfmsubcallinfo->minlength.defined = false;
  gfmsubcallinfo->maxlength.defined = false;
  gfmsubcallinfo->showmode = 0;
  gfmsubcallinfo->indexname = str_new();
  gfmsubcallinfo->queryfilenames = strarray_new();
  flagsoutputoption = strarray_new();

  op = option_parser_new("[options ...] -query queryfile [...]",
                         doms
                         ? "Compute matching statistics."
                         : "Compute length of minumum unique prefixes."
                         );
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");

  optionfmindex = option_new_string("fmi", "specify fmindex",
                                    gfmsubcallinfo->indexname,NULL);
  option_parser_add_option(op, optionfmindex);

  optionesaindex = option_new_string("esa", "specify suffix array",
                                     gfmsubcallinfo->indexname,NULL);
  option_parser_add_option(op, optionesaindex);

  optionpckindex = option_new_string("pck", "specify packed index",
                                     gfmsubcallinfo->indexname,NULL);
  option_parser_add_option(op, optionpckindex);

  option_exclude(optionfmindex,optionesaindex);
  option_exclude(optionpckindex,optionesaindex);
  option_exclude(optionpckindex,optionfmindex);

  optionquery = option_new_filenamearray("query", "specify queryfiles",
                                         gfmsubcallinfo->queryfilenames);
  option_is_mandatory(optionquery);
  option_parser_add_option(op, optionquery);

  optionmin = option_new_ulong_min("min",
                                   "only output length "
                                   "if >= given minimum length",
                                   &gfmsubcallinfo->minlength.
                                          valueunsignedlong,
                                   0,(unsigned long) 1);
  option_parser_add_option(op, optionmin);

  optionmax = option_new_ulong_min("max",
                                   "only output length "
                                   "if <= given maximum length",
                                   &gfmsubcallinfo->maxlength.
                                          valueunsignedlong,
                                   0,(unsigned long) 1);
  option_parser_add_option(op, optionmax);

  optionoutput = option_new_stringarray("output",
                   doms
                     ? "set output flags (sequence, querypos, subjectpos)"
                     : "set output flags (sequence, querypos)",
                   flagsoutputoption);
  option_parser_add_option(op, optionoutput);

  if (doms)
  {
    optionverify = option_new_bool("verify","verify witness positions",
                                   &gfmsubcallinfo->verifywitnesspos,
                                   false);
    option_is_development_option(optionverify);
    option_parser_add_option(op, optionverify);
  } else
  {
    gfmsubcallinfo->verifywitnesspos = false;
  }

  oprval = option_parser_parse(op, &parsed_args, argc, argv,
                               versionfunc,err);

  if (oprval == OPTIONPARSER_OK)
  {
    if (option_is_set(optionfmindex))
    {
      gfmsubcallinfo->indextype = Fmindextype;
    } else
    {
      if (option_is_set(optionesaindex))
      {
        gfmsubcallinfo->indextype = Esaindextype;
      } else
      {
        if (option_is_set(optionpckindex))
        {
          gfmsubcallinfo->indextype = Packedindextype;
        } else
        {
          error_set(err,"one of the options -esa, -pck must be used");
          oprval = OPTIONPARSER_ERROR;
        }
      }
    }
    if (oprval != OPTIONPARSER_ERROR)
    {
      if (option_is_set(optionmin))
      {
         gfmsubcallinfo->minlength.defined = true;
      }
      if (option_is_set(optionmax))
      {
         gfmsubcallinfo->maxlength.defined = true;
      }
      if (!option_is_set(optionmin) && !option_is_set(optionmax))
      {
        error_set(err,"one of the options -min or -max must be set");
        oprval = OPTIONPARSER_ERROR;
      }
    }
    if (oprval != OPTIONPARSER_ERROR)
    {
      if (gfmsubcallinfo->minlength.defined &&
          gfmsubcallinfo->maxlength.defined)
      {
        if (gfmsubcallinfo->maxlength.valueunsignedlong <
            gfmsubcallinfo->minlength.valueunsignedlong)
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
          if (doms)
          {
            if (optionaddbitmask(msgfmsubmodedesctable,
                                 sizeof (msgfmsubmodedesctable)/
                                 sizeof (msgfmsubmodedesctable[0]),
                                 &gfmsubcallinfo->showmode,
                                 "-output",
                                 strarray_get(flagsoutputoption,i),
                                 err) != 0)
            {
              oprval = OPTIONPARSER_ERROR;
              break;
            }
          } else
          {
            if (optionaddbitmask(gfmsubmodedesctable,
                                 sizeof (gfmsubmodedesctable)/
                                 sizeof (gfmsubmodedesctable[0]),
                                 &gfmsubcallinfo->showmode,
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

static bool dotestsequence(bool doms,const Gfmsubcallinfo *gfmsubcallinfo)
{
  if (gfmsubcallinfo->indextype == Packedindextype &&
      gfmsubcallinfo->verifywitnesspos &&
      doms)
  {
    return true;
  }
  return false;
}

static int gt_greedyfwdmat(bool doms,int argc, const char **argv,Error *err)
{
  Gfmsubcallinfo gfmsubcallinfo;
  Fmindex fmindex;
  Suffixarray suffixarray;
  void *packedindex = NULL;
  Verboseinfo *verboseinfo;
  bool haserr = false;
  Alphabet *alphabet = NULL;
  unsigned int prefixlength = 0;
  Seqpos totallength;

  error_check(err);
  switch (parsegfmsub(doms,&gfmsubcallinfo, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      str_delete(gfmsubcallinfo.indexname);
      strarray_delete(gfmsubcallinfo.queryfilenames);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      str_delete(gfmsubcallinfo.indexname);
      strarray_delete(gfmsubcallinfo.queryfilenames);
      return 0;
  }
  verboseinfo = newverboseinfo(false);
  if (gfmsubcallinfo.indextype == Fmindextype)
  {
    if (mapfmindex (&fmindex,gfmsubcallinfo.indexname,
                    verboseinfo,err) != 0)
    {
      haserr = true;
    } else
    {
      alphabet = fmindex.alphabet;
    }
    totallength = fmindex.bwtlength-1;
  } else
  {
    unsigned int mappedbits;

    if (gfmsubcallinfo.indextype == Esaindextype)
    {
      mappedbits = SARR_ESQTAB | SARR_SUFTAB
#undef WITHBCKTAB
#ifdef WITHBCKTAB
                   | SARR_BCKTAB
#endif
                   ;
    } else
    {
      if (dotestsequence(doms,&gfmsubcallinfo))
      {
        mappedbits = SARR_ESQTAB;
      } else
      {
        mappedbits = 0;
      }
    }
    if (mapsuffixarray(&suffixarray,
                       &totallength,
                       mappedbits,
                       gfmsubcallinfo.indexname,
                       verboseinfo,
                       err) != 0)
    {
      haserr = true;
    } else
    {
      alphabet = suffixarray.alpha;
      prefixlength = suffixarray.prefixlength;
    }
    if (!haserr)
    {
      if (gfmsubcallinfo.indextype == Packedindextype)
      {
        packedindex = loadvoidBWTSeqForSA(gfmsubcallinfo.indexname,
                                          &suffixarray,
                                          totallength, err);
        if (packedindex == NULL)
        {
          haserr = true;
        }
      }
    }
  }
  if (!haserr)
  {
    const void *theindex;
    Greedygmatchforwardfunction gmatchforwardfunction;

    if (gfmsubcallinfo.indextype == Fmindextype)
    {
      theindex = (const void *) &fmindex;
      if (doms)
      {
        gmatchforwardfunction = skfmmstats;
      } else
      {
        gmatchforwardfunction = skfmuniqueforward;
      }
    } else
    {
      if (gfmsubcallinfo.indextype == Esaindextype)
      {
        theindex = (const void *) &suffixarray;
        if (doms)
        {
          gmatchforwardfunction = suffixarraymstats;
        } else
        {
          gmatchforwardfunction = suffixarrayuniqueforward;
        }
      } else
      {
        assert(gfmsubcallinfo.indextype == Packedindextype);
        theindex = (const void *) packedindex;
        if (doms)
        {
          gmatchforwardfunction = voidpackedindexmstatsforward;
        } else
        {
          gmatchforwardfunction = voidpackedindexuniqueforward;
        }
      }
    }
    if (!haserr)
    {
#ifdef WITHBCKTAB
      if (prefixlength > 0 &&
          gfmsubcallinfo.indextype == Esaindextype &&
          runsubstringiteration(gmatchforwardfunction,
                                theindex,
                                totallength,
                                suffixarray.bcktab,
                                suffixarray.countspecialcodes,
                                alphabet,
                                prefixlength,
                                gfmsubcallinfo.queryfilenames,
                                err) != 0)

      {
        haserr = true;
      }
#endif
      if (!haserr &&
          findsubquerygmatchforward(dotestsequence(doms,&gfmsubcallinfo)
                                      ? suffixarray.encseq
                                      : NULL,
                                    theindex,
                                    totallength,
                                    gmatchforwardfunction,
                                    alphabet,
                                    gfmsubcallinfo.queryfilenames,
                                    gfmsubcallinfo.minlength,
                                    gfmsubcallinfo.maxlength,
                                    (gfmsubcallinfo.showmode & SHOWSEQUENCE)
                                           ? true : false,
                                    (gfmsubcallinfo.showmode & SHOWQUERYPOS)
                                           ? true : false,
                                    (gfmsubcallinfo.showmode & SHOWSUBJECTPOS)
                                           ? true : false,
                                    err) != 0)
      {
        haserr = true;
      }
    }
  }
  if (gfmsubcallinfo.indextype == Fmindextype)
  {
    freefmindex(&fmindex);
  } else
  {
    if (gfmsubcallinfo.indextype == Packedindextype && packedindex != NULL)
    {
      deletevoidBWTSeq(packedindex);
    }
    freesuffixarray(&suffixarray);
  }
  freeverboseinfo(&verboseinfo);
  str_delete(gfmsubcallinfo.indexname);
  strarray_delete(gfmsubcallinfo.queryfilenames);
  return haserr ? -1 : 0;
}

int gt_uniquesub(int argc, const char **argv, Error *err)
{
  return gt_greedyfwdmat(false,argc, argv, err);
}

int gt_matchingstatistics(int argc, const char **argv, Error *err)
{
  return gt_greedyfwdmat(true,argc, argv, err);
}
