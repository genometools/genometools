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

#include "core/defined-types.h"
#include "core/error.h"
#include "core/option_api.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "match/eis-voiditf.h"
#include "match/esa-map.h"
#include "match/fmindex.h"
#include "match/greedyfwdmat.h"
#include "match/optionargmode.h"
#include "match/sarr-def.h"
#include "match/stamp.h"
#include "tools/gt_matstat.h"
#include "tools/gt_uniquesub.h"
#include "match/esa-minunique.h"

#include "match/fmi-fwduni.pr"
#include "match/fmi-map.pr"

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
  GtStr *indexname;
  GtStrArray *queryfilenames;
  Indextype indextype;
} Gfmsubcallinfo;

static GtOPrval parsegfmsub(bool doms,
                            Gfmsubcallinfo *gfmsubcallinfo,
                            int argc,
                            const char **argv,
                            GtError *err)
{
  GtOptionParser *op;
  GtOption *optionmin, *optionmax, *optionoutput, *optionfmindex,
         *optionesaindex, *optionpckindex, *optionquery, *optionverify;
  GtOPrval oprval;
  GtStrArray *flagsoutputoption;
  int parsed_args;
  const Optionargmodedesc msgfmsubmodedesctable[] =
  {
    {"sequence","matching sequence",SHOWSEQUENCE},
    {"querypos","position in query sequence",SHOWQUERYPOS},
    {"subjectpos","position in subject sequence",SHOWSUBJECTPOS}
  };
  const Optionargmodedesc gfmsubmodedesctable[] =
  {
    {"sequence","matching sequence",SHOWSEQUENCE},
    {"querypos","position in query sequence",SHOWQUERYPOS}
  };

  gt_error_check(err);
  gfmsubcallinfo->minlength.defined = false;
  gfmsubcallinfo->maxlength.defined = false;
  gfmsubcallinfo->showmode = 0;
  gfmsubcallinfo->indexname = gt_str_new();
  gfmsubcallinfo->queryfilenames = gt_str_array_new();
  flagsoutputoption = gt_str_array_new();

  op = gt_option_parser_new("[options ...] -query queryfile [...]",
                         doms
                         ? "Compute matching statistics."
                         : "Compute length of minumum unique prefixes."
                         );
  gt_option_parser_set_mail_address(op,"<kurtz@zbh.uni-hamburg.de>");

  optionfmindex = gt_option_new_string("fmi", "specify fmindex",
                                    gfmsubcallinfo->indexname,NULL);
  gt_option_parser_add_option(op, optionfmindex);

  optionesaindex = gt_option_new_string("esa", "specify suffix array",
                                     gfmsubcallinfo->indexname,NULL);
  gt_option_parser_add_option(op, optionesaindex);

  optionpckindex = gt_option_new_string("pck", "specify packed index",
                                     gfmsubcallinfo->indexname,NULL);
  gt_option_parser_add_option(op, optionpckindex);

  gt_option_exclude(optionfmindex,optionesaindex);
  gt_option_exclude(optionpckindex,optionesaindex);
  gt_option_exclude(optionpckindex,optionfmindex);

  optionquery = gt_option_new_filename_array("query", "specify queryfiles",
                                             gfmsubcallinfo->queryfilenames);
  gt_option_is_mandatory(optionquery);
  gt_option_parser_add_option(op, optionquery);

  optionmin = gt_option_new_ulong_min("min",
                                   "only output length "
                                   "if >= given minimum length",
                                   &gfmsubcallinfo->minlength.
                                          valueunsignedlong,
                                   0,(unsigned long) 1);
  gt_option_parser_add_option(op, optionmin);

  optionmax = gt_option_new_ulong_min("max",
                                   "only output length "
                                   "if <= given maximum length",
                                   &gfmsubcallinfo->maxlength.
                                          valueunsignedlong,
                                   0,(unsigned long) 1);
  gt_option_parser_add_option(op, optionmax);

  optionoutput = gt_option_new_string_array("output",
                   doms
                     ? "set output flags (sequence, querypos, subjectpos)"
                     : "set output flags (sequence, querypos)",
                   flagsoutputoption);
  gt_option_parser_add_option(op, optionoutput);

  if (doms)
  {
    optionverify = gt_option_new_bool("verify","verify witness positions",
                                   &gfmsubcallinfo->verifywitnesspos,
                                   false);
    gt_option_is_development_option(optionverify);
    gt_option_parser_add_option(op, optionverify);
  } else
  {
    gfmsubcallinfo->verifywitnesspos = false;
  }

  gt_option_parser_refer_to_manual(op);
  oprval = gt_option_parser_parse(op, &parsed_args, argc, argv,
                               gt_versionfunc,err);

  if (oprval == GT_OPTION_PARSER_OK)
  {
    if (gt_option_is_set(optionfmindex))
    {
      gfmsubcallinfo->indextype = Fmindextype;
    } else
    {
      if (gt_option_is_set(optionesaindex))
      {
        gfmsubcallinfo->indextype = Esaindextype;
      } else
      {
        if (gt_option_is_set(optionpckindex))
        {
          gfmsubcallinfo->indextype = Packedindextype;
        } else
        {
          gt_error_set(err,"one of the options -esa, -pck must be used");
          oprval = GT_OPTION_PARSER_ERROR;
        }
      }
    }
    if (oprval != GT_OPTION_PARSER_ERROR)
    {
      if (gt_option_is_set(optionmin))
      {
         gfmsubcallinfo->minlength.defined = true;
      }
      if (gt_option_is_set(optionmax))
      {
         gfmsubcallinfo->maxlength.defined = true;
      }
      if (!gt_option_is_set(optionmin) && !gt_option_is_set(optionmax))
      {
        gt_error_set(err,"one of the options -min or -max must be set");
        oprval = GT_OPTION_PARSER_ERROR;
      }
    }
    if (oprval != GT_OPTION_PARSER_ERROR)
    {
      if (gfmsubcallinfo->minlength.defined &&
          gfmsubcallinfo->maxlength.defined)
      {
        if (gfmsubcallinfo->maxlength.valueunsignedlong <
            gfmsubcallinfo->minlength.valueunsignedlong)
        {
          gt_error_set(err,"minvalue must be smaller or equal than maxvalue");
          oprval = GT_OPTION_PARSER_ERROR;
        }
      }
    }
    if (oprval != GT_OPTION_PARSER_ERROR && gt_option_is_set(optionoutput))
    {
      if (gt_str_array_size(flagsoutputoption) == 0)
      {
        gt_error_set(err,"missing arguments to option -output");
        oprval = GT_OPTION_PARSER_ERROR;
      } else
      {
        unsigned long i;

        for (i=0; i<gt_str_array_size(flagsoutputoption); i++)
        {
          if (doms)
          {
            if (gt_optionargaddbitmask(msgfmsubmodedesctable,
                                 sizeof (msgfmsubmodedesctable)/
                                 sizeof (msgfmsubmodedesctable[0]),
                                 &gfmsubcallinfo->showmode,
                                 "-output",
                                 gt_str_array_get(flagsoutputoption,i),
                                 err) != 0)
            {
              oprval = GT_OPTION_PARSER_ERROR;
              break;
            }
          } else
          {
            if (gt_optionargaddbitmask(gfmsubmodedesctable,
                                    sizeof (gfmsubmodedesctable)/
                                    sizeof (gfmsubmodedesctable[0]),
                                    &gfmsubcallinfo->showmode,
                                    "-output",
                                    gt_str_array_get(flagsoutputoption,i),
                                    err) != 0)
            {
              oprval = GT_OPTION_PARSER_ERROR;
              break;
            }
          }
        }
      }
    }
  }
  gt_str_array_delete(flagsoutputoption);
  gt_option_parser_delete(op);
  if (oprval == GT_OPTION_PARSER_OK && parsed_args != argc)
  {
    gt_error_set(err,"superfluous program parameters");
    oprval = GT_OPTION_PARSER_ERROR;
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

static int gt_greedyfwdmat(bool doms,int argc, const char **argv,GtError *err)
{
  Gfmsubcallinfo gfmsubcallinfo;
  Fmindex fmindex;
  Suffixarray suffixarray;
  void *packedindex = NULL;
  GtLogger *logger = NULL;
  bool haserr = false;
  const GtAlphabet *alphabet = NULL;
#ifdef WITHBCKTAB
  unsigned int prefixlength = 0;
#endif
  unsigned long totallength;
  bool gt_mapfmindexfail = false;

  gt_error_check(err);
  switch (parsegfmsub(doms,&gfmsubcallinfo, argc, argv, err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR:
      gt_str_delete(gfmsubcallinfo.indexname);
      gt_str_array_delete(gfmsubcallinfo.queryfilenames);
      return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT:
      gt_str_delete(gfmsubcallinfo.indexname);
      gt_str_array_delete(gfmsubcallinfo.queryfilenames);
      return 0;
  }
  logger = gt_logger_new(false, GT_LOGGER_DEFLT_PREFIX, stdout);
  if (gfmsubcallinfo.indextype == Fmindextype)
  {
    if (gt_mapfmindex(&fmindex,gt_str_get(gfmsubcallinfo.indexname),
                      logger, err) != 0)
    {
      haserr = true;
      gt_mapfmindexfail = true;
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
    if (gt_mapsuffixarray(&suffixarray,
                       mappedbits,
                       gt_str_get(gfmsubcallinfo.indexname),
                       logger,
                       err) != 0)
    {
      haserr = true;
      totallength = 0;
    } else
    {
      alphabet = gt_encseq_alphabet(suffixarray.encseq);
#ifdef WITHBCKTAB
      prefixlength = suffixarray.prefixlength;
#endif
      totallength = gt_encseq_total_length(suffixarray.encseq);
    }
    if (!haserr)
    {
      if (gfmsubcallinfo.indextype == Packedindextype)
      {
        packedindex =
          gt_loadvoidBWTSeqForSA(gt_str_get(gfmsubcallinfo.indexname),
                                 false,
                                 err);
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
        gmatchforwardfunction = gt_skfmmstats;
      } else
      {
        gmatchforwardfunction = gt_skfmuniqueforward;
      }
    } else
    {
      if (gfmsubcallinfo.indextype == Esaindextype)
      {
        theindex = (const void *) &suffixarray;
        if (doms)
        {
          gmatchforwardfunction = gt_suffixarraymstats;
        } else
        {
          gmatchforwardfunction = gt_suffixarrayuniqueforward;
        }
      } else
      {
        gt_assert(gfmsubcallinfo.indextype == Packedindextype);
        theindex = (const void *) packedindex;
        if (doms)
        {
          gmatchforwardfunction = gt_voidpackedindexmstatsforward;
        } else
        {
          gmatchforwardfunction = gt_voidpackedindexuniqueforward;
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
          gt_findsubquerygmatchforward(dotestsequence(doms,&gfmsubcallinfo)
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
    if (!gt_mapfmindexfail)
    {
      gt_freefmindex(&fmindex);
    }
  } else
  {
    if (gfmsubcallinfo.indextype == Packedindextype && packedindex != NULL)
    {
      gt_deletevoidBWTSeq(packedindex);
    }
    gt_freesuffixarray(&suffixarray);
  }
  gt_logger_delete(logger);
  gt_str_delete(gfmsubcallinfo.indexname);
  gt_str_array_delete(gfmsubcallinfo.queryfilenames);
  return haserr ? -1 : 0;
}

int gt_uniquesub(int argc, const char **argv, GtError *err)
{
  return gt_greedyfwdmat(false,argc, argv, err);
}

int gt_matstat(int argc, const char **argv, GtError *err)
{
  return gt_greedyfwdmat(true,argc, argv, err);
}
