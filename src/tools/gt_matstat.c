/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2013 Center for Bioinformatics, University of Hamburg

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
#include "core/ma.h"
#include "core/option_api.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "match/eis-voiditf.h"
#include "match/esa-map.h"
#include "match/esa-minunique.h"
#include "match/fmindex.h"
#include "match/greedyfwdmat.h"
#include "match/optionargmode.h"
#include "match/sarr-def.h"
#include "tools/gt_matstat.h"
#include "tools/gt_uniquesub.h"

#include "match/fmi-fwduni.h"
#include "match/fmi-map.h"

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
  GtStrArray *queryfilenames, *flagsoutputoption;
  Indextype indextype;
  bool doms;
  GtOption *optionmin, *optionmax, *optionoutput, *optionfmindex,
           *optionesaindex, *optionpckindex, *optionquery, *optionverify;
} Gfmsubcallinfo;

static void* gt_matstat_arguments_new_generic(bool doms)
{
  Gfmsubcallinfo *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->minlength.defined = false;
  arguments->maxlength.defined = false;
  arguments->showmode = 0;
  arguments->indexname = gt_str_new();
  arguments->queryfilenames = gt_str_array_new();
  arguments->flagsoutputoption = gt_str_array_new();
  arguments->doms = doms;
  return arguments;
}

static void* gt_matstat_arguments_new_matstat(void)
{
  return gt_matstat_arguments_new_generic(true);
}

static void* gt_matstat_arguments_new_uniquesub(void)
{
  return gt_matstat_arguments_new_generic(false);
}

static void gt_matstat_arguments_delete(void *tool_arguments)
{
  Gfmsubcallinfo *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_array_delete(arguments->queryfilenames);
  gt_str_array_delete(arguments->flagsoutputoption);
  gt_str_delete(arguments->indexname);
  gt_free(arguments);
}

static GtOptionParser* gt_matstat_option_parser_new(void *tool_arguments)
{
  Gfmsubcallinfo *arguments = tool_arguments;
  GtOptionParser *op;
  gt_assert(arguments);

  op = gt_option_parser_new("[options ...] -query queryfile [...]",
                         arguments->doms
                         ? "Compute matching statistics."
                         : "Compute length of minimum unique prefixes.");
  gt_option_parser_set_mail_address(op,"<kurtz@zbh.uni-hamburg.de>");

  arguments->optionfmindex = gt_option_new_string("fmi", "specify fmindex",
                                                  arguments->indexname,NULL);
  gt_option_parser_add_option(op, arguments->optionfmindex);

  arguments->optionesaindex = gt_option_new_string("esa",
                                                   "specify suffix array",
                                                   arguments->indexname,NULL);
  gt_option_parser_add_option(op, arguments->optionesaindex);

  arguments->optionpckindex = gt_option_new_string("pck",
                                                   "specify packed index",
                                                   arguments->indexname,NULL);
  gt_option_parser_add_option(op, arguments->optionpckindex);

  gt_option_exclude(arguments->optionfmindex,arguments->optionesaindex);
  gt_option_exclude(arguments->optionpckindex,arguments->optionesaindex);
  gt_option_exclude(arguments->optionpckindex,arguments->optionfmindex);

  arguments->optionquery = gt_option_new_filename_array("query",
                                                     "specify queryfiles",
                                                     arguments->queryfilenames);
  gt_option_is_mandatory(arguments->optionquery);
  gt_option_parser_add_option(op, arguments->optionquery);

  arguments->optionmin = gt_option_new_uword_min("min",
                                   "only output length "
                                   "if >= given minimum length",
                                   &arguments->minlength.
                                          valueunsignedlong,
                                   0,(GtUword) 1);
  gt_option_parser_add_option(op, arguments->optionmin);

  arguments->optionmax = gt_option_new_uword_min("max",
                                   "only output length "
                                   "if <= given maximum length",
                                   &arguments->maxlength.
                                          valueunsignedlong,
                                   0,(GtUword) 1);
  gt_option_parser_add_option(op, arguments->optionmax);

  arguments->optionoutput = gt_option_new_string_array("output",
                   arguments->doms
                     ? "set output flags (sequence, querypos, subjectpos)"
                     : "set output flags (sequence, querypos)",
                   arguments->flagsoutputoption);
  gt_option_parser_add_option(op, arguments->optionoutput);

  if (arguments->doms) {
    arguments->optionverify = gt_option_new_bool("verify",
                                                 "verify witness positions",
                                                 &arguments->verifywitnesspos,
                                                 false);
    gt_option_is_development_option(arguments->optionverify);
    gt_option_parser_add_option(op, arguments->optionverify);
  } else
  {
    arguments->verifywitnesspos = false;
  }

  gt_option_parser_refer_to_manual(op);

  return op;
}

static int gt_matstat_arguments_check(GT_UNUSED int rest_argc,
                                      void *tool_arguments, GtError *err)
{
  Gfmsubcallinfo *arguments = tool_arguments;
  int had_err = 0;
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
  gt_assert(arguments);

  if (gt_option_is_set(arguments->optionfmindex))
    {
      arguments->indextype = Fmindextype;
    } else
    {
      if (gt_option_is_set(arguments->optionesaindex))
      {
        arguments->indextype = Esaindextype;
      } else
      {
        if (gt_option_is_set(arguments->optionpckindex))
        {
          arguments->indextype = Packedindextype;
        } else
        {
          gt_error_set(err,"one of the options -esa, -pck must be used");
          return -1;
        }
      }
    }
    if (gt_option_is_set(arguments->optionmin))
    {
       arguments->minlength.defined = true;
    }
    if (gt_option_is_set(arguments->optionmax))
    {
       arguments->maxlength.defined = true;
    }
    if (!gt_option_is_set(arguments->optionmin) &&
          !gt_option_is_set(arguments->optionmax))
    {
      gt_error_set(err,"one of the options -min or -max must be set");
      return -1;
    }
    if (arguments->minlength.defined &&
        arguments->maxlength.defined)
    {
      if (arguments->maxlength.valueunsignedlong <
          arguments->minlength.valueunsignedlong)
      {
        gt_error_set(err,"minvalue must be smaller or equal than maxvalue");
        return -1;
      }
    }
    if (gt_option_is_set(arguments->optionoutput))
    {
      if (gt_str_array_size(arguments->flagsoutputoption) == 0)
      {
        gt_error_set(err,"missing arguments to option -output");
        return -1;
      } else
      {
        GtUword i;
        for (i=0; i<gt_str_array_size(arguments->flagsoutputoption); i++)
        {
          if (arguments->doms)
          {
            if (gt_optionargaddbitmask(msgfmsubmodedesctable,
                                 sizeof (msgfmsubmodedesctable)/
                                 sizeof (msgfmsubmodedesctable[0]),
                                 &arguments->showmode,
                                 "-output",
                                 gt_str_array_get(arguments->flagsoutputoption,
                                                  i),
                                 err) != 0)
            {
              return -1;
            }
          } else
          {
            if (gt_optionargaddbitmask(gfmsubmodedesctable,
                                  sizeof (gfmsubmodedesctable)/
                                  sizeof (gfmsubmodedesctable[0]),
                                  &arguments->showmode,
                                  "-output",
                                  gt_str_array_get(arguments->flagsoutputoption,
                                                   i),
                                  err) != 0)
            {
              return -1;
            }
          }
        }
      }
    }

  return had_err;
}

static bool dotestsequence(const Gfmsubcallinfo *gfmsubcallinfo)
{
  if (gfmsubcallinfo->indextype == Packedindextype &&
      gfmsubcallinfo->verifywitnesspos &&
      gfmsubcallinfo->doms)
  {
    return true;
  }
  return false;
}

static int gt_matstat_runner(GT_UNUSED int argc, GT_UNUSED const char **argv,
                             GT_UNUSED int parsed_args,
                             void *tool_arguments, GtError *err)
{
  Gfmsubcallinfo *arguments = tool_arguments;
  Fmindex fmindex;
  Suffixarray suffixarray;
  void *packedindex = NULL;
  GtLogger *logger = NULL;
  bool haserr = false;
  const GtAlphabet *alphabet = NULL;
#ifdef WITHBCKTAB
  unsigned int prefixlength = 0;
#endif
  GtUword totallength;
  bool gt_mapfmindexfail = false;
  gt_error_check(err);
  gt_assert(arguments);

  logger = gt_logger_new(false, GT_LOGGER_DEFLT_PREFIX, stdout);
  if (arguments->indextype == Fmindextype)
  {
    if (gt_mapfmindex(&fmindex,gt_str_get(arguments->indexname),
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

    if (arguments->indextype == Esaindextype)
    {
      mappedbits = SARR_ESQTAB | SARR_SUFTAB
#undef WITHBCKTAB
#ifdef WITHBCKTAB
                   | SARR_BCKTAB
#endif
                   ;
    } else
    {
      if (dotestsequence(arguments))
      {
        mappedbits = SARR_ESQTAB;
      } else
      {
        mappedbits = 0;
      }
    }
    if (gt_mapsuffixarray(&suffixarray,
                       mappedbits,
                       gt_str_get(arguments->indexname),
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
      if (arguments->indextype == Packedindextype)
      {
        packedindex =
          gt_loadvoidBWTSeqForSA(gt_str_get(arguments->indexname),
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

    if (arguments->indextype == Fmindextype)
    {
      theindex = (const void *) &fmindex;
      if (arguments->doms)
      {
        gmatchforwardfunction = gt_skfmmstats;
      } else
      {
        gmatchforwardfunction = gt_skfmuniqueforward;
      }
    } else
    {
      if (arguments->indextype == Esaindextype)
      {
        theindex = (const void *) &suffixarray;
        if (arguments->doms)
        {
          gmatchforwardfunction = gt_suffixarraymstats;
        } else
        {
          gmatchforwardfunction = gt_suffixarrayuniqueforward;
        }
      } else
      {
        gt_assert(arguments->indextype == Packedindextype);
        theindex = (const void *) packedindex;
        if (arguments->doms)
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
          arguments->indextype == Esaindextype &&
          runsubstringiteration(gmatchforwardfunction,
                                theindex,
                                totallength,
                                suffixarray.bcktab,
                                suffixarray.countspecialcodes,
                                alphabet,
                                prefixlength,
                                arguments->queryfilenames,
                                err) != 0)

      {
        haserr = true;
      }
#endif
      if (!haserr &&
          gt_findsubquerygmatchforward(dotestsequence(arguments)
                                      ? suffixarray.encseq
                                      : NULL,
                                      theindex,
                                      totallength,
                                      gmatchforwardfunction,
                                      alphabet,
                                      arguments->queryfilenames,
                                      arguments->minlength,
                                      arguments->maxlength,
                                      (arguments->showmode & SHOWSEQUENCE)
                                             ? true : false,
                                      (arguments->showmode & SHOWQUERYPOS)
                                             ? true : false,
                                      (arguments->showmode & SHOWSUBJECTPOS)
                                             ? true : false,
                                      err) != 0)
      {
        haserr = true;
      }
    }
  }
  if (arguments->indextype == Fmindextype)
  {
    if (!gt_mapfmindexfail)
    {
      gt_freefmindex(&fmindex);
    }
  } else
  {
    if (arguments->indextype == Packedindextype && packedindex != NULL)
    {
      gt_deletevoidBWTSeq(packedindex);
    }
    gt_freesuffixarray(&suffixarray);
  }
  gt_logger_delete(logger);

  return haserr ? -1 : 0;;
}

GtTool* gt_matstat(void)
{
  return gt_tool_new(gt_matstat_arguments_new_matstat,
                     gt_matstat_arguments_delete,
                     gt_matstat_option_parser_new,
                     gt_matstat_arguments_check,
                     gt_matstat_runner);
}

GtTool* gt_uniquesub(void)
{
  return gt_tool_new(gt_matstat_arguments_new_uniquesub,
                     gt_matstat_arguments_delete,
                     gt_matstat_option_parser_new,
                     gt_matstat_arguments_check,
                     gt_matstat_runner);
}
