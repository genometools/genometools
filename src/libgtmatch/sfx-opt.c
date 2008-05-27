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

#include <string.h>
#include <inttypes.h>
#include "libgtcore/error.h"
#include "libgtcore/getbasename.h"
#include "libgtcore/ma.h"
#include "libgtcore/option.h"
#include "libgtcore/str.h"
#include "libgtcore/versionfunc.h"
#include "readmode-def.h"
#include "sfx-optdef.h"
#include "verbose-def.h"
#include "stamp.h"
#include "libgtmatch/eis-bwtseq-param.h"

static OPrval parse_options(int *parsed_args,
                            bool doesa,
                            Suffixeratoroptions *so,
                            int argc, const char **argv, Error *err)
{
  OptionParser *op;
  Option *option,
         *optionsmap,
         *optiondna,
         *optionprotein,
         *optionplain,
         *optionsuf,
         *optionlcp,
         *optionbwt,
         *optionpl,
         *optionmaxdepth,
         *optioncmpcharbychar,
         *optionindexname,
         *optiondb,
         *optionii,
         *optionsat,
         *optiondir,
         *optionstorespecialcodes,
         *optionmaxwidthrealmedian,
         *optionmaxbltriesort,
         *optiondes;
  OPrval oprval;
  Str *dirarg = str_new();

  error_check(err);
  op = option_parser_new("[option ...] (-db file [...] | -ii index)",
                         doesa ? "Compute enhanced suffix array."
                               : "Compute packed index.");
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  optiondb = option_new_filenamearray("db","specify database files",
                                      so->filenametab);
  option_parser_add_option(op, optiondb);

  optionii = option_new_filename("ii","specify sequence index created "
                                      " previously by -tis option",
                                 so->str_inputindex);
  option_is_mandatory_either(optiondb,optionii);
  option_parser_add_option(op, optionii);

  optionsmap = option_new_string("smap",
                                 "specify file containing a symbol mapping",
                                 so->str_smap, NULL);
  option_parser_add_option(op, optionsmap);

  optiondna = option_new_bool("dna","input is DNA sequence",
                              &so->isdna,false);
  option_parser_add_option(op, optiondna);

  optionprotein = option_new_bool("protein","input is Protein sequence",
                                  &so->isprotein,false);
  option_parser_add_option(op, optionprotein);

  optionplain = option_new_bool("plain","process as plain text",
                                &so->isplain,false);
  option_parser_add_option(op, optionplain);

  optiondir = option_new_string("dir",
                                "specify reading direction "
                                "(fwd, cpl, rev, rcl)",
                                 dirarg, "fwd");
  option_parser_add_option(op, optiondir);

  optionindexname = option_new_string("indexname",
                                      "specify name for index to be generated",
                                      so->str_indexname, NULL);
  option_parser_add_option(op, optionindexname);

  optionpl = option_new_uint_min("pl",
                                 "specify prefix length for bucket sort\n"
                                 "recommendation: use without argument;\n"
                                 "then a reasonable prefix length is "
                                 "automatically determined",
                                 &so->prefixlength,
                                 PREFIXLENGTH_AUTOMATIC,
                                 1U);
  option_argument_is_optional(optionpl);
  option_parser_add_option(op, optionpl);

  if (doesa)
  {
    optionmaxdepth = option_new_uint_min("maxdepth",
                                         "restrict suffix sorting to prefixes "
                                         "of the given length",
                                         &so->maxdepth.valueunsignedint,
                                         MAXDEPTH_AUTOMATIC,
                                         1U);
    option_is_development_option(optionmaxdepth);
    option_argument_is_optional(optionmaxdepth);
    option_parser_add_option(op, optionmaxdepth);

  } else
  {
    optionmaxdepth = NULL;
  }
  optioncmpcharbychar = option_new_bool("cmpcharbychar",
                         "compare suffixes by character wise comparisons",
                         &so->sfxstrategy.cmpcharbychar,false);
  option_is_development_option(optioncmpcharbychar);
  option_parser_add_option(op, optioncmpcharbychar);

  optionmaxwidthrealmedian = option_new_ulong(
                                     "maxwidthrealmedian",
                                     "compute real median for intervals of "
                                     "at most the given widthprefixes",
                                     &so->sfxstrategy.maxwidthrealmedian,
                                     1U);
  option_is_development_option(optionmaxwidthrealmedian);
  option_parser_add_option(op, optionmaxwidthrealmedian);

  optionmaxbltriesort = option_new_ulong("maxbltriesort",
                                         "all intervals of specified size and "
                                         "smaller are sorted by blind trie "
                                         "sorting algorithm",
                                         &so->sfxstrategy.maxbltriesort,
                                         10U);
  option_is_development_option(optionmaxbltriesort);
  option_parser_add_option(op, optionmaxbltriesort);

  optionstorespecialcodes
    = option_new_bool("storespecialcodes",
                      "store special codes (this may speed up the program)",
                      &so->sfxstrategy.storespecialcodes,false);
  option_is_development_option(optionstorespecialcodes);
  option_parser_add_option(op, optionstorespecialcodes);

  option = option_new_uint_min("parts",
                               "specify number of parts in which the "
                               "sequence is processed",
                               &so->numofparts,
                               1U,
                               1U);
  option_is_development_option(option);
  option_parser_add_option(op, option);

  optionsat = option_new_string("sat",
                                "specify kind of sequence representation",
                                so->str_sat, NULL);
  option_parser_add_option(op, optionsat);

  option = option_new_bool("tis",
                           "output transformed and encoded input "
                           "sequence to file",
                           &so->outtistab,
                           false);
  option_parser_add_option(op, option);

  optiondes = option_new_bool("des",
                              "output sequence descriptions to file ",
                              &so->outdestab,
                              false);
  option_parser_add_option(op, optiondes);

  if (doesa)
  {
    optionsuf = option_new_bool("suf",
                                "output suffix array (suftab) to file",
                                &so->outsuftab,
                                false);
    option_parser_add_option(op, optionsuf);

    optionlcp = option_new_bool("lcp",
                                "output lcp table (lcptab) to file",
                                &so->outlcptab,
                                false);
    option_parser_add_option(op, optionlcp);

    optionbwt = option_new_bool("bwt",
                                "output Burrows-Wheeler Transformation "
                                "(bwttab) to file",
                                &so->outbwttab,
                                false);
    option_parser_add_option(op, optionbwt);
  } else
  {
    optionsuf = optionlcp = optionbwt = NULL;
  }
  option = option_new_bool("bck",
                           "output bucket table to file",
                           &so->outbcktab,
                           false);
  option_parser_add_option(op, option);

  if (!doesa)
  {
    registerPackedIndexOptions(op, &so->bwtIdxParams, BWTDEFOPT_CONSTRUCTION,
                               so->str_indexname);
  }
  option = option_new_bool("showtime",
                           "show the time of the different computation phases",
                           &so->showtime,
                           false);
  option_parser_add_option(op, option);

  option = option_new_bool("v",
                           "be verbose ",
                           &so->beverbose,
                           false);
  option_parser_add_option(op, option);

  option_exclude(optionii, optiondb);
  option_exclude(optionii, optiondir);
  option_exclude(optionii, optionsmap);
  option_exclude(optionii, optiondna);
  option_exclude(optionii, optionprotein);
  option_exclude(optionii, optionplain);
  option_exclude(optionii, optionsat);
  option_exclude(optionsmap, optiondna);
  option_exclude(optionsmap, optionprotein);
  option_exclude(optiondna, optionprotein);
  if (doesa)
  {
    assert(optionmaxdepth != NULL);
    option_exclude(optionmaxdepth, optionlcp);
                   /* because lcp table may be incorrect. XXX change later */
    option_exclude(optionmaxdepth, optionbwt);
                   /* because bwt table may be incorrect. XXX change later */
  }
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc,
                               err);
  if (oprval == OPTIONPARSER_OK)
  {
    if (option_is_set(optiondb) && strarray_size(so->filenametab) == 0)
    {
      error_set(err,"missing argument to option -db");
      oprval = OPTIONPARSER_ERROR;
    } else
    {
      if (!option_is_set(optionindexname))
      {
        if (option_is_set(optiondb))
        {
          if (strarray_size(so->filenametab) > 1UL)
          {
            error_set(err,"if more than one input file is given, then "
                              "option -indexname is mandatory");
            oprval = OPTIONPARSER_ERROR;
          } else
          {
            char *basenameptr;

            basenameptr = getbasename(strarray_get(so->filenametab,0));
            str_set(so->str_indexname,basenameptr);
            ma_free(basenameptr);
          }
        } else
        {
          char *basenameptr;

          basenameptr = getbasename(str_get(so->str_inputindex));
          str_set(so->str_indexname,basenameptr);
          ma_free(basenameptr);
        }
      }
    }
  }
  if (oprval == OPTIONPARSER_OK)
  {
    if (option_is_set(optionplain))
    {
      if (!option_is_set(optiondna) &&
         !option_is_set(optionprotein) &&
         !option_is_set(optionsmap))
      {
        error_set(err,"if option -plain is used, then any of the options "
                          "-dna, -protein, or -smap is mandatory");
        oprval = OPTIONPARSER_ERROR;
      }
    }
  }
  if (oprval == OPTIONPARSER_OK && doesa)
  {
    assert(optionmaxdepth != NULL);
    if (option_is_set(optionmaxdepth))
    {
      so->maxdepth.defined = true;
    }
  }
  if (oprval == OPTIONPARSER_OK && !doesa)
  {
    computePackedIndexDefaults(&so->bwtIdxParams, BWTBaseFeatures);
  }
  option_parser_delete(op);
  if (oprval == OPTIONPARSER_OK && *parsed_args != argc)
  {
    error_set(err,"superfluous program parameters");
    oprval = OPTIONPARSER_ERROR;
  }
  if (oprval == OPTIONPARSER_OK)
  {
    int retval = parsereadmode(str_get(dirarg),err);

    if (retval < 0)
    {
      oprval = OPTIONPARSER_ERROR;
    } else
    {
      so->readmode = (Readmode) retval;
    }
  }
  str_delete(dirarg);
  return oprval;
}

static void showoptions(const Suffixeratoroptions *so)
{
  unsigned long i;

  if (str_length(so->str_smap) > 0)
  {
    showdefinitelyverbose("smap=\"%s\"",str_get(so->str_smap));
  }
  if (so->isdna)
  {
    showdefinitelyverbose("dna=yes");
  }
  if (so->isprotein)
  {
    showdefinitelyverbose("protein=yes");
  }
  if (so->isplain)
  {
    showdefinitelyverbose("plain=yes");
  }
  showdefinitelyverbose("indexname=\"%s\"",str_get(so->str_indexname));
  if (so->prefixlength == PREFIXLENGTH_AUTOMATIC)
  {
    showdefinitelyverbose("prefixlength=automatic");
  } else
  {
    showdefinitelyverbose("prefixlength=%u",so->prefixlength);
  }
  showdefinitelyverbose("storespecialcodes=%s",
                        so->sfxstrategy.storespecialcodes ? "true" : "false");
  showdefinitelyverbose("parts=%u",so->numofparts);
  for (i=0; i<strarray_size(so->filenametab); i++)
  {
    showdefinitelyverbose("inputfile[%lu]=%s",i,
                          strarray_get(so->filenametab,i));
  }
  if (str_length(so->str_inputindex) > 0)
  {
    showdefinitelyverbose("inputindex=%s",str_get(so->str_inputindex));
  }
  assert(str_length(so->str_indexname) > 0);
  showdefinitelyverbose("indexname=%s",str_get(so->str_indexname));
  showdefinitelyverbose("outtistab=%s,outsuftab=%s,outlcptab=%s,"
                        "outbwttab=%s,outbcktab=%s,outdestab=%s",
          so->outtistab ? "true" : "false",
          so->outsuftab ? "true" : "false",
          so->outlcptab ? "true" : "false",
          so->outbwttab ? "true" : "false",
          so->outbcktab ? "true" : "false",
          so->outdestab ? "true" : "false");
}

void wrapsfxoptions(Suffixeratoroptions *so)
{
  /* no checking if error occurs, since errors have been output before */
  str_delete(so->str_indexname);
  str_delete(so->str_inputindex);
  str_delete(so->str_smap);
  str_delete(so->str_sat);
  strarray_delete(so->filenametab);
}

int suffixeratoroptions(Suffixeratoroptions *so,
                        bool doesa,
                        int argc,
                        const char **argv,
                        Error *err)
{
  int parsed_args, retval = 0;
  OPrval rval;

  error_check(err);
  so->isdna = false;
  so->isprotein = false;
  so->str_inputindex = str_new();
  so->str_indexname = str_new();
  so->str_smap = str_new();
  so->str_sat = str_new();
  so->filenametab = strarray_new();
  so->prefixlength = PREFIXLENGTH_AUTOMATIC;
  so->maxdepth.defined = false;
  so->maxdepth.valueunsignedint = MAXDEPTH_AUTOMATIC;
  so->outsuftab = false;
  so->outlcptab = false;
  so->outbwttab = false;
  so->outbcktab = false;
  rval = parse_options(&parsed_args, doesa, so, argc, argv, err);
  if (rval == OPTIONPARSER_ERROR)
  {
    retval = -1;
  } else
  {
    if (rval == OPTIONPARSER_REQUESTS_EXIT)
    {
      retval = 2;
    } else
    {
      if (rval == OPTIONPARSER_OK)
      {
        if (so->beverbose)
        {
          showoptions(so);
        }
      }
    }
  }
  return retval;
}
