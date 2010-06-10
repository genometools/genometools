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
#include "core/error.h"
#include "core/basename_api.h"
#include "core/ma_api.h"
#include "core/option.h"
#include "core/str.h"
#include "core/versionfunc.h"
#include "core/readmode.h"
#include "sfx-optdef.h"
#include "core/logger.h"
#include "stamp.h"
#include "eis-bwtseq-param.h"

static GtOPrval parse_options(int *parsed_args,
                              bool doesa,
                              Suffixeratoroptions *so,
                              int argc, const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOption *option,
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
         *optionalgbounds,
         *optionparts,
         *optiondifferencecover,
         *optiondes,
         *optionsds,
         *optionkys,
         *optionshowprogress;
  GtOPrval oprval;
  const char *maxdepthmsg = "option of -maxdepth must the keyword abs, the "
                            "keyword he or an integer";
  GtStr *dirarg = gt_str_new();

  gt_error_check(err);
  op = gt_option_parser_new("[option ...] (-db file [...] | -ii index)",
                            doesa ? "Compute enhanced suffix array."
                                  : "Compute packed index.");
  gt_option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  optiondb = gt_option_new_filenamearray("db","specify database files",
                                         so->fn2encopt.filenametab);
  gt_option_parser_add_option(op, optiondb);

  optionii = gt_option_new_filename("ii","specify sequence index created "
                                    "previously by -tis option",
                                    so->inputindex);
  gt_option_is_mandatory_either(optiondb,optionii);
  gt_option_parser_add_option(op, optionii);

  optionsmap = gt_option_new_string("smap",
                                    "specify file containing a symbol mapping",
                                    so->fn2encopt.smap, NULL);
  gt_option_parser_add_option(op, optionsmap);

  optiondna = gt_option_new_bool("dna","input is DNA sequence",
                                 &so->fn2encopt.isdna,false);
  gt_option_parser_add_option(op, optiondna);

  optionprotein = gt_option_new_bool("protein","input is Protein sequence",
                                     &so->fn2encopt.isprotein,false);
  gt_option_parser_add_option(op, optionprotein);

  optionplain = gt_option_new_bool("plain","process as plain text",
                                   &so->fn2encopt.isplain,false);
  gt_option_parser_add_option(op, optionplain);

  optiondir = gt_option_new_string("dir",
                                   "specify reading direction "
                                   "(fwd, cpl, rev, rcl)",
                                   dirarg, "fwd");
  gt_option_parser_add_option(op, optiondir);

  optionindexname = gt_option_new_string("indexname",
                                         "specify name for index to "
                                         "be generated",
                                         so->fn2encopt.indexname, NULL);
  gt_option_parser_add_option(op, optionindexname);

  optionpl = gt_option_new_uint_min("pl",
                                    "specify prefix length for bucket sort\n"
                                    "recommendation: use without argument;\n"
                                    "then a reasonable prefix length is "
                                    "automatically determined",
                                    &so->prefixlength,
                                    PREFIXLENGTH_AUTOMATIC,
                                    1U);
  gt_option_argument_is_optional(optionpl);
  gt_option_parser_add_option(op, optionpl);

  optiondifferencecover = gt_option_new_uint_min("dc",
                                    "specify difference cover value",
                                    &so->sfxstrategy.differencecover,
                                    0,
                                    4U);
  gt_option_argument_is_optional(optiondifferencecover);
  gt_option_parser_add_option(op, optiondifferencecover);

  if (doesa)
  {
    optionmaxdepth = gt_option_new_string("maxdepth",
                                          "stop shallow suffix sorting "
                                          "at prefixes of the given length",
                                          so->maxdepth,
                                          NULL);
    gt_option_is_development_option(optionmaxdepth);
    gt_option_argument_is_optional(optionmaxdepth);
    gt_option_parser_add_option(op, optionmaxdepth);

  } else
  {
    optionmaxdepth = NULL;
  }
  optioncmpcharbychar = gt_option_new_bool("cmpcharbychar",
                                           "compare suffixes character "
                                           "by character",
                                           &so->sfxstrategy.cmpcharbychar,
                                           false);
  gt_option_is_development_option(optioncmpcharbychar);
  gt_option_parser_add_option(op, optioncmpcharbychar);

  optionmaxwidthrealmedian = gt_option_new_ulong("maxwidthrealmedian",
                                                 "compute real median for "
                                                 "intervals of at most the "
                                                 "given widthprefixes",
                                                 &so->sfxstrategy.
                                                      maxwidthrealmedian,
                                                 1U);
  gt_option_is_development_option(optionmaxwidthrealmedian);
  gt_option_parser_add_option(op, optionmaxwidthrealmedian);

  optionalgbounds
    = gt_option_new_stringarray("algbds",
                                "length boundaries for the different "
                                "algorithms to sort buckets of suffixes\n"
                                "first number: maxbound for insertion sort\n"
                                "second number: maxbound for blindtrie sort\n"
                                "third number: maxbound for counting sort\n",
                                so->algbounds);
  gt_option_is_development_option(optionalgbounds);
  gt_option_parser_add_option(op, optionalgbounds);
  so->optionalgboundsref = gt_option_ref(optionalgbounds);

  optionstorespecialcodes
    = gt_option_new_bool("storespecialcodes",
                         "store special codes (this may speed up the program)",
                         &so->sfxstrategy.storespecialcodes,false);
  gt_option_is_development_option(optionstorespecialcodes);
  gt_option_parser_add_option(op, optionstorespecialcodes);

  optionparts = gt_option_new_uint_min("parts",
                                       "specify number of parts in which the "
                                       "sequence is processed",
                                       &so->numofparts,
                                       1U,
                                       1U);
  gt_option_is_development_option(optionparts);
  gt_option_parser_add_option(op, optionparts);

  optionsat = gt_option_new_string("sat",
                                   "specify kind of sequence representation",
                                   so->fn2encopt.sat, NULL);
  gt_option_parser_add_option(op, optionsat);

  option = gt_option_new_bool("tis",
                              "output transformed and encoded input "
                              "sequence to file (deprecated, kept for "
                              "compatibility reasons)",
                              &so->fn2encopt.outtistab,
                              false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  option = gt_option_new_bool("ssp",
                              "output sequence separator positions to file",
                              &so->fn2encopt.outssptab,
                              false);
  gt_option_parser_add_option(op, option);

  optiondes = gt_option_new_bool("des",
                                 "output sequence descriptions to file",
                                 &so->fn2encopt.outdestab,
                                 false);
  gt_option_parser_add_option(op, optiondes);

  optionsds = gt_option_new_bool("sds",
                                 "output sequence description separator "
                                 "positions to file",
                                 &so->fn2encopt.outsdstab,
                                 false);
  gt_option_parser_add_option(op, optionsds);

  optionkys = gt_option_new_string("kys",
                                   "output/sort according to keys of the form "
                                   "|key| in fasta header",
                                   so->optionkysargumentstring,
                                   "nosort");
  gt_option_argument_is_optional(optionkys);
  gt_option_parser_add_option(op, optionkys);

  if (doesa)
  {
    optionsuf = gt_option_new_bool("suf",
                                   "output suffix array (suftab) to file",
                                   &so->outsuftab,
                                   false);
    gt_option_parser_add_option(op, optionsuf);

    optionlcp = gt_option_new_bool("lcp",
                                   "output lcp table (lcptab) to file",
                                   &so->outlcptab,
                                   false);
    gt_option_parser_add_option(op, optionlcp);

    optionbwt = gt_option_new_bool("bwt",
                                   "output Burrows-Wheeler Transformation "
                                   "(bwttab) to file",
                                   &so->outbwttab,
                                   false);
    gt_option_parser_add_option(op, optionbwt);

    option = gt_option_new_bool("bck",
                                "output bucket table to file",
                                &so->outbcktab,
                                false);
    gt_option_parser_add_option(op, option);
  } else
  {
    optionsuf = optionlcp = optionbwt = NULL;
    gt_registerPackedIndexOptions(op, &so->bwtIdxParams, BWTDEFOPT_CONSTRUCTION,
                                  so->fn2encopt.indexname);
  }

  optionshowprogress
    = gt_option_new_bool("showprogress",
                         "show a progress bar",
                         &so->showprogress,
                         false);
  gt_option_parser_add_option(op, optionshowprogress);

  option = gt_option_new_bool("iterscan",
                              "use iteratorbased-kmer scanning",
                              &so->iteratorbasedkmerscanning,
                              false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_verbose(&so->beverbose);
  gt_option_parser_add_option(op, option);

  gt_option_exclude(optionii, optiondb);
  gt_option_exclude(optionii, optionsmap);
  gt_option_exclude(optionii, optiondna);
  gt_option_exclude(optionii, optionprotein);
  gt_option_exclude(optionii, optionplain);
  gt_option_exclude(optionii, optionsat);
  gt_option_exclude(optionsmap, optiondna);
  gt_option_exclude(optionsmap, optionprotein);
  gt_option_exclude(optiondna, optionprotein);
  gt_option_imply(optionkys,optionsds);
  gt_option_imply(optionsds,optiondes);
  if (optionmaxdepth != NULL)
  {
    gt_option_exclude(optionmaxdepth, optiondifferencecover);
  }
  if (optionlcp != NULL)
  {
    gt_option_exclude(optionlcp, optiondifferencecover);
  }
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  if (oprval == GT_OPTION_PARSER_OK)
  {
    if (gt_option_is_set(optiondb) &&
        gt_str_array_size(so->fn2encopt.filenametab) == 0)
    {
      gt_error_set(err,"missing argument to option -db");
      oprval = GT_OPTION_PARSER_ERROR;
    } else
    {
      if (!gt_option_is_set(optionindexname))
      {
        if (gt_option_is_set(optiondb))
        {
          if (gt_str_array_size(so->fn2encopt.filenametab) > 1UL)
          {
            gt_error_set(err,"if more than one input file is given, then "
                             "option -indexname is mandatory");
            oprval = GT_OPTION_PARSER_ERROR;
          } else
          {
            char *basenameptr;

            basenameptr
              = gt_basename(gt_str_array_get(so->fn2encopt.filenametab,0));
            gt_str_set(so->fn2encopt.indexname,basenameptr);
            gt_free(basenameptr);
          }
        } else
        {
          char *basenameptr;

          basenameptr = gt_basename(gt_str_get(so->inputindex));
          gt_str_set(so->fn2encopt.indexname,basenameptr);
          gt_free(basenameptr);
        }
      }
    }
  }
  if (oprval == GT_OPTION_PARSER_OK)
  {
    if (gt_option_is_set(optionplain))
    {
      if (!gt_option_is_set(optiondna) &&
         !gt_option_is_set(optionprotein) &&
         !gt_option_is_set(optionsmap))
      {
        gt_error_set(err,"if option -plain is used, then any of the options "
                         "-dna, -protein, or -smap is mandatory");
        oprval = GT_OPTION_PARSER_ERROR;
      }
    }
  }
  so->sfxstrategy.absoluteinversesuftab = false;
  if (oprval == GT_OPTION_PARSER_OK && doesa)
  {
    gt_assert(optionmaxdepth != NULL);
    if (gt_option_is_set(optionmaxdepth))
    {
      if (gt_str_length(so->maxdepth) > 0)
      {
        if (strcmp(gt_str_get(so->maxdepth),"abs") == 0)
        {
          so->sfxstrategy.absoluteinversesuftab = true;
        } else
        {
          if (strcmp(gt_str_get(so->maxdepth),"he") == 0)
          {
            so->sfxstrategy.hashexceptions = true;
          } else
          {
            long readint;

            so->sfxstrategy.hashexceptions = true;
            if (sscanf(gt_str_get(so->maxdepth),"%ld",&readint) == 1 &&
                readint >= 1L)
            {
              so->sfxstrategy.ssortmaxdepth.defined = true;
              so->sfxstrategy.ssortmaxdepth.valueunsignedint
                = (unsigned int) readint;
            } else
            {
              gt_error_set(err, "%s", maxdepthmsg);
              oprval = GT_OPTION_PARSER_ERROR;
            }
          }
        }
      }
      if (oprval != GT_OPTION_PARSER_ERROR && so->numofparts > 1U)
      {
        if (so->sfxstrategy.ssortmaxdepth.valueunsignedint
            == MAXDEPTH_AUTOMATIC)
        {
          so->sfxstrategy.streamsuftab = true;
        } else
        {
          gt_error_set(err,"option -maxdepth with argument can only be used "
                           "either without -parts or with option -parts 1");
          oprval = GT_OPTION_PARSER_ERROR;
        }
      }
      if (oprval != GT_OPTION_PARSER_ERROR)
      {
        if (so->sfxstrategy.differencecover == 0)
        {
          so->sfxstrategy.ssortmaxdepth.defined = true;
          if (so->sfxstrategy.ssortmaxdepth.valueunsignedint !=
              MAXDEPTH_AUTOMATIC)
          {
            so->sfxstrategy.cmpcharbychar = true;
          }
        }
      }
    }
  }
  if (so->sfxstrategy.differencecover > 0 && optionlcp != NULL &&
      gt_option_is_set(optionlcp))
  {
    gt_error_set(err,"option -maxdepth with argument dc can not be combined "
                     "with option lcp");
    oprval = GT_OPTION_PARSER_ERROR;
  }
  if (oprval == GT_OPTION_PARSER_OK && !doesa)
  {
    gt_computePackedIndexDefaults(&so->bwtIdxParams, BWTBaseFeatures);
  }
  if (oprval == GT_OPTION_PARSER_OK && gt_option_is_set(optionkys))
  {
    so->outkystab = true;
    if (strcmp(gt_str_get(so->optionkysargumentstring),"sort") == 0)
    {
      so->outkyssort = true;
    } else
    {
      if (strcmp(gt_str_get(so->optionkysargumentstring),"nosort") != 0)
      {
        gt_error_set(err,"illegal argument to option -kys: either use no "
                         "argument or argument \"sort\"");
        oprval = GT_OPTION_PARSER_ERROR;
      }
    }
  }
  gt_option_parser_delete(op);
  if (oprval == GT_OPTION_PARSER_OK && *parsed_args != argc)
  {
    gt_error_set(err,"superfluous program parameters");
    oprval = GT_OPTION_PARSER_ERROR;
  }
  if (oprval == GT_OPTION_PARSER_OK)
  {
    int retval = gt_readmode_parse(gt_str_get(dirarg),err);

    if (retval < 0)
    {
      oprval = GT_OPTION_PARSER_ERROR;
    } else
    {
      so->readmode = (GtReadmode) retval;
    }
  }
  gt_str_delete(dirarg);
  return oprval;
}

static void showoptions(const Suffixeratoroptions *so)
{
  unsigned long i;
  GtLogger *logger = gt_logger_new(true, GT_LOGGER_DEFLT_PREFIX, stdout);

  if (gt_str_length(so->fn2encopt.smap) > 0)
  {
    gt_logger_log_force(logger, "smap=\"%s\"", gt_str_get(so->fn2encopt.smap));
  }
  if (so->fn2encopt.isdna)
  {
    gt_logger_log_force(logger, "dna=yes");
  }
  if (so->fn2encopt.isprotein)
  {
    gt_logger_log_force(logger, "protein=yes");
  }
  if (so->fn2encopt.isplain)
  {
    gt_logger_log_force(logger, "plain=yes");
  }
  gt_logger_log_force(logger, "indexname=\"%s\"",
                      gt_str_get(so->fn2encopt.indexname));
  if (so->prefixlength == PREFIXLENGTH_AUTOMATIC)
  {
    gt_logger_log_force(logger, "prefixlength=automatic");
  } else
  {
    gt_logger_log_force(logger, "prefixlength=%u",so->prefixlength);
  }
  gt_logger_log_force(logger, "storespecialcodes=%s",
                        so->sfxstrategy.storespecialcodes ? "true" : "false");
  gt_logger_log_force(logger, "parts=%u",so->numofparts);
  for (i=0; i<gt_str_array_size(so->fn2encopt.filenametab); i++)
  {
    gt_logger_log_force(logger, "inputfile[%lu]=%s",i,
                          gt_str_array_get(so->fn2encopt.filenametab,i));
  }
  if (gt_str_length(so->inputindex) > 0)
  {
    gt_logger_log_force(logger, "inputindex=%s",gt_str_get(so->inputindex));
  }
  gt_assert(gt_str_length(so->fn2encopt.indexname) > 0);
  gt_logger_log_force(logger, "indexname=%s",
                              gt_str_get(so->fn2encopt.indexname));
  gt_logger_log_force(logger, "outtistab=%s,outsuftab=%s,outlcptab=%s,"
                        "outbwttab=%s,outbcktab=%s,outdestab=%s,"
                        "outsdstab=%s,outssptab=%s,outkystab=%s",
          so->fn2encopt.outtistab ? "true" : "false",
          so->outsuftab ? "true" : "false",
          so->outlcptab ? "true" : "false",
          so->outbwttab ? "true" : "false",
          so->outbcktab ? "true" : "false",
          so->fn2encopt.outdestab ? "true" : "false",
          so->fn2encopt.outsdstab ? "true" : "false",
          so->fn2encopt.outssptab ? "true" : "false",
          so->outkystab ? (so->outkyssort ? "true with sort" : "true") :
                          "false"
          );
  gt_logger_delete(logger);
}

void gt_wrapsfxoptions(Suffixeratoroptions *so)
{
  /* no checking if error occurs, since errors have been output before */
  gt_str_delete(so->fn2encopt.indexname);
  gt_str_delete(so->inputindex);
  gt_str_delete(so->fn2encopt.smap);
  gt_str_delete(so->fn2encopt.sat);
  gt_str_delete(so->optionkysargumentstring);
  gt_str_delete(so->maxdepth);
  gt_str_array_delete(so->fn2encopt.filenametab);
  gt_str_array_delete(so->algbounds);
  gt_option_delete(so->optionalgboundsref);
}

#define READMAXBOUND(COMP,IDX)\
        if (retval == 0)\
        {\
          arg = gt_str_array_get(so->algbounds,IDX);\
          if (sscanf(arg,"%ld",&readint) != 1 || readint <= 0)\
          {\
            gt_error_set(err,"option -algbds: all arguments must be positive "\
                             "numbers");\
            retval = -1;\
          }\
          so->sfxstrategy.COMP = (unsigned long) readint;\
        }

int gt_suffixeratoroptions(Suffixeratoroptions *so,
                        bool doesa,
                        int argc,
                        const char **argv,
                        GtError *err)
{
  int parsed_args, retval = 0;
  GtOPrval rval;
  GtLogger *logger = gt_logger_new(true, GT_LOGGER_DEFLT_PREFIX, stdout);

  gt_error_check(err);
  so->fn2encopt.isdna = false;
  so->fn2encopt.isprotein = false;
  so->inputindex = gt_str_new();
  so->fn2encopt.indexname = gt_str_new();
  so->fn2encopt.smap = gt_str_new();
  so->fn2encopt.sat = gt_str_new();
  so->maxdepth = gt_str_new();
  so->outkystab = false;
  so->outkyssort = false;
  so->optionkysargumentstring = gt_str_new();
  so->fn2encopt.filenametab = gt_str_array_new();
  so->algbounds = gt_str_array_new();
  so->prefixlength = PREFIXLENGTH_AUTOMATIC;
  so->sfxstrategy.ssortmaxdepth.defined = false;
  so->sfxstrategy.ssortmaxdepth.valueunsignedint = MAXDEPTH_AUTOMATIC;
  so->sfxstrategy.streamsuftab = false;
  so->sfxstrategy.hashexceptions = false;
  so->outsuftab = false; /* if !doesa this is not defined */
  so->outlcptab = false;
  so->outbwttab = false;
  so->outbcktab = false;
  rval = parse_options(&parsed_args, doesa, so, argc, argv, err);
  if (rval == GT_OPTION_PARSER_ERROR)
  {
    retval = -1;
  } else
  {
    if (rval == GT_OPTION_PARSER_REQUESTS_EXIT)
    {
      retval = 2;
    } else
    {
      if (rval == GT_OPTION_PARSER_OK)
      {
        if (so->beverbose)
        {
          showoptions(so);
        }
      }
    }
  }
  if (gt_option_is_set(so->optionalgboundsref))
  {
    const char *arg;
    long readint;

    if (gt_str_array_size(so->algbounds) != 3)
    {
      gt_error_set(err,"option -algbds must have exactly 3 arguments");
      retval = -1;
    }
    READMAXBOUND(maxinsertionsort,0);
    READMAXBOUND(maxbltriesort,1);
    READMAXBOUND(maxcountingsort,2);
  } else
  {
    so->sfxstrategy.maxinsertionsort = MAXINSERTIONSORTDEFAULT;
    so->sfxstrategy.maxbltriesort = MAXBLTRIESORTDEFAULT;
    so->sfxstrategy.maxcountingsort = MAXCOUNTINGSORTDEFAULT;
  }
  if (retval != -1)
  {
    if (so->beverbose)
    {
      gt_logger_log_force(logger, "maxinsertionsort=%lu",
                            so->sfxstrategy.maxinsertionsort);
      gt_logger_log_force(logger, "maxbltriesort=%lu",
                            so->sfxstrategy.maxbltriesort);
      gt_logger_log_force(logger, "maxcountingsort=%lu",
                            so->sfxstrategy.maxcountingsort);
    }
    if (so->sfxstrategy.maxinsertionsort > so->sfxstrategy.maxbltriesort)
    {
      gt_error_set(err,"first argument of option -algbds must not be larger "
                       "than second argument");
      retval = -1;
    } else
    {
      if (so->sfxstrategy.maxbltriesort > so->sfxstrategy.maxcountingsort)
      {
        gt_error_set(err,"second argument of option -algbds must not be larger "
                         "than third argument");
        retval = -1;
      }
    }
  }
  gt_logger_delete(logger);
  return retval;
}
