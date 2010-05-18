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

#include "core/error.h"
#include "core/logger.h"
#include "core/option.h"
#include "core/versionfunc.h"
#include "match/sarr-def.h"
#include "match/stamp.h"
#include "match/esa-seqread.h"
#include "match/esa-map.h"
#include "match/echoseq.h"
#include "match/sfx-suftaborder.h"
#include "match/test-mappedstr.pr"
#include "tools/gt_sfxmap.h"

typedef struct
{
  bool usestream,
       verbose,
       inputtis,
       inputsuf,
       inputdes,
       inputsds,
       inputbwt,
       inputlcp,
       inputbck,
       inputssp;
  unsigned long scantrials,
                multicharcmptrials,
                delspranges;
} Sfxmapoptions;

static void deletethespranges(const GtEncseq *encseq,
                              unsigned long delspranges)
{
  GtSpecialrangeiterator *sri;
  GtRange range;
  unsigned long rangewidth, nextpos = 0, totallength;
  const unsigned long fastawidth = 70UL;

  sri = gt_specialrangeiterator_new(encseq,true);
  printf(">\n");
  while (gt_specialrangeiterator_next(sri,&range))
  {
    gt_assert(range.end > range.start);
    rangewidth = range.end - range.start;
    if (rangewidth > (unsigned long) delspranges)
    {
      if (range.start == 0)
      {
        nextpos = range.end;
      } else
      {
        if (range.start > nextpos)
        {
          gt_encseq2symbolstring(stdout,
                              encseq,
                              GT_READMODE_FORWARD,
                              nextpos,
                              range.start + delspranges - nextpos,
                              fastawidth);
          nextpos = range.end;
        }
      }
    }
  }
  totallength = gt_encseq_total_length(encseq);
  if (nextpos < totallength-1)
  {
    gt_encseq2symbolstring(stdout,
                        encseq,
                        GT_READMODE_FORWARD,
                        nextpos,
                        totallength - nextpos,
                        fastawidth);
  }
  gt_specialrangeiterator_delete(sri);
}

static GtOPrval parse_options(Sfxmapoptions *sfxmapoptions,
                              int *parsed_args,
                              int argc,
                              const char **argv,
                              GtError *err)
{
  GtOptionParser *op;
  GtOption *optionstream, *optionverbose, *optionscantrials,
         *optionmulticharcmptrials, *optionbck, *optionsuf,
         *optiondes, *optionsds, *optionbwt, *optionlcp, *optiontis, *optionssp,
         *optiondelspranges;
  GtOPrval oprval;

  gt_error_check(err);
  op = gt_option_parser_new("[options] indexname",
                            "Map or Stream <indexname> and check consistency.");
  gt_option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  optionstream = gt_option_new_bool("stream","stream the index",
                                 &sfxmapoptions->usestream,false);
  gt_option_parser_add_option(op, optionstream);

  optionscantrials = gt_option_new_ulong("scantrials",
                                         "specify number of scan trials",
                                         &sfxmapoptions->scantrials,0);
  gt_option_parser_add_option(op, optionscantrials);

  optionmulticharcmptrials
    = gt_option_new_ulong("multicharcmptrials",
                          "specify number of multichar cmp trials",
                          &sfxmapoptions->multicharcmptrials,0);
  gt_option_parser_add_option(op, optionmulticharcmptrials);

  optiondelspranges = gt_option_new_ulong("delspranges",
                                          "delete ranges of special values",
                                           &sfxmapoptions->delspranges,
                                           0);
  gt_option_parser_add_option(op, optiondelspranges);

  optiontis = gt_option_new_bool("tis","input the transformed input sequence",
                                 &sfxmapoptions->inputtis,
                                 false);
  gt_option_parser_add_option(op, optiontis);

  optiondes = gt_option_new_bool("des","input the descriptions",
                                 &sfxmapoptions->inputdes,
                                 false);
  gt_option_parser_add_option(op, optiondes);

  optionsds = gt_option_new_bool("sds","input the description end positions",
                                 &sfxmapoptions->inputsds,
                                 false);
  gt_option_parser_add_option(op, optionsds);

  optionsuf = gt_option_new_bool("suf","input the suffix array",
                                 &sfxmapoptions->inputsuf,
                                 false);
  gt_option_parser_add_option(op, optionsuf);

  optionlcp = gt_option_new_bool("lcp","input the lcp-table",
                                 &sfxmapoptions->inputlcp,
                                 false);
  gt_option_parser_add_option(op, optionlcp);

  optionbwt = gt_option_new_bool("bwt",
                                 "input the Burrows-Wheeler Transformation",
                                 &sfxmapoptions->inputbwt,
                                 false);
  gt_option_parser_add_option(op, optionbwt);

  optionbck = gt_option_new_bool("bck","input the bucket table",
                                 &sfxmapoptions->inputbck,
                                 false);
  gt_option_parser_add_option(op, optionbck);

  optionssp = gt_option_new_bool("ssp","input the sequence separator table",
                                 &sfxmapoptions->inputssp,
                                 false);
  gt_option_parser_add_option(op, optionssp);

  optionverbose = gt_option_new_bool("v","be verbose",&sfxmapoptions->verbose,
                                     false);
  gt_option_parser_add_option(op, optionverbose);

  gt_option_parser_set_min_max_args(op, 1U, 2U);
  gt_option_imply(optionlcp,optionsuf);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  gt_option_parser_delete(op);
  return oprval;
}

int gt_sfxmap(int argc, const char **argv, GtError *err)
{
  const char *indexname;
  bool haserr = false;
  Suffixarray suffixarray;
  int parsed_args;
  GtLogger *logger;
  Sfxmapoptions sfxmapoptions;
  unsigned int demand = 0;

  gt_error_check(err);

  switch (parse_options(&sfxmapoptions,&parsed_args, argc, argv, err))
  {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR: return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT: return 0;
  }
  if (argc < 2)
  {
    gt_error_set(err,"missing arguments");
    return -1;
  }
  if (parsed_args != argc - 1)
  {
    gt_error_set(err,"last argument must be indexname");
    return -1;
  }
  indexname = argv[parsed_args];
  logger = gt_logger_new(sfxmapoptions.verbose, GT_LOGGER_DEFLT_PREFIX, stdout);
  if (sfxmapoptions.inputtis || sfxmapoptions.delspranges > 0 ||
      sfxmapoptions.inputsuf)
  {
    demand |= SARR_ESQTAB;
  }
  if (sfxmapoptions.inputdes)
  {
    demand |= SARR_DESTAB;
  }
  if (sfxmapoptions.inputsds)
  {
    demand |= SARR_SDSTAB;
  }
  if (sfxmapoptions.inputsuf)
  {
    demand |= SARR_SUFTAB;
  }
  if (sfxmapoptions.inputlcp)
  {
    demand |= SARR_LCPTAB;
  }
  if (sfxmapoptions.inputbwt)
  {
    demand |= SARR_BWTTAB;
  }
  if (sfxmapoptions.inputbck)
  {
    demand |= SARR_BCKTAB;
  }
  if (sfxmapoptions.inputssp)
  {
    demand |= SARR_SSPTAB;
  }
  if ((sfxmapoptions.usestream ? streamsuffixarray
                               : gt_mapsuffixarray)(&suffixarray,
                                                   demand,
                                                   indexname,
                                                   logger,
                                                   err) != 0)
  {
    haserr = true;
  }
  if (!haserr && suffixarray.encseq != NULL)
  {
    if (sfxmapoptions.delspranges > 0)
    {
      deletethespranges(suffixarray.encseq,sfxmapoptions.delspranges);
    } else
    {
      if (!haserr && sfxmapoptions.inputtis)
      {
        int readmode;

        for (readmode = 0; readmode < 4; readmode++)
        {
          if (gt_alphabet_is_dna(
            gt_encseq_alphabet(suffixarray.encseq)) ||
               ((GtReadmode) readmode) == GT_READMODE_FORWARD ||
               ((GtReadmode) readmode) == GT_READMODE_REVERSE)
          {
            gt_logger_log(logger, "testencseq(readmode=%s)",
                                   gt_readmode_show((GtReadmode) readmode));
            if (gt_encseq_check_consistency(suffixarray.encseq,
                               gt_encseq_filenames(suffixarray.encseq),
                               (GtReadmode) readmode,
                               sfxmapoptions.scantrials,
                               sfxmapoptions.multicharcmptrials,
                               err) != 0)
            {
              haserr = true;
              break;
            }
          }
        }
      }
      if (!haserr && sfxmapoptions.inputtis)
      {
        gt_logger_log(logger, "checkspecialrangesfast");
        if (gt_encseq_check_specialranges(suffixarray.encseq) != 0)
        {
          haserr = true;
        }
      }
      if (!haserr && sfxmapoptions.inputtis)
      {
        gt_logger_log(logger, "gt_encseq_check_markpos");
        gt_encseq_check_markpos(suffixarray.encseq);
      }
      if (!haserr && sfxmapoptions.inputtis &&
          suffixarray.readmode == GT_READMODE_FORWARD &&
          suffixarray.prefixlength > 0)
      {
        gt_logger_log(logger, "verifymappedstr");
        if (gt_verifymappedstr(suffixarray.encseq,suffixarray.prefixlength,
                            err) != 0)
        {
          haserr = true;
        }
      }
      if (!haserr && sfxmapoptions.inputsuf && !sfxmapoptions.usestream)
      {
        Sequentialsuffixarrayreader *ssar;

        if (sfxmapoptions.inputlcp)
        {
          ssar = gt_newSequentialsuffixarrayreaderfromfile(indexname,
                                                        SARR_LCPTAB |
                                                        SARR_ESQTAB,
                                                        SEQ_scan,
                                                        err);
        } else
        {
          ssar = NULL;
        }
        gt_logger_log(logger, "checkentiresuftab");
        gt_checkentiresuftab(__FILE__,
                          __LINE__,
                          suffixarray.encseq,
                          suffixarray.readmode,
                          suffixarray.suftab,
                          gt_encseq_total_length(suffixarray.encseq)+1,
                          ssar,
                          false, /* specialsareequal  */
                          false,  /* specialsareequalatdepth0 */
                          0,
                          err);
        if (ssar != NULL)
        {
          gt_freeSequentialsuffixarrayreader(&ssar);
        }
        gt_logger_log(logger, "okay");
      }
      if (!haserr && sfxmapoptions.inputbwt)
      {
        unsigned long totallength, bwtdifferentconsecutive = 0, idx, longest;

        gt_assert(suffixarray.longest.defined);
        longest = suffixarray.longest.valueunsignedlong;
        printf("longest=%lu\n",(unsigned long) longest);
        totallength = gt_encseq_total_length(suffixarray.encseq);
        printf("totallength=%lu\n",(unsigned long) totallength);
        if (!sfxmapoptions.usestream)
        {
          for (idx = (unsigned long) 1; idx<totallength; idx++)
          {
            if (suffixarray.bwttab[idx-1] != suffixarray.bwttab[idx] ||
                ISSPECIAL(suffixarray.bwttab[idx]))
            {
              bwtdifferentconsecutive++;
            }
          }
        } else
        {
          GtUchar prevcc;

          if (readnextGtUcharfromstream(&prevcc,&suffixarray.bwttabstream) == 1)
          {
            GtUchar cc;
            while (readnextGtUcharfromstream(&cc,&suffixarray.bwttabstream)
                   == 1)
            {
              if (prevcc != cc || ISSPECIAL(cc))
              {
                bwtdifferentconsecutive++;
              }
            }
          }
        }
        printf("bwtdifferentconsecutive=%lu (%.4f)\n",
               bwtdifferentconsecutive,
               (double) bwtdifferentconsecutive/totallength);
      }
    }
  }
  if (!haserr && sfxmapoptions.inputdes)
  {
    gt_logger_log(logger, "checkallsequencedescriptions");
    gt_encseq_check_descriptions(suffixarray.encseq);
  }
  gt_freesuffixarray(&suffixarray);
  gt_logger_delete(logger);
  return haserr ? -1 : 0;
}
