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
#include "core/option.h"
#include "core/versionfunc.h"
#include "match/sarr-def.h"
#include "match/verbose-def.h"
#include "match/stamp.h"
#include "match/esa-seqread.h"
#include "match/esa-map.h"
#include "match/echoseq.h"
#include "tools/gt_sfxmap.h"

#include "match/test-encseq.pr"
#include "match/test-mappedstr.pr"
#include "match/sfx-suftaborder.pr"

typedef struct
{
  bool usestream,
       verbose,
       inputtis,
       inputsuf,
       inputdes,
       inputbwt,
       inputlcp,
       inputbck,
       inputssp;
  unsigned long scantrials,
                multicharcmptrials,
                delspranges;
} Sfxmapoptions;

DECLAREREADFUNCTION(GtUchar);

static void deletethespranges(const Encodedsequence *encseq,
                              unsigned long delspranges)
{
  Specialrangeiterator *sri;
  Sequencerange range;
  Seqpos rangewidth, nextpos = 0, totallength;
  const unsigned long fastawidth = 70UL;

  sri = newspecialrangeiterator(encseq,true);
  printf(">\n");
  while (nextspecialrangeiterator(&range,sri))
  {
    gt_assert(range.rightpos > range.leftpos);
    rangewidth = range.rightpos - range.leftpos;
    if (rangewidth > (Seqpos) delspranges)
    {
      if (range.leftpos == 0)
      {
        nextpos = range.rightpos;
      } else
      {
        if (range.leftpos > nextpos)
        {
          encseq2symbolstring(stdout,
                              encseq,
                              Forwardmode,
                              nextpos,
                              range.leftpos + delspranges - nextpos,
                              fastawidth);
          nextpos = range.rightpos;
        }
      }
    }
  }
  totallength = getencseqtotallength(encseq);
  if (nextpos < totallength-1)
  {
    encseq2symbolstring(stdout,
                        encseq,
                        Forwardmode,
                        nextpos,
                        totallength - nextpos,
                        fastawidth);
  }
  freespecialrangeiterator(&sri);
}

static OPrval parse_options(Sfxmapoptions *sfxmapoptions,
                            int *parsed_args,
                            int argc,
                            const char **argv,
                            GtError *err)
{
  GtOptionParser *op;
  GtOption *optionstream, *optionverbose, *optionscantrials,
         *optionmulticharcmptrials, *optionbck, *optionsuf,
         *optiondes, *optionbwt, *optionlcp, *optiontis, *optionssp,
         *optiondelspranges;
  OPrval oprval;

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
  GtStr *indexname = NULL;
  bool haserr = false;
  Suffixarray suffixarray;
  int parsed_args;
  Verboseinfo *verboseinfo;
  Sfxmapoptions sfxmapoptions;
  unsigned int demand = 0;

  gt_error_check(err);

  switch (parse_options(&sfxmapoptions,&parsed_args, argc, argv, err))
  {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
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
  indexname = gt_str_new_cstr(argv[parsed_args]);
  verboseinfo = newverboseinfo(sfxmapoptions.verbose);
  if (sfxmapoptions.inputtis || sfxmapoptions.delspranges > 0 ||
      sfxmapoptions.inputsuf)
  {
    demand |= SARR_ESQTAB;
  }
  if (sfxmapoptions.inputdes)
  {
    demand |= SARR_DESTAB;
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
                               : mapsuffixarray)(&suffixarray,
                                                 demand,
                                                 indexname,
                                                 verboseinfo,
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
          if (isdnaalphabet(getencseqAlphabet(suffixarray.encseq)) ||
             ((Readmode) readmode) == Forwardmode ||
             ((Readmode) readmode) == Reversemode)
          {
            showverbose(verboseinfo,"testencodedsequence(readmode=%s)",
                                    showreadmode((Readmode) readmode));
            if (testencodedsequence(getencseqfilenametab(suffixarray.encseq),
                                    suffixarray.encseq,
                                    (Readmode) readmode,
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
        showverbose(verboseinfo,"checkspecialrangesfast");
        if (checkspecialrangesfast(suffixarray.encseq) != 0)
        {
          haserr = true;
        }
      }
      if (!haserr && sfxmapoptions.inputtis)
      {
        showverbose(verboseinfo,"checkmarkpos");
        checkmarkpos(suffixarray.encseq);
      }
      if (!haserr && sfxmapoptions.inputtis &&
          suffixarray.readmode == Forwardmode &&
          suffixarray.prefixlength > 0)
      {
        showverbose(verboseinfo,"verifymappedstr");
        if (verifymappedstr(suffixarray.encseq,suffixarray.prefixlength,
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
          ssar = newSequentialsuffixarrayreaderfromfile(indexname,
                                                        SARR_LCPTAB,
                                                        SEQ_scan,
                                                        err);
        } else
        {
          ssar = NULL;
        }
        showverbose(verboseinfo,"checkentiresuftab");
        checkentiresuftab(suffixarray.encseq,
                          suffixarray.readmode,
                          suffixarray.suftab,
                          ssar,
                          false, /* specialsareequal  */
                          false,  /* specialsareequalatdepth0 */
                          0,
                          err);
        if (ssar != NULL)
        {
          freeSequentialsuffixarrayreader(&ssar);
        }
        showverbose(verboseinfo,"okay");
      }
      if (!haserr && sfxmapoptions.inputbwt)
      {
        Seqpos totallength, bwtdifferentconsecutive = 0, idx, longest;

        gt_assert(suffixarray.longest.defined);
        longest = suffixarray.longest.valueseqpos;
        printf("longest=%lu\n",(unsigned long) longest);
        totallength = getencseqtotallength(suffixarray.encseq);
        printf("totallength=%lu\n",(unsigned long) totallength);
        if (!sfxmapoptions.usestream)
        {
          for (idx = (Seqpos) 1; idx<totallength; idx++)
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
        printf("bwtdifferentconsecutive=" FormatSeqpos " (%.4f)\n",
               PRINTSeqposcast(bwtdifferentconsecutive),
               (double) bwtdifferentconsecutive/totallength);
      }
    }
  }
  if (!haserr && sfxmapoptions.inputdes)
  {
    showverbose(verboseinfo,"checkallsequencedescriptions");
    checkallsequencedescriptions(suffixarray.encseq);
  }
  gt_str_delete(indexname);
  freesuffixarray(&suffixarray);
  freeverboseinfo(&verboseinfo);
  return haserr ? -1 : 0;
}
