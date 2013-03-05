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
#include "core/encseq.h"
#include "core/error.h"
#include "core/option_api.h"
#include "core/str.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "match/cutendpfx.h"
#include "match/enum-patt.h"
#include "match/esa-map.h"
#include "match/esa-mmsearch.h"
#include "match/qgram2code.h"
#include "match/sarr-def.h"
#include "tools/gt_patternmatch.h"

typedef struct
{
  unsigned long minpatternlen, maxpatternlen, numofsamples;
  bool showpatt, usebcktab, immediate;
  GtStr *indexname;
} Pmatchoptions;

static void comparemmsis(const GtMMsearchiterator *mmsi1,
                         const GtMMsearchiterator *mmsi2)
{
  if (gt_mmsearchiterator_isempty(mmsi1))
  {
    if (!gt_mmsearchiterator_isempty(mmsi2))
    {
      fprintf(stderr,"mmsi1 is empty but mmsi2 not\n");
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  } else
  {
    if (gt_mmsearchiterator_isempty(mmsi2))
    {
      fprintf(stderr,"mmsi2 is empty but mmsi1 not\n");
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if (!gt_mmsearchiterator_identical(mmsi1,mmsi2))
    {
      fprintf(stderr,"mmsi1 and mmsi2 are different\n");
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

#define UNDEFREFSTART totallength

static int callpatternmatcher(const Pmatchoptions *pmopt, GtError *err)
{
  Suffixarray suffixarray;
  unsigned long totallength = 0;
  bool haserr = false;
  const GtUchar *pptr;
  unsigned long patternlen;
  unsigned int demand = SARR_SUFTAB | SARR_ESQTAB;

  if (pmopt->usebcktab)
  {
    demand |= SARR_BCKTAB;
  }
  if (gt_mapsuffixarray(&suffixarray,
                     demand,
                     gt_str_get(pmopt->indexname),
                     NULL,
                     err) != 0)
  {
    haserr = true;
  } else
  {
    totallength = gt_encseq_total_length(suffixarray.encseq);
  }
  if (!haserr)
  {
    unsigned long trial;
    unsigned long dbstart;
    Enumpatterniterator *epi;
    GT_UNUSED unsigned int firstspecial;
    GtMMsearchiterator *mmsibck = NULL, *mmsiimm = NULL;
    GtBucketspecification bucketspec;
    Bucketenumerator *bucketenumerator;
    Lcpinterval itv;
    unsigned long refstart;
    GtEncseqReader *esr1, *esr2;
    GT_UNUSED int retval;
    unsigned long idx, maxlcp;
    GtCodetype code = 0;
    const GtCodetype **multimappower;
    const GtAlphabet *alpha;

    if (pmopt->usebcktab)
    {
      multimappower = gt_bcktab_multimappower(suffixarray.bcktab);
    } else
    {
      multimappower = NULL;
    }
    epi = gt_newenumpatterniterator(pmopt->minpatternlen,
                                    pmopt->maxpatternlen,
                                    suffixarray.encseq,
                                    err);
    esr1 = gt_encseq_create_reader_with_readmode(suffixarray.encseq,
                                                 suffixarray.readmode, 0);
    esr2 = gt_encseq_create_reader_with_readmode(suffixarray.encseq,
                                                 suffixarray.readmode, 0);
    alpha = gt_encseq_alphabet(suffixarray.encseq);
    for (trial = 0; trial < pmopt->numofsamples; trial++)
    {
      pptr = gt_nextEnumpatterniterator(&patternlen,epi);
      if (pmopt->showpatt)
      {
        gt_alphabet_decode_seq_to_fp(alpha,stdout,pptr,patternlen);
        printf("\n");
      }
      if (pmopt->usebcktab)
      {
        if (patternlen < (unsigned long) suffixarray.prefixlength)
        {
          mmsibck = NULL;
          bucketenumerator
            = gt_newbucketenumerator(suffixarray.bcktab,
                                  suffixarray.prefixlength,
                                  pptr,
                                  (unsigned int) patternlen);
          refstart = UNDEFREFSTART;
          while (gt_nextbucketenumerator(&itv,bucketenumerator))
          {
            if (refstart == UNDEFREFSTART)
            {
              refstart = ESASUFFIXPTRGET(suffixarray.suftab,itv.left);
            } else
            {
              for (idx=itv.left; idx<=itv.right; idx++)
              {
                retval = gt_encseq_check_comparetwosuffixes(
                                        suffixarray.encseq,
                                        suffixarray.readmode,
                                        &maxlcp,
                                        false,
                                        false,
                                        patternlen,
                                        refstart,
                                        ESASUFFIXPTRGET(suffixarray.suftab,idx),
                                        esr1,
                                        esr2);
                gt_assert(retval == 0 && maxlcp == patternlen);
              }
            }
          }
          gt_freebucketenumerator(bucketenumerator);
        } else
        {
          firstspecial = qgram2code(&code,
                                    multimappower,
                                    suffixarray.prefixlength,
                                    pptr);
          gt_assert(firstspecial == suffixarray.prefixlength);
          gt_bcktab_calcboundaries(&bucketspec,
                                   suffixarray.bcktab,
                                   code);
          if (bucketspec.nonspecialsinbucket == 0)
          {
            mmsibck = NULL;
          } else
          {
            mmsibck
              = gt_mmsearchiterator_new_complete_olain(
                                       suffixarray.encseq,
                                       suffixarray.suftab,
                                       bucketspec.left,
                                       bucketspec.left +
                                         bucketspec.nonspecialsinbucket-1,
                                       (unsigned long) suffixarray.prefixlength,
                                       suffixarray.readmode,
                                       pptr,
                                       patternlen);
          }
        }
      }
      if (pmopt->immediate)
      {
        mmsiimm = gt_mmsearchiterator_new_complete_olain(
                                            suffixarray.encseq,
                                            suffixarray.suftab,
                                            0,  /* leftbound */
                                            totallength, /* rightbound */
                                            0, /* offset */
                                            suffixarray.readmode,
                                            pptr,
                                            patternlen);
      }
      if (pmopt->immediate && pmopt->usebcktab)
      {
        comparemmsis(mmsibck,mmsiimm);
      }
      if (pmopt->usebcktab && mmsibck != NULL)
      {
        while (gt_mmsearchiterator_next(&dbstart,mmsibck))
        {
          /* Nothing */;
        }
        gt_mmsearchiterator_delete(mmsibck);
        mmsibck = NULL;
      }
      if (pmopt->immediate)
      {
        while (gt_mmsearchiterator_next(&dbstart,mmsiimm))
        {
          /* Nothing */;
        }
        gt_mmsearchiterator_delete(mmsiimm);
        mmsiimm = NULL;
      }
    }
    gt_encseq_reader_delete(esr1);
    gt_encseq_reader_delete(esr2);
    if (pmopt->showpatt)
    {
      gt_showPatterndistribution(epi);
    }
    gt_freeEnumpatterniterator(epi);
  }
  gt_freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}

static GtOPrval parse_options(Pmatchoptions *pmopt,
                              int *parsed_args,
                              int argc, const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOption *option, *optionimm, *optionbck;
  GtOPrval oprval;

  gt_error_check(err);
  op = gt_option_parser_new("[options] -ii indexname",
                         "Perform pattern matches.");
  gt_option_parser_set_mail_address(op, "<kurtz@zbh.uni-hamburg.de>");

  option = gt_option_new_ulong("minpl","Specify minimum length of pattern",
                           &pmopt->minpatternlen,
                           (unsigned long) 20);
  gt_option_parser_add_option(op, option);
  option = gt_option_new_ulong("maxpl","Specify maximum length of pattern",
                            &pmopt->maxpatternlen,
                            (unsigned long) 30);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_ulong("samples","Specify number of samples",
                            &pmopt->numofsamples,
                           (unsigned long) 100000);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("s","Show generated pattern",
                            &pmopt->showpatt,
                            false);
  gt_option_parser_add_option(op, option);

  optionbck = gt_option_new_bool("bck","Use the bucket boundaries",
                              &pmopt->usebcktab,
                              false);
  gt_option_parser_add_option(op, optionbck);

  optionimm = gt_option_new_bool("imm","Start with offset 0",
                              &pmopt->immediate,
                              false);
  gt_option_parser_add_option(op, optionimm);

  option = gt_option_new_string("ii",
                             "Specify input index",
                             pmopt->indexname, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  oprval = gt_option_parser_parse(op, parsed_args, argc, argv,
                               gt_versionfunc, err);
  gt_option_parser_delete(op);

  return oprval;
}

int gt_patternmatch(int argc, const char **argv, GtError *err)
{
  bool haserr = false;
  int parsed_args;
  Pmatchoptions pmopt;
  GtOPrval oprval;

  gt_error_check(err);

  pmopt.indexname = gt_str_new();
  oprval = parse_options(&pmopt,&parsed_args, argc, argv, err);
  if (oprval == GT_OPTION_PARSER_OK)
  {
    gt_assert(parsed_args == argc);
    if (callpatternmatcher(&pmopt,err) != 0)
    {
      haserr = true;
    }
  }
  gt_str_delete(pmopt.indexname);
  if (oprval == GT_OPTION_PARSER_REQUESTS_EXIT)
  {
    return 0;
  }
  if (oprval == GT_OPTION_PARSER_ERROR)
  {
    return -1;
  }
  return haserr ? -1 : 0;
}
