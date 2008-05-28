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
#include "libgtcore/error.h"
#include "libgtcore/str.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/stamp.h"
#include "libgtmatch/enum-patt-def.h"
#include "libgtmatch/esa-mmsearch-def.h"
#include "libgtmatch/qgram2code.h"
#include "libgtmatch/bcktab.h"
#include "libgtmatch/spacedef.h"
#include "libgtmatch/cutendpfx.h"

#include "libgtmatch/esa-map.pr"
#include "libgtmatch/sfx-cmpsuf.pr"

#include "tools/gt_patternmatch.h"

typedef struct
{
  unsigned long minpatternlen, maxpatternlen, numofsamples;
  bool showpatt, usebcktab, immediate;
  Str *indexname;
} Pmatchoptions;

static void comparemmsis(const MMsearchiterator *mmsi1,
                         const MMsearchiterator *mmsi2)
{
  if (isemptymmsearchiterator(mmsi1))
  {
    if (!isemptymmsearchiterator(mmsi2))
    {
      fprintf(stderr,"mmsi1 is empty but mmsi2 not\n");
      exit(EXIT_FAILURE); /* programming error */
    }
  } else
  {
    if (isemptymmsearchiterator(mmsi2))
    {
      fprintf(stderr,"mmsi2 is empty but mmsi1 not\n");
      exit(EXIT_FAILURE); /* programming error */
    }
    if (!identicalmmsearchiterators(mmsi1,mmsi2))
    {
      fprintf(stderr,"mmsi1 and mmsi2 are different\n");
      exit(EXIT_FAILURE); /* programming error */
    }
  }
}

#define UNDEFREFSTART totallength

static int callpatternmatcher(const Pmatchoptions *pmopt, Error *err)
{
  Suffixarray suffixarray;
  Seqpos totallength;
  bool haserr = false;
  const Uchar *pptr;
  unsigned long patternlen;
  unsigned int demand = SARR_SUFTAB | SARR_ESQTAB;

  if (pmopt->usebcktab)
  {
    demand |= SARR_BCKTAB;
  }
  if (mapsuffixarray(&suffixarray,
                     &totallength,
                     demand,
                     pmopt->indexname,
                     NULL,
                     err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    unsigned long trial;
    Seqpos dbstart;
    Enumpatterniterator *epi;
    unsigned int firstspecial;
    MMsearchiterator *mmsibck, *mmsiimm;
    Bucketspecification bucketspec;
    Bucketenumerator *bucketenumerator;
    Lcpinterval itv;
    Seqpos refstart;
    Encodedsequencescanstate *esr1, *esr2;
    int retval;
    Seqpos idx, maxlcp;
    Codetype code = 0;
    const Codetype **multimappower;

    if (pmopt->usebcktab)
    {
      multimappower = bcktab_multimappower(suffixarray.bcktab);
    } else
    {
      multimappower = NULL;
    }
    epi = newenumpatterniterator(pmopt->minpatternlen,
                                 pmopt->maxpatternlen,
                                 suffixarray.encseq,
                                 getnumofcharsAlphabet(suffixarray.alpha),
                                 err);
    esr1 = newEncodedsequencescanstate();
    esr2 = newEncodedsequencescanstate();
    for (trial = 0; trial < pmopt->numofsamples; trial++)
    {
      pptr = nextEnumpatterniterator(&patternlen,epi);
      if (pmopt->showpatt)
      {
        showsymbolstring(suffixarray.alpha,pptr,patternlen);
        printf("\n");
      }
      if (pmopt->usebcktab)
      {
        if (patternlen < (unsigned long) suffixarray.prefixlength)
        {
          mmsibck = NULL;
          bucketenumerator
            = newbucketenumerator(suffixarray.bcktab,
                                  suffixarray.prefixlength,
                                  pptr,
                                  (unsigned int) patternlen);
          refstart = UNDEFREFSTART;
          while (nextbucketenumerator(&itv,bucketenumerator))
          {
            if (refstart == UNDEFREFSTART)
            {
              refstart = suffixarray.suftab[itv.left];
            } else
            {
              for (idx=itv.left; idx<=itv.right; idx++)
              {
                retval = comparetwosuffixes(suffixarray.encseq,
                                            suffixarray.readmode,
                                            &maxlcp,
                                            false,
                                            false,
                                            (Seqpos) patternlen,
                                            refstart,
                                            suffixarray.suftab[idx],
                                            esr1,
                                            esr2);
                assert(retval == 0 && maxlcp == (Seqpos) patternlen);
              }
            }
          }
          freebucketenumerator(&bucketenumerator);
        } else
        {
          firstspecial = qgram2code(&code,
                                    multimappower,
                                    suffixarray.prefixlength,
                                    pptr);
          assert(firstspecial == suffixarray.prefixlength);
          calcbucketboundaries(&bucketspec,
                               suffixarray.bcktab,
                               code);
          if (bucketspec.nonspecialsinbucket == 0)
          {
            mmsibck = NULL;
          } else
          {
            mmsibck
              = newmmsearchiterator(suffixarray.encseq,
                                    suffixarray.suftab,
                                    bucketspec.left,
                                    bucketspec.left +
                                      bucketspec.nonspecialsinbucket-1,
                                    (Seqpos) suffixarray.prefixlength,
                                    suffixarray.readmode,
                                    pptr,
                                    patternlen);
          }
        }
      }
      if (pmopt->immediate)
      {
        mmsiimm = newmmsearchiterator(suffixarray.encseq,
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
        while (nextmmsearchiterator(&dbstart,mmsibck))
        {
          /* Nothing */;
        }
        freemmsearchiterator(&mmsibck);
      }
      if (pmopt->immediate)
      {
        while (nextmmsearchiterator(&dbstart,mmsiimm))
        {
          /* Nothing */;
        }
        freemmsearchiterator(&mmsiimm);
      }
    }
    freeEncodedsequencescanstate(&esr1);
    freeEncodedsequencescanstate(&esr2);
    if (pmopt->showpatt)
    {
      showPatterndistribution(epi);
    }
    freeEnumpatterniterator(&epi);
  }
  freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}

static OPrval parse_options(Pmatchoptions *pmopt,
                            int *parsed_args,
                            int argc, const char **argv, Error *err)
{
  OptionParser *op;
  Option *option, *optionimm, *optionbck;
  OPrval oprval;

  error_check(err);
  op = option_parser_new("[options] -ii indexname",
                         "Perform pattern matches.");
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");

  option = option_new_ulong("minpl","Specify minimum length of pattern",
                           &pmopt->minpatternlen,
                           (unsigned long) 20);
  option_parser_add_option(op, option);
  option = option_new_ulong("maxpl","Specify maximum length of pattern",
                            &pmopt->maxpatternlen,
                            (unsigned long) 30);
  option_parser_add_option(op, option);

  option = option_new_ulong("samples","Specify number of samples",
                            &pmopt->numofsamples,
                           (unsigned long) 100000);
  option_parser_add_option(op, option);

  option = option_new_bool("s","Show generated pattern",
                            &pmopt->showpatt,
                            false);
  option_parser_add_option(op, option);

  optionbck = option_new_bool("bck","Use the bucket boundaries",
                              &pmopt->usebcktab,
                              false);
  option_parser_add_option(op, optionbck);

  optionimm = option_new_bool("imm","Start with offset 0",
                              &pmopt->immediate,
                              false);
  option_parser_add_option(op, optionimm);

  option = option_new_string("ii",
                             "Specify input index",
                             pmopt->indexname, NULL);
  option_parser_add_option(op, option);
  option_is_mandatory(option);

  oprval = option_parser_parse(op, parsed_args, argc, argv,
                               versionfunc, err);
  option_parser_delete(op);

  return oprval;
}

int gt_patternmatch(int argc, const char **argv, Error *err)
{
  bool haserr = false;
  int parsed_args;
  Pmatchoptions pmopt;
  OPrval oprval;

  error_check(err);

  pmopt.indexname = str_new();
  oprval = parse_options(&pmopt,&parsed_args, argc, argv, err);
  if (oprval == OPTIONPARSER_OK)
  {
    assert(parsed_args == argc);
    if (callpatternmatcher(&pmopt,err) != 0)
    {
      haserr = true;
    }
  }
  str_delete(pmopt.indexname);
  if (oprval == OPTIONPARSER_REQUESTS_EXIT)
  {
    return 0;
  }
  if (oprval == OPTIONPARSER_ERROR)
  {
    return -1;
  }
  return haserr ? -1 : 0;
}
