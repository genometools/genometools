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

#include "core/arraydef.h"
#include "core/disc_distri.h"
#include "core/fa.h"
#include "core/minmax.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "esa-fileend.h"
#include "sfx-lcpvalues.h"
#include "turnwheels.h"

GT_DECLAREARRAYSTRUCT(Largelcpvalue);

typedef struct
{
  FILE *outfplcptab,
       *outfpllvtab;
  GtArrayLargelcpvalue largelcpvalues;
  unsigned long maxbranchdepth,
                totalnumoflargelcpvalues,
                countoutputlcpvalues;
  void *reservoir;
  size_t sizereservoir;
  uint8_t *smalllcpvalues; /* pointer into reservoir */
} Lcpoutput2file;

typedef struct
{
  GtLcpvalues tableoflcpvalues;
  Lcpoutput2file *lcp2file;
  double lcptabsum;
  GtDiscDistri *distlcpvalues;
} Lcpsubtab;

typedef struct
{
  bool defined;
  GtCodetype code;
  unsigned int prefixindex;
#undef SKDEBUG
#ifdef SKDEBUG
  unsigned long startpos;
#endif
} Suffixwithcode;

struct Outlcpinfo
{
  Turningwheel *turnwheel;
  unsigned long numsuffixes2output;
  unsigned int minchanged;
  size_t sizeofinfo;
  Suffixwithcode previoussuffix;
  bool previousbucketwasempty;
  Lcpsubtab lcpsubtab;
};

/* Now some functions related to the computation of lcp values follows */

static unsigned long computelocallcpvalue(const Suffixwithcode *previoussuffix,
                                          const Suffixwithcode *currentsuffix,
                                          unsigned int minchanged)
{
  unsigned int lcpvalue;

  if (previoussuffix->code == currentsuffix->code)
  {
    lcpvalue = MIN(previoussuffix->prefixindex,
                   currentsuffix->prefixindex);
  } else
  {
    gt_assert(previoussuffix->code < currentsuffix->code);
    lcpvalue = MIN(minchanged,MIN(previoussuffix->prefixindex,
                                  currentsuffix->prefixindex));
  }
  return (unsigned long) lcpvalue;
}

static void outsmalllcpvalues(Lcpoutput2file *lcp2file,
                              unsigned long numoflcps)
{
  gt_assert (lcp2file != NULL);
  lcp2file->countoutputlcpvalues += numoflcps;
  gt_assert(lcp2file->outfplcptab != NULL);
  gt_xfwrite(lcp2file->smalllcpvalues,
             sizeof (*lcp2file->smalllcpvalues),
             (size_t) numoflcps,
             lcp2file->outfplcptab);
}

static unsigned int bucketends(Lcpsubtab *lcpsubtab,
                               Suffixwithcode *previoussuffix,
                               GT_UNUSED unsigned long firstspecialsuffix,
                               unsigned int minchanged,
                               unsigned long specialsinbucket,
                               GtCodetype code,
                               const GtBcktab *bcktab)
{
  unsigned long lcpvalue;
  unsigned int maxprefixindex, minprefixindex;
  Suffixwithcode firstspecialsuffixwithcode;

  /*
     there is at least one element in the bucket. if there is more than
     one element in the bucket, then we insert them using the
     information from the bcktab
  */
  if (specialsinbucket > 1UL)
  {
    if (lcpsubtab->lcp2file != NULL)
    {
      maxprefixindex = gt_bcktab_pfxidx2lcpvalues_uint8(
                          &minprefixindex,
                          lcpsubtab->lcp2file->smalllcpvalues,
                          specialsinbucket,
                          bcktab,
                          code);
      if (lcpsubtab->lcp2file->maxbranchdepth < (unsigned long) maxprefixindex)
      {
        lcpsubtab->lcp2file->maxbranchdepth = (unsigned long) maxprefixindex;
      }
    } else
    {
      maxprefixindex = gt_bcktab_pfxidx2lcpvalues_ulong(
                          &minprefixindex,
                          lcpsubtab->tableoflcpvalues.bucketoflcpvalues +
                          lcpsubtab->tableoflcpvalues.lcptaboffset,
                          specialsinbucket,
                          bcktab,
                          code);
#ifndef NDEBUG
      {
        unsigned long idx;
        for (idx=0; idx<specialsinbucket; idx++)
        {
          GT_SETIBIT(lcpsubtab->tableoflcpvalues.isset,
                     lcpsubtab->tableoflcpvalues.lcptaboffset+idx);
        }
      }
#endif
    }
  } else
  {
    minprefixindex = maxprefixindex = gt_bcktab_singletonmaxprefixindex(bcktab,
                                                                        code);
  }
  firstspecialsuffixwithcode.code = code;
  firstspecialsuffixwithcode.prefixindex = maxprefixindex;
#ifdef SKDEBUG
  firstspecialsuffixwithcode.startpos = firstspecialsuffix;
  /*
  gt_bcktab_consistencyofsuffix(__LINE__,
                                encseq,readmode,bcktab,numofchars,
                                &firstspecialsuffixwithcode);
  */
#endif
  lcpvalue = computelocallcpvalue(previoussuffix,
                                  &firstspecialsuffixwithcode,
                                  minchanged);
  if (lcpsubtab->lcp2file != NULL)
  {
    if (lcpsubtab->lcp2file->maxbranchdepth < lcpvalue)
    {
      lcpsubtab->lcp2file->maxbranchdepth = lcpvalue;
    }
    lcpsubtab->lcp2file->smalllcpvalues[0] = (uint8_t) lcpvalue;
  } else
  {
    lcpsubtab->tableoflcpvalues.bucketoflcpvalues
               [lcpsubtab->tableoflcpvalues.lcptaboffset] = lcpvalue;
#ifndef NDEBUG
    GT_SETIBIT(lcpsubtab->tableoflcpvalues.isset,
               lcpsubtab->tableoflcpvalues.lcptaboffset);
#endif
  }
  return minprefixindex;
}

Outlcpinfo *gt_Outlcpinfo_new(const char *indexname,
                              unsigned int numofchars,
                              unsigned int prefixlength,
                              bool withdistribution,
                              GtError *err)
{
  bool haserr = false;
  Outlcpinfo *outlcpinfo;

  outlcpinfo = gt_malloc(sizeof (*outlcpinfo));
  outlcpinfo->sizeofinfo = sizeof (*outlcpinfo);
  outlcpinfo->lcpsubtab.lcptabsum = 0.0;
  if (withdistribution)
  {
    outlcpinfo->lcpsubtab.distlcpvalues = gt_disc_distri_new();
  } else
  {
    outlcpinfo->lcpsubtab.distlcpvalues = NULL;
  }
  if (indexname == NULL)
  {
    outlcpinfo->lcpsubtab.lcp2file = NULL;
  } else
  {
    outlcpinfo->lcpsubtab.lcp2file
      = gt_malloc(sizeof (*outlcpinfo->lcpsubtab.lcp2file));
    outlcpinfo->sizeofinfo += sizeof (*outlcpinfo->lcpsubtab.lcp2file);
    outlcpinfo->lcpsubtab.lcp2file->countoutputlcpvalues = 0;
    outlcpinfo->lcpsubtab.lcp2file->maxbranchdepth = 0;
    outlcpinfo->lcpsubtab.lcp2file->totalnumoflargelcpvalues = 0;
    outlcpinfo->lcpsubtab.lcp2file->reservoir = NULL;
    outlcpinfo->lcpsubtab.lcp2file->sizereservoir = 0;
    outlcpinfo->lcpsubtab.lcp2file->smalllcpvalues = NULL;
    GT_INITARRAY(&outlcpinfo->lcpsubtab.lcp2file->largelcpvalues,
                 Largelcpvalue);
    outlcpinfo->lcpsubtab.lcp2file->outfplcptab
      = gt_fa_fopen_with_suffix(indexname,LCPTABSUFFIX,"wb",err);
    if (outlcpinfo->lcpsubtab.lcp2file->outfplcptab == NULL)
    {
      haserr = true;
    }
    if (!haserr)
    {
      outlcpinfo->lcpsubtab.lcp2file->outfpllvtab
        = gt_fa_fopen_with_suffix(indexname,LARGELCPTABSUFFIX,"wb",err);
      if (outlcpinfo->lcpsubtab.lcp2file->outfpllvtab == NULL)
      {
        haserr = true;
      }
    }
  }
  outlcpinfo->numsuffixes2output = 0;
  outlcpinfo->minchanged = 0;
  if (!haserr && prefixlength > 0)
  {
    outlcpinfo->turnwheel = gt_turningwheel_new(prefixlength,numofchars);
    outlcpinfo->sizeofinfo += gt_turningwheel_size();
  } else
  {
    outlcpinfo->turnwheel = NULL;
  }
#ifdef SKDEBUG
  outlcpinfo->previoussuffix.startpos = 0;
#endif
  outlcpinfo->previoussuffix.code = 0;
  outlcpinfo->previoussuffix.prefixindex = 0;
  outlcpinfo->previoussuffix.defined = false;
  outlcpinfo->previousbucketwasempty = false;
  if (haserr)
  {
    gt_free(outlcpinfo);
    return NULL;
  }
  return outlcpinfo;
}

size_t gt_Outlcpinfo_size(const Outlcpinfo *outlcpinfo)
{
  gt_assert(outlcpinfo != NULL);
  return outlcpinfo->sizeofinfo;
}

void gt_Outlcpinfo_reinit(Outlcpinfo *outlcpinfo,
                          unsigned int numofchars,
                          unsigned int prefixlength,
                          unsigned long numoflcpvalues)
{
  if (outlcpinfo != NULL)
  {
    size_t sizeofbucketoflcpvalues;
    if (prefixlength > 0)
    {
      outlcpinfo->turnwheel = gt_turningwheel_new(prefixlength,numofchars);
      outlcpinfo->sizeofinfo += gt_turningwheel_size();
    } else
    {
      outlcpinfo->turnwheel = NULL;
    }
    sizeofbucketoflcpvalues = sizeof (*outlcpinfo->lcpsubtab.tableoflcpvalues.
                                      bucketoflcpvalues) * numoflcpvalues;
    outlcpinfo->lcpsubtab.tableoflcpvalues.bucketoflcpvalues
      = gt_malloc(sizeofbucketoflcpvalues);
    outlcpinfo->sizeofinfo += sizeofbucketoflcpvalues;
#ifndef NDEBUG
    GT_INITBITTAB(outlcpinfo->lcpsubtab.tableoflcpvalues.isset,
                  numoflcpvalues);
#endif
    outlcpinfo->sizeofinfo += GT_NUMOFINTSFORBITS(numoflcpvalues) *
                              sizeof (GtBitsequence);
    outlcpinfo->lcpsubtab.tableoflcpvalues.numoflargelcpvalues = 0;
    outlcpinfo->lcpsubtab.tableoflcpvalues.numofentries = numoflcpvalues;
    outlcpinfo->lcpsubtab.tableoflcpvalues.lcptaboffset = 0;
  }
}

static void outlcpvalues(Lcpsubtab *lcpsubtab,
                         unsigned long bucketleft,
                         unsigned long bucketright,
                         unsigned long posoffset)
{
  unsigned long idx, lcpvalue;
  Largelcpvalue *largelcpvalueptr;

  if (lcpsubtab->lcp2file == NULL)
  {
    return;
  }
  lcpsubtab->lcp2file->largelcpvalues.nextfreeLargelcpvalue = 0;
  if (lcpsubtab->tableoflcpvalues.numoflargelcpvalues > 0 &&
      lcpsubtab->tableoflcpvalues.numoflargelcpvalues >=
      lcpsubtab->lcp2file->largelcpvalues.allocatedLargelcpvalue)
  {
    lcpsubtab->lcp2file->largelcpvalues.spaceLargelcpvalue
      = gt_realloc(lcpsubtab->lcp2file->largelcpvalues.spaceLargelcpvalue,
                   sizeof (*lcpsubtab->lcp2file->largelcpvalues.
                           spaceLargelcpvalue) *
                   lcpsubtab->tableoflcpvalues.numoflargelcpvalues);
    lcpsubtab->lcp2file->largelcpvalues.allocatedLargelcpvalue
      = lcpsubtab->tableoflcpvalues.numoflargelcpvalues;
  }
  for (idx=bucketleft; idx<=bucketright; idx++)
  {
    lcpvalue = lcpsubtab_getvalue(&lcpsubtab->tableoflcpvalues,0,idx);
    if (lcpsubtab->lcp2file->maxbranchdepth < lcpvalue)
    {
      lcpsubtab->lcp2file->maxbranchdepth = lcpvalue;
    }
    if (lcpvalue < (unsigned long) LCPOVERFLOW)
    {
      lcpsubtab->lcp2file->smalllcpvalues[idx-bucketleft]
        = (uint8_t) lcpvalue;
    } else
    {
      gt_assert(lcpsubtab->lcp2file->largelcpvalues.nextfreeLargelcpvalue
                < lcpsubtab->lcp2file->largelcpvalues.
                                             allocatedLargelcpvalue);
      largelcpvalueptr
        = lcpsubtab->lcp2file->largelcpvalues.spaceLargelcpvalue +
          lcpsubtab->lcp2file->largelcpvalues.nextfreeLargelcpvalue++;
      largelcpvalueptr->position = posoffset+idx;
      largelcpvalueptr->value = lcpvalue;
      lcpsubtab->lcp2file->smalllcpvalues[idx-bucketleft] = LCPOVERFLOW;
    }
    lcpsubtab->lcptabsum += (double) lcpvalue;
    if (lcpsubtab->distlcpvalues != NULL)
    {
      gt_disc_distri_add(lcpsubtab->distlcpvalues, lcpvalue);
    }
  }
  outsmalllcpvalues(lcpsubtab->lcp2file,
                    (unsigned long) (bucketright - bucketleft + 1));
  if (lcpsubtab->lcp2file->largelcpvalues.nextfreeLargelcpvalue > 0)
  {
    lcpsubtab->lcp2file->totalnumoflargelcpvalues
      += lcpsubtab->lcp2file->largelcpvalues.nextfreeLargelcpvalue;
    gt_assert(lcpsubtab->lcp2file->outfpllvtab != NULL);
    gt_xfwrite(lcpsubtab->lcp2file->largelcpvalues.spaceLargelcpvalue,
               sizeof (*lcpsubtab->lcp2file->largelcpvalues.
                                   spaceLargelcpvalue),
               (size_t) lcpsubtab->lcp2file->largelcpvalues.
                                   nextfreeLargelcpvalue,
               lcpsubtab->lcp2file->outfpllvtab);
  }
}

#define NUMBEROFZEROS 1024

static unsigned long outmany0lcpvalues(unsigned long many,
                                       FILE *outfplcptab)
{
  unsigned long i, countout;
  uint8_t outvalues[NUMBEROFZEROS] = {0};

  countout = many/NUMBEROFZEROS;
  for (i=0; i<countout; i++)
  {
    gt_xfwrite(outvalues,sizeof (uint8_t),(size_t) NUMBEROFZEROS,outfplcptab);
  }
  gt_xfwrite(outvalues,sizeof (uint8_t),(size_t) many % NUMBEROFZEROS,
             outfplcptab);
  return many;
}

void gt_Outlcpinfo_delete(Outlcpinfo *outlcpinfo)
{
  if (outlcpinfo == NULL)
  {
    return;
  }
  gt_turningwheel_delete(outlcpinfo->turnwheel);
  if (outlcpinfo->lcpsubtab.lcp2file != NULL)
  {
    if (outlcpinfo->lcpsubtab.lcp2file->countoutputlcpvalues <
        outlcpinfo->numsuffixes2output)
    {
      outlcpinfo->lcpsubtab.lcp2file->countoutputlcpvalues
        += outmany0lcpvalues(outlcpinfo->numsuffixes2output -
                             outlcpinfo->lcpsubtab.lcp2file
                                                  ->countoutputlcpvalues,
                             outlcpinfo->lcpsubtab.lcp2file->outfplcptab);
    }
    gt_assert(outlcpinfo->lcpsubtab.lcp2file->countoutputlcpvalues ==
              outlcpinfo->numsuffixes2output);
    GT_FREEARRAY(&outlcpinfo->lcpsubtab.lcp2file->largelcpvalues,
                 Largelcpvalue);
    gt_fa_fclose(outlcpinfo->lcpsubtab.lcp2file->outfplcptab);
    gt_fa_fclose(outlcpinfo->lcpsubtab.lcp2file->outfpllvtab);
    gt_free(outlcpinfo->lcpsubtab.lcp2file->reservoir);
    outlcpinfo->lcpsubtab.lcp2file->smalllcpvalues = NULL;
    outlcpinfo->lcpsubtab.lcp2file->reservoir = NULL;
    outlcpinfo->lcpsubtab.lcp2file->sizereservoir = 0;
    gt_free(outlcpinfo->lcpsubtab.lcp2file);
  } else
  {
    gt_free(outlcpinfo->lcpsubtab.tableoflcpvalues.bucketoflcpvalues);
#ifndef NDEBUG
    gt_free(outlcpinfo->lcpsubtab.tableoflcpvalues.isset);
#endif
  }
  outlcpinfo->lcpsubtab.tableoflcpvalues.bucketoflcpvalues = NULL;
#ifndef NDEBUG
  outlcpinfo->lcpsubtab.tableoflcpvalues.isset = NULL;
#endif
  outlcpinfo->lcpsubtab.tableoflcpvalues.numofentries = 0;
  if (outlcpinfo->lcpsubtab.distlcpvalues != NULL)
  {
    gt_disc_distri_show(outlcpinfo->lcpsubtab.distlcpvalues,NULL);
    gt_disc_distri_delete(outlcpinfo->lcpsubtab.distlcpvalues);
  }
  gt_free(outlcpinfo);
}

void gt_Outlcpinfo_check_lcpvalues(const GtEncseq *encseq,
                                   GtReadmode readmode,
                                   const GtSuffixsortspace *sortedsample,
                                   unsigned long effectivesamplesize,
                                   unsigned long maxdepth,
                                   const Outlcpinfo *outlcpinfosample)
{
  GT_UNUSED int cmp;
  unsigned long idx, reallcp, startpos1, startpos2, currentlcp;

  startpos1 = gt_suffixsortspace_getdirect(sortedsample,0);
  for (idx=1UL; idx<effectivesamplesize; idx++)
  {
    startpos2 = gt_suffixsortspace_getdirect(sortedsample,idx);
    cmp = gt_encseq_check_comparetwosuffixes(encseq,
                                             readmode,
                                             &reallcp,
                                             false,
                                             false,
                                             maxdepth,
                                             startpos1,
                                             startpos2,
                                             NULL,
                                             NULL);
    gt_assert(GT_ISIBITSET(outlcpinfosample->lcpsubtab.tableoflcpvalues
                                                      .isset,idx));
    currentlcp
      = outlcpinfosample->lcpsubtab.tableoflcpvalues.bucketoflcpvalues[idx];
    gt_assert(cmp <= 0);
    gt_assert(reallcp <= maxdepth);
    if (currentlcp > maxdepth || reallcp < currentlcp)
    {
      fprintf(stderr,"idx=%lu,suffixpair=%lu,%lu: "
                     "reallcp = %lu != %lu = currentlcp\n",
                      idx,startpos1,startpos2,reallcp,currentlcp);
      gt_encseq_showatstartposwithdepth(stderr,encseq,readmode,startpos1,20UL);
      fprintf(stderr,"\n");
      gt_encseq_showatstartposwithdepth(stderr,encseq,readmode,startpos2,20UL);
      fprintf(stderr,"\n");
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    startpos1 = startpos2;
  }
}

unsigned long gt_Outlcpinfo_numoflargelcpvalues(const Outlcpinfo *outlcpinfo)
{
  gt_assert(outlcpinfo->lcpsubtab.lcp2file != NULL);
  return outlcpinfo->lcpsubtab.lcp2file->totalnumoflargelcpvalues;
}

double gt_Outlcpinfo_lcptabsum(const Outlcpinfo *outlcpinfo)
{
  gt_assert(outlcpinfo != NULL);
  return outlcpinfo->lcpsubtab.lcptabsum;
}

void gt_Outlcpinfo_numsuffixes2output_set(Outlcpinfo *outlcpinfo,
                                          unsigned long numsuffixes2output)
{
  outlcpinfo->numsuffixes2output = numsuffixes2output;
}

unsigned long gt_Outlcpinfo_maxbranchdepth(const Outlcpinfo *outlcpinfo)
{
  gt_assert(outlcpinfo->lcpsubtab.lcp2file != NULL);
  return outlcpinfo->lcpsubtab.lcp2file->maxbranchdepth;
}

void gt_Outlcpinfo_prebucket(Outlcpinfo *outlcpinfo,
                             GtCodetype code,
                             unsigned long lcptaboffset)
{
  if (outlcpinfo != NULL)
  {
    outlcpinfo->lcpsubtab.tableoflcpvalues.numoflargelcpvalues = 0;
    if (outlcpinfo->lcpsubtab.lcp2file == NULL)
    {
      outlcpinfo->lcpsubtab.tableoflcpvalues.lcptaboffset = lcptaboffset;
    }
    if (code > 0)
    {
      (void) gt_turningwheel_next(outlcpinfo->turnwheel);
      if (outlcpinfo->previousbucketwasempty)
      {
        outlcpinfo->minchanged
          = MIN(outlcpinfo->minchanged,
                gt_turningwheel_minchanged(outlcpinfo->turnwheel));
      } else
      {
        outlcpinfo->minchanged
          = gt_turningwheel_minchanged(outlcpinfo->turnwheel);
      }
    }
  }
}

void gt_Outlcpinfo_nonspecialsbucket(Outlcpinfo *outlcpinfo,
                                     unsigned int prefixlength,
                                     GT_UNUSED GtSuffixsortspace *sssp,
                                     GtLcpvalues *tableoflcpvalues,
                                     const GtBucketspecification *bucketspec,
                                     GtCodetype code)
{
  if (outlcpinfo != NULL)
  {
    unsigned long lcpvalue;
    Suffixwithcode firstsuffixofbucket;

    if (outlcpinfo->previoussuffix.defined)
    {

      /* compute lcpvalue of first element of bucket with
         last element of previous bucket */
      firstsuffixofbucket.code = code;
      firstsuffixofbucket.prefixindex = prefixlength;
#ifdef SKDEBUG
      firstsuffixofbucket.startpos
        = gt_suffixsortspace_get(sssp,0,bucketspec->left);
      /*
      gt_bcktab_consistencyofsuffix(__LINE__,
                                    encseq,readmode,bcktab,numofchars,
                                    &firstsuffixofbucket);
      */
#endif
      lcpvalue = computelocallcpvalue(&outlcpinfo->previoussuffix,
                                      &firstsuffixofbucket,
                                      outlcpinfo->minchanged);
    } else
    {
      /* first part first code */
      lcpvalue = 0;
    }
    lcptab_update(tableoflcpvalues,0,0,lcpvalue);
    /* all other lcp-values are computed and they can be output */
    outlcpvalues(&outlcpinfo->lcpsubtab,
                 0,
                 bucketspec->nonspecialsinbucket-1,
                 bucketspec->left);
    /* previoussuffix becomes last nonspecial element in current bucket */
    outlcpinfo->previoussuffix.code = code;
    outlcpinfo->previoussuffix.prefixindex = prefixlength;
#ifdef SKDEBUG
    outlcpinfo->previoussuffix.startpos
      = gt_suffixsortspace_get(sssp,0,
                               bucketspec->left
                                 + bucketspec->nonspecialsinbucket - 1);
    /*
    gt_bcktab_consistencyofsuffix(__LINE__,
                                  encseq,readmode,bcktab,numofchars,
                                  &outlcpinfo->previoussuffix);
    */
#endif
  }
}

void gt_Outlcpinfo_postbucket(Outlcpinfo *outlcpinfo,
                              unsigned int prefixlength,
                              GtSuffixsortspace *sssp,
                              const GtBcktab *bcktab,
                              const GtBucketspecification *bucketspec,
                              GtCodetype code)
{
  if (outlcpinfo != NULL)
  {
    if (bucketspec->specialsinbucket > 0)
    {
      unsigned int minprefixindex;
      unsigned long suffixvalue
        = gt_suffixsortspace_get(sssp,
                                 0,
                                 bucketspec->left
                                   + bucketspec->nonspecialsinbucket);
      minprefixindex = bucketends(&outlcpinfo->lcpsubtab,
                                  &outlcpinfo->previoussuffix,
                                  /* first special element in bucket */
                                  suffixvalue,
                                  outlcpinfo->minchanged,
                                  bucketspec->specialsinbucket,
                                  code,
                                  bcktab);
      if (outlcpinfo->lcpsubtab.lcp2file != NULL)
      {
        outsmalllcpvalues(outlcpinfo->lcpsubtab.lcp2file,
                          bucketspec->specialsinbucket);
      }
      /* there is at least one special element: this is the last element
         in the bucket, and thus the previoussuffix for the next round */
      outlcpinfo->previoussuffix.defined = true;
      outlcpinfo->previoussuffix.code = code;
      outlcpinfo->previoussuffix.prefixindex = minprefixindex;
#ifdef SKDEBUG
      outlcpinfo->previoussuffix.startpos
        = gt_suffixsortspace_get(sssp,
                                 0,
                                 bucketspec->left
                                   + bucketspec->nonspecialsinbucket +
                                     bucketspec->specialsinbucket - 1);
      /*
        gt_bcktab_consistencyofsuffix(__LINE__,
                                      encseq,readmode,bcktab,numofchars,
                                      &outlcpinfo->previoussuffix);
      */
#endif
    } else
    {
      if (bucketspec->nonspecialsinbucket > 0)
      {
        /* if there is at least one element in the bucket, then the last
           one becomes the next previous suffix */
        outlcpinfo->previoussuffix.defined = true;
        outlcpinfo->previoussuffix.code = code;
        outlcpinfo->previoussuffix.prefixindex = prefixlength;
#ifdef SKDEBUG
        outlcpinfo->previoussuffix.startpos
          = gt_suffixsortspace_get(sssp,
                                   0,
                                   bucketspec.left
                                     + bucketspec.nonspecialsinbucket-1);
        /*
        gt_bcktab_consistencyofsuffix(__LINE__,
                                      encseq,readmode,bcktab,numofchars,
                                      &outlcpinfo->previoussuffix);
        */
#endif
      }
    }
    if (bucketspec->nonspecialsinbucket + bucketspec->specialsinbucket == 0)
    {
      outlcpinfo->previousbucketwasempty = true;
    } else
    {
      outlcpinfo->previousbucketwasempty = false;
    }
  }
}

GtLcpvalues *gt_Outlcpinfo_resizereservoir(Outlcpinfo *outlcpinfo,
                                           const GtBcktab *bcktab)
{
  Lcpsubtab *lcpsubtab = &outlcpinfo->lcpsubtab;

  if (lcpsubtab->lcp2file != NULL)
  {
    size_t sizeforlcpvalues; /* in bytes */

    gt_assert(bcktab != NULL);
    sizeforlcpvalues = gt_bcktab_sizeforlcpvalues(bcktab);
    if (lcpsubtab->lcp2file->sizereservoir < sizeforlcpvalues)
    {
      lcpsubtab->lcp2file->sizereservoir = sizeforlcpvalues;
      lcpsubtab->lcp2file->reservoir
        = gt_realloc(lcpsubtab->lcp2file->reservoir,
                     lcpsubtab->lcp2file->sizereservoir);
          /* point to the same area, since this is not used simultaneously */
          /* be careful for the parallel version */
      lcpsubtab->lcp2file->smalllcpvalues
        = (uint8_t *) lcpsubtab->lcp2file->reservoir;
#ifndef NDEBUG
      lcpsubtab->tableoflcpvalues.isset = NULL;
#endif
      lcpsubtab->tableoflcpvalues.bucketoflcpvalues
        = (unsigned long *) lcpsubtab->lcp2file->reservoir;
      lcpsubtab->tableoflcpvalues.lcptaboffset = 0;
      lcpsubtab->tableoflcpvalues.numofentries
        = (unsigned long) lcpsubtab->lcp2file->sizereservoir/
          sizeof (*lcpsubtab->tableoflcpvalues.bucketoflcpvalues);
    }
  }
  return &lcpsubtab->tableoflcpvalues;
}
