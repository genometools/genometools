/*
  Copyright (c) 2009-2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009-2012 Center for Bioinformatics, University of Hamburg

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

/* This module implements the difference cover based suffix sorting
   method from the following paper:

  @INPROCEEDINGS{BUR:KAER:2003,
  author = {Burkhardt, S. and K{\"a}rkk{\"a}inen, J.},
  title = {{Fast Lightweight Suffix Array Construction and Checking}},
  booktitle = {{Proceedings of the 14th Annual Symposium on Combinatorial
                Pattern Matching (CPM)}},
  year = {2003},
  editor = {{Baeza-Yates, R. and Ch{\'a}vez, E. and Crochemore, M.}},
  volume = {2676},
  series = {LNCS},
  pages = {200-210},
  publisher = {Springer-Verlag}
}

  The computation of the lcp follows the general concept of

  @INPROCEEDINGS{PUG:TUR:2008,
  author = {Puglisi, S.J. and Turpin, A.},
  title = {Space-Time Tradeoffs for Longest-Common-Prefix Array Computation},
  booktitle = {Proceedings of Algorithms and Computation, 19th International
               Symposium, {ISAAC} 2008, Gold Coast, Australia,
               December 15-17, 2008. Proceedings},
  year = {2008},
  editor = {Hong, S.-H. and Nagamochi, H. and Fukunaga, T.},
  volume = {5369},
  series = {Lecture Notes in Computer Science},
  pages = {124--135},
  publisher = {Springer},
  url = {http://dx.doi.org/10.1007/978-3-540-92182-0}
}

  based on the method of BUR:KAER:2003 to map sample positions to
  indexes in the sorted array of suffixes.

*/

#include <stdbool.h>
#include <math.h>
#include "core/arraydef.h"
#include "core/assert_api.h"
#include "core/error_api.h"
#include "core/divmodmul.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#include "core/mathsupport.h"
#include "core/unused_api.h"
#include "core/intbits.h"
#include "core/encseq.h"
#include "core/logger.h"
#include "core/encseq.h"
#include "core/stack-inlined.h"
#include "core/spacecalc.h"
#include "intcode-def.h"
#include "bcktab.h"
#include "initbasepower.h"
#include "sfx-strategy.h"
#include "sfx-diffcov.h"
#include "sfx-suffixgetset.h"
#include "sfx-bentsedg.h"
#include "sfx-apfxlen.h"
#include "sfx-suftaborder.h"
#include "sfx-enumcodes.h"

typedef unsigned char Diffrank;
#define Diffrankmax ((Diffrank) 255)
typedef unsigned short Diffvalue;
#define Diffvaluemax ((Diffvalue) 65535)

#define GT_MODV(VAL) ((VAL) & dcov->vmodmask)
#define GT_DIVV(VAL) ((VAL) >> dcov->logmod)

typedef struct
{
  GtUword key,
                suffixstart;
} GtDcItventry;

typedef struct
{
  GtUword blisbl, /* bucketleftindex + subbucketleft */
          width;
} GtDcPairsuffixptr;

GT_DECLAREARRAYSTRUCT(GtDcPairsuffixptr);

typedef GtDcPairsuffixptr GtInl_Queueelem;

#include "queue-inline.h"

typedef struct
{
  GtUword blisbl, /* bucketleftindex + subbucketleft */
          width,
          count,
          totalwidth,
          maxwidth,
          depth;
  bool defined;
} GtDcFirstwithnewdepth;

typedef struct
{
  unsigned int offset;
  GtUword idx1, idx2;
} GtLcptrace;

typedef uint32_t GtDifferencecover_Inversesuftabtype;
#define GT_DIFFERENCECOVER_INVERSESUFTAB_UNDEF UINT32_MAX
#define GT_DIFFERENCECOVER_MAX_SAMPLESIZE      (UINT32_MAX-1)

struct GtDifferencecover
{
  unsigned int vparam,
               logmod,
               size,
               vmodmask,
               hvalue,  /* not necessary */
               numofchars,
               prefixlength;
  GtUword *coverrank_evaluated;
  GtBitsequence *coverrank_bits;
  Diffvalue *diffvalues,  /* points to the difference cover */
            *diff2pos;    /* table d from BUR:KAER:2003 */
  size_t requiredspace;
  GtUword totallength;
  GtLeftborder *leftborder; /* points to bcktab->leftborder */
  GtBcktab *bcktab;
  const GtEncseq *encseq;
  GtReadmode readmode;
  GtUword samplesize, effectivesamplesize, samplesize_upperbound;
  const GtCodetype **multimappower;
  GtCodetype *filltable;
  GtDifferencecover_Inversesuftabtype *inversesuftab;
  GtUword allocateditvinfo,
                currentqueuesize,
                maxqueuesize,
                currentdepth;
  GtCodetype maxcode;
  GtDcFirstwithnewdepth firstwithnewdepth;
  GtInl_Queue *rangestobesorted;
  GtDcItventry *itvinfo;
  GtArrayGtDcPairsuffixptr firstgeneration;
  GtLcpvalues *samplelcpvalues;
  GtUword firstgenerationtotalwidth,
                firstgenerationcount;
  GtLogger *logger;
  GtSuffixsortspace *sortedsample;
  GtRMQ *rmq;
};

/* Compute difference cover on the fly */

#define UScast(X) ((Diffvalue) X)
#define UCcast(X) ((Diffrank) X)

#include "tab-diffcover.h"

static GtUword dc_suffixptrget(const GtDifferencecover *dcov,GtUword idx)
{
  gt_suffixsortspace_nooffsets(dcov->sortedsample);
  return gt_suffixsortspace_get(dcov->sortedsample,0,idx);
}

static void dc_suffixptrset(GtDifferencecover *dcov,GtUword idx,GtUword value)
{
  gt_suffixsortspace_nooffsets(dcov->sortedsample);
  gt_suffixsortspace_set(dcov->sortedsample,0,idx,value);
}

static void dc_fillcoverrank(GtDifferencecover *dcov)
{
  unsigned int i;
  Diffrank j;
  const GtUword step = (GtUword) GT_DIVV(dcov->totallength) + 1;
  GtUword sum;

  dcov->coverrank_evaluated
    = gt_malloc(sizeof (*dcov->coverrank_evaluated) * dcov->vparam);
  GT_INITBITTAB(dcov->coverrank_bits,dcov->vparam);
  dcov->requiredspace += sizeof (*dcov->coverrank_evaluated) * dcov->vparam;
  gt_assert(dcov->size <= Diffrankmax);
  for (i=0; i<dcov->vparam; i++)
  {
    dcov->coverrank_evaluated[i] = ULONG_MAX; /* initialize as undefined */
  }
  for (sum = 0, j=0; j<dcov->size; j++, sum += step)
  {
    Diffvalue d = dcov->diffvalues[j];
    dcov->coverrank_evaluated[d] = sum; /* jth value from difference cover
                                           gets rank j * step. */
    GT_SETIBIT(dcov->coverrank_bits,d);
  }
}

static bool dc_is_in_differencecover(const GtDifferencecover *dcov,
                                     GtUword modpos)
{
  return GT_ISIBITSET(dcov->coverrank_bits,modpos) ? true : false;
}

static void dc_filldiff2pos(GtDifferencecover *dcov)
{
  Diffvalue *iptr, *jptr;

  gt_assert(dcov->diff2pos == NULL);
  dcov->diff2pos = gt_malloc(sizeof (*dcov->diff2pos) * dcov->vparam);
  dcov->requiredspace += sizeof (*dcov->diff2pos) * dcov->vparam;
  for (iptr=dcov->diffvalues + dcov->size - 1; iptr>=dcov->diffvalues; iptr--)
  {
    for (jptr=dcov->diffvalues; jptr<dcov->diffvalues + dcov->size; jptr++)
    {
      dcov->diff2pos[GT_MODV(*jptr - *iptr)] = *iptr;
    }
  }
}

#ifdef WITHcomputehvalue

/* XXX: the following function is currently not used */

static unsigned int cd_computehvalue(const GtDifferencecover *dcov,
                                     GtUword totallength)
{
  Diffvalue next;
  unsigned int h, nmodv = GT_MODV(totallength);

  for (h = 0; h < dcov->size; h++)
  {
    if (dcov->diffvalues[h] <= (Diffvalue) nmodv)
    {
      if (h + 1 < dcov->size)
      {
        next = dcov->diffvalues[h+1];
      } else
      {
        next = (Diffvalue) dcov->vparam;
      }
      if ((Diffvalue) nmodv < next)
      {
        return h;
      }
    }
  }
  gt_assert(false);
  return 0;
}
#endif

static int gt_differencecover_vparamverify(const GtDifferencecover *dcov,
                                           GtError *err)
{
  if (dcov->vparam < dcov->prefixlength)
  {
    gt_error_set(err,"difference cover modulo %u is too small, use larger "
                     "parameter for option -dc",dcov->vparam);
    return -1;
  }
  return 0;
}

GtDifferencecover *gt_differencecover_new(unsigned int vparam,
                                          const GtEncseq *encseq,
                                          GtReadmode readmode,
                                          unsigned int outerprefixlength,
                                          GtLogger *logger)
{
  unsigned int offset = 0, v = 1U;
  GtDifferencecover *dcov;
  bool found = false;

  dcov = gt_malloc(sizeof (*dcov));
  dcov->requiredspace = sizeof (*dcov);
  dcov->numofchars = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
  dcov->totallength = gt_encseq_total_length(encseq);
  dcov->logger = logger;
  dcov->samplesize_upperbound = 0;
  dcov->rmq = NULL;
  for (dcov->logmod = 0;
       dcov->logmod < (unsigned int) (sizeof (differencecoversizes)/
                                      sizeof (differencecoversizes[0]));
       dcov->logmod++)
  {
    if (v == vparam)
    {
      dcov->size = differencecoversizes[dcov->logmod];
      dcov->samplesize_upperbound
        = (GtUword) (GT_DIVV(dcov->totallength) + 1) * dcov->size;
      dcov->diffvalues = differencecovertab + offset;
      if (dcov->samplesize_upperbound <= GT_DIFFERENCECOVER_MAX_SAMPLESIZE)
      {
        found = true;
        break;
      }
    }
    offset += differencecoversizes[dcov->logmod];
    v = GT_MULT2(v);
  }
  if (!found)
  {
    gt_free(dcov);
    return NULL;
  }
  if (outerprefixlength == 0)
  {
    dcov->prefixlength = 0;
  } else
  {
    dcov->prefixlength
      = gt_recommendedprefixlength(dcov->numofchars,
                                   dcov->samplesize_upperbound,
                                   GT_RECOMMENDED_MULTIPLIER_DEFAULT,
                                   true);
    if (outerprefixlength > 0 && dcov->prefixlength > outerprefixlength)
    {
      dcov->prefixlength = outerprefixlength;
    }
    gt_assert(dcov->prefixlength > 0);
  }
  dcov->vparam = 1U << (dcov->logmod);
  dcov->vmodmask = dcov->vparam-1;
#ifdef WITHcomputehvalue
  dcov->hvalue = dc_computehvalue(dcov,totallength);
#endif
  dcov->encseq = encseq;
  dcov->readmode = readmode;
  dcov->bcktab = NULL;
  dcov->sortedsample = NULL;
  dcov->filltable = NULL;
  dcov->multimappower = NULL;
  dc_fillcoverrank(dcov);
  dcov->diff2pos = NULL; /* this is later initialized */
  dcov->allocateditvinfo = 0;
  dcov->itvinfo = NULL;
  dcov->currentdepth = 0;
  dcov->firstwithnewdepth.defined = false;
  dcov->firstwithnewdepth.depth = 0;
  dcov->firstwithnewdepth.totalwidth = 0;
  dcov->firstwithnewdepth.count = 0;
  dcov->firstwithnewdepth.blisbl = 0;
  dcov->firstwithnewdepth.width = 0;
  dcov->firstwithnewdepth.maxwidth = 0;
  dcov->currentqueuesize = 0;
  dcov->maxqueuesize = 0;
  dcov->inversesuftab = NULL;
  dcov->samplelcpvalues = NULL;
  dcov->firstgenerationtotalwidth = 0;
  dcov->firstgenerationcount = 0;
  GT_INITARRAY(&dcov->firstgeneration,GtDcPairsuffixptr);
  return dcov;
}

GtUword gt_differencecover_samplesize(const GtDifferencecover *dcov)
{
  return dcov->samplesize;
}

size_t gt_differencecover_requiredspace(const GtDifferencecover *dcov)
{
  return dcov->requiredspace;
}

bool gt_differencecover_is_empty(const GtDifferencecover *dcov)
{
  gt_assert(dcov != NULL);
  return dcov->effectivesamplesize == 0 ? true : false;
}

void gt_differencecover_delete(GtDifferencecover *dcov)
{
  if (dcov != NULL)
  {
    gt_assert(dcov->bcktab == NULL);
    gt_assert(dcov->sortedsample == NULL);
    gt_assert(dcov->filltable == NULL);
    gt_assert(dcov->multimappower == NULL);
    gt_free(dcov->coverrank_evaluated);
    dcov->coverrank_evaluated = NULL;
    gt_free(dcov->coverrank_bits);
    dcov->coverrank_bits = NULL;
    gt_free(dcov->diff2pos);
    dcov->diff2pos = NULL;
    gt_free(dcov->inversesuftab);
    dcov->inversesuftab = NULL;
    gt_rmq_delete(dcov->rmq);
    dcov->rmq = NULL;
    gt_free(dcov);
  }
}

/* the following implements the \mu function from BAE:KAER:2003 */

static GtUword dc_differencecover_packsamplepos(const GtDifferencecover *dcov,
                                                GtUword pos)
{
  gt_assert(dc_is_in_differencecover(dcov,(GtUword) GT_MODV(pos)));
  return dcov->coverrank_evaluated[GT_MODV(pos)] + GT_DIVV(pos);
}

GT_DECLAREARRAYSTRUCT(Codeatposition);

static GtUword dc_derivespecialcodesonthefly(GtDifferencecover *dcov,
                                             const GtArrayCodeatposition
                                                           *codelist)
{
  unsigned int prefixindex, unitsnotspecial;
  Enumcodeatposition *ecp;
  Specialcontext specialcontext;
  GtUword countderived = 0, pos, sampleindex;
  GtCodetype code;
  GtEncseqReader *esr1 = gt_encseq_create_reader_with_readmode(dcov->encseq,
                                                               dcov->readmode,
                                                               0);

  for (prefixindex=1U; prefixindex < dcov->prefixlength; prefixindex++)
  {
    /* XXX use one structure and reinit it */
    ecp = gt_Enumcodeatposition_new(dcov->encseq,dcov->readmode,
                                    dcov->prefixlength,
                                    dcov->numofchars);
    while (gt_Enumcodeatposition_next(&specialcontext,ecp))
    {
      if (prefixindex <= specialcontext.maxprefixindex)
      {
        gt_assert(specialcontext.position >= (GtUword) prefixindex);
        pos = (GtUword) (specialcontext.position - prefixindex);
        if (dc_is_in_differencecover(dcov,(GtUword) GT_MODV(pos)))
        {
          if (codelist != NULL)
          {
            gt_assert(countderived < codelist->nextfreeCodeatposition);
            gt_assert(codelist->spaceCodeatposition[countderived].maxprefixindex
                      == prefixindex);
            gt_assert(codelist->spaceCodeatposition[countderived].position
                      == pos);
          }
          code = gt_encseq_extractprefixcode(&unitsnotspecial,
                                             dcov->encseq,
                                             dcov->filltable,
                                             dcov->readmode,
                                             esr1,
                                             dcov->multimappower,
                                             pos,
                                             dcov->prefixlength);
          if (codelist != NULL)
          {
            gt_assert((GtCodetype) codelist->spaceCodeatposition[
                                           countderived].code == code);
          }
          /*
          printf("%u "GT_WU"\n",prefixindex,
                            (GtUword)
                            (specialcontext.position-prefixindex));
          */
          countderived++;
          gt_bcktab_updatespecials(dcov->bcktab,code,prefixindex);
          gt_assert(code > 0);
          sampleindex = gt_bcktab_leftborder_insertionindex(dcov->leftborder,
                                                            code);
          gt_assert(sampleindex < dcov->effectivesamplesize);
          dc_suffixptrset(dcov,sampleindex,pos);
        }
      }
    }
    gt_Enumcodeatposition_delete(ecp);
    ecp = NULL;
  }
  if (codelist != NULL)
  {
    gt_assert(countderived == codelist->nextfreeCodeatposition);
  }
  gt_encseq_reader_delete(esr1);
  return countderived;
}

static int dc_compareCodeatpositon(const void *vala,const void *valb)
{
  const Codeatposition *a = (const Codeatposition *) vala;
  const Codeatposition *b = (const Codeatposition *) valb;

  if (a->maxprefixindex < b->maxprefixindex)
  {
    return -1;
  }
  if (a->maxprefixindex > b->maxprefixindex)
  {
    return 1;
  }
  if (a->position < b->position)
  {
    return 1;
  }
  if (a->position > b->position)
  {
    return -1;
  }
  gt_assert(false);
  return 0;
}

static void dc_validate_samplepositons(const GtDifferencecover *dcov)
{
  GtUword pos;
  unsigned int modvalue;
  Diffvalue *diffptr, *afterend;
  GtUword idx;
  GtBitsequence *sampleidxused = NULL;

  GT_INITBITTAB(sampleidxused,dcov->samplesize_upperbound);
  diffptr = dcov->diffvalues;
  afterend = dcov->diffvalues + dcov->size;
  for (pos = 0, modvalue = 0; pos <= dcov->totallength; pos++)
  {
    gt_assert((GtUword) modvalue == (GtUword) GT_MODV(pos));
    gt_assert(diffptr == afterend || *diffptr >= (Diffvalue) modvalue);
    if (diffptr < afterend && (Diffvalue) modvalue == *diffptr)
    {
      idx = dc_differencecover_packsamplepos(dcov,pos);
      gt_assert(sampleidxused != NULL);
      if (GT_ISIBITSET(sampleidxused,idx))
      {
        fprintf(stderr,
                "sample index "GT_WU" for pos "GT_WU" already used before\n",
                idx,(GtUword) pos);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
      GT_SETIBIT(sampleidxused,idx);
      diffptr++;
    }
    if (modvalue < dcov->vmodmask)
    {
      modvalue++;
    } else
    {
      modvalue = 0;
      diffptr = dcov->diffvalues;
    }
  }
  gt_free(sampleidxused);
}

static void dc_inversesuftab_set(GtDifferencecover *dcov,
                                 GtUword pos,
                                 GtUword sampleindex)
{
  gt_assert (dcov != NULL && sampleindex < dcov->samplesize);
  dcov->inversesuftab[dc_differencecover_packsamplepos(dcov,pos)]
    = (uint32_t) sampleindex;
}

static GtUword dc_inversesuftab_get(const GtDifferencecover *dcov,GtUword pos)
{
  gt_assert (dcov != NULL);
  return (GtUword)
         dcov->inversesuftab[dc_differencecover_packsamplepos(dcov,pos)];
}

static void dc_initinversesuftabnonspecials(GtDifferencecover *dcov)
{
  GtUword sampleindex;

  gt_assert(dcov != NULL);
  for (sampleindex=0; sampleindex < dcov->effectivesamplesize; sampleindex++)
  {
    GtUword pos = dc_suffixptrget(dcov,sampleindex);
    dc_inversesuftab_set(dcov,pos,sampleindex);
  }
}

static GtUword dc_insertfullspecialrangesample(GtDifferencecover *dcov,
                                               GtUword specialidx,
                                               GtUword leftpos,
                                               GtUword rightpos)
{
  GtUword pos;

  gt_assert(dcov != NULL && leftpos < rightpos);
  if (GT_ISDIRREVERSE(dcov->readmode))
  {
    pos = rightpos - 1;
  } else
  {
    pos = leftpos;
  }
  while (true)
  {
    if (GT_ISDIRREVERSE(dcov->readmode))
    {
      GtUword revpos;

      gt_assert(pos < dcov->totallength);
      revpos = GT_REVERSEPOS(dcov->totallength,pos);
      if (dc_is_in_differencecover(dcov,(GtUword) GT_MODV(revpos)))
      {
        dc_inversesuftab_set(dcov,revpos,specialidx);
        specialidx++;
      }
      if (pos == leftpos)
      {
        break;
      }
      pos--;
    } else
    {
      if (dc_is_in_differencecover(dcov,(GtUword) GT_MODV(pos)))
      {
        dc_inversesuftab_set(dcov,pos,specialidx);
        specialidx++;
      }
      if (pos == rightpos-1)
      {
        break;
      }
      pos++;
    }
  }
  return specialidx;
}

static void dc_initinversesuftabspecials(GtDifferencecover *dcov)
{
  GtUword idx;

  gt_assert(dcov != NULL);
  dcov->inversesuftab = gt_malloc(sizeof (*dcov->inversesuftab) *
                                  dcov->samplesize_upperbound);
  for (idx = 0; idx < dcov->samplesize_upperbound; idx++)
  {
    dcov->inversesuftab[idx] = GT_DIFFERENCECOVER_INVERSESUFTAB_UNDEF;
  }
  dcov->requiredspace += sizeof (*dcov->inversesuftab) *
                         dcov->samplesize_upperbound;
  if (gt_encseq_has_specialranges(dcov->encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;
    GtUword specialidx;

    sri = gt_specialrangeiterator_new(dcov->encseq,
                                      GT_ISDIRREVERSE(dcov->readmode)
                                      ? false : true);
    specialidx = dcov->effectivesamplesize;
    while (gt_specialrangeiterator_next(sri,&range))
    {
      /*
      printf("specialrange "GT_WU" "GT_WU"\n",(GtUword) range.leftpos,
                                      (GtUword) range.rightpos);
      */
      specialidx = dc_insertfullspecialrangesample(dcov,
                                                   specialidx,
                                                   range.start,
                                                   range.end);
    }
    gt_specialrangeiterator_delete(sri);
    sri = NULL;
  }
  if (dc_is_in_differencecover(dcov,(GtUword) GT_MODV(dcov->totallength)))
  {
    gt_assert(dcov->samplesize > 0);
    dc_inversesuftab_set(dcov,dcov->totallength,dcov->samplesize-1);
  }
}

static void dc_updatewidth (GtDifferencecover *dcov,GtUword width,
                            unsigned int depth)
{
  if (width > 1UL)
  {
    gt_assert(dcov != NULL);
    dcov->firstgenerationtotalwidth += width;
    dcov->firstgenerationcount++;
    if (dcov->allocateditvinfo < width)
    {
      dcov->allocateditvinfo = width;
    }
    if (dcov->currentdepth == 0)
    {
      dcov->currentdepth = (GtUword) depth;
    } else
    {
      gt_assert(dcov->currentdepth == (GtUword) depth);
    }
  }
}

static void dc_initinversesuftabnonspecialsadjust(GtDifferencecover *dcov)
{
  GtCodetype code;
  GtUword idx = 0;
  const GtCodetype mincode = 0;
  unsigned int rightchar = 0; /* as mincode is 0 */

  for (code = mincode; code <= dcov->maxcode; code++)
  {
    GtBucketspecification bucketspec;

    rightchar = gt_bcktab_calcboundsparts(&bucketspec,
                                          dcov->bcktab,
                                          code,
                                          dcov->maxcode,
                                          dcov->effectivesamplesize,
                                          rightchar);
    for (/* Nothing */; idx < bucketspec.left; idx++)
    {
      dc_inversesuftab_set(dcov,dc_suffixptrget(dcov,idx),idx);
    }
    dc_updatewidth (dcov,bucketspec.nonspecialsinbucket,dcov->prefixlength);
    for (/* Nothing */;
         idx < bucketspec.left + bucketspec.nonspecialsinbucket;
         idx++)
    {
      dc_inversesuftab_set(dcov,dc_suffixptrget(dcov,idx),
                           bucketspec.left);
    }
  }
  for (/* Nothing */; idx < dcov->effectivesamplesize; idx++)
  {
    dc_inversesuftab_set(dcov,dc_suffixptrget(dcov,idx),idx);
  }
}

static void dc_anchorleftmost(GtDifferencecover *dcov,
                              GtUword blisbl,
                              GtUword width)
{
  GtUword idx;

  for (idx = blisbl; idx < blisbl + width; idx++)
  {
    dc_inversesuftab_set(dcov,dc_suffixptrget(dcov,idx),blisbl);
  }
}

static void dc_showintervalsizes(GtUword count,
                                 GtUword totalwidth,
                                 GtUword effectivesamplesize,
                                 GtUword maxwidth,
                                 GtUword depth,
                                 GtLogger *logger)
{
  gt_logger_log(logger,
              "level "GT_WU""
              ": (intervals="GT_WU",total="GT_WU","
              "avg=%.2f,%.2f%% of all,maxwidth="GT_WU")",
              depth,
              count,
              totalwidth,
              (double) totalwidth/count,
              100.0 * (double) totalwidth/effectivesamplesize,
              maxwidth);
}

static void dc_processunsortedrange(GtDifferencecover *dcov,
                                    GtUword blisbl,
                                    GtUword width,
                                    GtUword depth)
{
  GtDcPairsuffixptr pairelem;

  gt_assert(width >= 2UL && depth > 0);
  gt_assert(!dcov->firstwithnewdepth.defined ||
            (dcov->firstwithnewdepth.depth > 0 &&
             dcov->firstwithnewdepth.depth <= depth));
  if (dcov->firstwithnewdepth.defined &&
      dcov->firstwithnewdepth.depth == depth)
  {
    dcov->firstwithnewdepth.count++;
    dcov->firstwithnewdepth.totalwidth += width;
    if (dcov->firstwithnewdepth.maxwidth < width)
    {
      dcov->firstwithnewdepth.maxwidth = width;
    }
  } else
  {
    if (dcov->firstwithnewdepth.defined)
    {
      dc_showintervalsizes(dcov->firstwithnewdepth.count,
                           dcov->firstwithnewdepth.totalwidth,
                           dcov->effectivesamplesize,
                           dcov->firstwithnewdepth.maxwidth,
                           dcov->firstwithnewdepth.depth,
                           dcov->logger);
    } else
    {
      dcov->firstwithnewdepth.defined = true;
    }
    gt_logger_log(dcov->logger,"enter new level "GT_WU"",depth);
    dcov->firstwithnewdepth.blisbl = blisbl;
    dcov->firstwithnewdepth.width = width;
    dcov->firstwithnewdepth.depth = depth;
    dcov->firstwithnewdepth.count = 1UL;
    dcov->firstwithnewdepth.totalwidth = width;
    dcov->firstwithnewdepth.maxwidth = width;
  }
  pairelem.blisbl = blisbl;
  pairelem.width = width;
  gt_inl_queue_add(dcov->rangestobesorted,pairelem,false);
  dcov->currentqueuesize++;
  if (dcov->maxqueuesize < dcov->currentqueuesize)
  {
    dcov->maxqueuesize = dcov->currentqueuesize;
  }
}

static int dc_compareitv(const void *a,const void *b)
{
  const GtDcItventry *itva = (const GtDcItventry *) a,
                     *itvb = (const GtDcItventry *) b;

  if (itva->key < itvb->key)
  {
    return -1;
  }
  if (itva->key > itvb->key)
  {
    return 1;
  }
  return 0;
}

static void dc_setlcpvaluesofrunsortedrange(GtLcpvalues *samplelcpvalues,
                                            GtUword blisbl,
                                            GtUword width,
                                            GtUword lcpvalue)
{
  GtUword idx;

  for (idx = blisbl+1; idx < blisbl + width; idx++)
  {
    gt_lcptab_update(samplelcpvalues,0,idx,lcpvalue);
  }
}

static void dc_sortsuffixesonthislevel(GtDifferencecover *dcov,
                                       GtUword blisbl,
                                       GtUword width)
{
  GtUword idx, rangestart, startpos;

  gt_assert(dcov != NULL);
  if (dcov->itvinfo == NULL)
  {
    dcov->itvinfo = gt_malloc(sizeof (*dcov->itvinfo) *
                              dcov->allocateditvinfo);
  }
  if (dcov->firstwithnewdepth.blisbl == blisbl &&
      dcov->firstwithnewdepth.width == width)
  {
    dcov->currentdepth = dcov->firstwithnewdepth.depth;
  }
  gt_assert(dcov->allocateditvinfo >= width);
  /*
  printf("new interval of width "GT_WU"\n",width);
  */
  for (idx=0; idx<width; idx++)
  {
    startpos = dc_suffixptrget(dcov,blisbl+idx);
    dcov->itvinfo[idx].suffixstart = startpos;
    dcov->itvinfo[idx].key
      = dc_inversesuftab_get(dcov,startpos + dcov->currentdepth);
  }
  qsort(dcov->itvinfo,(size_t) width,sizeof (*dcov->itvinfo),dc_compareitv);
  for (idx=0; idx<width; idx++)
  {
    dc_suffixptrset(dcov,blisbl+idx,dcov->itvinfo[idx].suffixstart);
  }
  rangestart = 0;
  for (idx=1UL; idx<width; idx++)
  {
    if (dcov->itvinfo[idx-1].key != dcov->itvinfo[idx].key)
    {
      if (rangestart + 1 < idx)
      {
        dc_processunsortedrange(dcov,
                                blisbl + rangestart,
                                idx - rangestart,
                                GT_MULT2(dcov->currentdepth));
        if (dcov->samplelcpvalues != NULL)
        {
          dc_setlcpvaluesofrunsortedrange(dcov->samplelcpvalues,
                                          blisbl + rangestart,
                                          idx - rangestart,
                                          GT_MULT2(dcov->currentdepth));
        }
        dc_anchorleftmost(dcov, blisbl + rangestart, idx - rangestart);
      } else
      {
        GtUword currentsuftabentry
          = dc_suffixptrget(dcov,blisbl+rangestart);

        gt_assert(rangestart + 1 == idx);
        dc_inversesuftab_set(dcov,currentsuftabentry,
                             blisbl+rangestart);
      }
      if (dcov->samplelcpvalues != NULL)
      {
        gt_lcptab_update(dcov->samplelcpvalues,
                         0,blisbl + idx,dcov->currentdepth);
      }
      rangestart = idx;
    }
  }
  if (rangestart + 1 < width)
  {
    dc_processunsortedrange(dcov,
                            blisbl + rangestart,
                            width - rangestart,
                            GT_MULT2(dcov->currentdepth));
    if (dcov->samplelcpvalues != NULL)
    {
      dc_setlcpvaluesofrunsortedrange(dcov->samplelcpvalues,
                                      blisbl + rangestart,
                                      width - rangestart,
                                      GT_MULT2(dcov->currentdepth));
    }
    dc_anchorleftmost(dcov, blisbl + rangestart, width - rangestart);
  } else
  {
    GtUword currentsuftabentry
      = dc_suffixptrget(dcov,blisbl+rangestart);
    dc_inversesuftab_set(dcov,currentsuftabentry,blisbl+rangestart);
  }
}

static void dc_bcktab2firstlevelintervals(GtDifferencecover *dcov)
{
  GtCodetype code;
  unsigned int rightchar;
  const GtCodetype mincode = 0;

  gt_assert(dcov != NULL);
  gt_logger_log(dcov->logger,"maxbucketsize="GT_WU"",dcov->allocateditvinfo);
  rightchar = (unsigned int) (mincode % dcov->numofchars);
  for (code = 0; code <= dcov->maxcode; code++)
  {
    GtBucketspecification bucketspec;

    rightchar = gt_bcktab_calcboundsparts(&bucketspec,
                                          dcov->bcktab,
                                          code,
                                          dcov->maxcode,
                                          dcov->effectivesamplesize,
                                          rightchar);
    if (bucketspec.nonspecialsinbucket > 1UL)
    {
      dc_sortsuffixesonthislevel(dcov,
                                 bucketspec.left,
                                 bucketspec.nonspecialsinbucket);
    }
  }
}

static void dc_addunsortedrange(void *voiddcov,
                                GT_UNUSED GtSuffixsortspace *sssp,
                                GtUword blisbl,
                                GtUword width,
                                GT_UNUSED GtUword depth)
{
  GtDifferencecover *dcov = (GtDifferencecover *) voiddcov;
  GtDcPairsuffixptr *ptr;

  gt_assert(depth >= (GtUword) dcov->vparam);
  dc_updatewidth (dcov,width,dcov->vparam);
  GT_GETNEXTFREEINARRAY(ptr,&dcov->firstgeneration,GtDcPairsuffixptr,1024);
  ptr->blisbl = blisbl;
  ptr->width = width;
}

static void dc_sortremainingsamples(GtDifferencecover *dcov)
{
  GtDcPairsuffixptr *pairptr;

  if (dcov->firstgenerationcount > 0)
  {
    dc_showintervalsizes(dcov->firstgenerationcount,
                         dcov->firstgenerationtotalwidth,
                         dcov->effectivesamplesize,
                         dcov->allocateditvinfo,
                         dcov->currentdepth,
                         dcov->logger);
  }
  if (dcov->inversesuftab == NULL)
  { /* now maxdepth > prefixlength */
    dc_initinversesuftabspecials(dcov);
    dc_initinversesuftabnonspecials(dcov);
  } else
  {
    gt_assert(dcov->firstgeneration.nextfreeGtDcPairsuffixptr == 0);
  }
  for (pairptr = dcov->firstgeneration.spaceGtDcPairsuffixptr;
       pairptr < dcov->firstgeneration.spaceGtDcPairsuffixptr +
                 dcov->firstgeneration.nextfreeGtDcPairsuffixptr;
       pairptr++)
  {
    dc_anchorleftmost(dcov,pairptr->blisbl,pairptr->width);
  }
  for (pairptr = dcov->firstgeneration.spaceGtDcPairsuffixptr;
       pairptr < dcov->firstgeneration.spaceGtDcPairsuffixptr +
                 dcov->firstgeneration.nextfreeGtDcPairsuffixptr;
       pairptr++)
  {
    dc_sortsuffixesonthislevel(dcov,pairptr->blisbl, pairptr->width);
  }
  GT_FREEARRAY(&dcov->firstgeneration,GtDcPairsuffixptr);
  while (!gt_inl_queue_isempty(dcov->rangestobesorted))
  {
    GtDcPairsuffixptr thispair;
    thispair = gt_inl_queue_get(dcov->rangestobesorted);
    gt_assert(dcov->currentqueuesize > 0);
    dcov->currentqueuesize--;
    dc_sortsuffixesonthislevel(dcov,thispair.blisbl,thispair.width);
  }
  gt_logger_log(dcov->logger,"maxqueuesize="GT_WU"",dcov->maxqueuesize);
  gt_free(dcov->itvinfo);
  dcov->itvinfo = NULL;
  gt_inl_queue_delete(dcov->rangestobesorted);
  dcov->rangestobesorted = NULL;
}

static void dc_init_sfxstrategy_for_sample(Sfxstrategy *sfxstrategy,
                                           const Sfxstrategy *mainsfxstrategy,
                                           bool bitwise_cmp_ok,
                                           GtUword effectivesamplesize,
                                           GtUword totallength,
                                           GtLogger *logger)
{
  if (mainsfxstrategy != NULL)
  {
    double sampledproportion = (double) effectivesamplesize/totallength;

    *sfxstrategy = *mainsfxstrategy;
#define SETMAXCOUNT(COMP)\
    if (mainsfxstrategy->COMP >= 1UL)\
    {\
      sfxstrategy->COMP = MAX(2UL,mainsfxstrategy->COMP * sampledproportion);\
    }
    SETMAXCOUNT(maxcountingsort);
    SETMAXCOUNT(maxbltriesort);
    SETMAXCOUNT(maxinsertionsort);
  } else
  {
    defaultsfxstrategy(sfxstrategy,bitwise_cmp_ok);
  }
  gt_logger_log(logger,"samplesort.maxinsertionsort="GT_WU"",
                sfxstrategy->maxinsertionsort);
  gt_logger_log(logger,"samplesort.maxbltriesort="GT_WU"",
                sfxstrategy->maxbltriesort);
  gt_logger_log(logger,"samplesort.maxcountingsort="GT_WU"",
                sfxstrategy->maxcountingsort);
}

static void dc_fill_samplelcpvalues(bool cmpcharbychar,GtDifferencecover *dcov)
{

  GtUword suffix, kvalue, lcpinherit, start0, start1, currentlcpvalue;
  GtDifferencecover_Inversesuftabtype *inversesuftabptr = dcov->inversesuftab;
  unsigned int svalue;
  GT_UNUSED int retval;
  GtCommonunits commonunits;
  GtUchar cc1, cc2;
  GtEncseqReader *esr1, *esr2;

  esr1 = gt_encseq_create_reader_with_readmode(dcov->encseq,dcov->readmode,0);
  esr2 = gt_encseq_create_reader_with_readmode(dcov->encseq,dcov->readmode,0);
  for (svalue = 0; svalue < dcov->size; svalue++)
  {
    lcpinherit = 0;
    for (kvalue = (GtUword) svalue; kvalue < dcov->samplesize;
         kvalue += dcov->size)
    {
      do
      {
        if (inversesuftabptr >= dcov->inversesuftab +
                                dcov->samplesize_upperbound)
        {
          suffix = (GtUword) GT_DIFFERENCECOVER_INVERSESUFTAB_UNDEF;
          break;
        }
        suffix = (GtUword) *inversesuftabptr;
        inversesuftabptr++;
      } while (suffix == (GtUword)GT_DIFFERENCECOVER_INVERSESUFTAB_UNDEF);
      if (suffix > 0 && suffix < dcov->effectivesamplesize)
      {
        currentlcpvalue = gt_lcptab_getvalue(dcov->samplelcpvalues,0,suffix);
        if (currentlcpvalue < (GtUword) dcov->vparam)
        {
          lcpinherit = 0;
        } else
        {
          if (lcpinherit < currentlcpvalue)
          {
            lcpinherit = currentlcpvalue;
          }
          start0 = gt_suffixsortspace_get(dcov->sortedsample,0,suffix-1);
          start1 = gt_suffixsortspace_get(dcov->sortedsample,0,suffix);
          if (cmpcharbychar)
          {
            while (start0 + lcpinherit < dcov->totallength &&
                   start1 + lcpinherit < dcov->totallength)
            {
              cc1 = gt_encseq_get_encoded_char(dcov->encseq,start0+lcpinherit,
                                               dcov->readmode);
              cc2 = gt_encseq_get_encoded_char(dcov->encseq,start1+lcpinherit,
                                               dcov->readmode);
              if (ISSPECIAL(cc1) || ISSPECIAL(cc2) || cc1 != cc2)
              {
                break;
              }
              lcpinherit++;
            }
          } else
          {
            retval = gt_encseq_compare_viatwobitencoding(&commonunits,
                                                         dcov->encseq,
                                                         dcov->encseq,
                                                         dcov->readmode,
                                                         esr1,
                                                         esr2,
                                                         start0,
                                                         start1,
                                                         lcpinherit,
                                                         0);
            gt_assert(retval <= 0);
            lcpinherit = commonunits.finaldepth;
          }
          gt_lcptab_update(dcov->samplelcpvalues,0,suffix,lcpinherit);
          lcpinherit = lcpinherit > (GtUword) dcov->vparam
                         ? (lcpinherit - (GtUword) dcov->vparam)
                         : 0;
        }
      }
    }
  }
  gt_encseq_reader_delete(esr1);
  gt_encseq_reader_delete(esr2);
}

static void dc_verify_inversesuftab(const GtDifferencecover *dcov)
{
  GtUword idx;

  for (idx=0; idx < dcov->effectivesamplesize; idx++)
  {
    GtUword idx2 = dc_inversesuftab_get(dcov,dc_suffixptrget(dcov,idx));
    if (idx != idx2)
    {
      fprintf(stderr,"idx = "GT_WU" != "GT_WU" = idx2\n",idx,idx2);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

static void dc_differencecover_sortsample(GtDifferencecover *dcov,
                                          GtOutlcpinfo *outlcpinfosample,
                                          const Sfxstrategy *mainsfxstrategy,
                                          GtTimer *sfxprogress,
                                          bool withcheck)
{
  GtUword pos, sampleindex, posinserted, fullspecials = 0, specials = 0;
  unsigned int modvalue, unitsnotspecial;
  Diffvalue *diffptr, *afterend;
  GtCodetype code;
  GtArrayCodeatposition codelist;
  Codeatposition *codeptr;
  GtEncseqReader *esr1;

  dcov->samplesize = 0;
  dcov->bcktab = gt_bcktab_new(dcov->numofchars,
                               dcov->prefixlength,
                               dcov->totallength+1,
                               true, /* storespecialcodes */
                               true, /* withspecialsuffixes */
                               NULL,
                               NULL);
  dcov->multimappower = gt_bcktab_multimappower(dcov->bcktab);
  dcov->maxcode = gt_bcktab_numofallcodes(dcov->bcktab) - 1;
  esr1 = gt_encseq_create_reader_with_readmode(dcov->encseq,dcov->readmode,0);
  dcov->rangestobesorted = gt_inl_queue_new(MAX(16UL,GT_DIV2(dcov->maxcode)));
  dcov->filltable = gt_filllargestchartable(dcov->numofchars,
                                            dcov->prefixlength);
  gt_assert(dcov->bcktab != NULL);
  dcov->leftborder = gt_bcktab_leftborder(dcov->bcktab);
  GT_INITARRAY(&codelist,Codeatposition);
  diffptr = dcov->diffvalues;
  afterend = dcov->diffvalues + dcov->size;
  for (pos = 0, modvalue = 0; pos <= dcov->totallength; pos++)
  {
    if (diffptr < afterend && (Diffvalue) modvalue == *diffptr)
    {
      if (pos < dcov->totallength)
      {
        code = gt_encseq_extractprefixcode(&unitsnotspecial,
                                           dcov->encseq,
                                           dcov->filltable,
                                           dcov->readmode,
                                           esr1,
                                           dcov->multimappower,
                                           pos,
                                           dcov->prefixlength);
      } else
      {
        code = 0;
        unitsnotspecial = 0;
      }
      dcov->samplesize++;
      if (unitsnotspecial > 0)
      {
        gt_bcktab_leftborder_addcode(dcov->leftborder,code);
        if (unitsnotspecial < dcov->prefixlength)
        {
          if (withcheck)
          {
            GT_GETNEXTFREEINARRAY(codeptr,&codelist,Codeatposition,128);
            gt_assert(codelist.spaceCodeatposition != NULL);
            codeptr->position = pos;
            gt_assert(code <= (GtCodetype) GT_MAXCODEVALUE);
            codeptr->code = (unsigned int) code;
            gt_assert(unitsnotspecial <= (unsigned int) GT_MAXPREFIXLENGTH);
            codeptr->maxprefixindex = unitsnotspecial;
          }
          specials++;
        }
      } else
      {
        fullspecials++;
      }
      diffptr++;
    }
    if (modvalue < dcov->vmodmask)
    {
      modvalue++;
    } else
    {
      modvalue = 0;
      diffptr = dcov->diffvalues;
    }
  }
  dcov->effectivesamplesize = dcov->samplesize - fullspecials;
  gt_assert(dcov->samplesize <= dcov->samplesize_upperbound);
  (void) gt_bcktab_leftborderpartialsums(NULL,NULL,dcov->bcktab);
  gt_logger_log(dcov->logger,
              ""GT_WU" positions are sampled (%.2f%%), pl=%u",
              dcov->samplesize,
              100.0 * (double) dcov->samplesize/(dcov->totallength+1),
              dcov->prefixlength);
  gt_logger_log(dcov->logger,"specials="GT_WU", fullspecials="GT_WU"",
                specials,fullspecials);
  if (withcheck && codelist.nextfreeCodeatposition > 1UL)
  {
    gt_assert(codelist.spaceCodeatposition != NULL);
    qsort(codelist.spaceCodeatposition,
          (size_t) codelist.nextfreeCodeatposition,
          sizeof (*codelist.spaceCodeatposition),dc_compareCodeatpositon);
  }
  if (dcov->effectivesamplesize > 0)
  {
    gt_assert(mainsfxstrategy != NULL);
    dcov->sortedsample = gt_suffixsortspace_new(dcov->effectivesamplesize,
                                                dcov->totallength,
                                                mainsfxstrategy->suftabuint,
                                                dcov->logger);
  } else
  {
    gt_assert(dcov->sortedsample == NULL);
  }
  posinserted = dc_derivespecialcodesonthefly(dcov,
                                              withcheck ? &codelist : NULL);
  GT_FREEARRAY(&codelist,Codeatposition);
  diffptr = dcov->diffvalues;
  afterend = dcov->diffvalues + dcov->size;
  for (pos = 0, modvalue = 0; pos < dcov->totallength; pos++)
  {
    if (diffptr < afterend && (Diffvalue) modvalue == *diffptr)
    {
      /* XXX: Use a function to extract the code in constant time for
         twobitencoding. */
      code = gt_encseq_extractprefixcode(&unitsnotspecial,
                                         dcov->encseq,
                                         dcov->filltable,
                                         dcov->readmode,
                                         esr1,
                                         dcov->multimappower,
                                         pos,
                                         dcov->prefixlength);
      if (unitsnotspecial == dcov->prefixlength)
      {
        sampleindex = gt_bcktab_leftborder_insertionindex(dcov->leftborder,
                                                          code);
        gt_assert(sampleindex < dcov->effectivesamplesize);
        dc_suffixptrset(dcov,sampleindex,pos);
        posinserted++;
      }
      diffptr++;
    }
    if (modvalue < dcov->vmodmask)
    {
      modvalue++;
    } else
    {
      modvalue = 0;
      diffptr = dcov->diffvalues;
    }
  }
  dcov->multimappower = NULL;
  gt_free(dcov->filltable);
  dcov->filltable = NULL;
  gt_assert(posinserted == dcov->effectivesamplesize);
  if (withcheck && dcov->effectivesamplesize > 0)
  {
    gt_checksortedsuffixes(__FILE__,
                           __LINE__,
                           dcov->encseq,
                           dcov->readmode,
                           dcov->sortedsample,
                           0,
                           dcov->effectivesamplesize,
                           false, /* specialsareequal  */
                           false,  /* specialsareequalatdepth0 */
                           (GtUword) dcov->prefixlength);
  }
  gt_Outlcpinfo_reinit(outlcpinfosample,dcov->numofchars,dcov->prefixlength,
                       dcov->effectivesamplesize);
  if (outlcpinfosample != NULL)
  {
    dcov->requiredspace += gt_Outlcpinfo_size(outlcpinfosample);
    dcov->samplelcpvalues = gt_Outlcpinfo_lcpvalues_ref(outlcpinfosample);
    gt_assert(dcov->samplelcpvalues != NULL);
  }
  if (dcov->vparam == dcov->prefixlength)
  {
    dc_initinversesuftabspecials(dcov);
    dc_initinversesuftabnonspecialsadjust(dcov);
    dc_bcktab2firstlevelintervals(dcov);
  } else
  {
    GtUint64 bucketiterstep = 0;
    Sfxstrategy sfxstrategy;

    gt_assert (dcov->vparam > dcov->prefixlength);
    dc_init_sfxstrategy_for_sample(&sfxstrategy,
                                   mainsfxstrategy,
                                   gt_encseq_bitwise_cmp_ok(dcov->encseq)
                                     ? false : true,
                                   dcov->effectivesamplesize,
                                   dcov->totallength,
                                   dcov->logger);
    /* now sort the suffix sample up to a prefix of length vparam */
    if (dcov->effectivesamplesize > 0)
    {
      gt_bcktab_determinemaxsize(dcov->bcktab, 0, dcov->maxcode,
                                 dcov->effectivesamplesize);
      gt_sortallbuckets(dcov->sortedsample,
                        dcov->effectivesamplesize,
                        NULL,
                        dcov->encseq,
                        dcov->readmode,
                        0, /* mincode */
                        dcov->maxcode,
                        dcov->bcktab,
                        dcov->numofchars,
                        dcov->prefixlength,
                        outlcpinfosample,
                        dcov->vparam,
                        &sfxstrategy,
                        dc_addunsortedrange,
                        (void *) dcov,
                        &bucketiterstep,
                        dcov->logger);
    }
    if (withcheck && dcov->effectivesamplesize > 0)
    {
      gt_checksortedsuffixes(__FILE__,
                             __LINE__,
                             dcov->encseq,
                             dcov->readmode,
                             dcov->sortedsample,
                             0,
                             dcov->effectivesamplesize,
                             false, /* specialsareequal  */
                             false,  /* specialsareequalatdepth0 */
                             (GtUword) dcov->vparam);
    }
  }
  gt_bcktab_delete(dcov->bcktab);
  dcov->bcktab = NULL;
  dc_sortremainingsamples(dcov);
  if (withcheck && dcov->effectivesamplesize > 0)
  {
    gt_checksortedsuffixes(__FILE__,
                           __LINE__,
                           dcov->encseq,
                           dcov->readmode,
                           dcov->sortedsample,
                           0,
                           dcov->effectivesamplesize,
                           false, /* specialsareequal  */
                           false,  /* specialsareequalatdepth0 */
                           0);
    dc_verify_inversesuftab(dcov);
    if (outlcpinfosample != NULL)
    {
      gt_Outlcpinfo_check_lcpvalues(dcov->encseq,
                                    dcov->readmode,
                                    dcov->sortedsample,
                                    dcov->effectivesamplesize,
                                    outlcpinfosample,
                                    false);
    }
  }
  if (outlcpinfosample != NULL)
  {
    gt_assert(mainsfxstrategy != NULL);
    dc_fill_samplelcpvalues(mainsfxstrategy->cmpcharbychar ||
                            !gt_encseq_bitwise_cmp_ok(dcov->encseq),dcov);
    if (withcheck)
    {
      gt_Outlcpinfo_check_lcpvalues(dcov->encseq,
                                    dcov->readmode,
                                    dcov->sortedsample,
                                    dcov->effectivesamplesize,
                                    outlcpinfosample,
                                    true);
    }
  }
  gt_suffixsortspace_delete(dcov->sortedsample,false);
  dcov->sortedsample = NULL;
  if (dcov->effectivesamplesize > 0 && outlcpinfosample != NULL)
  {
    size_t rmqsize;

    if (sfxprogress != NULL)
    {
      gt_timer_show_progress(sfxprogress,"preparing the RMQ",stdout);
    }
    dcov->rmq = gt_lcpvalues_rmq_new(dcov->samplelcpvalues);
    rmqsize = gt_rmq_size(dcov->rmq);
    gt_logger_log(dcov->logger,
                  "RMQ requires %.2f MB (%.2f bytes per sample position)",
                  GT_MEGABYTES(rmqsize),
                  (double) rmqsize/dcov->samplesize);
    dcov->requiredspace += rmqsize;
  }
  gt_encseq_reader_delete(esr1);
  dc_filldiff2pos(dcov);
}

static void dc_differencecover_sortsample0(GtDifferencecover *dcov,
                                           GtOutlcpinfo *outlcpinfosample,
                                           const Sfxstrategy *mainsfxstrategy,
                                           GT_UNUSED GtTimer *sfxprogress,
                                           bool withcheck)
{
  GtUword pos, posinserted, fullspecials = 0;
  unsigned int modvalue;
  Diffvalue *diffptr, *afterend;
  Sfxstrategy sfxstrategy;
  GtUchar cc;

  sfxstrategy.suftabuint = false;
  dcov->samplesize = 0;
  dcov->bcktab = NULL;
  dcov->multimappower = NULL;
  dcov->maxcode = 0;
  dcov->rangestobesorted = gt_inl_queue_new(MAX(16UL,GT_DIV2(dcov->maxcode)));
  dcov->filltable = NULL;
  dcov->leftborder = NULL;
  diffptr = dcov->diffvalues;
  afterend = dcov->diffvalues + dcov->size;
  for (pos = 0, modvalue = 0; pos <= dcov->totallength; pos++)
  {
    if (diffptr < afterend && (Diffvalue) modvalue == *diffptr)
    {
      dcov->samplesize++;
      if (pos < dcov->totallength)
      {
        cc = gt_encseq_get_encoded_char(dcov->encseq,pos,dcov->readmode);
        if (ISSPECIAL(cc))
        {
          fullspecials++;
        }
      } else
      {
        fullspecials++;
      }
      diffptr++;
    }
    if (modvalue < dcov->vmodmask)
    {
      modvalue++;
    } else
    {
      modvalue = 0;
      diffptr = dcov->diffvalues;
    }
  }
  dcov->effectivesamplesize = dcov->samplesize - fullspecials;
  gt_logger_log(dcov->logger,
              ""GT_WU" positions are sampled (%.2f) pl=%u",
              dcov->samplesize,
              100.0 * (double) dcov->samplesize/(dcov->totallength+1),
              dcov->prefixlength);
  gt_logger_log(dcov->logger,"fullspecials="GT_WU"",fullspecials);
  dcov->sortedsample = gt_suffixsortspace_new(dcov->effectivesamplesize,
                                              dcov->totallength,
                                              sfxstrategy.suftabuint,
                                              dcov->logger);
  posinserted = 0;
  diffptr = dcov->diffvalues;
  afterend = dcov->diffvalues + dcov->size;
  for (pos = 0, modvalue = 0; pos < dcov->totallength; pos++)
  {
    if (diffptr < afterend && (Diffvalue) modvalue == *diffptr)
    {
      if (pos < dcov->totallength)
      {
        cc = gt_encseq_get_encoded_char(dcov->encseq,pos,dcov->readmode);
        if (ISNOTSPECIAL(cc))
        {
          dc_suffixptrset(dcov,posinserted,pos);
          posinserted++;
        }
      }
      diffptr++;
    }
    if (modvalue < dcov->vmodmask)
    {
      modvalue++;
    } else
    {
      modvalue = 0;
      diffptr = dcov->diffvalues;
    }
  }
  gt_assert(posinserted == dcov->effectivesamplesize);

  dc_init_sfxstrategy_for_sample(&sfxstrategy,
                                 mainsfxstrategy,
                                 gt_encseq_bitwise_cmp_ok(dcov->encseq)
                                   ? false : true,
                                 dcov->effectivesamplesize,
                                 dcov->totallength,
                                 dcov->logger);
  gt_Outlcpinfo_reinit(outlcpinfosample,dcov->numofchars,dcov->prefixlength,
                       dcov->effectivesamplesize);
  gt_sortallsuffixesfromstart(dcov->sortedsample,
                              dcov->effectivesamplesize,
                              dcov->encseq,
                              dcov->readmode,
                              outlcpinfosample,
                              dcov->vparam,
                              &sfxstrategy,
                              dc_addunsortedrange,
                              (void *) dcov,
                              dcov->logger);
  if (withcheck && dcov->effectivesamplesize > 0)
  {
    gt_checksortedsuffixes(__FILE__,
                           __LINE__,
                           dcov->encseq,
                           dcov->readmode,
                           dcov->sortedsample,
                           0,
                           dcov->effectivesamplesize,
                           false, /* specialsareequal  */
                           false,  /* specialsareequalatdepth0 */
                           (GtUword) dcov->vparam);
  }
  dc_sortremainingsamples(dcov);
  if (withcheck && dcov->effectivesamplesize > 0)
  {
    GtUword idx;

    gt_checksortedsuffixes(__FILE__,
                           __LINE__,
                           dcov->encseq,
                           dcov->readmode,
                           dcov->sortedsample,
                           0,
                           dcov->effectivesamplesize,
                           false, /* specialsareequal  */
                           false,  /* specialsareequalatdepth0 */
                           0);
    for (idx=0; idx < dcov->effectivesamplesize; idx++)
    {
      GtUword idx2 = dc_inversesuftab_get(dcov,dc_suffixptrget(dcov,idx));
      if (idx != idx2)
      {
        fprintf(stderr,"idx = "GT_WU" != "GT_WU" = idx2\n",idx,idx2);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
    }
  }
  gt_suffixsortspace_delete(dcov->sortedsample,false);
  dcov->sortedsample = NULL;
  dc_filldiff2pos(dcov);
}

GtDifferencecover *gt_differencecover_prepare_sample(
                                        unsigned int vparam,
                                        const GtEncseq *encseq,
                                        GtReadmode readmode,
                                        unsigned int prefixlength,
                                        const Sfxstrategy *sfxstrategy,
                                        GtOutlcpinfo *outlcpinfosample,
                                        GtLogger *logger,
                                        GtTimer *sfxprogress,
                                        GtError *err)
{
  GtDifferencecover *dcov = NULL;

  gt_assert(vparam > 0);
  if (sfxprogress != NULL)
  {
    gt_timer_show_progress(sfxprogress,
                           (outlcpinfosample == NULL)
                           ? "sorting difference cover sample"
                           : ("sorting difference cover sample "
                             "and determine their lcp values"),
                           stdout);
  }
  dcov = gt_differencecover_new(vparam,encseq,readmode,
                                prefixlength,logger);
  if (dcov == NULL)
  {
    gt_error_set(err,"no difference cover modulo %u found",vparam);
  } else
  {
    if (gt_differencecover_vparamverify(dcov,err) != 0)
    {
      gt_differencecover_delete(dcov);
      dcov = NULL;
    } else
    {
      gt_assert(sfxstrategy != NULL);
      gt_logger_log(logger,"presorting sample suffixes according to "
                           "difference cover modulo %u",vparam);
      if (prefixlength > 0)
      {
        dc_differencecover_sortsample(dcov,outlcpinfosample,sfxstrategy,
                                      sfxprogress,sfxstrategy->dccheck);
      } else
      {
        dc_differencecover_sortsample0(dcov,outlcpinfosample,sfxstrategy,
                                       sfxprogress,sfxstrategy->dccheck);
      }
    }
  }
  return dcov;
}

void gt_differencecover_check(const GtEncseq *encseq,GtReadmode readmode)
{
  GtDifferencecover *dcov;
  size_t logmod;
  const size_t startlogmod = (size_t) 4;
  unsigned int vparam;
  bool withcheck = true;

  printf("sizeof (differencecovertab)="GT_WU"\n",
          (GtUword) sizeof (differencecovertab));
  for (logmod = startlogmod;
       logmod < sizeof (differencecoversizes)/sizeof (differencecoversizes[0]);
       logmod++)
  {
    vparam = 1U << logmod;
    dcov = gt_differencecover_new(vparam,encseq,readmode,0,NULL);
    if (dcov == NULL)
    {
      fprintf(stderr,"no difference cover for v=%u\n",vparam);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    printf("v=%u (size=%u)\n",dcov->vparam,dcov->size);
    if (withcheck)
    {
      dc_validate_samplepositons(dcov);
    }
    dc_differencecover_sortsample(dcov,NULL,NULL,NULL,withcheck);
    gt_differencecover_delete(dcov);
  }
  printf("# %u difference covers checked\n",
          (unsigned int) (logmod - startlogmod));
}

/* The remaining code only reads the difference cover, so it can be
   given as a const pointer */

/* The following is the \delta function from BUR:KAER:2003. */

static unsigned int dc_differencecover_delta(const GtDifferencecover *dcov,
                                             GtUword pos1,
                                             GtUword pos2)
{
  gt_assert(dcov != NULL);
  return (unsigned int) GT_MODV(dcov->diff2pos[GT_MODV(pos2 - pos1)] - pos1);
}

void gt_differencecover_completelargelcpvalues(const GtDifferencecover *dcov,
                                               const GtSuffixsortspace *sssp,
                                               GtLcpvalues *tableoflcpvalues,
                                               GtUword width,
                                               GtUword posoffset)
{
  GtUword idx, lcpvalue, s0, s1, r0, r1;

  gt_assert(width > 0 && sssp != NULL && tableoflcpvalues != NULL &&
            dcov->rmq != NULL);
  for (idx = 1UL; idx < width; idx++)
  {
    lcpvalue = gt_lcptab_getvalue(tableoflcpvalues,0,idx);
    if (lcpvalue >= (GtUword) dcov->vparam)
    {
      s0 = gt_suffixsortspace_get(sssp, posoffset, idx-1);
      s1 = gt_suffixsortspace_get(sssp, posoffset, idx);
      lcpvalue = (GtUword) dc_differencecover_delta(dcov,s0,s1);
      r0 = dc_inversesuftab_get(dcov,s0 + lcpvalue);
      r1 = dc_inversesuftab_get(dcov,s1 + lcpvalue);
      if (r0 < dcov->effectivesamplesize && r1 < dcov->effectivesamplesize)
      {
        gt_assert(r0 < r1);
        lcpvalue += gt_rmq_find_min_value(dcov->rmq, r0+1, r1);
      }
      gt_lcptab_update(tableoflcpvalues,0,idx,lcpvalue);
    }
  }
}

static GtUword gt_differencecover_eval_lcp(const GtLcptrace *lcptrace,
                                           const GtDifferencecover *dcov)
{
  if (lcptrace->idx1 < dcov->effectivesamplesize &&
      lcptrace->idx2 < dcov->effectivesamplesize)
  {
    gt_assert(dcov->rmq != NULL);
    return (GtUword)
            gt_rmq_find_min_value(dcov->rmq,
                                  MIN(lcptrace->idx1,lcptrace->idx2)+1,
                                  MAX(lcptrace->idx1,lcptrace->idx2))
             + lcptrace->offset;
  }
  return (GtUword) lcptrace->offset;
}

static int gt_differencecover_compare (const GtDifferencecover *dcov,
                                       GtLcptrace *lcptrace,
                                       GtUword suffixpos1,
                                       GtUword suffixpos2)
{
  gt_assert(suffixpos1 < dcov->totallength &&
            suffixpos2 < dcov->totallength);
  lcptrace->offset = dc_differencecover_delta(dcov,suffixpos1,suffixpos2);
  lcptrace->idx1 = dc_inversesuftab_get(dcov,suffixpos1 + lcptrace->offset);
  lcptrace->idx2 = dc_inversesuftab_get(dcov,suffixpos2 + lcptrace->offset);
  gt_assert(lcptrace->idx1 != lcptrace->idx2);
  return lcptrace->idx1 < lcptrace->idx2 ? -1 : 1;
}

#ifdef QSORTNAME
#undef QSORTNAME
#endif

#define QSORTNAME(NAME) dc_##NAME

#define dc_ARRAY_GET(ARR,RELIDX)\
        gt_suffixsortspace_get(ARR,subbucketleft,RELIDX)

#define dc_ARRAY_SET(ARR,RELIDX,VALUE)\
        gt_suffixsortspace_set(ARR,subbucketleft,RELIDX,VALUE)

typedef GtSuffixsortspace * QSORTNAME(Datatype);

static int QSORTNAME(qsortcmparr) (GtUword a,
                                   GtUword b,
                                   const GtDifferencecover *dcov,
                                   GtLcptrace *lcptrace,
                                   GtUword subbucketleft,
                                   const QSORTNAME(Datatype) sssp)
{
  GtUword suffixpos1, suffixpos2;

  gt_assert(sssp != NULL);
  suffixpos1 = dc_ARRAY_GET(sssp,a);
  suffixpos2 = dc_ARRAY_GET(sssp,b);
  return gt_differencecover_compare (dcov, lcptrace, suffixpos1, suffixpos2);
}

typedef void * QSORTNAME(Sorttype);

/*
 * Qsort routine from Bentley & McIlroy's ``Engineering a Sort Function''.
 */

#ifndef GT_QSORT_ARR_SWAP
#define GT_QSORT_ARR_SWAP(ARR,A,B)\
        if ((A) != (B))\
        {\
          tmp = QSORTNAME(ARRAY_GET)(ARR,A);\
          QSORTNAME(ARRAY_SET)(ARR,A,QSORTNAME(ARRAY_GET)(ARR,B));\
          QSORTNAME(ARRAY_SET)(ARR,B,tmp);\
        }
#endif

#ifndef GT_QSORT_ARR_VECSWAP
#define GT_QSORT_ARR_VECSWAP(ARR,A,B,N)\
        aidx = A;\
        bidx = B;\
        while ((N)-- > 0)\
        {\
          tmp = QSORTNAME(ARRAY_GET)(ARR,aidx);\
          QSORTNAME(ARRAY_SET)(ARR,aidx,QSORTNAME(ARRAY_GET)(ARR,bidx));\
          QSORTNAME(ARRAY_SET)(ARR,bidx,tmp);\
          aidx++;\
          bidx++;\
        }
#endif

static inline GtUword QSORTNAME(gt_inlined_qsort_arr_r_med3)
                     (GtUword a, GtUword b, GtUword c,
                      const GtDifferencecover *dcov,
                      GtLcptrace *lcptrace,
                      GtUword subbucketleft,
                      const QSORTNAME(Datatype) data)
{
  return
    QSORTNAME(qsortcmparr) (a, b, dcov, lcptrace, subbucketleft, data) < 0
      ? (QSORTNAME(qsortcmparr) (b, c, dcov, lcptrace,
                                 subbucketleft, data) < 0
        ? b
        : (QSORTNAME(qsortcmparr) (a, c, dcov, lcptrace,
                                   subbucketleft, data) < 0
          ? c : a))
          : (QSORTNAME(qsortcmparr) (b, c, dcov, lcptrace,
                                     subbucketleft, data) > 0
            ? b
            : (QSORTNAME(qsortcmparr) (a, c, dcov, lcptrace,
                                       subbucketleft, data) < 0
              ? a
              : c));
}

#ifndef GT_STACK_INTERVALARRAYTOBESORTED_DEFINED
typedef struct
{
  GtUword startindex,
                len;
} Intervalarrtobesorted;

GT_STACK_DECLARESTRUCT(Intervalarrtobesorted,32UL);
#define GT_STACK_INTERVALARRAYTOBESORTED_DEFINED
#endif

static void QSORTNAME(gt_inlinedarr_qsort_r) (QSORTNAME(Datatype) data,
                                              GtLcpvalues *sssplcpvalues,
                                              const GtDifferencecover *dcov,
                                              GtUword insertionsortthreshold,
                                              bool handlenotswapped,
                                              GtUword subbucketleft,
                                              GtUword len)
{
  GtUword tmp, pa, pb, pc, pd, pl, pm, pn, aidx, bidx, s,
          smallermaxlcp, greatermaxlcp;
  int r;
  bool swapped;
  GtStackIntervalarrtobesorted intervalstack;
  Intervalarrtobesorted current;
  GtLcptrace lcptrace;

  GT_STACK_INIT(&intervalstack,32UL);
  current.startindex = 0;
  current.len = len;
  GT_STACK_PUSH(&intervalstack,current);
  if (insertionsortthreshold <= 2UL)
  {
    insertionsortthreshold = 6UL;
  }
  while (!GT_STACK_ISEMPTY(&intervalstack))
  {
    swapped = false;
    current = GT_STACK_POP(&intervalstack);
    if (current.len <= insertionsortthreshold)
    {
      for (pm = current.startindex + 1;
           pm < current.startindex + current.len; pm++)
      {
        for (pl = pm; pl > current.startindex; pl--)
        {
          r = QSORTNAME(qsortcmparr) (pl - 1, pl, dcov, &lcptrace,
                                      subbucketleft, data);
          if (sssplcpvalues != NULL)
          {
            if (pl < pm && r > 0)
            {
              gt_lcptab_update(sssplcpvalues,subbucketleft,pl+1,
                               gt_lcptab_getvalue(sssplcpvalues,
                                                  subbucketleft,pl));
            }
            gt_lcptab_update(sssplcpvalues,subbucketleft,pl,
                             gt_differencecover_eval_lcp(&lcptrace, dcov));
          }
          if (r <= 0)
          {
            break;
          }
          GT_QSORT_ARR_SWAP (data, pl, pl - 1);
        }
      }
      continue;
    }
    pm = current.startindex + GT_DIV2 (current.len);
    if (current.len > 7UL)
    {
      pl = current.startindex;
      pn = current.startindex + current.len - 1;
      if (current.len > 40UL)
      {
        s = GT_DIV8 (current.len);
        pl = QSORTNAME(gt_inlined_qsort_arr_r_med3) (pl, pl + s,
                                                     pl + GT_MULT2 (s), dcov,
                                                     &lcptrace, subbucketleft,
                                                     data);
        gt_assert(pm >= s);
        pm = QSORTNAME(gt_inlined_qsort_arr_r_med3) (pm - s, pm,
                                                     pm + s, dcov,
                                                     &lcptrace, subbucketleft,
                                                     data);
        gt_assert(pn >= GT_MULT2(s));
        pn = QSORTNAME(gt_inlined_qsort_arr_r_med3) (pn - GT_MULT2 (s),
                                                     pn - s, pn, dcov,
                                                     &lcptrace, subbucketleft,
                                                     data);
      }
      pm = QSORTNAME(gt_inlined_qsort_arr_r_med3) (pl, pm, pn, dcov, &lcptrace,
                                                   subbucketleft,  data);
    }
    GT_QSORT_ARR_SWAP (data, current.startindex, pm);
    pa = pb = current.startindex + 1;
    pc = pd = current.startindex + current.len - 1;
    smallermaxlcp = greatermaxlcp = 0;
    while (1)
    {
      while (pb <= pc)
      {
        r = QSORTNAME(qsortcmparr) (pb, current.startindex, dcov, &lcptrace,
                                    subbucketleft, data);
        if (r > 0)
        {
          if (sssplcpvalues != NULL)
          {
            GtUword tmplcplen = gt_differencecover_eval_lcp(&lcptrace,dcov);
            GT_UPDATE_MAX(greatermaxlcp,tmplcplen);
          }
          break;
        }
        if (r == 0)
        {
          swapped = true;
          GT_QSORT_ARR_SWAP (data, pa, pb);
          pa++;
        } else
        {
          if (sssplcpvalues != NULL)
          {
            GtUword tmplcplen = gt_differencecover_eval_lcp(&lcptrace,dcov);
            GT_UPDATE_MAX(smallermaxlcp,tmplcplen);
          }
        }
        pb++;
      }
      while (pb <= pc)
      {
        r = QSORTNAME(qsortcmparr) (pc, current.startindex, dcov, &lcptrace,
                                    subbucketleft, data);
        if (r < 0)
        {
          if (sssplcpvalues != NULL)
          {
            GtUword tmplcplen = gt_differencecover_eval_lcp(&lcptrace,dcov);
            GT_UPDATE_MAX(smallermaxlcp,tmplcplen);
          }
          break;
        }
        if (r == 0)
        {
          swapped = true;
          GT_QSORT_ARR_SWAP (data, pc, pd);
          gt_assert(pd > 0);
          pd--;
        } else
        {
          if (sssplcpvalues != NULL)
          {
            GtUword tmplcplen = gt_differencecover_eval_lcp(&lcptrace,dcov);
            GT_UPDATE_MAX(greatermaxlcp,tmplcplen);
          }
        }
        gt_assert(pc > 0);
        pc--;
      }
      if (pb > pc)
      {
        break;
      }
      GT_QSORT_ARR_SWAP (data, pb, pc);
      swapped = true;
      pb++;
      gt_assert(pc > 0);
      pc--;
    }
    /* The following switch is not explained in the above mentioned
       paper and therefore we ignore it. */
    if (handlenotswapped && !swapped)
    {                                  /* Switch to insertion sort */
      gt_assert(current.len <= 40UL);
      for (pm = current.startindex + 1;
           pm < current.startindex + current.len; pm++)
      {
        for (pl = pm; pl > current.startindex; pl--)
        {
          r = QSORTNAME(qsortcmparr) (pl - 1, pl, dcov, &lcptrace,
                                      subbucketleft, data);
          if (r <= 0)
          {
            break;
          }
          GT_QSORT_ARR_SWAP (data, pl, pl - 1);
        }
      }
      continue;
    }
    pn = current.startindex + current.len;
    gt_assert(pa >= current.startindex);
    gt_assert(pb >= pa);
    s = MIN ((GtUword) (pa - current.startindex),
             (GtUword) (pb - pa));
    gt_assert(pb >= s);
    GT_QSORT_ARR_VECSWAP (data, current.startindex, pb - s, s);
    gt_assert(pd >= pc);
    gt_assert(pn > pd);
    s = MIN ((GtUword) (pd - pc), (GtUword) (pn - pd - 1));
    gt_assert(pn > s);
    GT_QSORT_ARR_VECSWAP (data, pb, pn - s, s);
    gt_assert(pb >= pa);
    if ((s = (GtUword) (pb - pa)) > 0)
    {
      if (sssplcpvalues != NULL)
      {
        /*
          left part has suffix with lcp up to length smallermaxlcp w.r.t.
          to the pivot. This lcp belongs to a suffix on the left
          which is at a minimum distance to the pivot and thus to an
          element in the final part of the left side.
        */
        gt_lcptab_update(sssplcpvalues,
                         subbucketleft,current.startindex + s,
                         smallermaxlcp);
      }
      if (s > 1UL)
      {
        current.len = s;
        GT_STACK_PUSH(&intervalstack,current);
      }
    }
    gt_assert(pd >= pc);
    if ((s = (GtUword) (pd - pc)) > 0)
    {
      if (sssplcpvalues != NULL)
      {
        /*
          right part has suffix with lcp up to length largermaxlcp w.r.t.
          to the pivot. This lcp belongs to a suffix on the right
          which is at a minimum distance to the pivot and thus to an
          element in the first part of the right side.
        */
        gt_assert(pn >= s);
        gt_lcptab_update(sssplcpvalues,subbucketleft,pn - s, greatermaxlcp);
      }
      if (s > 1UL)
      {
        gt_assert(pn >= s);
        current.startindex = pn - s;
        current.len = s;
        GT_STACK_PUSH(&intervalstack,current);
      }
    }
  }
  GT_STACK_DELETE(&intervalstack);
}

void gt_differencecover_sortunsortedbucket(GtSuffixsortspace *sssp,
                                           GtLcpvalues *sssplcpvalues,
                                           const GtDifferencecover *dcov,
                                           GtUword blisbl,
                                           GtUword width,
                                           GT_UNUSED GtUword depth)
{
  const GtUword bucketleftidx = gt_suffixsortspace_bucketleftidx_get(sssp);

  gt_assert(dcov != NULL &&
            depth >= (GtUword) dcov->vparam &&
            dcov->diff2pos != NULL &&
            width >= 2UL &&
            sssp != NULL &&
            blisbl >= bucketleftidx);
  /* blisbl = bucketleftindex + subbucketleft already contains
     bucketleftindex, therefore we cannot use
     gt_suffixsortspace_set or
     gt_suffixsortspace_get as these use the index
     bucketleftindex + subbucketleft + idx - partoffset
     = blisbl - partoffset + idx. Thus, instead we use the functions
     to directly access the suffix sortspace. */
  QSORTNAME(gt_inlinedarr_qsort_r) (sssp, sssplcpvalues, dcov, 6UL, false,
                                    blisbl - bucketleftidx, width);
}
