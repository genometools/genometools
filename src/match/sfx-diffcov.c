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
  unsigned long key,
                suffixstart;
} GtDcItventry;

typedef struct
{
  unsigned long blisbl, /* bucketleftindex + subbucketleft */
                width;
} GtDcPairsuffixptr;

GT_DECLAREARRAYSTRUCT(GtDcPairsuffixptr);

typedef GtDcPairsuffixptr Inl_Queueelem;

#include "queue-inline.h"

typedef struct
{
  unsigned long blisbl, /* bucketleftindex + subbucketleft */
                width,
                count,
                totalwidth,
                maxwidth,
                depth;
  bool defined;
} GtDcFirstwithnewdepth;

struct GtDifferencecover
{
  unsigned int vparam,
               logmod,
               size,
               vmodmask,
               hvalue,  /* not necessary */
               numofchars,
               prefixlength;
  unsigned long *coverrank_evaluated;
  Diffrank *coverrank;
  GtBitsequence *coverrank_bits;
  Diffvalue *diffvalues,  /* points to the difference cover */
            *diff2pos;    /* table d from BUR:KAER:2003 */
  size_t requiredspace;
  unsigned long totallength;
  GtLeftborder *leftborder; /* points to bcktab->leftborder */
  GtBcktab *bcktab;
  const GtEncseq *encseq;
  GtReadmode readmode;
  unsigned long samplesize, effectivesamplesize, maxsamplesize;
  const GtCodetype **multimappower;
  GtCodetype *filltable;
  GtEncseqReader *esr;
  unsigned long *shat,
                *inversesuftab,
                allocateditvinfo,
                currentqueuesize,
                maxqueuesize,
                currentdepth;
  GtCodetype maxcode;
  GtDcFirstwithnewdepth firstwithnewdepth;
  Inl_Queue *rangestobesorted;
  GtDcItventry *itvinfo;
  GtArrayGtDcPairsuffixptr firstgeneration;
  GtLcpvalues *samplelcpvalues;
  unsigned long firstgenerationtotalwidth,
                firstgenerationcount;
  GtLogger *logger;
  GtSuffixsortspace *sssp,
                    *sortedsample;
  unsigned long sortoffset;
};

/* Compute difference cover on the fly */

#define UScast(X) ((Diffvalue) X)
#define UCcast(X) ((Diffrank) X)

#include "tab-diffcover.h"

static unsigned long dc_suffixptrget(const GtDifferencecover *dcov,
                                     unsigned long idx)
{
  gt_suffixsortspace_nooffsets(dcov->sortedsample);
  return gt_suffixsortspace_get(dcov->sortedsample,0,idx);
}

static void dc_suffixptrset(const GtDifferencecover *dcov,
                            unsigned long idx,unsigned long value)
{
  gt_suffixsortspace_nooffsets(dcov->sortedsample);
  gt_suffixsortspace_set(dcov->sortedsample,0,idx,value);
}

static void dc_fillcoverrank(GtDifferencecover *dcov)
{
  unsigned int i;
  Diffrank j;
  const unsigned long step = GT_DIVV(dcov->totallength) + 1;
  unsigned long sum;

  dcov->coverrank_evaluated
    = gt_malloc(sizeof (*dcov->coverrank_evaluated) * dcov->vparam);
  GT_INITBITTAB(dcov->coverrank_bits,dcov->vparam);
  dcov->coverrank = gt_malloc(sizeof (*dcov->coverrank) * dcov->vparam);
  dcov->requiredspace += sizeof (*dcov->coverrank_evaluated) * dcov->vparam;
  gt_assert(dcov->size <= Diffrankmax);
  for (i=0; i<dcov->vparam; i++)
  {
    dcov->coverrank_evaluated[i] = ULONG_MAX; /* initialize as undefined */
    dcov->coverrank[i] = dcov->size;          /* initialize as undefined */
  }
  for (sum = 0, j=0; j<dcov->size; j++, sum += step)
  {
    Diffvalue d = dcov->diffvalues[j];
    dcov->coverrank_evaluated[d] = sum; /* jth value from difference cover
                                           gets rank j * step. */
    GT_SETIBIT(dcov->coverrank_bits,d);
    dcov->coverrank[d] = j;
  }
}

static bool dc_is_in_differencecover(const GtDifferencecover *dcov,
                                     unsigned long modpos)
{
  return GT_ISIBITSET(dcov->coverrank_bits,modpos) ? true : false;
}

static void dc_filldiff2pos(GtDifferencecover *dcov)
{
  Diffvalue *iptr, *jptr;

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
                                     unsigned long totallength)
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

void gt_differencecoversetsuffixsortspace(GtDifferencecover *dcov,
                                          GtSuffixsortspace *sssp)
{
  dcov->sssp = sssp;
}

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
  dcov->sssp = NULL;
  dcov->shat = NULL;
  for (dcov->logmod = 0;
       dcov->logmod < (unsigned int) (sizeof (differencecoversizes)/
                                      sizeof (differencecoversizes[0]));
       dcov->logmod++)
  {
    if (v == vparam)
    {
      dcov->size = differencecoversizes[dcov->logmod];
      dcov->diffvalues = differencecovertab + offset;
      found = true;
      break;
    }
    offset += differencecoversizes[dcov->logmod];
    v = GT_MULT2(v);
  }
  if (!found)
  {
    gt_free(dcov);
    return NULL;
  }
  dcov->maxsamplesize = (unsigned long) (GT_DIVV(dcov->totallength) + 1) *
                                         dcov->size;
  if (outerprefixlength == 0)
  {
    dcov->prefixlength = 0;
  } else
  {
    dcov->prefixlength
      = gt_recommendedprefixlength(dcov->numofchars,
                                   dcov->maxsamplesize,
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
  dcov->esr = NULL;
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

size_t gt_differencecover_requiredspace(const GtDifferencecover *dcov)
{
  return dcov->requiredspace;
}

/* The following is the \delta function from BUR:KAER:2003. */

static unsigned int dc_differencecover_offset(const GtDifferencecover *dcov,
                                              unsigned long pos1,
                                              unsigned long pos2)
{
  return (unsigned int) GT_MODV(dcov->diff2pos[GT_MODV(pos2 - pos1)] - pos1);
}

void gt_differencecover_delete(GtDifferencecover *dcov)
{
  if (dcov != NULL)
  {
    gt_assert(dcov->bcktab == NULL);
    gt_assert(dcov->sortedsample == NULL);
    gt_assert(dcov->filltable == NULL);
    gt_assert(dcov->multimappower == NULL);
    gt_assert(dcov->esr == NULL);

    gt_free(dcov->coverrank_evaluated);
    dcov->coverrank_evaluated = NULL;
    gt_free(dcov->coverrank_bits);
    dcov->coverrank_bits = NULL;
    gt_free(dcov->coverrank);
    dcov->coverrank = NULL;
    gt_free(dcov->diff2pos);
    dcov->diff2pos = NULL;
    gt_free(dcov->inversesuftab);
    dcov->inversesuftab = NULL;
    gt_free(dcov->shat);
    dcov->shat = NULL;
    gt_free(dcov);
  }
}

/* the following implements the \mu function from BAE:KAER:2003 */

static unsigned long dc_differencecover_packsamplepos(
                                                 const GtDifferencecover *dcov,
                                                 unsigned long pos)
{
  gt_assert(dc_is_in_differencecover(dcov,GT_MODV(pos)));
  return dcov->coverrank_evaluated[GT_MODV(pos)] + GT_DIVV(pos);
}

GT_DECLAREARRAYSTRUCT(Codeatposition);

static unsigned long dc_derivespecialcodesonthefly(GtDifferencecover *dcov,
                                                   const GtArrayCodeatposition
                                                           *codelist)
{
  unsigned int prefixindex, unitsnotspecial;
  Enumcodeatposition *ecp;
  Specialcontext specialcontext;
  unsigned long countderived = 0, pos, sampleindex;
  GtCodetype code;

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
        gt_assert(specialcontext.position >= (unsigned long) prefixindex);
        pos = (unsigned long) (specialcontext.position - prefixindex);
        if (dc_is_in_differencecover(dcov,GT_MODV(pos)))
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
                                             dcov->esr,
                                             dcov->multimappower,
                                             pos,
                                             dcov->prefixlength);
          if (codelist != NULL)
          {
            gt_assert((GtCodetype) codelist->spaceCodeatposition[
                                           countderived].code == code);
          }
          /*
          printf("%u %lu\n",prefixindex,
                            (unsigned long)
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
  unsigned long pos;
  unsigned int modvalue;
  Diffvalue *diffptr, *afterend;
  unsigned long idx;
  GtBitsequence *sampleidxused = NULL;

  GT_INITBITTAB(sampleidxused,dcov->maxsamplesize);
  diffptr = dcov->diffvalues;
  afterend = dcov->diffvalues + dcov->size;
  for (pos = 0, modvalue = 0; pos <= dcov->totallength; pos++)
  {
    gt_assert((unsigned long) modvalue == GT_MODV(pos));
    gt_assert(diffptr == afterend || *diffptr >= (Diffvalue) modvalue);
    if (diffptr < afterend && (Diffvalue) modvalue == *diffptr)
    {
      idx = dc_differencecover_packsamplepos(dcov,pos);
      gt_assert(sampleidxused != NULL);
      if (GT_ISIBITSET(sampleidxused,idx))
      {
        fprintf(stderr,"sample index %lu for pos %lu already used before\n",
                       idx,(unsigned long) pos);
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
                                 unsigned long pos,
                                 unsigned long sampleindex)
{
  gt_assert (sampleindex < dcov->samplesize);
  dcov->inversesuftab[dc_differencecover_packsamplepos(dcov,pos)] = sampleindex;
}

static unsigned long dc_inversesuftab_get(const GtDifferencecover *dcov,
                                          unsigned long pos)
{
  return dcov->inversesuftab[dc_differencecover_packsamplepos(dcov,pos)];
}

static void dc_initinversesuftabnonspecials(GtDifferencecover *dcov)
{
  unsigned long sampleindex, pos;

  for (sampleindex=0; sampleindex < dcov->effectivesamplesize; sampleindex++)
  {
    pos = dc_suffixptrget(dcov,sampleindex);
    dc_inversesuftab_set(dcov,pos,sampleindex);
  }
}

static unsigned long dc_insertfullspecialrangesample(GtDifferencecover *dcov,
                                                     unsigned long specialidx,
                                                     unsigned long leftpos,
                                                     unsigned long rightpos)
{
  unsigned long pos;

  gt_assert(leftpos < rightpos);
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
      unsigned long revpos;

      gt_assert(pos < dcov->totallength);
      revpos = GT_REVERSEPOS(dcov->totallength,pos);
      if (dc_is_in_differencecover(dcov,GT_MODV(revpos)))
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
      if (dc_is_in_differencecover(dcov,GT_MODV(pos)))
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
  dcov->inversesuftab = gt_malloc(sizeof (*dcov->inversesuftab) *
                                  dcov->maxsamplesize);
  dcov->requiredspace += sizeof (*dcov->inversesuftab) * dcov->maxsamplesize;
  if (gt_encseq_has_specialranges(dcov->encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;
    unsigned long specialidx;

    sri = gt_specialrangeiterator_new(dcov->encseq,
                                      GT_ISDIRREVERSE(dcov->readmode)
                                      ? false : true);
    specialidx = dcov->effectivesamplesize;
    while (gt_specialrangeiterator_next(sri,&range))
    {
      /*
      printf("specialrange %lu %lu\n",(unsigned long) range.leftpos,
                                      (unsigned long) range.rightpos);
      */
      specialidx = dc_insertfullspecialrangesample(dcov,
                                                   specialidx,
                                                   range.start,
                                                   range.end);
    }
    gt_specialrangeiterator_delete(sri);
    sri = NULL;
  }
  if (dc_is_in_differencecover(dcov,GT_MODV(dcov->totallength)))
  {
    gt_assert(dcov->samplesize > 0);
    dc_inversesuftab_set(dcov,dcov->totallength,dcov->samplesize-1);
  }
}

static void dc_updatewidth (GtDifferencecover *dcov,unsigned long width,
                            unsigned int depth)
{
  if (width > 1UL)
  {
    dcov->firstgenerationtotalwidth += width;
    dcov->firstgenerationcount++;
    if (dcov->allocateditvinfo < width)
    {
      dcov->allocateditvinfo = width;
    }
    if (dcov->currentdepth == 0)
    {
      dcov->currentdepth = (unsigned long) depth;
    } else
    {
      gt_assert(dcov->currentdepth == (unsigned long) depth);
    }
  }
}

static void dc_initinversesuftabnonspecialsadjust(GtDifferencecover *dcov)
{
  GtCodetype code;
  GtBucketspecification bucketspec;
  unsigned long idx = 0;
  const GtCodetype mincode = 0;
  unsigned int rightchar = 0; /* as mincode is 0 */

  for (code = mincode; code <= dcov->maxcode; code++)
  {
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
      dc_inversesuftab_set(dcov,dc_suffixptrget(dcov,idx),bucketspec.left);
    }
  }
  for (/* Nothing */; idx < dcov->effectivesamplesize; idx++)
  {
    dc_inversesuftab_set(dcov,dc_suffixptrget(dcov,idx),idx);
  }
}

static void dc_anchorleftmost(GtDifferencecover *dcov,
                              unsigned long blisbl,
                              unsigned long width)
{
  unsigned long idx;

  for (idx = blisbl; idx < blisbl + width; idx++)
  {
    dc_inversesuftab_set(dcov,dc_suffixptrget(dcov,idx),blisbl);
  }
}

static void dc_showintervalsizes(unsigned long count,
                                 unsigned long totalwidth,
                                 unsigned long effectivesamplesize,
                                 unsigned long maxwidth,
                                 unsigned long depth,
                                 GtLogger *logger)
{
  gt_logger_log(logger,
              "level %lu"
              ": (intervals=%lu,total=%lu,avg=%.2f,%.2f%% of all,maxwidth=%lu)",
              depth,
              count,
              totalwidth,
              (double) totalwidth/count,
              100.0 * (double) totalwidth/effectivesamplesize,
              maxwidth);
}

static void dc_processunsortedrange(GtDifferencecover *dcov,
                                    unsigned long blisbl,
                                    unsigned long width,
                                    unsigned long depth)
{
  GtDcPairsuffixptr *pairelem;

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
    gt_logger_log(dcov->logger,"enter new level %lu",depth);
    dcov->firstwithnewdepth.blisbl = blisbl;
    dcov->firstwithnewdepth.width = width;
    dcov->firstwithnewdepth.depth = depth;
    dcov->firstwithnewdepth.count = 1UL;
    dcov->firstwithnewdepth.totalwidth = width;
    dcov->firstwithnewdepth.maxwidth = width;
  }
  pairelem = gt_malloc(sizeof (*pairelem));
  pairelem->blisbl = blisbl;
  pairelem->width = width;
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
                                            unsigned long blisbl,
                                            unsigned long width,
                                            unsigned long lcpvalue)
{
  unsigned long idx;

  for (idx = blisbl+1; idx < blisbl + width; idx++)
  {
    gt_lcptab_update(samplelcpvalues,0,idx,lcpvalue);
  }
}

static void dc_sortsuffixesonthislevel(GtDifferencecover *dcov,
                                       unsigned long blisbl,
                                       unsigned long width)
{
  unsigned long idx, rangestart, startpos;

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
  printf("new interval of width %lu\n",width);
  */
  for (idx=0; idx<width; idx++)
  {
    startpos = dc_suffixptrget(dcov,blisbl+idx);
    dcov->itvinfo[idx].suffixstart = startpos;
    dcov->itvinfo[idx].key
      = dc_inversesuftab_get(dcov,startpos + dcov->currentdepth);
    /*
    printf("key=%lu,startpos=%lu,currentdepth=%lu\n",
                    dcov->itvinfo[idx].key,
                    startpos,
                    dcov->currentdepth);
    */
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
        unsigned long currentsuftabentry
          = dc_suffixptrget(dcov,blisbl+rangestart);
        dc_inversesuftab_set(dcov,currentsuftabentry,blisbl+rangestart);
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
    unsigned long currentsuftabentry
      = dc_suffixptrget(dcov,blisbl+rangestart);
    dc_inversesuftab_set(dcov,currentsuftabentry,blisbl+rangestart);
  }
}

static void dc_bcktab2firstlevelintervals(GtDifferencecover *dcov)
{
  GtCodetype code;
  unsigned int rightchar;
  GtBucketspecification bucketspec;
  const GtCodetype mincode = 0;

  printf("# maxbucketsize=%lu\n",dcov->allocateditvinfo);
  rightchar = (unsigned int) (mincode % dcov->numofchars);
  for (code = 0; code <= dcov->maxcode; code++)
  {
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
                                unsigned long blisbl,
                                unsigned long width,
                                GT_UNUSED unsigned long depth)
{
  GtDifferencecover *dcov = (GtDifferencecover *) voiddcov;
  GtDcPairsuffixptr *ptr;

  gt_assert(dcov->sssp == NULL);
  gt_assert(depth >= (unsigned long) dcov->vparam);
  dc_updatewidth (dcov,width,dcov->vparam);
  GT_GETNEXTFREEINARRAY(ptr,&dcov->firstgeneration,GtDcPairsuffixptr,1024);
  ptr->blisbl = blisbl;
  ptr->width = width;
}

#ifdef  QSORTNAME
#undef  QSORTNAME
#endif

#define QSORTNAME(NAME) dc_##NAME

#define dc_ARRAY_GET(ARR,RELIDX)\
        gt_suffixsortspace_get(data->sssp,data->sortoffset,RELIDX)

#define dc_ARRAY_SET(ARR,RELIDX,VALUE)\
        gt_suffixsortspace_set(data->sssp,data->sortoffset,RELIDX,VALUE)

typedef GtDifferencecover * QSORTNAME(Datatype);

int gt_differencecover_compare (const GtDifferencecover *dcov,
                                unsigned long suffixpos1,
                                unsigned long suffixpos2)
{
  unsigned int offset;
  unsigned long idx1, idx2;

  gt_assert(suffixpos1 < dcov->totallength);
  gt_assert(suffixpos2 < dcov->totallength);
  offset = dc_differencecover_offset(dcov,suffixpos1,suffixpos2);
  idx1 = dc_inversesuftab_get(dcov,suffixpos1 + offset);
  idx2 = dc_inversesuftab_get(dcov,suffixpos2 + offset);
  if (idx1 < idx2)
  {
    return -1;
  }
  if (idx1 > idx2)
  {
    return 1;
  }
  gt_assert(false);
  return 0;
}

static int QSORTNAME(qsortcmparr) (
                  GT_UNUSED const void *subbucket,
                  unsigned long a,
                  unsigned long b,
                  const QSORTNAME(Datatype) data)
{
  unsigned long suffixpos1, suffixpos2;

  gt_assert(data->sssp != NULL);
  suffixpos1 = dc_ARRAY_GET(NULL,a);
  suffixpos2 = dc_ARRAY_GET(NULL,b);
  return gt_differencecover_compare (data, suffixpos1, suffixpos2);
}

typedef void * QSORTNAME(Sorttype);

#include "qsort-array.gen"

void gt_differencecover_sortunsortedbucket(void *data,
                                           unsigned long blisbl,
                                           unsigned long width,
                                           GT_UNUSED unsigned long depth)
{
  GtDifferencecover *dcov = (GtDifferencecover *) data;

  gt_assert(depth >= (unsigned long) dcov->vparam);
  gt_assert(dcov->diff2pos != NULL);
  gt_assert(width >= 2UL);
  gt_assert(dcov->sssp != NULL);
  gt_assert(blisbl >= gt_suffixsortspace_bucketleftidx_get(dcov->sssp));
  /* blisbl = bucketleftindex + subbucketleft already contains
     bucketleftindex, therefore we cannot use
     gt_suffixsortspace_set or
     gt_suffixsortspace_get as these use the index
     bucketleftindex + subbucketleft + idx - partoffset
     = blisbl - partoffset + idx. Thus, instead we use the functions
     to directly access the suffix sortspace. */
  dcov->sortoffset = blisbl - gt_suffixsortspace_bucketleftidx_get(dcov->sssp);
  QSORTNAME(gt_inlinedarr_qsort_r) (6UL, false, NULL,width,data);
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
    GtDcPairsuffixptr *thispairptr;
    thispairptr = (GtDcPairsuffixptr*) gt_inl_queue_get(dcov->rangestobesorted);
    gt_assert(dcov->currentqueuesize > 0);
    dcov->currentqueuesize--;
    dc_sortsuffixesonthislevel(dcov,thispairptr->blisbl,
                               thispairptr->width);
    gt_free(thispairptr);
  }
  gt_logger_log(dcov->logger,"maxqueuesize=%lu",dcov->maxqueuesize);
  gt_free(dcov->itvinfo);
  dcov->itvinfo = NULL;
  gt_inl_queue_delete(dcov->rangestobesorted);
  dcov->rangestobesorted = NULL;
}

static void dc_init_sfxstrategy_for_sample(Sfxstrategy *sfxstrategy,
                                           const Sfxstrategy *mainsfxstrategy,
                                           bool bitwise_cmp_ok,
                                           unsigned long effectivesamplesize,
                                           unsigned long totallength,
                                           GtLogger *logger)
{
  if (mainsfxstrategy != NULL)
  {
    double sampledproportion = (double) effectivesamplesize/totallength;

    *sfxstrategy = *mainsfxstrategy;
#define SETMAXCOUNT(COMP)\
    if (mainsfxstrategy->COMP >= 1UL)\
    {\
      sfxstrategy->COMP = MAX(1UL,mainsfxstrategy->COMP * sampledproportion);\
    }
    SETMAXCOUNT(maxcountingsort);
    SETMAXCOUNT(maxbltriesort);
    SETMAXCOUNT(maxinsertionsort);
  } else
  {
    defaultsfxstrategy(sfxstrategy,bitwise_cmp_ok);
  }
  gt_logger_log(logger,"samplesort.maxinsertionsort=%lu",
                sfxstrategy->maxinsertionsort);
  gt_logger_log(logger,"samplesort.maxbltriesort=%lu",
                sfxstrategy->maxbltriesort);
  gt_logger_log(logger,"samplesort.maxcountingsort=%lu",
                sfxstrategy->maxcountingsort);
}

static void dc_fill_lcpvalues(GtDifferencecover *dcov)
{

  unsigned long suffix, idx;
  unsigned long j, kvalue, lcpvalue, start0, start1;
  unsigned int svalue;
  GtUchar cc1, cc2;

  dcov->shat = gt_malloc(sizeof (*dcov->shat) * dcov->samplesize);
  for (j = 0; j < dcov->samplesize; j++)
  {
    dcov->shat[idx] = ULONG_MAX;
  }
  for (j = 0; j < dcov->effectivesamplesize; j++)
  {
    suffix = gt_suffixsortspace_get(dcov->sortedsample,0,j);
    gt_assert(dc_is_in_differencecover(dcov,GT_MODV(suffix)));
    idx = dcov->size * GT_DIVV(suffix) + dcov->coverrank[GT_MODV(suffix)];
    gt_assert(idx < dcov->samplesize);
    dcov->shat[idx] = j;
  }
  for (svalue = 0; svalue < dcov->size; svalue++)
  {
    kvalue = (unsigned long) svalue;
    lcpvalue = 0;
    while (kvalue < dcov->samplesize)
    {
      suffix = dcov->shat[kvalue];
      if (suffix > 0 && suffix < dcov->effectivesamplesize)
      {
        start0 = gt_suffixsortspace_get(dcov->sortedsample,0,suffix-1);
        start1 = gt_suffixsortspace_get(dcov->sortedsample,0,suffix);
        while (start0 + lcpvalue < dcov->totallength &&
               start1 + lcpvalue < dcov->totallength)
        {
          cc1 = gt_encseq_get_encoded_char(dcov->encseq,start0+lcpvalue,
                                           dcov->readmode);
          cc2 = gt_encseq_get_encoded_char(dcov->encseq,start1+lcpvalue,
                                           dcov->readmode);
          if (ISSPECIAL(cc1) || ISSPECIAL(cc2) || cc1 != cc2)
          {
            break;
          }
          lcpvalue++;
        }
        gt_lcptab_update(dcov->samplelcpvalues,0,suffix,lcpvalue);
        lcpvalue = lcpvalue > (unsigned long) dcov->vparam
                   ? (lcpvalue - (unsigned long) dcov->vparam) : 0;
      }
      kvalue += dcov->size;
    }
  }
}

static void dc_differencecover_sortsample(GtDifferencecover *dcov,
                                          GtOutlcpinfo *outlcpinfosample,
                                          const Sfxstrategy *mainsfxstrategy,
                                          bool withcheck)
{
  unsigned long pos, sampleindex, posinserted, fullspecials = 0, specials = 0;
  unsigned int modvalue, unitsnotspecial;
  Diffvalue *diffptr, *afterend;
  GtCodetype code;
  GtArrayCodeatposition codelist;
  Codeatposition *codeptr;

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
  dcov->esr = gt_encseq_create_reader_with_readmode(dcov->encseq,
                                                    dcov->readmode,
                                                    0);
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
                                           dcov->esr,
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
            gt_assert(code <= (GtCodetype) MAXCODEVALUE);
            codeptr->code = (unsigned int) code;
            gt_assert(unitsnotspecial <= (unsigned int) MAXPREFIXLENGTH);
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
  (void) gt_bcktab_leftborderpartialsums(NULL,NULL,dcov->bcktab);
  gt_logger_log(dcov->logger,
              "%lu positions are sampled (%.2f%%) pl=%u",
              dcov->samplesize,
              100.0 * (double) dcov->samplesize/(dcov->totallength+1),
              dcov->prefixlength);
  gt_logger_log(dcov->logger,"specials=%lu, fullspecials=%lu",
                specials,fullspecials);
  if (withcheck)
  {
    qsort(codelist.spaceCodeatposition,
          (size_t) codelist.nextfreeCodeatposition,
          sizeof (*codelist.spaceCodeatposition),dc_compareCodeatpositon);
  }
  dcov->sortedsample = gt_suffixsortspace_new(dcov->effectivesamplesize,
                                              dcov->totallength,
                                              false,
                                              dcov->logger);
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
                                         dcov->esr,
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
  gt_encseq_reader_delete(dcov->esr);
  dcov->esr = NULL;
  gt_assert(posinserted == dcov->effectivesamplesize);
  if (withcheck)
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
                           (unsigned long) dcov->prefixlength);
  }
  gt_Outlcpinfo_reinit(outlcpinfosample,dcov->numofchars,dcov->prefixlength,
                       dcov->effectivesamplesize);
  if (outlcpinfosample != NULL)
  {
    dcov->requiredspace += gt_Outlcpinfo_size(outlcpinfosample);
    dcov->samplelcpvalues = gt_Outlcpinfo_lcpvalues_ref(outlcpinfosample);
  }
  if (dcov->vparam == dcov->prefixlength)
  {
    dc_initinversesuftabspecials(dcov);
    dc_initinversesuftabnonspecialsadjust(dcov);
    dc_bcktab2firstlevelintervals(dcov);
  } else
  {
    unsigned long long bucketiterstep = 0;
    Sfxstrategy sfxstrategy;

    gt_assert (dcov->vparam > dcov->prefixlength);
    dc_init_sfxstrategy_for_sample(&sfxstrategy,
                                   mainsfxstrategy,
                                   gt_encseq_bitwise_cmp_ok(dcov->encseq)
                                     ? false : true,
                                   dcov->effectivesamplesize,
                                   dcov->totallength,
                                   dcov->logger);
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
    if (withcheck)
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
                             (unsigned long) dcov->vparam);
    }
  }
  gt_bcktab_delete(dcov->bcktab);
  dcov->bcktab = NULL;
  dc_sortremainingsamples(dcov);
  if (withcheck)
  {
    unsigned long idx;

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
      unsigned long idx2 = dc_inversesuftab_get(dcov,dc_suffixptrget(dcov,idx));
      if (idx != idx2)
      {
        fprintf(stderr,"idx = %lu != %lu = idx2\n",idx,idx2);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
    }
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
    dc_fill_lcpvalues(dcov);
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
  gt_assert(dcov->diff2pos == NULL);
  dc_filldiff2pos(dcov);
}

static void dc_differencecover_sortsample0(GtDifferencecover *dcov,
                                           GtOutlcpinfo *outlcpinfosample,
                                           const Sfxstrategy *mainsfxstrategy,
                                           bool withcheck)
{
  unsigned long pos, posinserted, fullspecials = 0;
  unsigned int modvalue;
  Diffvalue *diffptr, *afterend;
  Sfxstrategy sfxstrategy;
  GtUchar cc;

  dcov->samplesize = 0;
  dcov->bcktab = NULL;
  dcov->multimappower = NULL;
  dcov->maxcode = 0;
  dcov->esr = NULL;
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
              "%lu positions are sampled (%.2f) pl=%u",
              dcov->samplesize,
              100.0 * (double) dcov->samplesize/(dcov->totallength+1),
              dcov->prefixlength);
  gt_logger_log(dcov->logger,"fullspecials=%lu",fullspecials);
  dcov->sortedsample = gt_suffixsortspace_new(dcov->effectivesamplesize,
                                              dcov->totallength,
                                              false,
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
  if (withcheck)
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
                           (unsigned long) dcov->vparam);
  }
  dc_sortremainingsamples(dcov);
  if (withcheck)
  {
    unsigned long idx;

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
      unsigned long idx2 = dc_inversesuftab_get(dcov,dc_suffixptrget(dcov,idx));
      if (idx != idx2)
      {
        fprintf(stderr,"idx = %lu != %lu = idx2\n",idx,idx2);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
    }
  }
  gt_suffixsortspace_delete(dcov->sortedsample,false);
  dcov->sortedsample = NULL;
  gt_assert(dcov->diff2pos == NULL);
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
      gt_logger_log(logger,"presorting sample suffixes according to "
                           "difference cover modulo %u",vparam);
      (prefixlength > 0 ? dc_differencecover_sortsample
                        : dc_differencecover_sortsample0)
                          (dcov,outlcpinfosample,sfxstrategy,
                           sfxstrategy->dccheck);
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

  printf("sizeof (differencecovertab)=%lu\n",
          (unsigned long) sizeof (differencecovertab));
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
    dc_differencecover_sortsample(dcov,NULL,NULL,withcheck);
    gt_differencecover_delete(dcov);
  }
  printf("# %u difference covers checked\n",
          (unsigned int) (logmod - startlogmod));
}
