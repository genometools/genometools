/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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
#include "stamp.h"

typedef unsigned char Diffrank;
#define Diffrankmax ((Diffrank) 255)
typedef unsigned short Diffvalue;
#define Diffvaluemax ((Diffvalue) 65535)

#define MODV(VAL) ((VAL) & dcov->vmodmask)
#define DIVV(VAL) ((VAL) >> dcov->logmod)

typedef struct
{
  unsigned long key,
                suffixstart;
} Itventry;

typedef struct
{
  unsigned long subbucketleft,
                width;
} Pairsuffixptr;

GT_DECLAREARRAYSTRUCT(Pairsuffixptr);

typedef Pairsuffixptr Inl_Queueelem;

#include "queue-inline.h"

typedef struct
{
  unsigned long subbucketleft,
                width,
                count,
                totalwidth,
                maxwidth,
                depth;
  bool defined;
} Firstwithnewdepth;

struct Differencecover
{
  unsigned int vparam,
               logmod,
               size,
               vmodmask,
               hvalue,  /* not necessary */
               numofchars,
               prefixlength;
  Diffrank *coverrank;
  Diffvalue *diffvalues, *diff2pos;
  size_t requiredspace;
  unsigned long totallength,
                *leftborder; /* points to bcktab->leftborder */
  Bcktab *bcktab;
  const GtEncseq *encseq;
  GtReadmode readmode;
  unsigned long samplesize, effectivesamplesize, maxsamplesize;
  const GtCodetype **multimappower;
  GtCodetype *filltable;
  GtEncseqReader *esr;
  unsigned long *inversesuftab,
                allocateditvinfo,
                currentqueuesize,
                maxqueuesize,
                currentdepth;
  GtCodetype maxcode;
  Firstwithnewdepth firstwithnewdepth;
  Inl_Queue *rangestobesorted;
  Itventry *itvinfo;
  GtArrayPairsuffixptr firstgeneration;
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

static unsigned long suffixptrgetdcov(const Differencecover *dcov,
                                      unsigned long idx)
{
  return gt_suffixsortspace_getdirect(dcov->sortedsample,idx);
}

static void suffixptrsetdcov(const Differencecover *dcov,
                             unsigned long idx,unsigned long value)
{
  gt_suffixsortspace_setdirect(dcov->sortedsample,idx,value);
}

static void fillcoverrank(Differencecover *dcov)
{
  unsigned int i;
  Diffrank j;

  dcov->coverrank = gt_malloc(sizeof (*dcov->coverrank) * dcov->vparam);
  dcov->requiredspace += sizeof (*dcov->coverrank) * dcov->vparam;
  gt_assert(dcov->size <= Diffrankmax);
  for (i=0; i<dcov->vparam; i++)
  {
    dcov->coverrank[i] = dcov->size;
  }
  for (j=0; j<dcov->size; j++)
  {
    dcov->coverrank[dcov->diffvalues[j]] = j;
  }
}

static bool checkifindifferencecover(const Differencecover *dcov,
                                     unsigned long modpos)
{
  return dcov->coverrank[modpos] == dcov->size ? false : true;
}

static void filldiff2pos(Differencecover *dcov)
{
  Diffvalue *iptr, *jptr;

  dcov->diff2pos = gt_malloc(sizeof (*dcov->diff2pos) * dcov->vparam);
  dcov->requiredspace += sizeof (*dcov->diff2pos) * dcov->vparam;
  for (iptr=dcov->diffvalues + dcov->size - 1; iptr>=dcov->diffvalues; iptr--)
  {
    for (jptr=dcov->diffvalues; jptr<dcov->diffvalues + dcov->size; jptr++)
    {
      dcov->diff2pos[MODV(*jptr - *iptr)] = *iptr;
    }
  }
}

#ifdef WITHcomputehvalue

/* XXX: the following function is currently not used */

static unsigned int computehvalue(const Differencecover *dcov,
                                  unsigned long totallength)
{
  Diffvalue next;
  unsigned int h, nmodv = MODV(totallength);

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

void gt_differencecoversetsuffixsortspace(Differencecover *dcov,
                                          GtSuffixsortspace *sssp)
{
  dcov->sssp = sssp;
}

static int gt_differencecover_vparamverify(const Differencecover *dcov,
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

Differencecover *gt_differencecover_new(unsigned int vparam,
                                        const GtEncseq *encseq,
                                        GtReadmode readmode,
                                        unsigned int outerprefixlength,
                                        GtLogger *logger)
{
  unsigned int offset = 0, v = 1U;
  Differencecover *dcov;
  bool found = false;

  dcov = gt_malloc(sizeof (*dcov));
  dcov->requiredspace = sizeof (*dcov);
  dcov->numofchars = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
  dcov->totallength = gt_encseq_total_length(encseq);
  dcov->logger = logger;
  dcov->sssp = NULL;
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
  dcov->maxsamplesize = (unsigned long) (DIVV(dcov->totallength) + 1) *
                                         dcov->size;
  if (outerprefixlength == 0)
  {
    dcov->prefixlength = 0;
  } else
  {
    dcov->prefixlength = gt_recommendedprefixlength(dcov->numofchars,
                                                    dcov->maxsamplesize);
    if (outerprefixlength > 0 && dcov->prefixlength > outerprefixlength)
    {
      dcov->prefixlength = outerprefixlength;
    }
    gt_assert(dcov->prefixlength > 0);
  }
  dcov->vparam = 1U << (dcov->logmod);
  dcov->vmodmask = dcov->vparam-1;
#ifdef WITHcomputehvalue
  dcov->hvalue = computehvalue(dcov,totallength);
#endif
  dcov->encseq = encseq;
  dcov->readmode = readmode;
  dcov->bcktab = NULL;
  dcov->sortedsample = NULL;
  dcov->filltable = NULL;
  dcov->multimappower = NULL;
  fillcoverrank(dcov);
  dcov->diff2pos = NULL; /* this is later initialized */
  dcov->esr = NULL;
  dcov->allocateditvinfo = 0;
  dcov->itvinfo = NULL;
  dcov->currentdepth = 0;
  dcov->firstwithnewdepth.defined = false;
  dcov->firstwithnewdepth.depth = 0;
  dcov->firstwithnewdepth.totalwidth = 0;
  dcov->firstwithnewdepth.count = 0;
  dcov->firstwithnewdepth.subbucketleft = 0;
  dcov->firstwithnewdepth.width = 0;
  dcov->firstwithnewdepth.maxwidth = 0;
  dcov->currentqueuesize = 0;
  dcov->maxqueuesize = 0;
  dcov->inversesuftab = NULL;
  dcov->firstgenerationtotalwidth = 0;
  dcov->firstgenerationcount = 0;
  GT_INITARRAY(&dcov->firstgeneration,Pairsuffixptr);
  return dcov;
}

size_t gt_differencecover_requiredspace(const Differencecover *dcov)
{
  return dcov->requiredspace;
}

static unsigned int differencecover_offset(const Differencecover *dcov,
                                           unsigned long pos1,
                                           unsigned long pos2)
{
  return (unsigned int) MODV(dcov->diff2pos[MODV(pos2-pos1)] - pos1);
}

void gt_differencecover_delete(Differencecover *dcov)
{
  if (dcov != NULL)
  {
    gt_assert(dcov->bcktab == NULL);
    gt_assert(dcov->sortedsample == NULL);
    gt_assert(dcov->filltable == NULL);
    gt_assert(dcov->multimappower == NULL);
    gt_assert(dcov->esr == NULL);

    gt_free(dcov->coverrank);
    dcov->coverrank = NULL;
    gt_free(dcov->diff2pos);
    dcov->diff2pos = NULL;
    gt_free(dcov->inversesuftab);
    dcov->inversesuftab = NULL;
    gt_free(dcov);
  }
}

static unsigned long differencecover_packsamplepos(const Differencecover *dcov,
                                                   unsigned long pos)
{
  unsigned long result;

  result =  dcov->coverrank[MODV(pos)] * (DIVV(dcov->totallength) + 1) +
            DIVV(pos);
  gt_assert(result < dcov->maxsamplesize);
  return result;
}

GT_DECLAREARRAYSTRUCT(Codeatposition);

static unsigned long dcov_derivespecialcodesonthefly(Differencecover *dcov,
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
        if (checkifindifferencecover(dcov,MODV(pos)))
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
          gt_updatebckspecials(dcov->bcktab,code,dcov->numofchars,prefixindex);
          gt_assert(code > 0);
          sampleindex = --dcov->leftborder[code];
          gt_assert(sampleindex < dcov->effectivesamplesize);
          suffixptrsetdcov(dcov,sampleindex,pos);
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

static int compareCodeatpositon(const void *vala,const void *valb)
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

static void validate_samplepositons(const Differencecover *dcov)
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
    gt_assert((unsigned long) modvalue == MODV(pos));
    gt_assert(diffptr == afterend || *diffptr >= (Diffvalue) modvalue);
    if (diffptr < afterend && (Diffvalue) modvalue == *diffptr)
    {
      idx = differencecover_packsamplepos(dcov,pos);
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

static void inversesuftab_set(Differencecover *dcov,
                              unsigned long pos,
                              unsigned long sampleindex)
{
  gt_assert (sampleindex < dcov->samplesize);
  dcov->inversesuftab[differencecover_packsamplepos(dcov,pos)] = sampleindex;
}

static unsigned long inversesuftab_get(const Differencecover *dcov,
                                       unsigned long pos)
{
  return dcov->inversesuftab[differencecover_packsamplepos(dcov,pos)];
}

static void dc_initinversesuftabnonspecials(Differencecover *dcov)
{
  unsigned long sampleindex, pos;

  for (sampleindex=0; sampleindex < dcov->effectivesamplesize; sampleindex++)
  {
    pos = suffixptrgetdcov(dcov,sampleindex);
    inversesuftab_set(dcov,pos,sampleindex);
  }
}

static unsigned long insertfullspecialrangesample(Differencecover *dcov,
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
      if (checkifindifferencecover(dcov,MODV(revpos)))
      {
        inversesuftab_set(dcov,revpos,specialidx);
        specialidx++;
      }
      if (pos == leftpos)
      {
        break;
      }
      pos--;
    } else
    {
      if (checkifindifferencecover(dcov,MODV(pos)))
      {
        inversesuftab_set(dcov,pos,specialidx);
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

static void dc_initinversesuftabspecials(Differencecover *dcov)
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
      specialidx = insertfullspecialrangesample(dcov,
                                                specialidx,
                                                range.start,
                                                range.end);
    }
    gt_specialrangeiterator_delete(sri);
    sri = NULL;
  }
  if (checkifindifferencecover(dcov,MODV(dcov->totallength)))
  {
    gt_assert(dcov->samplesize > 0);
    inversesuftab_set(dcov,dcov->totallength,dcov->samplesize-1);
  }
}

static void dc_updatewidth (Differencecover *dcov,unsigned long width,
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

static void dc_initinversesuftabnonspecialsadjust(Differencecover *dcov)
{
  GtCodetype code;
  unsigned int rightchar = 0;
  Bucketspecification bucketspec;
  unsigned long idx = 0;
  const GtCodetype mincode = 0;

  for (code = mincode; code <= dcov->maxcode; code++)
  {
    rightchar = gt_calcbucketboundsparts(&bucketspec,
                                         dcov->bcktab,
                                         code,
                                         dcov->maxcode,
                                         dcov->effectivesamplesize,
                                         rightchar,
                                         dcov->numofchars);
    for (/* Nothing */; idx < bucketspec.left; idx++)
    {
      inversesuftab_set(dcov,suffixptrgetdcov(dcov,idx),idx);
    }
    dc_updatewidth (dcov,bucketspec.nonspecialsinbucket,dcov->prefixlength);
    for (/* Nothing */;
         idx < bucketspec.left + bucketspec.nonspecialsinbucket;
         idx++)
    {
      inversesuftab_set(dcov,suffixptrgetdcov(dcov,idx),bucketspec.left);
    }
  }
  for (/* Nothing */; idx < dcov->effectivesamplesize; idx++)
  {
    inversesuftab_set(dcov,suffixptrgetdcov(dcov,idx),idx);
  }
}

static void dc_anchorleftmost(Differencecover *dcov,
                              unsigned long subbucketleft,
                              unsigned long width)
{
  unsigned long idx;

  for (idx = subbucketleft; idx < subbucketleft + width; idx++)
  {
    inversesuftab_set(dcov,suffixptrgetdcov(dcov,idx),subbucketleft);
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

static void dc_processunsortedrange(Differencecover *dcov,
                                    unsigned long subbucketleft,
                                    unsigned long width,
                                    unsigned long depth)
{
  Pairsuffixptr *pairelem;

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
    dcov->firstwithnewdepth.subbucketleft = subbucketleft;
    dcov->firstwithnewdepth.width = width;
    dcov->firstwithnewdepth.depth = depth;
    dcov->firstwithnewdepth.count = 1UL;
    dcov->firstwithnewdepth.totalwidth = width;
    dcov->firstwithnewdepth.maxwidth = width;
  }
  pairelem = gt_malloc(sizeof (*pairelem));
  pairelem->subbucketleft = subbucketleft;
  pairelem->width = width;
  gt_inl_queue_add(dcov->rangestobesorted,pairelem,false);
  dcov->currentqueuesize++;
  if (dcov->maxqueuesize < dcov->currentqueuesize)
  {
    dcov->maxqueuesize = dcov->currentqueuesize;
  }
}

static int dcov_compareitv(const void *a,const void *b)
{
  const Itventry *itva = (const Itventry *) a,
                 *itvb = (const Itventry *) b;

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

static void dc_sortsuffixesonthislevel(Differencecover *dcov,
                                       unsigned long subbucketleft,
                                       unsigned long width)
{
  unsigned long idx, rangestart, startpos;

  if (dcov->itvinfo == NULL)
  {
    dcov->itvinfo = gt_malloc(sizeof (*dcov->itvinfo) *
                              dcov->allocateditvinfo);
  }
  if (dcov->firstwithnewdepth.subbucketleft == subbucketleft &&
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
    startpos = suffixptrgetdcov(dcov,subbucketleft+idx);
    dcov->itvinfo[idx].suffixstart = startpos;
    dcov->itvinfo[idx].key
      = inversesuftab_get(dcov,startpos + dcov->currentdepth);
    /*
    printf("key=%lu,startpos=%lu,currentdepth=%lu\n",
                    dcov->itvinfo[idx].key,
                    startpos,
                    dcov->currentdepth);
    */
  }
  qsort(dcov->itvinfo,(size_t) width,sizeof (*dcov->itvinfo),dcov_compareitv);
  for (idx=0; idx<width; idx++)
  {
    suffixptrsetdcov(dcov,subbucketleft+idx,dcov->itvinfo[idx].suffixstart);
  }
  rangestart = 0;
  for (idx=1UL; idx<width; idx++)
  {
    if (dcov->itvinfo[idx-1].key != dcov->itvinfo[idx].key)
    {
      if (rangestart + 1 < idx)
      {
        dc_processunsortedrange(dcov,
                                subbucketleft + rangestart,
                                idx - rangestart,
                                GT_MULT2(dcov->currentdepth));
        dc_anchorleftmost(dcov, subbucketleft + rangestart, idx - rangestart);
      } else
      {
        unsigned long currentsuftabentry
          = suffixptrgetdcov(dcov,subbucketleft+rangestart);
        inversesuftab_set(dcov,currentsuftabentry,subbucketleft+rangestart);
      }
      rangestart = idx;
    }
  }
  if (rangestart + 1 < width)
  {
    dc_processunsortedrange(dcov,
                            subbucketleft + rangestart,
                            width - rangestart,
                            GT_MULT2(dcov->currentdepth));
    dc_anchorleftmost(dcov, subbucketleft + rangestart, width - rangestart);
  } else
  {
    unsigned long currentsuftabentry
      = suffixptrgetdcov(dcov,subbucketleft+rangestart);
    inversesuftab_set(dcov,currentsuftabentry,subbucketleft+rangestart);
  }
}

static void dc_bcktab2firstlevelintervals(Differencecover *dcov)
{
  GtCodetype code;
  unsigned int rightchar;
  Bucketspecification bucketspec;
  const GtCodetype mincode = 0;

  printf("# maxbucketsize=%lu\n",dcov->allocateditvinfo);
  rightchar = (unsigned int) (mincode % dcov->numofchars);
  for (code = 0; code <= dcov->maxcode; code++)
  {
    rightchar = gt_calcbucketboundsparts(&bucketspec,
                                         dcov->bcktab,
                                         code,
                                         dcov->maxcode,
                                         dcov->effectivesamplesize,
                                         rightchar,
                                         dcov->numofchars);
    if (bucketspec.nonspecialsinbucket > 1UL)
    {
      dc_sortsuffixesonthislevel(dcov,
                                 bucketspec.left,
                                 bucketspec.nonspecialsinbucket);
    }
  }
}

static void dc_addunsortedrange(void *voiddcov,
                                unsigned long subbucketleft,
                                unsigned long width,
                                GT_UNUSED unsigned long depth)
{
  Differencecover *dcov = (Differencecover *) voiddcov;
  Pairsuffixptr *ptr;

  gt_assert(dcov->sssp == NULL);
  gt_assert(depth >= (unsigned long) dcov->vparam);
  dc_updatewidth (dcov,width,dcov->vparam);
  GT_GETNEXTFREEINARRAY(ptr,&dcov->firstgeneration,Pairsuffixptr,1024);
  ptr->subbucketleft = subbucketleft;
  ptr->width = width;
}

#ifdef  QSORTNAME
#undef  QSORTNAME
#endif

#define QSORTNAME(NAME) dc_##NAME

#ifdef QSORT_ARRAY_DECLARE
#undef QSORT_ARRAY_DECLARE
#endif

#define QSORT_ARRAY_DECLARE\
        Differencecover *dcov = (Differencecover *) data

#ifdef QSORT_ARRAY_GET
#undef QSORT_ARRAY_GET
#endif

#define QSORT_ARRAY_GET(ARR,RELIDX)\
        gt_suffixsortspace_getdirect(dcov->sssp,dcov->sortoffset+(RELIDX))

#ifdef QSORT_ARRAY_SET
#undef QSORT_ARRAY_SET
#endif

#define QSORT_ARRAY_SET(ARR,RELIDX,VALUE)\
        gt_suffixsortspace_setdirect(dcov->sssp,dcov->sortoffset+(RELIDX),VALUE)

int gt_differencecover_compare (const Differencecover *dcov,
                                unsigned long suffixpos1,
                                unsigned long suffixpos2)
{
  unsigned int offset;
  unsigned long idx1, idx2;

  gt_assert(suffixpos1 < dcov->totallength);
  gt_assert(suffixpos2 < dcov->totallength);
  offset = differencecover_offset(dcov,suffixpos1,suffixpos2);
  idx1 = inversesuftab_get(dcov,suffixpos1 + offset);
  idx2 = inversesuftab_get(dcov,suffixpos2 + offset);
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
                  const void *data)
{
  const Differencecover *dcov = (const Differencecover *) data;
  unsigned long suffixpos1, suffixpos2;

  gt_assert(dcov->sssp != NULL);
  suffixpos1 = QSORT_ARRAY_GET(NULL,a);
  suffixpos2 = QSORT_ARRAY_GET(NULL,b);
  return gt_differencecover_compare (dcov, suffixpos1, suffixpos2);
}

typedef void * QSORTNAME(Sorttype);

#include "qsort-array.gen"

void gt_differencecover_sortunsortedbucket(void *data,
                                           unsigned long subbucketleft,
                                           unsigned long width,
                                           GT_UNUSED unsigned long depth)
{
  Differencecover *dcov = (Differencecover *) data;

  gt_assert(depth >= (unsigned long) dcov->vparam);
  gt_assert(dcov->diff2pos != NULL);
  gt_assert(width >= 2UL);
  gt_assert(dcov->sssp != NULL);
  gt_assert(subbucketleft >= gt_suffixsortspace_offset_get(dcov->sssp));
  dcov->sortoffset = subbucketleft - gt_suffixsortspace_offset_get(dcov->sssp);
  QSORTNAME(gt_inlinedarr_qsort_r) (NULL,width,data);
}

static void dc_sortremainingsamples(Differencecover *dcov)
{
  Pairsuffixptr *pairptr;

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
    gt_assert(dcov->firstgeneration.nextfreePairsuffixptr == 0);
  }
  for (pairptr = dcov->firstgeneration.spacePairsuffixptr;
       pairptr < dcov->firstgeneration.spacePairsuffixptr +
                 dcov->firstgeneration.nextfreePairsuffixptr;
       pairptr++)
  {
    dc_anchorleftmost(dcov,pairptr->subbucketleft,pairptr->width);
  }
  for (pairptr = dcov->firstgeneration.spacePairsuffixptr;
       pairptr < dcov->firstgeneration.spacePairsuffixptr +
                 dcov->firstgeneration.nextfreePairsuffixptr;
       pairptr++)
  {
    dc_sortsuffixesonthislevel(dcov,pairptr->subbucketleft, pairptr->width);
  }
  GT_FREEARRAY(&dcov->firstgeneration,Pairsuffixptr);
  while (!gt_inl_queue_isempty(dcov->rangestobesorted))
  {
    Pairsuffixptr *thispairptr;
    thispairptr = (Pairsuffixptr*) gt_inl_queue_get(dcov->rangestobesorted);
    gt_assert(dcov->currentqueuesize > 0);
    dcov->currentqueuesize--;
    dc_sortsuffixesonthislevel(dcov,thispairptr->subbucketleft,
                               thispairptr->width);
    gt_free(thispairptr);
  }
  gt_logger_log(dcov->logger,"maxqueuesize=%lu",dcov->maxqueuesize);
  gt_free(dcov->itvinfo);
  dcov->itvinfo = NULL;
  gt_inl_queue_delete(dcov->rangestobesorted);
  dcov->rangestobesorted = NULL;
}

static void init_sfxstrategy_for_sample(Sfxstrategy *sfxstrategy,
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

static void gt_differencecover_sortsample(Differencecover *dcov,
                                          Outlcpinfo *outlcpinfosample,
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
  dcov->bcktab = gt_allocBcktab(dcov->numofchars, dcov->prefixlength, true,
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
      /* printf("pos mod %u in difference cover\n",dcov->vparam); */
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
        dcov->leftborder[code]++;
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
  (void) gt_bcktab_leftborderpartialsums(dcov->bcktab);
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
          sizeof (*codelist.spaceCodeatposition),compareCodeatpositon);
  }
  dcov->sortedsample = gt_suffixsortspace_new(dcov->effectivesamplesize,
                                              dcov->totallength,
                                              false);
  posinserted = dcov_derivespecialcodesonthefly(dcov,
                                                withcheck ? &codelist : NULL);
  GT_FREEARRAY(&codelist,Codeatposition);
  diffptr = dcov->diffvalues;
  afterend = dcov->diffvalues + dcov->size;
  for (pos = 0, modvalue = 0; pos < dcov->totallength; pos++)
  {
    if (diffptr < afterend && (Diffvalue) modvalue == *diffptr)
    {
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
        sampleindex = --dcov->leftborder[code];
        gt_assert(sampleindex < dcov->effectivesamplesize);
        suffixptrsetdcov(dcov,sampleindex,pos);
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
    init_sfxstrategy_for_sample(&sfxstrategy,
                                mainsfxstrategy,
                                gt_encseq_bitwise_cmp_ok(dcov->encseq)
                                  ? false : true,
                                dcov->effectivesamplesize,
                                dcov->totallength,
                                dcov->logger);
    gt_Outlcpinfo_reinit(outlcpinfosample,dcov->numofchars,dcov->prefixlength,
                         dcov->effectivesamplesize);
    if (outlcpinfosample != NULL)
    {
      dcov->requiredspace += gt_Outlcpinfo_size(outlcpinfosample);
    }
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
      unsigned long idx2 = inversesuftab_get(dcov,suffixptrgetdcov(dcov,idx));
      if (idx != idx2)
      {
        fprintf(stderr,"idx = %lu != %lu = idx2\n",idx,idx2);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
    }
    /*
    if (outlcpinfosample != NULL)
    {
      gt_Outlcpinfo_check_lcpvalues(dcov->encseq,
                                    dcov->readmode,
                                    dcov->sortedsample,
                                    dcov->effectivesamplesize,
                                    (unsigned long) dcov->vparam,
                                    outlcpinfosample);
    }
    */
  }
  gt_suffixsortspace_delete(dcov->sortedsample,false);
  dcov->sortedsample = NULL;
  gt_assert(dcov->diff2pos == NULL);
  filldiff2pos(dcov);
}

static void gt_differencecover_sortsample0(Differencecover *dcov,
                                           Outlcpinfo *outlcpinfosample,
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
      /* printf("pos mod %u in difference cover\n",dcov->vparam); */
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
                                              false);
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
          suffixptrsetdcov(dcov,posinserted,pos);
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

  init_sfxstrategy_for_sample(&sfxstrategy,
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
      unsigned long idx2 = inversesuftab_get(dcov,suffixptrgetdcov(dcov,idx));
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
  filldiff2pos(dcov);
}

Differencecover *gt_differencecover_prepare_sample(
                                        unsigned int vparam,
                                        const GtEncseq *encseq,
                                        GtReadmode readmode,
                                        unsigned int prefixlength,
                                        const Sfxstrategy *sfxstrategy,
                                        Outlcpinfo *outlcpinfosample,
                                        GtLogger *logger,
                                        GtTimer *sfxprogress,
                                        GtError *err)
{
  Differencecover *dcov = NULL;

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
      (prefixlength > 0 ? gt_differencecover_sortsample
                        : gt_differencecover_sortsample0)
                          (dcov,outlcpinfosample,sfxstrategy,
                           sfxstrategy->dccheck);
    }
  }
  return dcov;
}

void gt_differencecover_check(const GtEncseq *encseq,
                               GtReadmode readmode)
{
  Differencecover *dcov;
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
      validate_samplepositons(dcov);
    }
    gt_differencecover_sortsample(dcov,NULL,NULL,withcheck);
    gt_differencecover_delete(dcov);
  }
  printf("# %u difference covers checked\n",
          (unsigned int) (logmod - startlogmod));
}
