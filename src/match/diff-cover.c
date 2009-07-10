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
#include "core/ma.h"
#include "core/assert_api.h"
#include "core/unused_api.h"
#include "core/error_api.h"
#include "core/minmax.h"
#include "core/arraydef.h"
#include "core/mathsupport.h"
#include "core/qsort_r.h"
#include "divmodmul.h"
#include "intbits-tab.h"
#include "diff-cover.h"
#include "sfx-apfxlen.h"
#include "sfx-enumcodes.h"
#include "bcktab.h"
#include "initbasepower.h"
#include "encseq-def.h"
#include "sfx-suftaborder.h"
#include "sfx-bentsedg.h"
#include "verbose-def.h"
#include "stamp.h"

typedef unsigned char Diffrank;
#define Diffrankmax ((Diffrank) 255)
typedef unsigned short Diffvalue;
#define Diffvaluemax ((Diffvalue) 65535)

#define MODV(VAL) ((VAL) & dcov->vmodmask)
#define DIVV(VAL) ((VAL) >> dcov->logmod)

typedef struct
{
  unsigned long key;
  Seqpos suffixstart;
} Itventry;

typedef struct
{
  Seqpos *left,
         *right;
} Pairsuffixptr;

GT_DECLAREARRAYSTRUCT(Pairsuffixptr);

typedef Pairsuffixptr Inl_Queueelem;

#include "queue-inline.h"

typedef struct
{
  Seqpos *left,
         *right;
  unsigned long count,
                totalwidth,
                maxwidth;
  Seqpos depth;
  bool defined;
} Firstwithnewdepth;

struct Differencecover
{
  unsigned int vparam, logmod, size, vmodmask,
               hvalue,  /* not necessary */
               numofchars,
               prefixlength;
  Diffrank *coverrank;
  Diffvalue *diffvalues, *diff2pos;
  Seqpos totallength, *sortedsample,
         *leftborder; /* points to bcktab->leftborder */
  Bcktab *bcktab;
  const Encodedsequence *encseq;
  Readmode readmode;
  unsigned long samplesize, effectivesamplesize, maxsamplesize;
  const Codetype **multimappower;
  Codetype *filltable;
  Encodedsequencescanstate *esr;
  unsigned long *inversesuftab;
  unsigned long allocateditvinfo,
                currentqueuesize,
                maxqueuesize;
  Seqpos currentdepth;
  Codetype maxcode;
  Firstwithnewdepth firstwithnewdepth;
  Inl_Queue *rangestobesorted;
  Itventry *itvinfo;
  GtArrayPairsuffixptr firstgeneration;
  unsigned long firstgenerationtotalwidth,
                firstgenerationcount;
  Verboseinfo *verboseinfo;
};

/* Compute difference cover on the fly */

#define UScast(X) ((Diffvalue) X)
#define UCcast(X) ((Diffrank) X)

#include "tab-diffcover.h"

#undef QSORT_INTEGER
#ifdef QSORT_INTEGER
typedef unsigned long Sorttype;

static int qsortcmp (const Sorttype *a,const Sorttype *b,
                     const GT_UNUSED void *data)
{
  if ((*a) < (*b))
  {
    return -1;
  }
  if ((*a) > (*b))
  {
    return 1;
  }
  return 0;
}

#include "qsort-inplace.gen"

static void checkqsort(void)
{
#define MAXSIZE 100000

  int times;
  unsigned long idx, n = (unsigned long) MAXSIZE, a[MAXSIZE];

  for (times = 0; times < 100; times++)
  {
    for (idx = 0; idx < n; idx++)
    {
      a[idx] = gt_rand_max(1000UL);
    }
    gt_inlined_qsort_r (a,n,NULL);
    for (idx = 1UL; idx < n; idx++)
    {
      gt_assert(a[idx-1] <= a[idx]);
    }
  }
}
#endif

static void fillcoverrank(Differencecover *dcov)
{
  unsigned int i;
  Diffrank j;

  dcov->coverrank = gt_malloc(sizeof(*dcov->coverrank) * dcov->vparam);
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
                                     unsigned int modpos)
{
  return dcov->coverrank[modpos] == dcov->size ? false : true;
}

static void filldiff2pos(Differencecover *dcov)
{
  Diffvalue *iptr, *jptr;

  dcov->diff2pos = gt_malloc(sizeof(*dcov->diff2pos) * dcov->vparam);
  for (iptr=dcov->diffvalues + dcov->size - 1; iptr>=dcov->diffvalues; iptr--)
  {
    for (jptr=dcov->diffvalues; jptr<dcov->diffvalues + dcov->size; jptr++)
    {
      dcov->diff2pos[MODV(*jptr - *iptr)] = *iptr;
    }
  }
}

#ifdef WITHcomputehvalue

/* XXX: following function is currently not used */

static unsigned int computehvalue(const Differencecover *dcov,
                                  Seqpos totallength)
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

int differencecover_vparamverify(const Differencecover *dcov,GtError *err)
{
  if (dcov->vparam < dcov->prefixlength)
  {
    gt_error_set(err,"difference cover modulo %u is too small, use larger "
                     "parameter for option -dc",dcov->vparam);
    return -1;
  }
  return 0;
}

Differencecover *differencecover_new(unsigned int vparam,
                                     const Encodedsequence *encseq,
                                     Readmode readmode,
                                     Verboseinfo *verboseinfo)
{
  unsigned int offset = 0, v = 1U;
  Differencecover *dcov;
  bool found = false;

#ifdef QSORT_INTEGER
  checkqsort();
#endif
  dcov = gt_malloc(sizeof (*dcov));
  dcov->numofchars = getencseqAlphabetnumofchars(encseq);
  dcov->totallength = getencseqtotallength(encseq);
  dcov->verboseinfo = verboseinfo;
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
    v = MULT2(v);
  }
  if (!found)
  {
    gt_free(dcov);
    return NULL;
  }
  dcov->maxsamplesize = (unsigned long) (DIVV(dcov->totallength) + 1) *
                                         dcov->size;
  dcov->prefixlength = recommendedprefixlength(dcov->numofchars,
                                               (Seqpos) dcov->maxsamplesize);
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
  dcov->firstwithnewdepth.left = 0;
  dcov->firstwithnewdepth.right = 0;
  dcov->firstwithnewdepth.maxwidth = 0;
  dcov->currentqueuesize = 0;
  dcov->maxqueuesize = 0;
  dcov->inversesuftab = NULL;
  dcov->firstgenerationtotalwidth = 0;
  dcov->firstgenerationcount = 0;
  GT_INITARRAY(&dcov->firstgeneration,Pairsuffixptr);
  return dcov;
}

static unsigned int differencecover_offset(const Differencecover *dcov,
                                           Seqpos pos1,Seqpos pos2)
{
  return (unsigned int) MODV(dcov->diff2pos[MODV(pos2-pos1)] - pos1);
}

void differencecover_delete(Differencecover *dcov)
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

static unsigned long differencecover_packsamplepos(const Differencecover *dcov,
                                                   Seqpos pos)
{
  unsigned long result;

  result =  dcov->coverrank[MODV(pos)] * (DIVV(dcov->totallength) + 1) +
            (unsigned long) DIVV(pos);
  gt_assert(result < dcov->maxsamplesize);
  return result;
}

GT_DECLAREARRAYSTRUCT(Codeatposition);

static unsigned long derivespecialcodesonthefly(Differencecover *dcov,
                                                const GtArrayCodeatposition
                                                       *codelist)
{
  unsigned int prefixindex, unitsnotspecial;
  Enumcodeatposition *ecp;
  Specialcontext specialcontext;
  unsigned long countderived = 0;
  Seqpos pos, sampleindex;
  Codetype code;

  for (prefixindex=1U; prefixindex < dcov->prefixlength; prefixindex++)
  {
    /* XXX use one structure and reinit it */
    ecp = newEnumcodeatposition(dcov->encseq,dcov->readmode,
                                dcov->prefixlength,
                                dcov->numofchars);
    while (nextEnumcodeatposition(&specialcontext,ecp))
    {
      if (prefixindex <= specialcontext.maxprefixindex)
      {
        gt_assert(specialcontext.position >= (Seqpos) prefixindex);
        pos = (Seqpos) (specialcontext.position - prefixindex);
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
          code = extractprefixcode(&unitsnotspecial,
                                   dcov->encseq,
                                   dcov->filltable,
                                   dcov->readmode,
                                   dcov->esr,
                                   dcov->multimappower,
                                   pos,
                                   dcov->prefixlength);
          if (codelist != NULL)
          {
            gt_assert((Codetype) codelist->spaceCodeatposition[
                                           countderived].code == code);
          }
          /*
          printf("%u %lu\n",prefixindex,
                            (unsigned long)
                            (specialcontext.position-prefixindex));
          */
          countderived++;
          updatebckspecials(dcov->bcktab,code,dcov->numofchars,prefixindex);
          gt_assert(code > 0);
          sampleindex = --dcov->leftborder[code];
          gt_assert(sampleindex < (Seqpos) dcov->effectivesamplesize);
          dcov->sortedsample[sampleindex] = pos;
        }
      }
    }
    freeEnumcodeatposition(&ecp);
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
  Seqpos pos;
  unsigned int modvalue;
  Diffvalue *diffptr, *afterend;
  unsigned long idx;
  Bitsequence *sampleidxused = NULL;

  INITBITTAB(sampleidxused,dcov->maxsamplesize);
  diffptr = dcov->diffvalues;
  afterend = dcov->diffvalues + dcov->size;
  for (pos = 0, modvalue = 0; pos <= dcov->totallength; pos++)
  {
    gt_assert(modvalue == MODV(pos));
    gt_assert(diffptr == afterend || *diffptr >= (Diffvalue) modvalue);
    if (diffptr < afterend && (Diffvalue) modvalue == *diffptr)
    {
      idx = differencecover_packsamplepos(dcov,pos);
      gt_assert(sampleidxused != NULL);
      if (ISIBITSET(sampleidxused,idx))
      {
        fprintf(stderr,"sample index %lu for pos %lu already used before\n",
                       idx,(unsigned long) pos);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
      SETIBIT(sampleidxused,idx);
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

static void inversesuftab_set(Differencecover *dcov,Seqpos pos,
                              unsigned long sampleindex)
{
  gt_assert (sampleindex < dcov->samplesize);
  dcov->inversesuftab[differencecover_packsamplepos(dcov,pos)] = sampleindex;
}

static unsigned long inversesuftab_get(const Differencecover *dcov,Seqpos pos)
{
  return dcov->inversesuftab[differencecover_packsamplepos(dcov,pos)];
}

static void initinversesuftabnonspecials(Differencecover *dcov)
{
  unsigned long sampleindex;
  Seqpos pos;

  for (sampleindex=0; sampleindex<dcov->effectivesamplesize; sampleindex++)
  {
    pos = dcov->sortedsample[sampleindex];
    inversesuftab_set(dcov,pos,sampleindex);
  }
}

static unsigned long insertfullspecialrangesample(Differencecover *dcov,
                                                  unsigned long specialidx,
                                                  Seqpos leftpos,
                                                  Seqpos rightpos)
{
  Seqpos pos;

  gt_assert(leftpos < rightpos);
  if (ISDIRREVERSE(dcov->readmode))
  {
    pos = rightpos - 1;
  } else
  {
    pos = leftpos;
  }
  while (true)
  {
    if (ISDIRREVERSE(dcov->readmode))
    {
      Seqpos revpos = REVERSEPOS(dcov->totallength,pos);
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

static void initinversesuftabspecials(Differencecover *dcov)
{
  dcov->inversesuftab = gt_malloc(sizeof(*dcov->inversesuftab) *
                                  dcov->maxsamplesize);
  if (hasspecialranges(dcov->encseq))
  {
    Specialrangeiterator *sri;
    Sequencerange range;
    unsigned long specialidx;

    sri = newspecialrangeiterator(dcov->encseq,
                                  ISDIRREVERSE(dcov->readmode)
                                  ? false : true);
    specialidx = dcov->effectivesamplesize;
    while (nextspecialrangeiterator(&range,sri))
    {
      /*
      printf("specialrange %lu %lu\n",(unsigned long) range.leftpos,
                                      (unsigned long) range.rightpos);
      */
      specialidx = insertfullspecialrangesample(dcov,
                                                specialidx,
                                                range.leftpos,
                                                range.rightpos);
    }
    freespecialrangeiterator(&sri);
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
      dcov->currentdepth = (Seqpos) depth;
    } else
    {
      gt_assert(dcov->currentdepth == (Seqpos) depth);
    }
  }
}

static void dc_initinversesuftabnonspecialsadjust(Differencecover *dcov)
{
  Codetype code;
  unsigned int rightchar = 0;
  Bucketspecification bucketspec;
  unsigned long idx = 0;
  const Codetype mincode = 0;

  for (code = mincode; code <= dcov->maxcode; code++)
  {
    rightchar = calcbucketboundsparts(&bucketspec,
                                      dcov->bcktab,
                                      code,
                                      dcov->maxcode,
                                      (Seqpos) dcov->effectivesamplesize,
                                      rightchar,
                                      dcov->numofchars);
    for (/* Nothing */; idx < (unsigned long) bucketspec.left; idx++)
    {
      inversesuftab_set(dcov,dcov->sortedsample[idx],idx);
    }
    dc_updatewidth (dcov,bucketspec.nonspecialsinbucket,dcov->prefixlength);
    for (/* Nothing */;
         idx < (unsigned long) (bucketspec.left +
                                bucketspec.nonspecialsinbucket);
         idx++)
    {
      inversesuftab_set(dcov,dcov->sortedsample[idx],
                        (unsigned long) bucketspec.left);
    }
  }
  for (/* Nothing */; idx < dcov->effectivesamplesize; idx++)
  {
    inversesuftab_set(dcov,dcov->sortedsample[idx],idx);
  }
}

static void dc_anchorleftmost(Differencecover *dcov,Seqpos *left,
                              Seqpos *right)
{
  Seqpos *ptr;
  unsigned long baseindex = (unsigned long) (left - dcov->sortedsample);

  for (ptr = left; ptr <= right; ptr++)
  {
    inversesuftab_set(dcov,*ptr,baseindex);
  }
}

static void dc_showintervalsizes(unsigned long count,unsigned long totalwidth,
                                 unsigned long effectivesamplesize,
                                 unsigned long maxwidth,
                                 Seqpos depth,
                                 Verboseinfo *verboseinfo)
{
  showverbose(verboseinfo,
              "level " FormatSeqpos
              ": (intervals=%lu,total=%lu,avg=%.2f,%.2f%% of all,maxwidth=%lu)",
              PRINTSeqposcast(depth),
              count,
              totalwidth,
              (double) totalwidth/count,
              100.0 * (double) totalwidth/effectivesamplesize,
              maxwidth);
}

static void dc_processunsortedrange(Differencecover *dcov,
                                    Seqpos *left,Seqpos *right,Seqpos depth)
{
  Pairsuffixptr pairelem;
  unsigned long width;

  gt_assert(left < right && depth > 0);
  gt_assert(!dcov->firstwithnewdepth.defined ||
            (dcov->firstwithnewdepth.depth > 0 &&
             dcov->firstwithnewdepth.depth <= depth));
  width = (unsigned long) (right - left + 1);
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
                           dcov->verboseinfo);
    } else
    {
      dcov->firstwithnewdepth.defined = true;
    }
    showverbose(dcov->verboseinfo,
                "enter new level " FormatSeqpos,PRINTSeqposcast(depth));
    dcov->firstwithnewdepth.left = left;
    dcov->firstwithnewdepth.right = right;
    dcov->firstwithnewdepth.depth = depth;
    dcov->firstwithnewdepth.count = 1UL;
    dcov->firstwithnewdepth.totalwidth = width;
    dcov->firstwithnewdepth.maxwidth = width;
  }
  pairelem.left = left;
  pairelem.right = right;
  gt_inl_queue_add(dcov->rangestobesorted,pairelem,false);
  dcov->currentqueuesize++;
  if (dcov->maxqueuesize < dcov->currentqueuesize)
  {
    dcov->maxqueuesize = dcov->currentqueuesize;
  }
}

static int compareitv(const void *a,const void *b)
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
                                       Seqpos *left,
                                       Seqpos *right)
{
  unsigned long idx, rangestart;
  Seqpos startpos;
  const unsigned long width = (unsigned long) (right - left + 1);

  if (dcov->itvinfo == NULL)
  {
    dcov->itvinfo = gt_malloc(sizeof (*dcov->itvinfo) *
                              dcov->allocateditvinfo);
  }
  if (dcov->firstwithnewdepth.left == left &&
      dcov->firstwithnewdepth.right == right)
  {
    dcov->currentdepth = dcov->firstwithnewdepth.depth;
  }
  gt_assert(dcov->allocateditvinfo >= width);
  /*
  printf("new interval of width %lu\n",width);
  */
  for (idx=0; idx<width; idx++)
  {
    startpos = left[idx];
    dcov->itvinfo[idx].suffixstart = startpos;
    dcov->itvinfo[idx].key
      = inversesuftab_get(dcov,startpos + dcov->currentdepth);
    /*
    printf("key=%lu,startpos=%lu,currentdepth=%lu\n",
                    (unsigned long) dcov->itvinfo[idx].key,
                    (unsigned long) startpos,
                    (unsigned long) dcov->currentdepth);
    */
  }
  qsort(dcov->itvinfo,(size_t) width,sizeof(*dcov->itvinfo),compareitv);
  for (idx=0; idx<width; idx++)
  {
    left[idx] = dcov->itvinfo[idx].suffixstart;
  }
  rangestart = 0;
  for (idx=1UL; idx<width; idx++)
  {
    if (dcov->itvinfo[idx-1].key != dcov->itvinfo[idx].key)
    {
      if (rangestart + 1 < idx)
      {
        dc_processunsortedrange(dcov,
                                left + rangestart,
                                left + idx - 1,
                                MULT2(dcov->currentdepth));
        dc_anchorleftmost(dcov,
                          left + rangestart,
                          left + idx - 1);
      } else
      {
        Seqpos currentsuftabentry = left[rangestart];
        inversesuftab_set(dcov,currentsuftabentry,
                          (unsigned long) (left+rangestart-dcov->sortedsample));
      }
      rangestart = idx;
    }
  }
  if (rangestart + 1 < width)
  {
    dc_processunsortedrange(dcov,
                         left + rangestart,
                         left + width - 1,
                         MULT2(dcov->currentdepth));
    dc_anchorleftmost(dcov,
                      left + rangestart,
                      left + width - 1);
  } else
  {
    Seqpos currentsuftabentry = left[rangestart];
    inversesuftab_set(dcov,currentsuftabentry,
                      (unsigned long) (left+rangestart-dcov->sortedsample));
  }
}

static void dc_bcktab2firstlevelintervals(Differencecover *dcov)
{
  Codetype code;
  unsigned int rightchar;
  Bucketspecification bucketspec;
  const Codetype mincode = 0;

  printf("# maxbucketsize=%lu\n",dcov->allocateditvinfo);
  rightchar = (unsigned int) (mincode % dcov->numofchars);
  for (code = 0; code <= dcov->maxcode; code++)
  {
    rightchar = calcbucketboundsparts(&bucketspec,
                                      dcov->bcktab,
                                      code,
                                      dcov->maxcode,
                                      (Seqpos) dcov->effectivesamplesize,
                                      rightchar,
                                      dcov->numofchars);
    if (bucketspec.nonspecialsinbucket > 1UL)
    {
      dc_sortsuffixesonthislevel(dcov,
                                 dcov->sortedsample + bucketspec.left,
                                 dcov->sortedsample + bucketspec.left +
                                                      bucketspec.
                                                      nonspecialsinbucket-1);
    }
  }
}

static void dc_addunsortedrange(void *voiddcov,
                                Seqpos *left,
                                Seqpos *right,
                                GT_UNUSED Seqpos depth)
{
  Differencecover *dcov = (Differencecover *) voiddcov;
  Pairsuffixptr *ptr;

  gt_assert(depth >= (Seqpos) dcov->vparam);
  dc_updatewidth (dcov,(unsigned long) (right - left + 1),dcov->vparam);
  GT_GETNEXTFREEINARRAY(ptr,&dcov->firstgeneration,Pairsuffixptr,1024);
  ptr->left = left;
  ptr->right = right;
}

#ifdef QSORT_INTEGER
static int comparedcov_presortedsuffixes(const void *a,const void *b,
                                         void *data)
{
  const Differencecover *dcov = (const Differencecover *) data;
  const Seqpos suffixpos1 = *(const Seqpos *) a;
  const Seqpos suffixpos2 = *(const Seqpos *) b;
  unsigned long idx1, idx2;
  unsigned int offset;

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

#else

typedef Seqpos Sorttype;

static int qsortcmp (const Sorttype *a,const Sorttype *b,
                     const void *data)
{
  const Differencecover *dcov = (const Differencecover *) data;
  const Seqpos suffixpos1 = *a;
  const Seqpos suffixpos2 = *b;
  unsigned long idx1, idx2;
  unsigned int offset;

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

#include "qsort-inplace.gen"
#endif

void dc_sortunsortedbucket(void *data,
                           Seqpos *left,
                           Seqpos *right,
                           GT_UNUSED Seqpos depth)
{
#ifdef WITHCHECK
  const Differencecover *dcov = (const Differencecover *) data;

  gt_assert(depth >= (Seqpos) dcov->vparam);
  gt_assert(dcov->diff2pos != NULL);
#endif
  gt_assert(left < right);
#ifdef WITHCHECK
  checksortedsuffixes(__FILE__,
                      __LINE__,
                      dcov->encseq,
                      dcov->readmode,
                      left,
                      (Seqpos) (right - left + 1),
                      false, /* specialsareequal  */
                      false,  /* specialsareequalatdepth0 */
                      (Seqpos) dcov->vparam);
#endif
#ifdef QSORT_INTEGER
  gt_qsort_r(left,(size_t) (right - left + 1),sizeof(Seqpos),data,
             comparedcov_presortedsuffixes);
#else
  gt_inlined_qsort_r (left,(unsigned long) (right - left + 1),data);
#endif
}

static void dc_sortremainingsamples(Differencecover *dcov)
{
  Pairsuffixptr *pairptr, pair;

  if (dcov->firstgenerationcount > 0)
  {
    dc_showintervalsizes(dcov->firstgenerationcount,
                         dcov->firstgenerationtotalwidth,
                         dcov->effectivesamplesize,
                         dcov->allocateditvinfo,
                         dcov->currentdepth,
                         dcov->verboseinfo);
  }
  if (dcov->inversesuftab == NULL)
  { /* now maxdepth > prefixlength */
    initinversesuftabspecials(dcov);
    initinversesuftabnonspecials(dcov);
  } else
  {
    gt_assert(dcov->firstgeneration.nextfreePairsuffixptr == 0);
  }
  for (pairptr = dcov->firstgeneration.spacePairsuffixptr;
       pairptr < dcov->firstgeneration.spacePairsuffixptr +
                 dcov->firstgeneration.nextfreePairsuffixptr;
       pairptr++)
  {
    dc_anchorleftmost(dcov,pairptr->left,pairptr->right);
  }
  for (pairptr = dcov->firstgeneration.spacePairsuffixptr;
       pairptr < dcov->firstgeneration.spacePairsuffixptr +
                 dcov->firstgeneration.nextfreePairsuffixptr;
       pairptr++)
  {
    dc_sortsuffixesonthislevel(dcov,pairptr->left, pairptr->right);
  }
  GT_FREEARRAY(&dcov->firstgeneration,Pairsuffixptr);
  while (!gt_inl_queue_isempty(dcov->rangestobesorted))
  {
    pair = gt_inl_queue_get(dcov->rangestobesorted);
    gt_assert(dcov->currentqueuesize > 0);
    dcov->currentqueuesize--;
    dc_sortsuffixesonthislevel(dcov,pair.left,pair.right);
  }
  showverbose(dcov->verboseinfo,"maxqueuesize = %lu",dcov->maxqueuesize);
  gt_free(dcov->itvinfo);
  dcov->itvinfo = NULL;
  gt_inl_queue_delete(dcov->rangestobesorted);
  dcov->rangestobesorted = NULL;
}

void differencecover_sortsample(Differencecover *dcov,bool cmpcharbychar,
                                bool withcheck)
{
  Seqpos pos;
  unsigned int modvalue;
  Diffvalue *diffptr, *afterend;
  unsigned long fullspecials = 0, specials = 0;
  unsigned int unitsnotspecial;
  Codetype code;
  GtArrayCodeatposition codelist;
  Codeatposition *codeptr;
  Seqpos sampleindex;
  unsigned long posinserted;

  dcov->samplesize = 0;
  dcov->bcktab = allocBcktab(dcov->numofchars,
                             dcov->prefixlength,
                             true,
                             NULL,
                             NULL);
  dcov->multimappower = bcktab_multimappower(dcov->bcktab);
  dcov->esr = newEncodedsequencescanstate();
  dcov->maxcode = bcktab_numofallcodes(dcov->bcktab) - 1;
  dcov->rangestobesorted = gt_inl_queue_new(MAX(16UL,DIV2(dcov->maxcode)));
  gt_assert(dcov->bcktab != NULL);
  dcov->filltable = filllargestchartable(dcov->numofchars,dcov->prefixlength);
  dcov->leftborder = bcktab_leftborder(dcov->bcktab);
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
        code = extractprefixcode(&unitsnotspecial,
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
            gt_assert(code <= (Codetype) MAXCODEVALUE);
            codeptr->code = (unsigned int) code;
            gt_assert(unitsnotspecial <= MAXPREFIXLENGTH);
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
  bcktab_leftborderpartialsums(dcov->bcktab,(Seqpos) dcov->effectivesamplesize);
  showverbose(dcov->verboseinfo,
              "%lu positions are sampled (%.2f) pl=%u",
              dcov->samplesize,
              100.0 * (double) dcov->samplesize/(dcov->totallength+1),
              dcov->prefixlength);
  showverbose(dcov->verboseinfo,"specials = %lu, fullspecials=%lu",
              specials,fullspecials);
  if (withcheck)
  {
    qsort(codelist.spaceCodeatposition,
          (size_t) codelist.nextfreeCodeatposition,
          sizeof (*codelist.spaceCodeatposition),compareCodeatpositon);
  }
  dcov->sortedsample = gt_malloc(sizeof(*dcov->sortedsample) *
                                 dcov->effectivesamplesize);
  posinserted = derivespecialcodesonthefly(dcov,withcheck ? &codelist : NULL);
  GT_FREEARRAY(&codelist,Codeatposition);
  diffptr = dcov->diffvalues;
  afterend = dcov->diffvalues + dcov->size;
  for (pos = 0, modvalue = 0; pos < dcov->totallength; pos++)
  {
    if (diffptr < afterend && (Diffvalue) modvalue == *diffptr)
    {
      code = extractprefixcode(&unitsnotspecial,
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
        gt_assert(sampleindex < (Seqpos) dcov->effectivesamplesize);
        dcov->sortedsample[sampleindex] = pos;
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
  if (dcov->esr != NULL)
  {
    freeEncodedsequencescanstate(&dcov->esr);
  }
  gt_assert(posinserted == dcov->effectivesamplesize);
  if (withcheck)
  {
    checksortedsuffixes(__FILE__,
                        __LINE__,
                        dcov->encseq,
                        dcov->readmode,
                        dcov->sortedsample,
                        (Seqpos) dcov->effectivesamplesize,
                        false, /* specialsareequal  */
                        false,  /* specialsareequalatdepth0 */
                        (Seqpos) dcov->prefixlength);
  }
  if (dcov->vparam == dcov->prefixlength)
  {
    initinversesuftabspecials(dcov);
    dc_initinversesuftabnonspecialsadjust(dcov);
    dc_bcktab2firstlevelintervals(dcov);
  } else
  {
    Sfxstrategy sfxstrategy;

    gt_assert (dcov->vparam > dcov->prefixlength);
    if (cmpcharbychar)
    {
      defaultsfxstrategy(&sfxstrategy,true);
    } else
    {
      defaultsfxstrategy(&sfxstrategy,
                         possibletocmpbitwise(dcov->encseq) ? false : true);
    }
    sfxstrategy.differencecover = dcov->vparam;
    sortbucketofsuffixes(dcov->sortedsample,
                         NULL,
                         dcov->effectivesamplesize,
                         dcov->encseq,
                         dcov->readmode,
                         0, /* mincode */
                         dcov->maxcode,
                         dcov->bcktab,
                         dcov->numofchars,
                         dcov->prefixlength,
                         &sfxstrategy,
                         (void *) dcov,
                         dc_addunsortedrange,
                         NULL);
    if (withcheck)
    {
      checksortedsuffixes(__FILE__,
                          __LINE__,
                          dcov->encseq,
                          dcov->readmode,
                          dcov->sortedsample,
                          (Seqpos) dcov->effectivesamplesize,
                          false, /* specialsareequal  */
                          false,  /* specialsareequalatdepth0 */
                          (Seqpos) dcov->vparam);
    }
  }
  if (dcov->bcktab != NULL)
  {
    bcktab_delete(&dcov->bcktab);
  }
  dc_sortremainingsamples(dcov);
  if (withcheck)
  {
#ifndef NDEBUG
    unsigned long idx;
#endif

    checksortedsuffixes(__FILE__,
                        __LINE__,
                        dcov->encseq,
                        dcov->readmode,
                        dcov->sortedsample,
                        (Seqpos) dcov->effectivesamplesize,
                        false, /* specialsareequal  */
                        false,  /* specialsareequalatdepth0 */
                        0);
#ifndef NDEBUG
    for (idx=0; idx<dcov->effectivesamplesize; idx++)
    {
      unsigned long idx2 = inversesuftab_get(dcov,dcov->sortedsample[idx]);
      gt_assert(idx == idx2);
    }
#endif
  }
  gt_free(dcov->sortedsample);
  dcov->sortedsample = NULL;
  gt_assert(dcov->diff2pos == NULL);
  filldiff2pos(dcov);
}

void differencecovers_check(const Encodedsequence *encseq,Readmode readmode)
{
  Differencecover *dcov;
  size_t logmod;
  unsigned int vparam;
  bool withcheck = true;

  printf("sizeof(differencecovertab)=%lu\n",
          (unsigned long) sizeof (differencecovertab));
  for (logmod = (size_t) 4;
       logmod < sizeof (differencecoversizes)/sizeof (differencecoversizes[0]);
       logmod++)
  {
    vparam = 1U << logmod;
    dcov = differencecover_new(vparam,encseq,readmode,NULL);
    if (dcov == NULL)
    {
      fprintf(stderr,"no difference cover for v=%u\n",vparam);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    printf("v=%u (size=%u): ",dcov->vparam,dcov->size);
    if (withcheck)
    {
      validate_samplepositons(dcov);
    }
    differencecover_sortsample(dcov,false,withcheck);
    differencecover_delete(dcov);
  }
  printf("# %u difference covers checked\n",(unsigned int) logmod);
}
