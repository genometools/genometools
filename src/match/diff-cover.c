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
#include "divmodmul.h"
#include "intbits-tab.h"
#include "diff-cover.h"
#include "sfx-apfxlen.h"
#include "sfx-enumcodes.h"
#include "bcktab.h"
#include "initbasepower.h"
#include "encseq-def.h"
#include "sfx-suftaborder.h"
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
  unsigned long left,
                right;
} Pairsuffixptr;

typedef Pairsuffixptr Inl_Queueelem;

#include "queue-inline.h"

typedef struct
{
  unsigned long left,
         right;
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
  Bitsequence *isindifferencecover;
  unsigned long samplesize, effectivesamplesize, maxsamplesize;
  const Codetype **multimappower;
  Codetype *filltable;
  Encodedsequencescanstate *esr;
  unsigned long *inversesuftab;
  unsigned long allocateditvinfo,
                currentqueuesize,
                maxqueuesize;
  /* XXX check if this is really necessary */
  Seqpos currentdepth;
  Codetype maxcode;
  Firstwithnewdepth firstwithnewdepth;
  Inl_Queue *rangestobesorted;
  Itventry *itvinfo;
};

/* Compute difference cover on the fly */

#define UScast(X) ((Diffvalue) X)
#define UCcast(X) ((Diffrank) X)

#include "tab-diffcover.h"

static void fillcoverrank(Differencecover *dcov)
{
  unsigned int i;
  Diffrank j;

  dcov->coverrank = gt_malloc(sizeof(*dcov->coverrank) * dcov->vparam);
  for (i=0, j=0; i<dcov->vparam; i++)
  {
    dcov->coverrank[i] = j;
    if (j < dcov->size && dcov->diffvalues[j] <= (Diffvalue) i)
    {
      gt_assert(j < Diffrankmax);
      j++;
    }
  }
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

static Differencecover *differencecover_new(unsigned int vparam,
                                            const Encodedsequence *encseq,
                                            Readmode readmode)
{
  size_t logmod;
  unsigned int offset = 0, v = 1U;
  Differencecover *dcov;
  bool found = false;

  dcov = gt_malloc(sizeof (*dcov));
  for (logmod = 0;
       logmod < sizeof (differencecoversizes)/sizeof (differencecoversizes[0]);
       logmod++)
  {
    if (v == vparam)
    {
      dcov->size = differencecoversizes[logmod];
      dcov->diffvalues = differencecovertab + offset;
      found = true;
      break;
    }
    offset += differencecoversizes[logmod];
    v = MULT2(v);
  }
  if (!found)
  {
    gt_free(dcov);
    return NULL;
  }
  dcov->logmod = (unsigned int) logmod;
  dcov->vparam = 1U << logmod;
  dcov->vmodmask = dcov->vparam-1;
#ifdef WITHcomputehvalue
  dcov->hvalue = computehvalue(dcov,totallength);
#endif
  dcov->totallength = getencseqtotallength(encseq);
  dcov->encseq = encseq;
  dcov->readmode = readmode;
  dcov->numofchars = getencseqAlphabetnumofchars(encseq);
  INITBITTAB(dcov->isindifferencecover,dcov->vparam);
  dcov->maxsamplesize = (unsigned long) (DIVV(dcov->totallength) + 1) *
                                        dcov->size;
  fillcoverrank(dcov);
  filldiff2pos(dcov);
  dcov->esr = newEncodedsequencescanstate();
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
  return dcov;
}

/*
static unsigned int differencecover_offset(const Differencecover *dcov,
                                           Seqpos pos1,Seqpos pos2)
{
  return (unsigned int) MODV(dcov->diff2pos[MODV(pos2-pos1)] - pos1);
}
*/

static void differencecover_delete(Differencecover *dcov)
{
  gt_free(dcov->coverrank);
  dcov->coverrank = NULL;
  gt_free(dcov->diff2pos);
  dcov->diff2pos = NULL;
  gt_free(dcov->sortedsample);
  dcov->sortedsample = NULL;
  gt_assert(dcov->bcktab != NULL);
  gt_free(dcov->inversesuftab);
  dcov->inversesuftab = NULL;
  bcktab_delete(&dcov->bcktab);
  /* XXX the following tables can be deleted before the inversesuftab is
     computed */
  gt_free(dcov->isindifferencecover);
  dcov->isindifferencecover = NULL;
  gt_free(dcov->filltable);
  dcov->filltable = NULL;
  dcov->multimappower = NULL;
  if (dcov->esr != NULL)
  {
    freeEncodedsequencescanstate(&dcov->esr);
  }
  gt_free(dcov);
}

static unsigned long differencecover_packsamplepos(const Differencecover *dcov,
                                            Seqpos pos)
{
  return dcov->coverrank[MODV(pos)] * (DIVV(dcov->totallength) + 1) +
         (unsigned long) DIVV(pos);
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
        if (ISIBITSET(dcov->isindifferencecover,MODV(pos)))
        {
          if (codelist != NULL)
          {
            gt_assert(countderived < codelist->nextfreeCodeatposition);
            gt_assert(codelist->spaceCodeatposition[countderived].maxprefixindex
                      == prefixindex);
            gt_assert(codelist->spaceCodeatposition[countderived].position
                      == pos);
          /* XXX if prefixindex is small then directly extract characters
                 and compute code */
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
      gt_assert(idx < dcov->maxsamplesize && sampleidxused != NULL);
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
  unsigned long idx = differencecover_packsamplepos(dcov,pos);
  gt_assert(idx < dcov->maxsamplesize);
  dcov->inversesuftab[idx] = sampleindex;
}

static unsigned long inversesuftab_get(Differencecover *dcov,Seqpos pos)
{
  unsigned long idx = differencecover_packsamplepos(dcov,pos);
  gt_assert(idx < dcov->maxsamplesize);
  return dcov->inversesuftab[idx];
}

static void initinversesuftab(Differencecover *dcov)
{
  unsigned long sampleindex;
  Seqpos pos;

  dcov->inversesuftab = gt_malloc(sizeof(*dcov->inversesuftab) *
                                  dcov->maxsamplesize);
  /* XXX improve this as in sfx-remainsort */
  for (sampleindex=0; sampleindex<dcov->effectivesamplesize; sampleindex++)
  {
    pos = dcov->sortedsample[sampleindex];
    inversesuftab_set(dcov,pos,sampleindex);
  }
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
      for (pos = range.leftpos; pos < range.rightpos; pos++)
      {
        if (ISIBITSET(dcov->isindifferencecover,MODV(pos)))
        {
          inversesuftab_set(dcov,pos,specialidx);
          specialidx++;
        }
      }
    }
    freespecialrangeiterator(&sri);
  }
  if (ISIBITSET(dcov->isindifferencecover,MODV(dcov->totallength)))
  {
    gt_assert(dcov->samplesize > 0);
    inversesuftab_set(dcov,dcov->totallength,dcov->samplesize-1);
  }
}

static void dc_updatewidth (Differencecover *dcov,unsigned long width,
                            Seqpos depth)
{
  if (width > 1UL)
  {
    if (dcov->allocateditvinfo < width)
    {
      dcov->allocateditvinfo = width;
    }
    if (dcov->currentdepth == 0)
    {
      dcov->currentdepth = depth;
    } else
    {
      gt_assert(dcov->currentdepth == depth);
    }
  }
}

static void dc_initinversesuftabnonspecialsadjust(Differencecover *dcov)
{
  Codetype code;
  unsigned int rightchar;
  Bucketspecification bucketspec;
  Seqpos startpos;
  const Codetype mincode = 0;
  unsigned long idx;

  rightchar = (unsigned int) (mincode % dcov->numofchars);
  idx = 0;
  for (code = 0; code <= dcov->maxcode; code++)
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
      startpos = dcov->sortedsample[idx];
      inversesuftab_set(dcov,startpos,idx);
    }
    dc_updatewidth (dcov,bucketspec.nonspecialsinbucket,
                    (Seqpos) dcov->prefixlength);
    for (/* Nothing */;
         idx < (unsigned long) bucketspec.left+bucketspec.nonspecialsinbucket;
         idx++)
    {
      startpos = dcov->sortedsample[idx];
      inversesuftab_set(dcov,startpos,(unsigned long) bucketspec.left);
    }
  }
  for (/* Nothing */; idx < dcov->effectivesamplesize; idx++)
  {
    startpos = dcov->sortedsample[idx];
    inversesuftab_set(dcov,startpos,idx);
  }
}

static void dc_anchorleftmost(Differencecover *dcov,unsigned long left,
                              unsigned long right)
{
  unsigned long idx;

  for (idx = left; idx <= right; idx++)
  {
    inversesuftab_set(dcov,dcov->sortedsample[idx],left);
  }
}

static void dc_showintervalsizes(unsigned long count,unsigned long totalwidth,
                              Seqpos totallength,unsigned long maxwidth)
{
  printf("%lu\n(total=%lu,avg=%.2f,%.2f%% of all, maxwidth=%lu)\n",
          count,
          totalwidth,
          (double) totalwidth/count,
          100.0 * (double) totalwidth/totallength,
          maxwidth);
}

static void dc_processunsortedrange(Differencecover *dcov,
                                    unsigned long left,unsigned long right,
                                    Seqpos depth)
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
      printf("intervals in level " FormatSeqpos "=",
             PRINTSeqposcast(dcov->firstwithnewdepth.depth));
      dc_showintervalsizes(dcov->firstwithnewdepth.count,
                        dcov->firstwithnewdepth.totalwidth,
                        dcov->totallength,
                        dcov->firstwithnewdepth.maxwidth);
    } else
    {
      dcov->firstwithnewdepth.defined = true;
    }
    printf("enter new level with depth=" FormatSeqpos "\n",
            PRINTSeqposcast(depth));
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

static void dc_sortsuffixesonthislevel(Differencecover *dcov,unsigned long left,
                                       unsigned long right)
{
  unsigned long idx, rangestart;
  Seqpos startpos;
  const unsigned long width = right - left + 1;

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
  for (idx=0; idx<width; idx++)
  {
    startpos = dcov->sortedsample[left+idx];
    dcov->itvinfo[idx].suffixstart = startpos;
    dcov->itvinfo[idx].key
      = inversesuftab_get(dcov,startpos + dcov->currentdepth);
  }
  qsort(dcov->itvinfo,(size_t) width,sizeof(*dcov->itvinfo),compareitv);
  for (idx=0; idx<width; idx++)
  {
    dcov->sortedsample[left+idx] = dcov->itvinfo[idx].suffixstart;
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
        Seqpos currentsuftabentry = dcov->sortedsample[left+rangestart];
        inversesuftab_set(dcov,currentsuftabentry,left+rangestart);
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
    Seqpos currentsuftabentry = dcov->sortedsample[left+rangestart];
    inversesuftab_set(dcov,currentsuftabentry,left+rangestart);
  }
}

static void dc_bcktab2firstlevelintervals(Differencecover *dcov)
{
  Codetype code;
  unsigned int rightchar;
  Bucketspecification bucketspec;
  const Codetype mincode = 0;

  dc_initinversesuftabnonspecialsadjust(dcov);
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
                                 (unsigned long) bucketspec.left,
                                 (unsigned long) (bucketspec.left +
                                                  bucketspec.
                                                  nonspecialsinbucket-1));
    }
  }
}

static void dc_sortremainingsuffixes(Differencecover *dcov)
{
  Pairsuffixptr pairptr;

  while (!gt_inl_queue_isempty(dcov->rangestobesorted))
  {
    pairptr = gt_inl_queue_get(dcov->rangestobesorted);
    gt_assert(dcov->currentqueuesize > 0);
    dcov->currentqueuesize--;
    dc_sortsuffixesonthislevel(dcov,pairptr.left,pairptr.right);
  }
  printf("maxqueuesize = %lu\n",dcov->maxqueuesize);
  gt_free(dcov->itvinfo);
  dcov->itvinfo = NULL;
  gt_inl_queue_delete(dcov->rangestobesorted);
  dcov->rangestobesorted = NULL;
}

static void differencecover_sample(Differencecover *dcov,bool withcheck)
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
  dcov->prefixlength = recommendedprefixlength(dcov->numofchars,
                                               (Seqpos) dcov->maxsamplesize);
  dcov->bcktab = allocBcktab(dcov->numofchars,
                             dcov->prefixlength,
                             true,
                             NULL,
                             NULL);
  if (possibletocmpbitwise(dcov->encseq))
  {
    dcov->multimappower = NULL;
  } else
  {
    dcov->multimappower = bcktab_multimappower(dcov->bcktab);
  }
  STAMP;
  dcov->maxcode = bcktab_numofallcodes(dcov->bcktab) - 1;
  dcov->rangestobesorted = gt_inl_queue_new(MAX(16UL,DIV2(dcov->maxcode)));
  STAMP;
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
      SETIBIT(dcov->isindifferencecover,modvalue);
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
  printf("%lu positions are sampled (%.2f) wasted=%lu, pl=%u\n",
                                  dcov->samplesize,
                                  100.0 *
                                  (double) dcov->samplesize/
                                           (dcov->totallength+1),
                                  dcov->maxsamplesize - dcov->samplesize,
                                  dcov->prefixlength);
  printf("specials = %lu, fullspecials=%lu\n",specials,fullspecials);
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
  gt_assert(posinserted == dcov->effectivesamplesize);
  if (withcheck)
  {
    checksortedsuffixes(dcov->encseq,
                        dcov->readmode,
                        dcov->sortedsample,
                        (Seqpos) dcov->effectivesamplesize,
                        false, /* specialsareequal  */
                        false,  /* specialsareequalatdepth0 */
                        (Seqpos) dcov->prefixlength);
  }
  initinversesuftab(dcov);
  dc_bcktab2firstlevelintervals(dcov);
  dc_sortremainingsuffixes(dcov);
}

void differencecovers_check(Seqpos maxcheck,const Encodedsequence *encseq,
                            Readmode readmode)
{
  Differencecover *dcov;
  size_t logmod, next = 0;
  unsigned int j, vparam;
  bool withcheck = true;

  printf("sizeof(differencecovertab)=%lu\n",
          (unsigned long) sizeof (differencecovertab));
  if (maxcheck > getencseqtotallength(encseq))
  {
    maxcheck = getencseqtotallength(encseq);
  }
  for (logmod = 0;
       logmod < sizeof (differencecoversizes)/sizeof (differencecoversizes[0]);
       logmod++)
  {
    vparam = 1U << logmod;
    dcov = differencecover_new(vparam,encseq,readmode);
    if (dcov == NULL)
    {
      fprintf(stderr,"no difference cover for v=%u\n",vparam);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    for (j = 0; j<dcov->size; j++)
    {
      gt_assert(dcov->diffvalues[j] == differencecovertab[next]);
      next++;
    }
    printf("v=%u (size=%u): ",dcov->vparam,dcov->size);
    if (withcheck)
    {
      validate_samplepositons(dcov);
    }
    differencecover_sample(dcov,withcheck);
    differencecover_delete(dcov);
  }
  printf("# %u difference covers checked\n",(unsigned int) logmod);
  gt_assert(next == sizeof (differencecovertab)/sizeof (differencecovertab[0]));
}
