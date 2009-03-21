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

#include <stdio.h>
#include <string.h>
#include "core/chardef.h"
#include "core/queue.h"
#include "core/chardef.h"
#include "core/unused_api.h"
#include "core/ma_api.h"
#include "core/arraydef.h"
#include "seqpos-def.h"
#include "encseq-def.h"
#include "sfx-remainsort.h"
#include "bcktab.h"
#include "compressedtab.h"

typedef struct
{
  Seqpos key,
         suffixstart;
} Itventry;

typedef struct
{
  Seqpos *left,
         *right,
         *base;
} Pairsuffixptr;

DECLAREARRAYSTRUCT(Pairsuffixptr);

typedef struct
{
  Seqpos *left,
         *right,
         depth;
  unsigned long count, totalwidth, maxwidth;
} Firstwithnewdepth;

struct Rmnsufinfo
{
  Compressedtable *inversesuftab;
  Seqpos *suftab;
  GtQueue *rangestobesorted;
  Seqpos partwidth,
         totallength,
         currentdepth;
  const Encodedsequence *encseq;
  const Bcktab *bcktab;
  Readmode readmode;
  unsigned long allocateditvinfo,
                currentqueuesize,
                maxqueuesize;
  Itventry *itvinfo;
  ArrayPairsuffixptr firstgeneration;
  Firstwithnewdepth firstwithnewdepth;
};

static void initinversesuftab(Rmnsufinfo *rmnsufinfo)
{
  Seqpos idx;

  rmnsufinfo->inversesuftab = compressedtable_new(rmnsufinfo->totallength+1,
                                                  rmnsufinfo->totallength);
  for (idx=0; idx < rmnsufinfo->partwidth; idx++)
  {
    compressedtable_update(rmnsufinfo->inversesuftab,rmnsufinfo->suftab[idx],
                           idx);
  }
  if (hasspecialranges(rmnsufinfo->encseq))
  {
    Specialrangeiterator *sri;
    Sequencerange range;
    Seqpos specialidx;

    sri = newspecialrangeiterator(rmnsufinfo->encseq,
                                  ISDIRREVERSE(rmnsufinfo->readmode)
                                  ? false : true);
    specialidx = rmnsufinfo->partwidth;
    while (nextspecialrangeiterator(&range,sri))
    {
      for (idx = range.leftpos; idx < range.rightpos; idx++)
      {
        compressedtable_update(rmnsufinfo->inversesuftab,idx,specialidx);
        specialidx++;
      }
    }
    gt_assert(specialidx == rmnsufinfo->totallength);
    freespecialrangeiterator(&sri);
  }
  compressedtable_update(rmnsufinfo->inversesuftab,rmnsufinfo->totallength,
                         rmnsufinfo->totallength);
}

Rmnsufinfo *newRmnsufinfo(Seqpos *presortedsuffixes,
                          const Encodedsequence *encseq,
                          const Bcktab *bcktab,
                          Readmode readmode,
                          Seqpos partwidth,
                          unsigned long nonspecialsmaxbucketsize)
{
  Rmnsufinfo *rmnsufinfo;

  rmnsufinfo = gt_malloc(sizeof(Rmnsufinfo));
  rmnsufinfo->suftab = presortedsuffixes;
  rmnsufinfo->rangestobesorted = gt_queue_new();
  rmnsufinfo->partwidth = partwidth;
  rmnsufinfo->totallength = getencseqtotallength(encseq);
  rmnsufinfo->encseq = encseq;
  rmnsufinfo->readmode = readmode;
  rmnsufinfo->allocateditvinfo = nonspecialsmaxbucketsize;
  rmnsufinfo->itvinfo
    = gt_malloc(sizeof(Itventry) * rmnsufinfo->allocateditvinfo);
  rmnsufinfo->bcktab = bcktab;
  rmnsufinfo->currentqueuesize = 0;
  rmnsufinfo->maxqueuesize = 0;
  rmnsufinfo->firstwithnewdepth.depth = 0;
  rmnsufinfo->firstwithnewdepth.totalwidth = 0;
  rmnsufinfo->firstwithnewdepth.count = 0;
  rmnsufinfo->firstwithnewdepth.left = NULL;
  rmnsufinfo->firstwithnewdepth.right = NULL;
  rmnsufinfo->firstwithnewdepth.maxwidth = 0;
  rmnsufinfo->currentdepth = 0;
  INITARRAY(&rmnsufinfo->firstgeneration,Pairsuffixptr);
  initinversesuftab(rmnsufinfo);
  return rmnsufinfo;
}

static void showintervalsizes(unsigned long count,unsigned long totalwidth,
                              Seqpos totallength,unsigned long maxwidth)
{
  printf("%lu (total=%lu,avg=%.2f,%.2f%% of all, maxwidth=%lu)\n",
          count,
          totalwidth,
          (double) totalwidth/count,
          100.0 * (double) totalwidth/totallength,
          maxwidth);
}

void processunsortedrange(Rmnsufinfo *rmnsufinfo,
                          Seqpos *left,Seqpos *right,
                          Seqpos *base,Seqpos depth)
{
  Pairsuffixptr *pairptr;
  unsigned long width;

  gt_assert(left < right && depth > 0);
  gt_assert(rmnsufinfo->firstwithnewdepth.left == NULL ||
            (rmnsufinfo->firstwithnewdepth.depth > 0 &&
             rmnsufinfo->firstwithnewdepth.depth <= depth));
  width = (unsigned long) (right - left + 1);
  if (rmnsufinfo->firstwithnewdepth.left == NULL ||
      rmnsufinfo->firstwithnewdepth.depth < depth)
  {
    if (rmnsufinfo->firstwithnewdepth.left != NULL)
    {
      printf("intervals in previous level=");
      showintervalsizes(rmnsufinfo->firstwithnewdepth.count,
                        rmnsufinfo->firstwithnewdepth.totalwidth,
                        rmnsufinfo->totallength,
                        rmnsufinfo->firstwithnewdepth.maxwidth);
    }
    printf("enter new level with depth=" FormatSeqpos "\n",
            PRINTSeqposcast(depth));
    rmnsufinfo->firstwithnewdepth.left = left;
    rmnsufinfo->firstwithnewdepth.right = right;
    rmnsufinfo->firstwithnewdepth.depth = depth;
    rmnsufinfo->firstwithnewdepth.count = 1UL;
    rmnsufinfo->firstwithnewdepth.totalwidth = width;
    rmnsufinfo->firstwithnewdepth.maxwidth = width;
  } else
  {
    rmnsufinfo->firstwithnewdepth.count++;
    rmnsufinfo->firstwithnewdepth.totalwidth += width;
    if (rmnsufinfo->firstwithnewdepth.maxwidth < width)
    {
      rmnsufinfo->firstwithnewdepth.maxwidth = width;
    }
  }
  pairptr = gt_malloc(sizeof(Pairsuffixptr));
  pairptr->left = left;
  pairptr->right = right;
  pairptr->base = base;
  gt_queue_add(rmnsufinfo->rangestobesorted,pairptr);
  rmnsufinfo->currentqueuesize++;
  if (rmnsufinfo->maxqueuesize < rmnsufinfo->currentqueuesize)
  {
    rmnsufinfo->maxqueuesize = rmnsufinfo->currentqueuesize;
  }
}

static void inversesuftabrange(Rmnsufinfo *rmnsufinfo,Seqpos *left,
                               Seqpos *right,
                               GT_UNUSED Seqpos *base)
{
  Seqpos *ptr, startindex;

  startindex = (Seqpos) (left - rmnsufinfo->suftab);
  for (ptr = left; ptr <= right; ptr++)
  {
    compressedtable_update(rmnsufinfo->inversesuftab,*ptr,startindex);
  }
}

void addunsortedrange(Rmnsufinfo *rmnsufinfo,
                      Seqpos *left,Seqpos *right, Seqpos depth)
{
  Pairsuffixptr *ptr;

  inversesuftabrange(rmnsufinfo,left,right,left);
  GETNEXTFREEINARRAY(ptr,&rmnsufinfo->firstgeneration,Pairsuffixptr,1024);
  if (rmnsufinfo->currentdepth == 0)
  {
    rmnsufinfo->currentdepth = depth;
  } else
  {
    gt_assert(rmnsufinfo->currentdepth == depth);
  }
  ptr->left = left;
  ptr->right = right;
  ptr->base = left;
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

static void sortitv(Rmnsufinfo *rmnsufinfo,Seqpos *left,Seqpos *right,
                    Seqpos *base)
{
  Seqpos startindex;
  unsigned long idx, rangestart;
  const unsigned long width = (unsigned long) (right - left + 1);

  if (rmnsufinfo->firstwithnewdepth.left == left &&
      rmnsufinfo->firstwithnewdepth.right == right)
  {
    rmnsufinfo->currentdepth = rmnsufinfo->firstwithnewdepth.depth;
  }
  gt_assert(rmnsufinfo->allocateditvinfo >= width);
  for (idx=0; idx<width; idx++)
  {
    rmnsufinfo->itvinfo[idx].suffixstart = left[idx];
    gt_assert(left[idx]+rmnsufinfo->currentdepth <= rmnsufinfo->totallength);
    rmnsufinfo->itvinfo[idx].key
      = compressedtable_get(rmnsufinfo->inversesuftab,
                            left[idx]+rmnsufinfo->currentdepth);
  }
  qsort(rmnsufinfo->itvinfo,(size_t) width,sizeof(Itventry),compareitv);
  for (idx=0; idx<width; idx++)
  {
    left[idx] = rmnsufinfo->itvinfo[idx].suffixstart;
  }
  rangestart = 0;
  startindex = (Seqpos) (left - rmnsufinfo->suftab);
  for (idx=1UL; idx<width; idx++)
  {
    if (rmnsufinfo->itvinfo[idx-1].key != rmnsufinfo->itvinfo[idx].key)
    {
      if (rangestart + 1 < idx)
      {
        processunsortedrange(rmnsufinfo,
                             left + rangestart,
                             left + idx - 1,
                             base,
                             MULT2(rmnsufinfo->currentdepth));
        inversesuftabrange(rmnsufinfo,
                           left + rangestart,
                           left + idx - 1,
                           base);
      } else
      {
        compressedtable_update(rmnsufinfo->inversesuftab,left[rangestart],
                               startindex+rangestart);
      }
      rangestart = idx;
    }
  }
  if (rangestart + 1 < width)
  {
    processunsortedrange(rmnsufinfo,
                         left + rangestart,
                         left + width - 1,
                         base,
                         MULT2(rmnsufinfo->currentdepth));
    inversesuftabrange(rmnsufinfo,
                       left + rangestart,
                       left + width - 1,
                       base);
  } else
  {
    compressedtable_update(rmnsufinfo->inversesuftab,left[rangestart],
                           startindex+rangestart);
  }
}

static void sortremainingsuffixes(Rmnsufinfo *rmnsufinfo)
{
  Pairsuffixptr *pairptr;
  unsigned long totalwidth = 0, maxwidth = 0;

  for (pairptr = rmnsufinfo->firstgeneration.spacePairsuffixptr;
       pairptr < rmnsufinfo->firstgeneration.spacePairsuffixptr +
                 rmnsufinfo->firstgeneration.nextfreePairsuffixptr;
       pairptr++)
  {
    unsigned long width = (unsigned long) (pairptr->right - pairptr->left + 1);
    totalwidth += width;
    if (maxwidth < width)
    {
      maxwidth = width;
    }
    sortitv(rmnsufinfo,
            pairptr->left,
            pairptr->right,
            pairptr->base);
  }
  printf("number of intervals at base level " FormatSeqpos " is ",
          PRINTSeqposcast(rmnsufinfo->currentdepth));
  showintervalsizes(rmnsufinfo->firstgeneration.nextfreePairsuffixptr,
                    totalwidth,
                    rmnsufinfo->totallength,
                    maxwidth);
  FREEARRAY(&rmnsufinfo->firstgeneration,Pairsuffixptr);
  while (gt_queue_size(rmnsufinfo->rangestobesorted) > 0)
  {
    pairptr = gt_queue_get(rmnsufinfo->rangestobesorted);
    gt_assert(rmnsufinfo->currentqueuesize > 0);
    rmnsufinfo->currentqueuesize--;
    sortitv(rmnsufinfo,
            pairptr->left,
            pairptr->right,
            pairptr->base);
    gt_free(pairptr);
  }
  printf("maxqueuesize = %lu\n",rmnsufinfo->maxqueuesize);
  gt_free(rmnsufinfo->itvinfo);
  rmnsufinfo->itvinfo = NULL;
  gt_queue_delete(rmnsufinfo->rangestobesorted);
  rmnsufinfo->rangestobesorted = NULL;
}

Seqpos *lcp13_manzini(const Rmnsufinfo *rmnsufinfo)
{
  Seqpos pos, lcpvalue = 0, *lcptab;

  lcptab = gt_malloc(sizeof(Seqpos) * rmnsufinfo->partwidth);
  lcptab[0] = 0;
  for (pos=0; pos <= rmnsufinfo->totallength; pos++)
  {
    Seqpos fillpos = compressedtable_get(rmnsufinfo->inversesuftab,pos);
    if (fillpos > 0 && fillpos < rmnsufinfo->partwidth)
    {
      Seqpos previousstart = rmnsufinfo->suftab[fillpos-1];
      while (pos+lcpvalue < rmnsufinfo->totallength &&
             previousstart+lcpvalue < rmnsufinfo->totallength)
      {
        Uchar cc1, cc2;

        cc1 = getencodedchar(rmnsufinfo->encseq,pos+lcpvalue,
                             rmnsufinfo->readmode);
        cc2 = getencodedchar(rmnsufinfo->encseq,previousstart+lcpvalue,
                             rmnsufinfo->readmode);
        if (cc1 == cc2 && ISNOTSPECIAL(cc1))
        {
          lcpvalue++;
        } else
        {
          break;
        }
      }
      lcptab[fillpos] = lcpvalue;
    }
    if (lcpvalue > 0)
    {
      lcpvalue--;
    }
  }
  return lcptab;
}

static unsigned long *computeocclesstab(const Rmnsufinfo *rmnsufinfo)
{
  unsigned long *occless, numofchars, idx;

  numofchars = (unsigned long) getencseqAlphabetnumofchars(rmnsufinfo->encseq);
  occless = gt_malloc(sizeof (unsigned long) * numofchars);
  occless[0] = 0;
  for (idx = 1UL; idx < numofchars; idx++)
  {
    occless[idx] = occless[idx-1] +
                   getencseqcharactercount(rmnsufinfo->encseq,(Uchar) (idx-1));
    /*
    printf("occless[%lu]=%lu\n",(unsigned long) idx,
                                (unsigned long) occless[idx]);
    */
  }
  return occless;
}

/* for computing the ranknext-values of special positions, we only
   need the values inversesuftab[range.rightpos] in this order,
   where range is a special range
   Now, if range.rightpos = suffixarray[i] for some i, then
   inversesuftab[range.rightpos] = inversesuftab[suffixarray[i]] = i.
   Thus, in case where the inversesuftab is not available,
   we obtain these values by the following function:
*/

static void setrelevantfrominversetab(Compressedtable *rightposinverse,
                                      const Rmnsufinfo *rmnsufinfo)
{
  if (hasspecialranges(rmnsufinfo->encseq))
  {
    Seqpos idx;

    for (idx = 0; idx < rmnsufinfo->partwidth; idx++)
    {
      if (rmnsufinfo->suftab[idx] > 0)
      {
        Seqpos pos = rmnsufinfo->suftab[idx];
        Uchar cc = getencodedchar(rmnsufinfo->encseq,pos-1,
                                  rmnsufinfo->readmode);
        if (ISSPECIAL(cc))
        {
          compressedtable_update(rightposinverse,pos,idx);
          /*
          printf("(1) store rightposinverse[%lu]=%lu\n",
               (unsigned long) pos,
               (unsigned long) idx);
          */
        }
      }
    }
  }
}

static Seqpos *fillrightofpartwidth(const Compressedtable *rightposinverse,
                                    const Rmnsufinfo *rmnsufinfo)
{
  Specialrangeiterator *sri;
  Sequencerange range;
  Seqpos realspecialranges, *rightofpartwidth = NULL;
  unsigned long countranges = 0, nextrightofpartwidth = 0;

  realspecialranges = getencseqrealspecialranges(rmnsufinfo->encseq);
  sri = newspecialrangeiterator(rmnsufinfo->encseq,
                                ISDIRREVERSE(rmnsufinfo->readmode)
                                ? false : true);
  while (nextspecialrangeiterator(&range,sri))
  {
    if (range.rightpos < rmnsufinfo->partwidth)
    {
      countranges++;
    } else
    {
      if (nextrightofpartwidth == 0)
      {
        size_t allocsize = sizeof (Seqpos) * (realspecialranges - countranges);
        rightofpartwidth = gt_malloc(allocsize);
        printf("allocated %lu bytes for rightofpartwidth (%.2f)\n",
                 (unsigned long) allocsize,
                 (double) allocsize/rmnsufinfo->totallength);
      }
      gt_assert(rightofpartwidth != NULL && (Seqpos) nextrightofpartwidth <
                (realspecialranges - countranges));
      rightofpartwidth[nextrightofpartwidth++]
        = compressedtable_get(rightposinverse,range.rightpos);
      /*
      printf("(1) access rightposinverse[%lu]=%lu\n",
             (unsigned long) range.rightpos,
             (unsigned long) rightposinverse[range.rightpos]);
      */
    }
  }
  freespecialrangeiterator(&sri);
  return rightofpartwidth;
}

static void inversesuffixarray2specialranknext(
                         const Compressedtable *rightposinverse,
                         Compressedtable *ranknext,
                         const Rmnsufinfo *rmnsufinfo)
{
  if (hasspecialranges(rmnsufinfo->encseq))
  {
    Specialrangeiterator *sri;
    Sequencerange range;
    Seqpos specialcharacters, idx, *rightofpartwidth = NULL;
    Seqpos specialranklistindex, nextrightofpartwidth = 0;

    rightofpartwidth = fillrightofpartwidth(rightposinverse,rmnsufinfo);
    specialcharacters = getencseqspecialcharacters(rmnsufinfo->encseq);
    specialranklistindex = rmnsufinfo->partwidth;
    sri = newspecialrangeiterator(rmnsufinfo->encseq,
                                  ISDIRREVERSE(rmnsufinfo->readmode)
                                  ? false : true);
    nextrightofpartwidth = 0;
    while (nextspecialrangeiterator(&range,sri))
    {
      gt_assert(range.rightpos<=rmnsufinfo->totallength);
      for (idx = range.leftpos; idx < range.rightpos-1; idx++)
      {
        gt_assert(specialranklistindex < rmnsufinfo->totallength);
        compressedtable_update(ranknext,specialranklistindex,
                               specialranklistindex + 1);
        /*
        printf("(1) set ranknext[%lu] = %lu\n",
                  (unsigned long) specialranklistindex,
                  (unsigned long) ranknext[specialranklistindex]);
        fflush(stdout);
        */
        specialranklistindex++;
      }
      gt_assert(specialranklistindex < rmnsufinfo->totallength);
      if (range.rightpos < rmnsufinfo->partwidth)
      {
        compressedtable_update(ranknext,specialranklistindex,
                               compressedtable_get(rightposinverse,
                                                   range.rightpos));
        /*
        printf("(2) set ranknext[%lu] = %lu = rightposinverse[%lu]\n",
                  (unsigned long) specialranklistindex,
                  (unsigned long) ranknext[specialranklistindex],
                  (unsigned long) range.rightpos);
        fflush(stdout);
        */
      } else
      {
        compressedtable_update(ranknext,specialranklistindex,
                               rightofpartwidth[nextrightofpartwidth]);
        nextrightofpartwidth++;
      }
      specialranklistindex++;
    }
    gt_free(rightofpartwidth);
    gt_assert(specialranklistindex == rmnsufinfo->totallength);
    freespecialrangeiterator(&sri);
  }
}

static Seqpos sa2ranknext(Compressedtable *ranknext,
                          const Rmnsufinfo *rmnsufinfo)
{
  Seqpos idx, longest = 0;
  unsigned long *occless;

  gt_assert(rmnsufinfo->partwidth > 0);

  occless = computeocclesstab(rmnsufinfo);
  /* now inveresuftab is not used any more, and thus the
     ranknext array (which points to ranknext can savely be stored */
  for (idx=0; idx < rmnsufinfo->partwidth; idx++)
  {
    if (rmnsufinfo->suftab[idx] > 0)
    {
      Uchar cc = getencodedchar(rmnsufinfo->encseq,
                                rmnsufinfo->suftab[idx]-1,
                                rmnsufinfo->readmode);
      if (ISNOTSPECIAL(cc))
      {
        gt_assert(occless[cc] < (unsigned long) rmnsufinfo->partwidth);
        /*
        printf("(2) cc=%u: set ranknext[%lu]=%lu\n",
                        (unsigned int) cc,
                        (unsigned long) occless[cc],
                        (unsigned long) idx);
        */
        compressedtable_update(ranknext,(Seqpos) occless[cc],idx);
        occless[cc]++;
      }
    } else
    {
      longest = idx;
    }
  }
  if (hasspecialranges(rmnsufinfo->encseq))
  {
    Specialrangeiterator *sri;
    Sequencerange range;
    Seqpos specialidx, specialpos;

    sri = newspecialrangeiterator(rmnsufinfo->encseq,
                                  ISDIRREVERSE(rmnsufinfo->readmode)
                                  ? false : true);
    gt_assert(rmnsufinfo->partwidth > 0); /* otherwise all lcps would be 0 */
    specialidx = rmnsufinfo->partwidth;
    while (nextspecialrangeiterator(&range,sri))
    {
      for (specialpos = range.leftpos; specialpos < range.rightpos;
           specialpos++)
      {
        if (specialpos > 0)
        {
          Uchar cc = getencodedchar(rmnsufinfo->encseq,
                                    specialpos-1,
                                    rmnsufinfo->readmode);
          if (ISNOTSPECIAL(cc))
          {
            gt_assert(occless[cc] < (unsigned long) rmnsufinfo->partwidth);
            compressedtable_update(ranknext,(Seqpos) occless[cc],specialidx);
            occless[cc]++;
          }
        } else
        {
          longest = rmnsufinfo->partwidth;
        }
        specialidx++;
      }
    }
    if (getencseqlengthofspecialsuffix(rmnsufinfo->encseq))
    {
      compressedtable_update(ranknext,rmnsufinfo->totallength-1,
                             rmnsufinfo->totallength);
    }
    freespecialrangeiterator(&sri);
  }
  gt_free(occless);
  return longest;
}

static Compressedtable *lcp9_manzini(Compressedtable *spacefortab,
                                     Rmnsufinfo *rmnsufinfo)
{
  Seqpos pos, previousstart, nextfillpos = 0, fillpos, lcpvalue = 0;
  Compressedtable *lcptab, *ranknext, *rightposinverse;
  Seqpos previouscc1pos, previouscc2pos;
  Uchar cc1, cc2;
  Encodedsequencescanstate *esr1, *esr2;

  if (spacefortab == NULL)
  {
    rightposinverse = ranknext
                    = compressedtable_new(rmnsufinfo->totallength+1,
                                          rmnsufinfo->totallength);
    compressedtable_update(ranknext,rmnsufinfo->totallength,
                           rmnsufinfo->totallength);
    setrelevantfrominversetab(rightposinverse,rmnsufinfo);
  } else
  {
    rightposinverse = ranknext = spacefortab;
  }
  inversesuffixarray2specialranknext(rightposinverse,ranknext,rmnsufinfo);
  fillpos = sa2ranknext(ranknext,rmnsufinfo);
  printf("longest=" FormatSeqpos "\n",PRINTSeqposcast(fillpos));
  lcptab = ranknext;
  /* now ranknext and lcptab point to the same memory area. After reading
     ranknext at position fillpos, the same cell is used for storing
     the determined lcp-value */
  /* exploit the fact, that pos + lcpvalue is monotone */
  esr1 = newEncodedsequencescanstate();
  initEncodedsequencescanstate(esr1,rmnsufinfo->encseq,rmnsufinfo->readmode,0);
  cc1 = sequentialgetencodedchar(rmnsufinfo->encseq,
                                 esr1,0,rmnsufinfo->readmode);
  previouscc1pos = 0;
  esr2 = newEncodedsequencescanstate();
  previouscc2pos = rmnsufinfo->totallength;
  cc2 = 0;
  for (pos = 0; pos < rmnsufinfo->totallength; pos++)
  {
    if (pos < rmnsufinfo->totallength - 1)
    {
      nextfillpos = compressedtable_get(ranknext,fillpos);
    }
    if (fillpos > 0 && fillpos - 1 < rmnsufinfo->partwidth)
    {
      previousstart = rmnsufinfo->suftab[fillpos-1];
      while (pos+lcpvalue < rmnsufinfo->totallength &&
             previousstart+lcpvalue < rmnsufinfo->totallength)
      {
        gt_assert(pos + lcpvalue >= previouscc1pos);
        while (previouscc1pos < pos + lcpvalue)
        {
          previouscc1pos++;
          cc1 = sequentialgetencodedchar(rmnsufinfo->encseq,
                                         esr1,
                                         previouscc1pos,
                                         rmnsufinfo->readmode);
        }
        if (ISSPECIAL(cc1))
        {
          break;
        }
        if (previousstart+lcpvalue < previouscc2pos ||
            previousstart+lcpvalue > previouscc2pos+1)
        {
          previouscc2pos = previousstart+lcpvalue;
          initEncodedsequencescanstate(esr2,rmnsufinfo->encseq,
                                       rmnsufinfo->readmode,
                                       previouscc2pos);
          cc2 = sequentialgetencodedchar(rmnsufinfo->encseq,
                                         esr2,previouscc2pos,
                                         rmnsufinfo->readmode);
        } else
        {
          if (previousstart+lcpvalue == previouscc2pos+1)
          {
            previouscc2pos++;
            cc2 = sequentialgetencodedchar(rmnsufinfo->encseq,
                                           esr2,
                                           previouscc2pos,
                                           rmnsufinfo->readmode);
          } else
          {
            gt_assert(previousstart+lcpvalue == previouscc2pos);
          }
        }
        if (cc1 != cc2)
        {
          break;
        }
        lcpvalue++;
      }
      compressedtable_update(lcptab,fillpos,lcpvalue);
      if (lcpvalue > 0)
      {
        lcpvalue--;
      }
    }
    fillpos = nextfillpos;
  }
  freeEncodedsequencescanstate(&esr1);
  freeEncodedsequencescanstate(&esr2);
  return lcptab;
}

Compressedtable *wrapRmnsufinfo(Rmnsufinfo **rmnsufinfoptr,bool withlcptab)
{
  Rmnsufinfo *rmnsufinfo = *rmnsufinfoptr;
  Compressedtable *lcptab;

  sortremainingsuffixes(rmnsufinfo);
  if (withlcptab)
  {
#ifdef NOINVERSESUFTAB
    compressedtable_free(rmnsufinfo->inversesuftab);
    rmnsufinfo->inversesuftab = NULL;
    lcptab = lcp9_manzini(NULL,rmnsufinfo);
#else
    gt_assert(rmnsufinfo->inversesuftab != NULL);
    lcptab = lcp9_manzini(rmnsufinfo->inversesuftab,rmnsufinfo);
    rmnsufinfo->inversesuftab = NULL;
#endif
  } else
  {
    compressedtable_free(rmnsufinfo->inversesuftab,true);
    rmnsufinfo->inversesuftab = NULL;
    lcptab = NULL;
  }
  gt_free(rmnsufinfo);
  rmnsufinfoptr = NULL;
  return lcptab;
}
