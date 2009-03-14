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
#include "core/minmax.h"
#include "core/unused_api.h"
#include "core/ma_api.h"
#include "core/qsort_r.h"
#include "divmodmul.h"
#include "seqpos-def.h"
#include "encseq-def.h"
#include "sfx-remainsort.h"
#include "stamp.h"

typedef struct
{
  Seqpos key,
         suffixstart;
} Itventry;

typedef struct
{
  Seqpos *left,
         *right;
} Pairsuffixptr;

typedef struct
{
  Seqpos *left,
         *right,
         depth;
  unsigned long count, totalwidth;
} Firstwithnewdepth;

struct Rmnsufinfo
{
  Seqpos *inversesuftab,
         *presortedsuffixes;
  GtQueue *rangestobesorted;
  Seqpos partwidth,
         totallength,
         currentdepth;
  const Encodedsequence *encseq;
  Readmode readmode;
  unsigned long allocateditvinfo,
                currentqueuesize,
                maxqueuesize;
  Itventry *itvinfo;
  Firstwithnewdepth firstwithnewdepth;
};

Rmnsufinfo *newRmnsufinfo(Seqpos *presortedsuffixes,
                          const Encodedsequence *encseq,
                          Readmode readmode,
                          Seqpos partwidth)
{
  Rmnsufinfo *rmnsufinfo;

  rmnsufinfo = gt_malloc(sizeof(Rmnsufinfo));
  rmnsufinfo->presortedsuffixes = presortedsuffixes;
  rmnsufinfo->rangestobesorted = gt_queue_new();
  rmnsufinfo->partwidth = partwidth;
  rmnsufinfo->totallength = getencseqtotallength(encseq);
  rmnsufinfo->encseq = encseq;
  rmnsufinfo->readmode = readmode;
  rmnsufinfo->allocateditvinfo = 0;
  rmnsufinfo->itvinfo = NULL;
  rmnsufinfo->currentqueuesize = 0;
  rmnsufinfo->maxqueuesize = 0;
  rmnsufinfo->firstwithnewdepth.depth = 0;
  rmnsufinfo->firstwithnewdepth.totalwidth = 0;
  rmnsufinfo->firstwithnewdepth.count = 0;
  rmnsufinfo->firstwithnewdepth.left = NULL;
  rmnsufinfo->firstwithnewdepth.right = NULL;
  rmnsufinfo->currentdepth = 0;
  return rmnsufinfo;
}

void addunsortedrange(Rmnsufinfo *rmnsufinfo,
                      Seqpos *left,Seqpos *right,Seqpos depth)
{
  Pairsuffixptr *pairptr;
  unsigned long width;

  gt_assert(left < right && depth > 0);
  gt_assert(rmnsufinfo->firstwithnewdepth.left == NULL ||
            (rmnsufinfo->firstwithnewdepth.depth > 0 &&
             rmnsufinfo->firstwithnewdepth.depth <= depth));
  if (rmnsufinfo->firstwithnewdepth.left == NULL ||
      rmnsufinfo->firstwithnewdepth.depth < depth)
  {
    if (rmnsufinfo->firstwithnewdepth.left != NULL)
    {
      printf("intervals in previous level=%lu (total=%lu,avg=%.2f,"
             "%.2f%% of all)\n",
              rmnsufinfo->firstwithnewdepth.count,
              rmnsufinfo->firstwithnewdepth.totalwidth,
              (double) rmnsufinfo->firstwithnewdepth.totalwidth/
                       rmnsufinfo->firstwithnewdepth.count,
              100.0 * (double) rmnsufinfo->firstwithnewdepth.totalwidth/
                               rmnsufinfo->totallength);
    }
    printf("enter new level with depth=" FormatSeqpos "\n",
            PRINTSeqposcast(depth));
    rmnsufinfo->firstwithnewdepth.left = left;
    rmnsufinfo->firstwithnewdepth.right = right;
    rmnsufinfo->firstwithnewdepth.depth = depth;
    rmnsufinfo->firstwithnewdepth.count = 1UL;
    rmnsufinfo->firstwithnewdepth.totalwidth
      = (unsigned long) (right - left + 1);
  } else
  {
    rmnsufinfo->firstwithnewdepth.count++;
    rmnsufinfo->firstwithnewdepth.totalwidth
      += (unsigned long) (right - left + 1);
  }
  width = (unsigned long) (right - left + 1);
  if (rmnsufinfo->allocateditvinfo < width)
  {
    gt_assert(rmnsufinfo->itvinfo == NULL);
    rmnsufinfo->allocateditvinfo = width;
  }
  pairptr = gt_malloc(sizeof(Pairsuffixptr));
  pairptr->left = left;
  pairptr->right = right;
  gt_queue_add(rmnsufinfo->rangestobesorted,pairptr);
  rmnsufinfo->currentqueuesize++;
  if (rmnsufinfo->maxqueuesize < rmnsufinfo->currentqueuesize)
  {
    rmnsufinfo->maxqueuesize = rmnsufinfo->currentqueuesize;
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

static void inverserange(Rmnsufinfo *rmnsufinfo,Seqpos *left,Seqpos *right)
{
  Seqpos *ptr, startindex;

  startindex = (Seqpos) (left - rmnsufinfo->presortedsuffixes);
  for (ptr = left; ptr <= right; ptr++)
  {
    rmnsufinfo->inversesuftab[*ptr] = startindex;
  }
}

static void sortitv(Rmnsufinfo *rmnsufinfo,Seqpos *left,Seqpos *right)
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
    if (left[idx]+rmnsufinfo->currentdepth > rmnsufinfo->totallength)
    {
      fprintf(stderr,"left[%lu]+depth=%lu+%lu=%lu>%lu\n",
              idx,
              (unsigned long) left[idx],
              (unsigned long) rmnsufinfo->currentdepth,
              (unsigned long) (left[idx]+rmnsufinfo->currentdepth),
              (unsigned long) rmnsufinfo->totallength);
      exit(EXIT_FAILURE); /* programming error */
    }
    rmnsufinfo->itvinfo[idx].key
      = rmnsufinfo->inversesuftab[left[idx]+rmnsufinfo->currentdepth];
  }
  qsort(rmnsufinfo->itvinfo,(size_t) width,sizeof(Itventry),compareitv);
  for (idx=0; idx<width; idx++)
  {
    left[idx] = rmnsufinfo->itvinfo[idx].suffixstart;
  }
  rangestart = 0;
  startindex = (Seqpos) (left - rmnsufinfo->presortedsuffixes);
  for (idx=1UL; idx<width; idx++)
  {
    if (rmnsufinfo->itvinfo[idx-1].key != rmnsufinfo->itvinfo[idx].key)
    {
      if (rangestart + 1 < idx)
      {
        addunsortedrange(rmnsufinfo,
                         left + rangestart,
                         left + idx - 1,
                         MULT2(rmnsufinfo->currentdepth));
        inverserange(rmnsufinfo,
                     left + rangestart,
                     left + idx - 1);
      } else
      {
        rmnsufinfo->inversesuftab[left[rangestart]] = startindex+rangestart;
      }
      rangestart = idx;
    }
  }
  if (rangestart + 1 < width)
  {
    addunsortedrange(rmnsufinfo,
                     left + rangestart,
                     left + width - 1,
                     MULT2(rmnsufinfo->currentdepth));
    inverserange(rmnsufinfo,
                 left + rangestart,
                 left + width - 1);
  } else
  {
    rmnsufinfo->inversesuftab[left[rangestart]] = startindex+rangestart;
  }
}

static int putleftbound(void **elem,void *info, GT_UNUSED GtError *err)
{
  Pairsuffixptr *pairptr = *(Pairsuffixptr**) elem;

  inverserange((Rmnsufinfo *) info,pairptr->left,pairptr->right);
  return 0;
}

static void sortremainingsuffixes(Rmnsufinfo *rmnsufinfo)
{
  Pairsuffixptr *pairptr;
  Seqpos idx;

  printf("# countovermaxdepth=%lu\n",
           gt_queue_size(rmnsufinfo->rangestobesorted));
  rmnsufinfo->inversesuftab
    = gt_malloc(sizeof(Seqpos) * (rmnsufinfo->totallength+1));
  for (idx=0; idx < rmnsufinfo->partwidth; idx++)
  {
    rmnsufinfo->inversesuftab[rmnsufinfo->presortedsuffixes[idx]] = idx;
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
        rmnsufinfo->inversesuftab[idx] = specialidx++;
      }
    }
    gt_assert(specialidx == rmnsufinfo->totallength);
    freespecialrangeiterator(&sri);
  }
  /* now we can unmap the encoded sequence */
  rmnsufinfo->inversesuftab[rmnsufinfo->totallength] = rmnsufinfo->totallength;
  (void) gt_queue_iterate(rmnsufinfo->rangestobesorted,
                          putleftbound,
                          rmnsufinfo,
                          NULL);
  rmnsufinfo->itvinfo
    = gt_malloc(sizeof(Itventry) * (rmnsufinfo->allocateditvinfo));
  gt_assert(rmnsufinfo->currentqueuesize  ==
            gt_queue_size(rmnsufinfo->rangestobesorted));
  while (gt_queue_size(rmnsufinfo->rangestobesorted) > 0)
  {
    pairptr = gt_queue_get(rmnsufinfo->rangestobesorted);
    gt_assert(rmnsufinfo->currentqueuesize > 0);
    rmnsufinfo->currentqueuesize--;
    sortitv(rmnsufinfo,
            pairptr->left,
            pairptr->right);
    gt_free(pairptr);
  }
  printf("maxqueuesize = %lu\n",rmnsufinfo->maxqueuesize);
}

Seqpos *lcp13_manzini(const Rmnsufinfo *rmnsufinfo)
{
  Seqpos pos, lcpvalue = 0, *lcptab;

  lcptab = gt_malloc(sizeof(Seqpos) * rmnsufinfo->partwidth);
  lcptab[0] = 0;
  for (pos=0; pos <= rmnsufinfo->totallength; pos++)
  {
    Seqpos fillpos = rmnsufinfo->inversesuftab[pos];
    if (fillpos > 0 && fillpos < rmnsufinfo->partwidth)
    {
      Seqpos previousstart = rmnsufinfo->presortedsuffixes[fillpos-1];
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
  }
  return occless;
}

#undef WITHrankedbounds
#ifdef WITHrankedbounds
typedef struct
{
  Seqpos lowerbound,
         upperbound,
         rank;
} Rankedbounds;

static Rankedbounds *fillrankbounds(const Rmnsufinfo *rmnsufinfo)
{
  if (hasspecialranges(rmnsufinfo->encseq))
  {
    Specialrangeiterator *sri;
    Sequencerange range;
    Seqpos currentrank = 0, realspecialranges;
    Rankedbounds *rankedbounds, *rbptr;

    realspecialranges = getencseqrealspecialranges(rmnsufinfo->encseq);
    rankedbounds = gt_malloc(sizeof(Rankedbounds) * realspecialranges);
    sri = newspecialrangeiterator(rmnsufinfo->encseq,
                                  ISDIRREVERSE(rmnsufinfo->readmode)
                                  ? false : true);
    for (rbptr = rankedbounds; nextspecialrangeiterator(&range,sri); rbptr++)
    {
      rbptr->lowerbound = range.leftpos;
      rbptr->upperbound = range.rightpos;
      rbptr->rank = currentrank;
      currentrank += rbptr->upperbound - rbptr->lowerbound;
    }
    gt_assert(rbptr == rankedbounds + realspecialranges);
    freespecialrangeiterator(&sri);
    return rankedbounds;
  }
  return NULL;
}

static Seqpos frompos2rank(const Rankedbounds *leftptr,
                           const Rankedbounds *rightptr,
                           Seqpos specialpos)
{
  const Rankedbounds *midptr;

  while (leftptr <= rightptr)
  {
    midptr = leftptr + DIV2((unsigned long) (rightptr-leftptr));
    if (specialpos < midptr->lowerbound)
    {
      rightptr = midptr-1;
    } else
    {
      if (specialpos >= midptr->upperbound)
      {
        leftptr = midptr + 1;
      } else
      {
        return midptr->rank + specialpos - midptr->lowerbound;
      }
    }
  }
  fprintf(stderr,"frompos2rank: cannot find pos " FormatSeqpos
                 " in ranges",PRINTSeqposcast(specialpos));
  exit(EXIT_FAILURE);
  return 0;
}
#endif

#ifdef WIHTRANK2POS
static Seqpos fromrank2pos(const Rankedbounds *leftptr,
                           const Rankedbounds *rightptr,
                           Seqpos rank)
{
  const Rankedbounds *midptr;

  while (leftptr <= rightptr)
  {
    midptr = leftptr + DIV2((unsigned long) (rightptr-leftptr));
    if (rank < midptr->rank)
    {
      rightptr = midptr-1;
    } else
    {
      if (rank >= midptr->rank + (midptr->upperbound - midptr->lowerbound))
      {
        leftptr = midptr + 1;
      } else
      {
        return midptr->lowerbound + (rank - midptr->rank);
      }
    }
  }
  fprintf(stderr,"fromrank2rank: cannot find rank " FormatSeqpos
                 " in ranges",PRINTSeqposcast(rank));
  exit(EXIT_FAILURE); /* programming error */
  /*@ignore@*/
  return 0;
  /*@end@*/
}
#endif

typedef struct
{
  Seqpos specialrank,
         key;
} Specialrank;

static int compareSpecialrank(const void *a,const void *b)
{
  const Specialrank *aptr = (const Specialrank *) a,
                    *bptr = (const Specialrank *) b;

  if (aptr->key < bptr->key)
  {
    return -1;
  }
  if (aptr->key > bptr->key)
  {
    return 1;
  }
  gt_assert(false);
  return 0;
}

static Specialrank *fillspecialranklist(const Rmnsufinfo *rmnsufinfo)
{
  if (hasspecialranges(rmnsufinfo->encseq))
  {
    Specialrangeiterator *sri;
    Sequencerange range;
    Seqpos realspecialranges, specialrank;
    Specialrank *specialranklist, *rbptr;

    realspecialranges = getencseqrealspecialranges(rmnsufinfo->encseq);
    specialranklist = gt_malloc(sizeof(Specialrank) * realspecialranges);
    sri = newspecialrangeiterator(rmnsufinfo->encseq,
                                  ISDIRREVERSE(rmnsufinfo->readmode)
                                  ? false : true);
    rbptr = specialranklist;
    specialrank = 0;
    while (nextspecialrangeiterator(&range,sri))
    {
      gt_assert(rbptr < specialranklist + realspecialranges);
      gt_assert(range.rightpos<=rmnsufinfo->totallength);
      specialrank += range.rightpos - range.leftpos;
      rbptr->specialrank = specialrank - 1;
      rbptr->key = rmnsufinfo->inversesuftab[range.rightpos];
      /*
      printf("fill at %lu: key = %lu, value = %lu\n",
             (unsigned long) (rbptr - specialranklist),
             (unsigned long) rbptr->key,
             (unsigned long) rbptr->specialrank);
      */
      rbptr++;
    }
    gt_assert(rbptr == specialranklist + realspecialranges);
    freespecialrangeiterator(&sri);
    qsort(specialranklist,(size_t) realspecialranges,
          sizeof (Specialrank),compareSpecialrank);
    return specialranklist;
  }
  return NULL;
}

static Seqpos *fillnextinverselist(const Rmnsufinfo *rmnsufinfo)
{
  if (hasspecialranges(rmnsufinfo->encseq))
  {
    Specialrangeiterator *sri;
    Sequencerange range;
    Seqpos specialcharacters, specialrank;
    Seqpos *specialranklist, specialranklistindex = 0;

    specialcharacters = getencseqspecialcharacters(rmnsufinfo->encseq);
    specialranklist = gt_malloc(sizeof(Specialrank) * specialcharacters);
    sri = newspecialrangeiterator(rmnsufinfo->encseq,
                                  ISDIRREVERSE(rmnsufinfo->readmode)
                                  ? false : true);
    while (nextspecialrangeiterator(&range,sri))
    {
      gt_assert(range.rightpos<=rmnsufinfo->totallength);
      for (specialrank = range.leftpos; specialrank < range.rightpos;
           specialrank++)
      {
        gt_assert(specialranklistindex < specialcharacters);
        specialranklist[specialranklistindex++]
          = rmnsufinfo->inversesuftab[specialrank+1];
      }
      /*
      printf("fill at %lu: key = %lu, value = %lu\n",
             (unsigned long) (rbptr - specialranklist),
             (unsigned long) rbptr->key,
             (unsigned long) rbptr->specialrank);
      */
    }
    gt_assert(specialranklistindex == specialcharacters);
    freespecialrangeiterator(&sri);
    return specialranklist;
  }
  return NULL;
}

#ifdef WITHrankedbounds
static void updateranknext(bool insiderange,
                           Rmnsufinfo *rmnsufinfo,
                           Seqpos *ranknext,
                           Specialrank *specialranklist,
                           unsigned long *specialranklistindex,
                           unsigned long *occless,
                           const Rankedbounds *rankedbounds,
                           Seqpos realspecialranges,
                           Seqpos pos,
                           Seqpos idx)
{
  Uchar cc = getencodedchar(rmnsufinfo->encseq,pos,rmnsufinfo->readmode);

  if (ISNOTSPECIAL(cc))
  {
    gt_assert(occless[cc] < (unsigned long) rmnsufinfo->partwidth);
    ranknext[occless[cc]++] = idx;
  } else
  {
    Seqpos rank;

    if (insiderange)
    {
      gt_assert(idx > 0);
      rank = idx-1;
    } else
    {
      rank = rmnsufinfo->partwidth +
             specialranklist[*specialranklistindex].specialrank;
      (*specialranklistindex)++;
    }
    gt_assert(rank <= rmnsufinfo->totallength);
    ranknext[rank] = idx;
    {
      Seqpos rank2 = rmnsufinfo->partwidth +
                     frompos2rank(rankedbounds,
                                  rankedbounds + realspecialranges,
                                  pos);
      if (rank != rank2)
      {
        fprintf(stderr,"failure: rank = %lu != %lu rank2\n",
               (unsigned long) rank,
               (unsigned long) rank2);
        exit(EXIT_FAILURE); /* programming error */
      }
    }
  }
}
#endif

#ifdef WITHcomparewithsuffixarray
static int comparewithsuffixarray(const void *a,const void *b,void *data)
{
  const Rmnsufinfo *rmnsufinfo = (const Rmnsufinfo *) data;
  Seqpos aval = *((Seqpos *) a), bval = *((Seqpos *) b);

  if (aval < rmnsufinfo->partwidth)
  {
    if (bval < rmnsufinfo->partwidth)
    {
      if (rmnsufinfo->presortedsuffixes[aval] <
          rmnsufinfo->presortedsuffixes[bval])
      {
        return -1; /* a < b */
      }
      if (rmnsufinfo->presortedsuffixes[aval] >
          rmnsufinfo->presortedsuffixes[bval])
      {
        return 1; /* a > b */
      }
      gt_assert(false);
    } else
    {
      return -1; /* a < b */
    }
  } else
  {
    if (bval < rmnsufinfo->partwidth)
    {
      return 1; /* a > b */
    }
    if (aval < bval)
    {
      return -1; /* a < b */
    }
    if (aval > bval)
    {
      return 1; /* a > b */
    }
    gt_assert(false);
  }
  gt_assert(false);
  return 0;
}
#endif

static Seqpos sa2ranknext(Seqpos *ranknext,Rmnsufinfo *rmnsufinfo)
{
  Seqpos idx, longest = 0;
  unsigned long *occless,
                specialranklistindex;
  Specialrank *specialranklist;
  Seqpos *specialranklist2;
  Seqpos realspecialranges;
#ifdef WITHrankedbounds
  Rankedbounds *rankedbounds;

  rankedbounds = fillrankbounds(rmnsufinfo);
#endif
  realspecialranges = getencseqrealspecialranges(rmnsufinfo->encseq);
  gt_assert(rmnsufinfo->partwidth > 0);
  occless = computeocclesstab(rmnsufinfo);
  specialranklist = fillspecialranklist(rmnsufinfo);
  specialranklist2 = fillnextinverselist(rmnsufinfo);
  /* now inveresuftab is not used any more, and thus the
     ranknext array (which points to ranknext can savely be stored */
  specialranklistindex = 0;
  for (idx=0; idx < rmnsufinfo->partwidth; idx++)
  {
    if (rmnsufinfo->presortedsuffixes[idx] > 0)
    {
      /*
      updateranknext(false,rmnsufinfo,ranknext,specialranklist,
                     &specialranklistindex,
                     occless,rankedbounds,realspecialranges,
                     rmnsufinfo->presortedsuffixes[idx]-1,idx);
      */
      Uchar cc = getencodedchar(rmnsufinfo->encseq,
                                rmnsufinfo->presortedsuffixes[idx]-1,
                                rmnsufinfo->readmode);
      if (ISNOTSPECIAL(cc))
      {
        gt_assert(occless[cc] < (unsigned long) rmnsufinfo->partwidth);
        /*
        printf("(1) set ranknext[%lu]=%lu\n",
                             (unsigned long) occless[cc],
                             (unsigned long) idx);
        */
        ranknext[occless[cc]++] = idx;
      } else
      {
        Seqpos rank = rmnsufinfo->partwidth +
                      specialranklist[specialranklistindex].specialrank;
        ranknext[rank] = idx;
        /*
        printf("(2) set ranknext[%lu]=%lu with key %lu\n",
                             (unsigned long) rank,
                             (unsigned long) idx,
                             (unsigned long)
                             specialranklist[specialranklistindex].key);
        */
        specialranklistindex++;
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
      /*
      printf("range %lu %lu\n",(unsigned long) range.leftpos,
                               (unsigned long) range.rightpos);
      */
      for (specialpos = range.leftpos; specialpos < range.rightpos;
           specialpos++)
      {
        if (specialpos > 0)
        {
          /*
          updateranknext(true,rmnsufinfo,ranknext,specialranklist,
                         &specialranklistindex,occless,rankedbounds,
                         realspecialranges,specialpos-1,specialidx);
          */
          Uchar cc = getencodedchar(rmnsufinfo->encseq,
                                    specialpos-1,
                                    rmnsufinfo->readmode);
          if (ISNOTSPECIAL(cc))
          {
            gt_assert(occless[cc] < (unsigned long) rmnsufinfo->partwidth);
            /*
            printf("(3) set ranknext[%lu]=%lu\n",
                             (unsigned long) occless[cc],
                             (unsigned long) specialidx);
            */
            ranknext[occless[cc]++] = specialidx;
          } else
          {
            gt_assert(specialidx > 0);
            ranknext[specialidx-1] = specialidx;
            /*
            printf("(4) set ranknext[%lu]=%lu\n",
                             (unsigned long) specialidx-1,
                             (unsigned long) specialidx);
            */
          }
        } else
        {
          longest = rmnsufinfo->partwidth;
        }
        specialidx++;
      }
      if (range.rightpos == rmnsufinfo->totallength)
      {
        ranknext[rmnsufinfo->totallength-1] = rmnsufinfo->totallength;
      }
    }
    freespecialrangeiterator(&sri);
  }
#ifdef WITHrankedbounds
  gt_free(rankedbounds);
#endif
  {
    Seqpos fillpos = longest, pos, specialcharacters;
    Seqpos nextfillpos;

    specialranklistindex = 0;
    specialcharacters = getencseqspecialcharacters(rmnsufinfo->encseq);
    for (pos = 0; pos < rmnsufinfo->totallength; pos++)
    {
      if (fillpos < rmnsufinfo->partwidth)
      {
        nextfillpos = ranknext[fillpos];
      } else
      {
        gt_assert((Seqpos) specialranklistindex < specialcharacters);
        nextfillpos = specialranklist2[specialranklistindex];
        if (ranknext[fillpos] != specialranklist2[specialranklistindex])
        {
          fprintf(stderr,
                  "ranknext[%lu] = %lu != %lu = specialranklist[%lu]\n",
                  (unsigned long) fillpos,
                  (unsigned long) ranknext[fillpos],
                  (unsigned long) specialranklist2[specialranklistindex],
                  (unsigned long) specialranklistindex);
          exit(EXIT_FAILURE);
        }
        specialranklistindex++;
      }
      fillpos = nextfillpos;
    }
  }
  gt_free(specialranklist);
  gt_free(specialranklist2);
  gt_free(occless);
  return longest;
}

static Seqpos *lcp9_manzini(Rmnsufinfo *rmnsufinfo)
{
  Seqpos pos, previousstart, nextfillpos, fillpos, lcpvalue = 0, *lcptab,
         *ranknext, *specialranklist2;

  specialranklist2 = fillnextinverselist(rmnsufinfo);
  ranknext = rmnsufinfo->inversesuftab;
  fillpos = sa2ranknext(ranknext,rmnsufinfo);
  printf("longest=%lu\n",(unsigned long) fillpos);
  lcptab = rmnsufinfo->inversesuftab;
  /* now ranknext and lcptab point to the same memory area. After reading
     ranknext at position fillpos, the same cell is used for storing
     the determined lcp-value */
  for (pos = 0; pos < rmnsufinfo->totallength; pos++)
  {
    if (fillpos < rmnsufinfo->partwidth)
    {
      nextfillpos = ranknext[fillpos];
    } else
    {
      nextfillpos = specialranklist2[fillpos - rmnsufinfo->partwidth];
    }
    if (fillpos > 0 && fillpos - 1 < rmnsufinfo->partwidth)
    {
      previousstart = rmnsufinfo->presortedsuffixes[fillpos-1];
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
      if (lcpvalue > 0)
      {
        lcpvalue--;
      }
    }
    fillpos = nextfillpos;
  }
  rmnsufinfo->inversesuftab = NULL;
  gt_free(specialranklist2);
  return lcptab;
}

Seqpos *wrapRmnsufinfo(Rmnsufinfo **rmnsufinfoptr,bool withlcptab)
{
  Rmnsufinfo *rmnsufinfo = *rmnsufinfoptr;
  Seqpos *lcptab;

  sortremainingsuffixes(rmnsufinfo);
  if (withlcptab)
  {
    lcptab = lcp9_manzini(rmnsufinfo);
  } else
  {
    lcptab = NULL;
  }
  gt_free(rmnsufinfo->itvinfo);
  rmnsufinfo->itvinfo = NULL;
  gt_queue_delete(rmnsufinfo->rangestobesorted);
  rmnsufinfo->rangestobesorted = NULL;
  gt_free(rmnsufinfo->inversesuftab);
  rmnsufinfo->inversesuftab = NULL;
  gt_free(rmnsufinfo);
  rmnsufinfoptr = NULL;
  return lcptab;
}
