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

#include <math.h>
#include "core/ma.h"
#include "core/qsort_r.h"
#include "core/array2dim_api.h"
#include "core/unused_api.h"
#include "bcktab.h"
#include "kmer2string.h"
#include "sfx-copysort.h"
#include "verbose-def.h"
#include "stamp.h"

typedef struct
{
  bool hardworktodo,
       sorted;
  Seqpos bucketend;
} Bucketinfo;

struct GtBucketspec2
{
  Seqpos partwidth;
  const Encodedsequence *encseq;
  Readmode readmode;
  unsigned int numofchars, numofcharssquared, prefixlength, *order;
  Codetype expandfactor, expandfillsum;
  Bucketinfo *superbuckettab, **subbuckettab;
};

static Seqpos superbucketsize(const GtBucketspec2 *bucketspec2,
                              unsigned int bucketnum)
{
  if (bucketnum == 0)
  {
    return bucketspec2->superbuckettab[0].bucketend;
  }
  return bucketspec2->superbuckettab[bucketnum].bucketend -
         bucketspec2->superbuckettab[bucketnum-1].bucketend;
}

static int comparesuperbucketsizes(const void *a,const void *b,void *data)
{
  const GtBucketspec2 *bucketspec2 = (const GtBucketspec2 *) data;
  Seqpos size1 = superbucketsize(bucketspec2, *(const unsigned int *) a);
  Seqpos size2 = superbucketsize(bucketspec2, *(const unsigned int *) b);
  if (size1 < size2)
  {
    return -1;
  }
  if (size1 > size2)
  {
    return 1;
  }
  return 0;
}

static Seqpos getstartidx(const GtBucketspec2 *bucketspec2,
                          unsigned int first,
                          unsigned int second)
{
  gt_assert(first < bucketspec2->numofchars);
  gt_assert(second <= bucketspec2->numofchars);
  if (second > 0)
  {
    return bucketspec2->subbuckettab[first][second-1].bucketend;
  }
  if (first > 0)
  {
    return bucketspec2->superbuckettab[first-1].bucketend;
  }
  return 0;
}

static Seqpos getendidx(const GtBucketspec2 *bucketspec2,
                        unsigned int first,
                        unsigned int second)
{
  gt_assert(first < bucketspec2->numofchars);
  gt_assert(second <= bucketspec2->numofchars);
  if (second < bucketspec2->numofchars)
  {
    return bucketspec2->subbuckettab[first][second].bucketend;
  }
  return bucketspec2->superbuckettab[first].bucketend;
}

static void resetsorted(GtBucketspec2 *bucketspec2)
{
  unsigned int idx, idx2;

  for (idx = 0; idx<bucketspec2->numofchars; idx++)
  {
    bucketspec2->superbuckettab[idx].sorted = false;
    for (idx2 = 0; idx2<bucketspec2->numofchars; idx2++)
    {
      Seqpos startidx = getstartidx(bucketspec2,idx,idx2),
             endidx = getendidx(bucketspec2,idx,idx2);
      bucketspec2->subbuckettab[idx][idx2].sorted =
        (startidx < endidx) ? false : true;
    }
  }
}

static void determinehardwork(GtBucketspec2 *bucketspec2)
{
  unsigned int idx, idxsource, source, second;

  for (idxsource = 0; idxsource<bucketspec2->numofchars; idxsource++)
  {
    source = bucketspec2->order[idxsource];
    for (second = 0; second < bucketspec2->numofchars; second++)
    {
      if (!bucketspec2->subbuckettab[source][second].sorted && source != second)
      {
        bucketspec2->subbuckettab[source][second].hardworktodo = true;
        bucketspec2->subbuckettab[source][second].sorted = true;
      } else
      {
        bucketspec2->subbuckettab[source][second].hardworktodo = false;
      }
    }
    bucketspec2->superbuckettab[source].sorted = true;
    for (idx = 0; idx < bucketspec2->numofchars; idx++)
    {
      bucketspec2->subbuckettab[idx][source].sorted = true;
    }
  }
}

static Codetype expandtwocharcode(Codetype twocharcode,
                                  const GtBucketspec2 *bucketspec2)
{
  gt_assert(twocharcode < (Codetype) bucketspec2->numofcharssquared);
  return twocharcode * bucketspec2->expandfactor + bucketspec2->expandfillsum;
}

static void leftcontextofspecialchardist(Seqpos *dist,
                                         const Encodedsequence *encseq,
                                         Readmode readmode)
{
  Specialrangeiterator *sri;
  Sequencerange range;
  GtUchar cc;
  Seqpos totallength = getencseqtotallength(encseq);

  if (hasspecialranges(encseq))
  {
    sri = newspecialrangeiterator(encseq,ISDIRREVERSE(readmode) ? false : true);
    while (nextspecialrangeiterator(&range,sri))
    {
      gt_assert(range.leftpos < totallength);
      if (range.leftpos > 0)
      {
        cc = getencodedchar(encseq,range.leftpos-1,readmode);
        if (ISNOTSPECIAL(cc))
        {
          dist[cc]++;
        }
      }
    }
  }
  freespecialrangeiterator(&sri);
  if (getencseqlengthofspecialsuffix(encseq) == 0)
  {
    cc = getencodedchar(encseq,totallength-1,readmode);
    gt_assert(ISNOTSPECIAL(cc));
    dist[cc]++;
  }
}

#undef SHOWBUCKETSPEC2
#ifdef SHOWBUCKETSPEC2
static void showbucketspec2(const GtBucketspec2 *bucketspec2)
{
  unsigned int idx1, idx2;

  for (idx1 = 0; idx1 < bucketspec2->numofchars; idx1++)
  {
    for (idx2 = 0; idx2 < bucketspec2->numofchars; idx2++)
    {
      printf("subbucket[%u][%u]=" FormatSeqpos,idx1,idx2,
              PRINTSeqposcast(bucketspec2->subbuckettab[idx1][idx2].bucketend));
      if (bucketspec2->subbuckettab[idx1][idx2].sorted)
      {
        printf(" sorted\n");
      } else
      {
        printf("\n");
      }
    }
    printf("superbucket[%u]=" FormatSeqpos "\n",idx1,
              PRINTSeqposcast(bucketspec2->superbuckettab[idx1].bucketend));
  }
}

static void showexpandcode(const GtBucketspec2 *bucketspec2,
                           unsigned int prefixlength)
{
  Codetype ecode, code2;
  const GtUchar *characters = getencseqAlphabetcharacters(bucketspec2->encseq);

  for (code2 = 0; code2 < (Codetype) bucketspec2->numofcharssquared; code2++)
  {
    char buffer[100];

    ecode = expandtwocharcode(code2,bucketspec2);
    fromkmercode2string(buffer,
                        ecode,
                        bucketspec2->numofchars,
                        prefixlength,
                        (const char *) characters);
    printf("code2=%u = %lu %s\n",(unsigned int) code2,ecode,buffer);
  }
}
#endif

static void fill2subbuckets(GtBucketspec2 *bucketspec2,const Bcktab *bcktab)
{
  Codetype code, maxcode;
  unsigned int rightchar = 0, currentchar = 0;
  Bucketspecification bucketspec;
  Seqpos accubucketsize = 0;

  maxcode = bcktab_numofallcodes(bcktab) - 1;
  for (code = 0; code <= maxcode; code++)
  {
    rightchar = calcbucketboundsparts(&bucketspec,
                                      bcktab,
                                      code,
                                      maxcode,
                                      bucketspec2->partwidth,
                                      rightchar,
                                      bucketspec2->numofchars);
    accubucketsize += bucketspec.nonspecialsinbucket;
    if (rightchar == 0)
    {
      bucketspec2->subbuckettab[currentchar]
                               [bucketspec2->numofchars-1].bucketend
        = accubucketsize;
      accubucketsize += bucketspec.specialsinbucket;
      bucketspec2->superbuckettab[currentchar].bucketend = accubucketsize;
      currentchar++;
    } else
    {
      gt_assert(bucketspec.specialsinbucket == 0);
      bucketspec2->subbuckettab[currentchar]
                               [rightchar-1].bucketend
        = accubucketsize;
    }
  }
}

static void fillanysubbuckets(GtBucketspec2 *bucketspec2,
                              const Bcktab *bcktab)
{
  Codetype code2, maxcode;
  unsigned int idx, rightchar = 0, currentchar = 0;
  Seqpos rightbound, *specialchardist;

  maxcode = bcktab_numofallcodes(bcktab) - 1;
  bucketspec2->expandfactor
    = (Codetype) pow((double) bucketspec2->numofchars,
                     (double) (bucketspec2->prefixlength-2));
  bucketspec2->expandfillsum = bcktab_filltable(bcktab,2U);
#ifdef SHOWBUCKETSPEC2
  showexpandcode(bucketspec2,bucketspec2->prefixlength);
#endif
  specialchardist = gt_malloc(sizeof(*specialchardist) *
                              bucketspec2->numofchars);
  for (idx = 0; idx<bucketspec2->numofchars; idx++)
  {
    specialchardist[idx] = 0;
  }
  leftcontextofspecialchardist(specialchardist,bucketspec2->encseq,
                               bucketspec2->readmode);
  for (code2 = 0; code2 < (Codetype) bucketspec2->numofcharssquared; code2++)
  {
    Codetype ecode = expandtwocharcode(code2,bucketspec2);
    gt_assert(ecode / bucketspec2->expandfactor == code2);
    rightbound = calcbucketrightbounds(bcktab,
                                       ecode,
                                       maxcode,
                                       bucketspec2->partwidth);
    rightchar = (unsigned int) ((code2+1) % bucketspec2->numofchars);
    gt_assert((Codetype) currentchar == code2 / bucketspec2->numofchars);
    if (rightchar == 0)
    {
      gt_assert(rightbound >= specialchardist[currentchar]);
      gt_assert((Codetype) (bucketspec2->numofchars-1) ==
                code2 % bucketspec2->numofchars);
      bucketspec2->subbuckettab[currentchar]
                               [bucketspec2->numofchars-1].bucketend
        = rightbound - specialchardist[currentchar];
      bucketspec2->superbuckettab[currentchar].bucketend = rightbound;
      currentchar++;
    } else
    {
      gt_assert((Codetype) (rightchar-1) == code2 % bucketspec2->numofchars);
      bucketspec2->subbuckettab[currentchar][rightchar-1].bucketend
        = rightbound;
    }
  }
  gt_free(specialchardist);
}

GtBucketspec2 *gt_bucketspec2_new(const Bcktab *bcktab,
                                  const Encodedsequence *encseq,
                                  Readmode readmode,
                                  Seqpos partwidth,
                                  unsigned int numofchars)
{
  GtBucketspec2 *bucketspec2;
  unsigned int idx;

  gt_assert(numofchars > 0);
  bucketspec2 = gt_malloc(sizeof(*bucketspec2));
  bucketspec2->partwidth = partwidth;
  bucketspec2->prefixlength = bcktab_prefixlength(bcktab);
  bucketspec2->numofchars = numofchars;
  bucketspec2->numofcharssquared = numofchars * numofchars;
  bucketspec2->encseq = encseq;
  bucketspec2->readmode = readmode;
  bucketspec2->order = gt_malloc(sizeof(*bucketspec2->order) * numofchars);
  bucketspec2->superbuckettab
    = gt_malloc(sizeof(*bucketspec2->superbuckettab) * numofchars);
  gt_array2dim_malloc(bucketspec2->subbuckettab,(unsigned long) numofchars,
                      (unsigned long) numofchars);
  if (bucketspec2->prefixlength == 2U)
  {
    fill2subbuckets(bucketspec2,bcktab);
  } else
  {
    fillanysubbuckets(bucketspec2,bcktab);
  }
  for (idx = 0; idx<numofchars; idx++)
  {
    bucketspec2->order[idx] = idx;
  }
  gt_qsort_r(bucketspec2->order,(size_t) numofchars,
             sizeof (*bucketspec2->order),bucketspec2,
             comparesuperbucketsizes);
  resetsorted(bucketspec2);
#ifdef SHOWBUCKETSPEC2
  showbucketspec2(bucketspec2);
#endif
  determinehardwork(bucketspec2);
  resetsorted(bucketspec2);
  return bucketspec2;
}

static void forwardderive(bool check,
                          const GtBucketspec2 *bucketspec2,
                          GT_UNUSED Seqpos *suftab,
                          Seqpos **targetptr,
                          unsigned int source,
                          Seqpos *idx)
{
  Seqpos startpos;
  GtUchar cc;

  gt_assert (idx < targetptr[source]);
  for (; idx < targetptr[source]; idx++)
  {
    startpos = *idx;
    if (startpos > 0)
    {
      cc = getencodedchar(bucketspec2->encseq,startpos-1,bucketspec2->readmode);
      /*printf("fwd: superbucket[%u].sorted = %s\n",(unsigned int) cc,
                        bucketspec2->superbuckettab[cc].sorted ? "true" :
                                                                 "false"); */
      if (ISNOTSPECIAL(cc) && !bucketspec2->superbuckettab[cc].sorted)
      {
        if (check)
        {
          gt_assert(*(targetptr[cc]) == startpos - 1);
        } else
        {
          /*printf("fwd: suftab[%lu]=%lu from idx=%lu\n",
                        (unsigned long) (targetptr[cc] - suftab),
                        (unsigned long) (startpos-1),
                        (unsigned long) (idx - suftab)); */
          *(targetptr[cc]) = startpos - 1;
        }
        targetptr[cc]++;
      }
    }
  }
}

static void backwardderive(bool check,
                           const GtBucketspec2 *bucketspec2,
                           GT_UNUSED Seqpos *suftab,
                           Seqpos **targetptr,
                           unsigned int source,
                           Seqpos *idx)
{
  Seqpos startpos;
  GtUchar cc;

  gt_assert (idx > targetptr[source]);
  for (; idx > targetptr[source]; idx--)
  {
    startpos = *idx;
    if (startpos > 0)
    {
      cc = getencodedchar(bucketspec2->encseq,startpos-1,bucketspec2->readmode);
      /*printf("back: superbucket[%u].sorted = %s\n",(unsigned int) cc,
                        bucketspec2->superbuckettab[cc].sorted ? "true" :
                                                                 "false");*/
      if (ISNOTSPECIAL(cc) && !bucketspec2->superbuckettab[cc].sorted)
      {
        /*
        if (suftab[targetptr[cc]] != startpos - 1)
        {
          fprintf(stderr,"targetptr[%u]=%lu: suftab = %lu != "
                         "%lu = startpos - 1\n",
                         cc,
                         (unsigned long) targetptr[cc],
                         (unsigned long) suftab[targetptr[cc]],
                         (unsigned long) (startpos-1));
          exit(EXIT_FAILURE);
        }
        */
        if (check)
        {
          gt_assert(*(targetptr[cc]) == startpos - 1);
        } else
        {
          /*printf("back: suftab[%lu]=%lu from idx=%lu\n",
                            (unsigned long) (targetptr[cc] - suftab),
                            (unsigned long) (startpos-1),
                            (unsigned long) (idx - suftab));*/
          *(targetptr[cc]) = startpos - 1;
        }
        targetptr[cc]--;
      }
    }
  }
}

bool gt_hardworkbeforecopysort(const GtBucketspec2 *bucketspec2,
                               Codetype code)
{
  Codetype code2;
  unsigned int idx1, idx2;

  if (bucketspec2->prefixlength > 2U)
  {
    code2 = code / bucketspec2->expandfactor;
  } else
  {
    code2 = code;
  }
  idx1 = (unsigned int) (code2 / bucketspec2->numofchars);
  idx2 = (unsigned int) (code2 % bucketspec2->numofchars);
  return bucketspec2->subbuckettab[idx1][idx2].hardworktodo;
}

void gt_copysortsuffixes(bool check,
                         const GtBucketspec2 *bucketspec2,
                         Seqpos *suftab,
                         Verboseinfo *verboseinfo)
{
  Seqpos hardwork = 0, **targetptr;
  unsigned int idx, idxsource, source, second;

#ifdef WITHSUFFIXES
  {
    const Seqpos *ptr;
    for (ptr = suftab; ptr < suftab + bucketspec2->partwidth; ptr++)
    {
      showsequenceatstartpos(stdout,
                             ISDIRREVERSE(readmode) ? false : true,
                             ISDIRCOMPLEMENT(readmode) ? true : false,
                             encseq,
                             *ptr);
    }
  }
#endif
  targetptr = gt_malloc(sizeof(*targetptr) * bucketspec2->numofchars);
  for (idxsource = 0; idxsource<bucketspec2->numofchars; idxsource++)
  {
    source = bucketspec2->order[idxsource];
    for (second = 0; second < bucketspec2->numofchars; second++)
    {
      if (!bucketspec2->subbuckettab[source][second].sorted && source != second)
      {
        gt_assert(bucketspec2->subbuckettab[source][second].hardworktodo);
        showverbose(verboseinfo,"hard work for %u %u",source,second);
        hardwork += getendidx(bucketspec2,source,second) -
                    getstartidx(bucketspec2,source,second);
        bucketspec2->subbuckettab[source][second].sorted = true;
      } else
      {
        gt_assert(!bucketspec2->subbuckettab[source][second].hardworktodo);
      }
    }
    if (getstartidx(bucketspec2,source,0) <
        getstartidx(bucketspec2,source,source))
    {
      for (idx = 0; idx < bucketspec2->numofchars; idx++)
      {
        targetptr[idx] = suftab + getstartidx(bucketspec2,idx,source);
      }
      forwardderive(check,
                    bucketspec2,
                    suftab,
                    targetptr,
                    source,
                    suftab + getstartidx(bucketspec2,source,0));
    }
    if (getendidx(bucketspec2,source,source) <
        getendidx(bucketspec2,source,bucketspec2->numofchars))
    {
      for (idx = 0; idx < bucketspec2->numofchars; idx++)
      {
        targetptr[idx] = suftab + getendidx(bucketspec2,idx,source) - 1;
      }
      backwardderive(check,
                     bucketspec2,
                     suftab,
                     targetptr,
                     source,
                     suftab +
                     getendidx(bucketspec2,source,bucketspec2->numofchars) - 1);
    }
    for (idx = 0; idx < bucketspec2->numofchars; idx++)
    {
      bucketspec2->subbuckettab[idx][source].sorted = true;
    }
    bucketspec2->superbuckettab[source].sorted = true;
  }
  gt_free(targetptr);
  showverbose(verboseinfo,"hardwork = " FormatSeqpos " (%.2f)",
            PRINTSeqposcast(hardwork),
            (double) hardwork/getencseqtotallength(bucketspec2->encseq));
}

void gt_bucketspec2_delete(GtBucketspec2 *bucketspec2)
{
  gt_assert(bucketspec2 != NULL);
  gt_array2dim_delete(bucketspec2->subbuckettab);
  gt_free(bucketspec2->superbuckettab);
  gt_free(bucketspec2->order);
  gt_free(bucketspec2);
}
