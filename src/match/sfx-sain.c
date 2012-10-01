/*
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include <limits.h>
#include "core/intbits.h"
#include "core/minmax.h"
#include "core/unused_api.h"
#include "sfx-linlcp.h"
#include "sfx-sain.h"

#define GT_SSTARLENGTH_MAX 50

typedef struct
{
  unsigned long totallength,
                specialcharacters,
                numofchars,
                startoccupied,
                *bucketsize,
                *leftborder;
  union
  {
    const GtEncseq *encseq;
    const unsigned long *array;
  } seq;
  bool leftborderpoints2suftab,
       bucketsizepoints2suftab,
       hasencseq;
} GtSainseq;

static GtSainseq *gt_sain_seq_new_from_encseq(const GtEncseq *encseq)
{
  unsigned long idx;
  GtSainseq *sainseq = gt_malloc(sizeof *sainseq);

  sainseq->hasencseq = true;
  sainseq->seq.encseq = encseq;
  sainseq->totallength = gt_encseq_total_length(encseq);
  sainseq->specialcharacters = gt_encseq_specialcharacters(encseq);
  gt_assert(sainseq->totallength >= sainseq->specialcharacters);
  sainseq->numofchars = (unsigned long) gt_encseq_alphabetnumofchars(encseq);
  sainseq->bucketsize = gt_malloc(sizeof (*sainseq->bucketsize) *
                                  sainseq->numofchars);
  sainseq->leftborder = gt_malloc(sizeof (*sainseq->leftborder) *
                                  sainseq->numofchars);
  sainseq->leftborderpoints2suftab = false;
  sainseq->bucketsizepoints2suftab = false;
  for (idx = 0; idx<sainseq->numofchars; idx++)
  {
    sainseq->bucketsize[idx] = gt_encseq_charcount(encseq,(GtUchar) idx);
  }
  return sainseq;
}

static void gt_sain_determinebucketsize(GtSainseq *sainseq)
{
  unsigned long idx;

  for (idx = 0; idx < sainseq->numofchars; idx++)
  {
    sainseq->bucketsize[idx] = 0;
  }
  for (idx = 0; idx<sainseq->totallength; idx++)
  {
    gt_assert(sainseq->seq.array[idx] < sainseq->numofchars);
    sainseq->bucketsize[sainseq->seq.array[idx]]++;
  }
}

static GtSainseq *gt_sain_seq_new_from_array(unsigned long *arr,
                                             unsigned long len,
                                             unsigned long numofchars,
                                             unsigned long *suftab,
                                             unsigned long suftabentries)
{
  unsigned long maxused = GT_MULT2(len);
  GtSainseq *sainseq = gt_malloc(sizeof *sainseq);

  sainseq->hasencseq = false;
  sainseq->seq.array = arr;
  sainseq->totallength = len;
  sainseq->specialcharacters = 0;
  sainseq->numofchars = numofchars;
  sainseq->leftborderpoints2suftab = false;
  sainseq->bucketsizepoints2suftab = false;
  sainseq->startoccupied = suftabentries;
  if (maxused < suftabentries)
  {
    if (suftabentries - maxused >= numofchars)
    {
      printf("bucketsize: reclaim suftab at %lu for %lu elements\n",
              suftabentries - numofchars,numofchars);
      sainseq->bucketsize = suftab + suftabentries - numofchars;
      sainseq->bucketsizepoints2suftab = true;
      sainseq->startoccupied = suftabentries - numofchars;
    }
    if (suftabentries - maxused >= GT_MULT2(numofchars))
    {
      printf("leftborder: reclaim suftab at %lu for %lu elements\n",
              suftabentries - GT_MULT2(numofchars),GT_MULT2(numofchars));
      sainseq->leftborder = suftab + suftabentries - GT_MULT2(numofchars);
      sainseq->leftborderpoints2suftab = true;
      sainseq->startoccupied = suftabentries - GT_MULT2(numofchars);
    }
  }
  if (!sainseq->leftborderpoints2suftab)
  {
    sainseq->leftborder = gt_malloc(sizeof (*sainseq->leftborder) * numofchars);
    printf("leftborder requires %lu bytes\n",
           (unsigned long) sizeof (*sainseq->leftborder) * numofchars);
  }
  if (!sainseq->bucketsizepoints2suftab)
  {
    sainseq->bucketsize = gt_malloc(sizeof (*sainseq->bucketsize) * numofchars);
    printf("bucketsize requires %lu bytes\n",
           (unsigned long) sizeof (*sainseq->bucketsize) * numofchars);
  }
  gt_sain_determinebucketsize(sainseq);
  return sainseq;
}

static unsigned long gt_sain_seq_delete(GtSainseq *sainseq)
{
  unsigned long ret = ULONG_MAX;

  if (sainseq != NULL)
  {
    if (!sainseq->leftborderpoints2suftab)
    {
      gt_free(sainseq->leftborder);
    }
    if (!sainseq->bucketsizepoints2suftab)
    {
      gt_free(sainseq->bucketsize);
    }
    if (sainseq->bucketsizepoints2suftab ||
        sainseq->bucketsizepoints2suftab)
    {
      ret = sainseq->startoccupied;
    }
    gt_free(sainseq);
  }
  return ret;
}

static unsigned long gt_sain_seq_getchar(const GtSainseq *sainseq,
                                         unsigned long position)
{
  gt_assert(position < sainseq->totallength);
  if (sainseq->hasencseq)
  {
    GtUchar cc = gt_encseq_get_encoded_char(sainseq->seq.encseq,
                                            position,
                                            GT_READMODE_FORWARD);
    return ISSPECIAL(cc) ? GT_UNIQUEINT(position) : (unsigned long) cc;
  } else
  {
    return sainseq->seq.array[position];
  }
}

static unsigned long gt_sain_seq_getchar_nospecial(const GtSainseq *sainseq,
                                                   unsigned long position)
{
  gt_assert(position < sainseq->totallength);
  if (sainseq->hasencseq)
  {
    return (unsigned long)
           gt_encseq_get_encoded_char_nospecial(sainseq->seq.encseq,
                                                position,
                                                GT_READMODE_FORWARD);
  } else
  {
    return sainseq->seq.array[position];
  }
}

typedef struct
{
  GtBitsequence *isStype;
  unsigned long countStype,
                countSstartype,
                totalSstarlength,
                longerthanmax,
                lendist[GT_SSTARLENGTH_MAX+1];
  GtSainseq *sainseq;
} GtSaininfo;

static GtSaininfo *gt_sain_info_new(GtSainseq *sainseq)
{
  unsigned long position,
                idx,
                nextcc,
                nextSstartypepos;
  bool nextisStype = true;
  GtSaininfo *saininfo;

  saininfo = gt_malloc(sizeof *saininfo);
  saininfo->sainseq = sainseq;
  saininfo->totalSstarlength = 0;
  saininfo->countStype = 1UL;
  saininfo->countSstartype = 0;
  saininfo->longerthanmax = 0;
  GT_INITBITTAB(saininfo->isStype,saininfo->sainseq->totallength+1);
  GT_SETIBIT(saininfo->isStype,saininfo->sainseq->totallength);
  nextSstartypepos = saininfo->sainseq->totallength;
  nextcc = GT_UNIQUEINT(saininfo->sainseq->totallength);
#undef SAINSHOWSTATE
#ifdef SAINSHOWSTATE
  printf("%lu: S\n",saininfo->sainseq->totallength);
#endif
  for (idx = 0; idx <= (unsigned long) GT_SSTARLENGTH_MAX; idx++)
  {
    saininfo->lendist[idx] = 0;
  }
  for (position = saininfo->sainseq->totallength-1; /* Nothing */; position--)
  {
    bool currentisStype;
    unsigned long currentcc;

    currentcc = gt_sain_seq_getchar(saininfo->sainseq,position);
    if (currentcc < nextcc || (currentcc == nextcc && nextisStype))
    {
      saininfo->countStype++;
      currentisStype = true;
      GT_SETIBIT(saininfo->isStype,position);
#ifdef SAINSHOWSTATE
      printf("%lu: S\n",position);
#endif
    } else
    {
      currentisStype = false;
#ifdef SAINSHOWSTATE
      printf("%lu: L\n",position);
#endif
    }
    if (!currentisStype && nextisStype)
    {
      unsigned long currentlen;

      saininfo->countSstartype++;
      gt_assert(position < nextSstartypepos);
      currentlen = nextSstartypepos - position;
      saininfo->totalSstarlength += currentlen;
#ifdef SAINSHOWSTATE
      printf("Sstar: %lu\n",position+1);
#endif
      if (currentlen <= (unsigned long) GT_SSTARLENGTH_MAX)
      {
        saininfo->lendist[currentlen]++;
      } else
      {
        saininfo->longerthanmax++;
      }
      nextSstartypepos = position + 1;
    }
    nextisStype = currentisStype;
    nextcc = currentcc;
    if (position == 0)
    {
      break;
    }
  }
  gt_assert(GT_MULT2(saininfo->countSstartype) <=
            saininfo->sainseq->totallength);
  return saininfo;
}

static void gt_sain_info_delete(GtSaininfo *saininfo)
{
  if (saininfo != NULL)
  {
    gt_free(saininfo->isStype);
    gt_free(saininfo);
  }
}

static bool gt_sain_info_isSstartype(const GtSaininfo *saininfo,
                                       unsigned long position)
{
  gt_assert(position <= saininfo->sainseq->totallength);
  return position == saininfo->sainseq->totallength ||
         (position > 0 &&
         GT_ISIBITSET(saininfo->isStype,position) &&
         !GT_ISIBITSET(saininfo->isStype,position-1)) ? true : false;
}

#ifdef CRITICAL
unsigned long gt_sain_countDcriticalsubstrings (const GtSaininfo *saininfo,
                                                unsigned long d)
{
  unsigned long int i = 0, j = 0;

  gt_assert(d >= 2UL);
  while (i < saininfo->sainseq->totallength)
  {
    unsigned long h;
    bool isLMS = false;

    for (h = 1UL; h <= d; h++)
    {
      if (gt_sain_info_isSstartype(saininfo,i+h))
      {
        isLMS = true;
        break;
      }
    }
    if (j == 0 && !isLMS)
    {
      i += d;
      continue;
    }
    i = isLMS ? i + h : i + d;
    gt_assert(i>0);
    /*printf("crititical %lu\n",i-1);*/
    j++;
  }
  return j;
}
#endif

static void gt_sain_info_show(const GtSaininfo *saininfo)
{
  unsigned long idx;

  printf("S-type: %lu (%.2f)\n",saininfo->countStype,
                (double) saininfo->countStype/saininfo->sainseq->totallength);
  printf("Sstar-type: %lu (%.2f)\n",saininfo->countSstartype,
              (double) saininfo->countSstartype/saininfo->sainseq->totallength);
  printf("Sstar-type.length: %lu (%.2f)\n",saininfo->totalSstarlength,
              (double) saininfo->totalSstarlength/saininfo->countSstartype);
  for (idx = 0; idx <= (unsigned long) GT_SSTARLENGTH_MAX; idx++)
  {
    if (saininfo->lendist[idx] > 0)
    {
      printf("%lu %lu (%.2f)\n",idx,saininfo->lendist[idx],
               (double) saininfo->lendist[idx]/saininfo->countSstartype);
    }
  }
  if (saininfo->longerthanmax)
  {
    printf(">%d %lu (%.2f)\n",GT_SSTARLENGTH_MAX,saininfo->longerthanmax,
                 (double) saininfo->longerthanmax/saininfo->countSstartype);
  }
#ifdef CRITICAL
  {
    unsigned long d;
    for (d=2UL; d<10UL; d++)
    {
      unsigned long critical = gt_sain_countDcriticalsubstrings (saininfo,d);
      printf("d=%lu,critical=%lu (%.2f)\n",d,critical,
                        (double) critical/saininfo->sainseq->totallength);
    }
  }
#endif
}

static bool gt_sain_info_isStype(const GtSaininfo *saininfo,
                                 unsigned long position)
{
  gt_assert(position <= saininfo->sainseq->totallength);
  return GT_ISIBITSET(saininfo->isStype,position) ? true : false;
}

static void gt_sain_endbuckets(GtSainseq *sainseq)
{
  unsigned long charidx;

  sainseq->leftborder[0] = sainseq->bucketsize[0];
  for (charidx = 1UL; charidx < sainseq->numofchars; charidx++)
  {
    sainseq->leftborder[charidx] = sainseq->leftborder[charidx-1] +
                                   sainseq->bucketsize[charidx];
  }
}

static void gt_sain_startbuckets(GtSainseq *sainseq)
{
  unsigned long charidx;

  sainseq->leftborder[0] = 0;
  for (charidx = 1UL; charidx < sainseq->numofchars; charidx++)
  {
    sainseq->leftborder[charidx] = sainseq->leftborder[charidx-1] +
                                   sainseq->bucketsize[charidx-1];
  }
}

static void gt_sain_insertSstarsuffixes(const GtSaininfo *saininfo,
                                        unsigned long *suftab,
                                        GT_UNUSED unsigned long
                                                  nonspecialentries)
{
  unsigned long position;

  for (position = 0; position < saininfo->sainseq->totallength; position++)
  {
    if (gt_sain_info_isSstartype(saininfo,position))
    {
      unsigned long putidx,
                    cc = gt_sain_seq_getchar_nospecial(saininfo->sainseq,
                                                       position);
      gt_assert(saininfo->sainseq->leftborder[cc] > 0);
      putidx = --saininfo->sainseq->leftborder[cc];
      gt_assert(putidx < nonspecialentries);
      suftab[putidx] = position;
#ifdef SAINSHOWSTATE
      printf("Sstar.suftab[%lu]=%lu\n",putidx,position);
#endif
    }
  }
}

static void gt_sain_induceLtypesuffixes(const GtSaininfo *saininfo,
                                        unsigned long *suftab,
                                        unsigned long nonspecialentries)
{
  unsigned long idx;

  for (idx = 0; idx < nonspecialentries; idx++)
  {
    unsigned long position = suftab[idx];

    if (position != ULONG_MAX && position > 0)
    {
      gt_assert(position < saininfo->sainseq->totallength);
      if (!gt_sain_info_isStype(saininfo,position-1))
      {
        unsigned long cc = gt_sain_seq_getchar(saininfo->sainseq,position-1);
        if (cc < saininfo->sainseq->numofchars)
        {
          unsigned long putidx = saininfo->sainseq->leftborder[cc]++;
          gt_assert(putidx < nonspecialentries);
          suftab[putidx] = position-1;
#ifdef SAINSHOWSTATE
          printf("L-induce: suftab[%lu]=%lu\n",putidx,position-1);
#endif
        }
      }
    }
  }
}

static void gt_sain_induceStypesfromspecialranges(GT_UNUSED const GtSaininfo
                                                            *saininfo,
                                                  const GtEncseq *encseq,
                                                  unsigned long *suftab,
                                                  GT_UNUSED unsigned long
                                                            nonspecialentries)
{
  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;
    sri = gt_specialrangeiterator_new(encseq,false);
    while (gt_specialrangeiterator_next(sri,&range))
    {
      if (range.start > 0)
      {
        unsigned long putidx;
        GtUchar cc;

        gt_assert(gt_sain_info_isStype(saininfo,range.start-1));
        cc = gt_encseq_get_encoded_char(encseq,range.start-1,
                                        GT_READMODE_FORWARD);
        gt_assert(ISNOTSPECIAL(cc) && saininfo->sainseq->leftborder[cc] > 0);
        putidx = --saininfo->sainseq->leftborder[cc];
        gt_assert(putidx < nonspecialentries);
        suftab[putidx] = range.start-1;
#ifdef SAINSHOWSTATE
        printf("Srange-induce: suftab[%lu]=%lu in %d-bucket\n",putidx,
                          range.start-1,(int) cc);
#endif
      }
    }
    gt_specialrangeiterator_delete(sri);
  }
}

static void gt_sain_induceStypesuffixes(const GtSaininfo *saininfo,
                                        unsigned long *suftab,
                                        unsigned long nonspecialentries)
{
  unsigned long idx, cc;

  cc = gt_sain_seq_getchar(saininfo->sainseq,saininfo->sainseq->totallength-1);
  if (cc < saininfo->sainseq->numofchars)
  {
    unsigned long putidx = --saininfo->sainseq->leftborder[cc];
    gt_assert(putidx < nonspecialentries);
    suftab[putidx] = saininfo->sainseq->totallength-1;
#ifdef SAINSHOWSTATE
    printf("end S-induce: suftab[%lu]=%lu\n",putidx,
                                         saininfo->sainseq->totallength-1);
#endif
  }
  if (saininfo->sainseq->hasencseq)
  {
    gt_sain_induceStypesfromspecialranges(saininfo,
                                          saininfo->sainseq->seq.encseq,
                                          suftab,
                                          nonspecialentries);
  }
  if (nonspecialentries == 0)
  {
    return;
  }
  for (idx = nonspecialentries - 1; /* Nothing */; idx--)
  {
    unsigned long position = suftab[idx];

    if (position != ULONG_MAX && position > 0)
    {
      if (gt_sain_info_isStype(saininfo,position-1))
      {
        cc = gt_sain_seq_getchar(saininfo->sainseq,position-1);
        if (cc < saininfo->sainseq->numofchars)
        {
          unsigned long putidx = --saininfo->sainseq->leftborder[cc];
          gt_assert(putidx < nonspecialentries);
          suftab[putidx] = position-1;
#ifdef SAINSHOWSTATE
          printf("S-induce: suftab[%lu]=%lu\n",putidx,position-1);
#endif
        }
      }
    }
    if (idx == 0)
    {
      break;
    }
  }
}

static void gt_sain_moveSstar2front(const GtSaininfo *saininfo,
                                    unsigned long *suftab,
                                    unsigned long nonspecialentries)
{
  unsigned long ridx, position;

  for (ridx = 0; ridx < nonspecialentries; ridx++)
  {
    position = suftab[ridx];
    if (!gt_sain_info_isSstartype(saininfo,position))
    {
      break;
    }
  }
  if (ridx < saininfo->countSstartype)
  {
    unsigned long widx;

    for (widx = ridx, ridx++; /* Nothing */; ridx++)
    {
      gt_assert(widx < ridx && ridx < nonspecialentries);
      position = suftab[ridx];
      gt_assert(position != ULONG_MAX);
      if (gt_sain_info_isSstartype(saininfo,position))
      {
        suftab[widx++] = position;
        suftab[ridx] = ULONG_MAX;
        if (widx == saininfo->countSstartype)
        {
          break;
        }
      }
    }
  }
}

static int gt_sain_compare_substrings(const GtSaininfo *saininfo,
                                      bool withtype,
                                      unsigned long start1,
                                      unsigned long start2)
{
  bool firstcmp = true, previousisS = true;

  gt_assert(start1 <= saininfo->sainseq->totallength &&
            start2 <= saininfo->sainseq->totallength &&
            start1 != start2);

  while (true)
  {
    unsigned long cc1, cc2;

    if (start1 == saininfo->sainseq->totallength)
    {
      gt_assert(start1 > start2);
      return 1;
    }
    if (start2 == saininfo->sainseq->totallength)
    {
      gt_assert(start1 < start2);
      return -1;
    }
    cc1 = gt_sain_seq_getchar(saininfo->sainseq,start1);
    cc2 = gt_sain_seq_getchar(saininfo->sainseq,start2);
    if (cc1 < cc2)
    {
      return -1;
    }
    if (cc1 > cc2)
    {
      return 1;
    }
    if (withtype)
    {
      if (gt_sain_info_isStype(saininfo,start1))
      {
        if (gt_sain_info_isStype(saininfo,start2))
        {
          if (!firstcmp && !previousisS)
          {
            return 0; /* previous is L => Sstar in both stop with equality */
          }
        } else
        {
          /* S > L */
          return 1;
        }
      } else
      {
        if (gt_sain_info_isStype(saininfo,start2))
        {
          /* L < S */
          return -1;
        } else
        {
          /* L == L */
          previousisS = false;
        }
      }
      firstcmp = false;
    }
    start1++;
    start2++;
  }
}

static void gt_sain_setundefined(unsigned long *suftab,
                                 unsigned long start, unsigned long end)
{
  unsigned long idx;

  for (idx = start; idx <= end; idx++)
  {
    suftab[idx] = ULONG_MAX;
  }
}

static unsigned long gt_sain_assignSstarnames(const GtSaininfo *saininfo,
                                              unsigned long *suftab,
                                              unsigned long availableentries)
{
  unsigned long idx, previouspos, currentname = 0;

  gt_assert(availableentries > 0);
  gt_sain_setundefined(suftab,saininfo->countSstartype,availableentries-1);
  previouspos = suftab[0];
  suftab[saininfo->countSstartype + GT_DIV2(previouspos)] = 0;
  gt_assert(gt_sain_info_isSstartype(saininfo,previouspos));
  for (idx = 1UL; idx < saininfo->countSstartype; idx++)
  {
    int cmp;
    unsigned long position = suftab[idx];

    gt_assert(gt_sain_info_isSstartype(saininfo,position));
    cmp = gt_sain_compare_substrings(saininfo,true,previouspos,position);
    gt_assert(cmp != 1);
    if (cmp == -1)
    {
      currentname++;
    }
    gt_assert(saininfo->countSstartype + GT_DIV2(position) < availableentries);
    suftab[saininfo->countSstartype + GT_DIV2(position)] = currentname;
    previouspos = position;
  }
  return currentname+1;
}

static void gt_sain_movenames2front(const GtSaininfo *saininfo,
                                    unsigned long *suftab,
                                    GT_UNUSED unsigned long availableentries)
{
  unsigned long ridx, widx;
  const unsigned long maxridx = saininfo->countSstartype +
                                GT_DIV2(saininfo->sainseq->totallength);

  for (ridx = widx = saininfo->countSstartype; ridx <= maxridx; ridx++)
  {
    if (suftab[ridx] != ULONG_MAX)
    {
      if (widx < ridx)
      {
        gt_assert(widx < availableentries);
        suftab[widx++] = suftab[ridx];
      } else
      {
        gt_assert(widx == ridx);
        widx++;
      }
    }
  }
  gt_assert(widx == GT_MULT2(saininfo->countSstartype));
}

static void gt_sain_checkorder(const GtSaininfo *saininfo,
                               const unsigned long *suftab,
                               unsigned long start,
                               unsigned long end)
{
  unsigned long idx;

  for (idx = start+1; idx <= end; idx++)
  {
    int cmp
      = gt_sain_compare_substrings(saininfo,false,suftab[idx-1],suftab[idx]);

    if (cmp != -1)
    {
      fprintf(stderr,"%s: check interval [%lu,%lu] at idx=%lu: suffix "
                     "%lu >= %lu\n",__func__,start,end,idx,suftab[idx-1],
                     suftab[idx]);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

static void gt_sain_insertsortedSstarsuffixes(const GtSaininfo *saininfo,
                                              unsigned long *suftab,
                                              GT_UNUSED unsigned long
                                                        nonspecialsuffixes)
{
  unsigned long idx;

  if (saininfo->countSstartype == 0)
  {
    return;
  }
  for (idx = saininfo->countSstartype - 1; /* Nothing */; idx--)
  {
    unsigned long position = suftab[idx], putidx, cc;

    cc = gt_sain_seq_getchar_nospecial(saininfo->sainseq,position);
    gt_assert(saininfo->sainseq->leftborder[cc] > 0);
    putidx = --saininfo->sainseq->leftborder[cc];
    gt_assert(idx <= putidx);
    if (idx < putidx)
    {
      gt_assert(putidx < nonspecialsuffixes);
      suftab[putidx] = position;
      suftab[idx] = ULONG_MAX;
#ifdef SAINSHOWSTATE
      printf("insertsorted: suftab[%lu]=%lu\n",putidx,position);
      printf("insertsorted: suftab[%lu]=undef\n",idx);
#endif
    }
    if (idx == 0)
    {
      break;
    }
  }
}

static void gt_sain_filltailsuffixes(unsigned long *suftabtail,
                                     const GtEncseq *encseq)
{
  unsigned long specialcharacters = gt_encseq_specialcharacters(encseq),
                totallength = gt_encseq_total_length(encseq);

  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;
    unsigned long countspecial = 0;

    sri = gt_specialrangeiterator_new(encseq,true);
    while (gt_specialrangeiterator_next(sri,&range))
    {
      unsigned long idx;

      for (idx = range.start; idx < range.end; idx++)
      {
        gt_assert(countspecial < specialcharacters && idx < totallength);
        suftabtail[countspecial++] = idx;
      }
    }
    gt_assert(countspecial == specialcharacters);
    gt_specialrangeiterator_delete(sri);
  }
  suftabtail[specialcharacters] = totallength;
}

static void gt_sain_rec_sortsuffixes(GtSaininfo *saininfo,
                                     unsigned long *suftab,
                                     unsigned long nonspecialentries,
                                     unsigned long availableentries,
                                     unsigned long suftabentries,
                                     unsigned int sainmode)
{
  if (saininfo->countSstartype > 0)
  {
    unsigned long idx, startoccupied, numberofnames;

    gt_sain_endbuckets(saininfo->sainseq);
    gt_sain_insertSstarsuffixes(saininfo, suftab, nonspecialentries);
    gt_sain_startbuckets(saininfo->sainseq);
    gt_sain_induceLtypesuffixes(saininfo, suftab, nonspecialentries);
    gt_sain_endbuckets(saininfo->sainseq);
    gt_sain_induceStypesuffixes(saininfo, suftab, nonspecialentries);
    gt_sain_moveSstar2front(saininfo,suftab,nonspecialentries);
    numberofnames = gt_sain_assignSstarnames(saininfo,suftab,availableentries);
    gt_sain_movenames2front(saininfo,suftab,availableentries);
    gt_assert(numberofnames <= saininfo->countSstartype);
    if (numberofnames < saininfo->countSstartype)
    {
    /* Now the name sequence is in the range from
       saininfo->countSstartype .. 2 * saininfo->countSstartype - 1 */
      unsigned long position, *subseq = suftab + saininfo->countSstartype;
      GtSainseq *sainseq_rec;
      GtSaininfo *saininfo_rec;

      sainseq_rec = gt_sain_seq_new_from_array(subseq,
                                               saininfo->countSstartype,
                                               numberofnames,
                                               suftab,
                                               suftabentries);
      saininfo_rec = gt_sain_info_new(sainseq_rec);
      gt_sain_info_show(saininfo_rec);
      printf("recursively sort the named sequence of length %lu over %lu "
             "symbols (%.2f)\n",saininfo->countSstartype,numberofnames,
                         (double) numberofnames/saininfo->countSstartype);
      gt_sain_setundefined(suftab,0,saininfo->countSstartype-1);
      gt_sain_rec_sortsuffixes(saininfo_rec,suftab,
                               saininfo->countSstartype,
                               saininfo->countSstartype,
                               suftabentries,
                               sainmode);
      gt_sain_info_delete(saininfo_rec);
      startoccupied = gt_sain_seq_delete(sainseq_rec);
      gt_sain_setundefined(suftab,startoccupied,suftabentries-1);
      for (idx = 0; idx < saininfo->countSstartype; idx++)
      {
        gt_assert(saininfo->countSstartype + suftab[idx] < nonspecialentries);
        suftab[saininfo->countSstartype + suftab[idx]] = idx;
      }
      for (idx = saininfo->countSstartype, position = 0;
           position < saininfo->sainseq->totallength;
           position++)
      {
        if (gt_sain_info_isSstartype(saininfo,position))
        {
          unsigned long putidx = suftab[idx];

          gt_assert(putidx < nonspecialentries);
          suftab[putidx] = position;
          idx++;
        }
      }
    }
  }
  if (sainmode == 1U && saininfo->countSstartype > 0)
  {
    gt_sain_checkorder(saininfo,suftab,0,saininfo->countSstartype-1);
  }
  if (!saininfo->sainseq->hasencseq)
  {
    gt_sain_determinebucketsize(saininfo->sainseq);
  }
  gt_sain_setundefined(suftab,saininfo->countSstartype,availableentries-1);
  gt_sain_endbuckets(saininfo->sainseq);
  gt_sain_insertsortedSstarsuffixes(saininfo,suftab,nonspecialentries);
  gt_sain_startbuckets(saininfo->sainseq);
  gt_sain_induceLtypesuffixes(saininfo, suftab, nonspecialentries);
  gt_sain_endbuckets(saininfo->sainseq);
  gt_sain_induceStypesuffixes(saininfo, suftab, nonspecialentries);
  if (nonspecialentries > 0)
  {
    if (sainmode == 1U)
    {
      gt_sain_checkorder(saininfo,suftab,0,nonspecialentries-1);
    } else
    {
      if (saininfo->sainseq->hasencseq && sainmode == 2U)
      {
        gt_sain_filltailsuffixes(suftab + nonspecialentries,
                                 saininfo->sainseq->seq.encseq);
        gt_suftab_lightweightcheck(saininfo->sainseq->seq.encseq,
                                   GT_READMODE_FORWARD,
                                   saininfo->sainseq->totallength,
                                   suftab,
                                   NULL);
      }
    }
  }
}

/* sainmode = 1: check with own gt_sain_checkorder
   sainmode = 2: check with gt_suftab_lightweightcheck
   sainmode = 3: no check
*/

void gt_sain_sortsuffixes(const GtEncseq *encseq,unsigned int sainmode)
{
  unsigned long nonspecialentries, requiredentries, suftabentries, *suftab;
  GtSainseq *sainseq = gt_sain_seq_new_from_encseq(encseq);
  GtSaininfo *saininfo;

  saininfo = gt_sain_info_new(sainseq);
  gt_sain_info_show(saininfo);
  nonspecialentries = sainseq->totallength - sainseq->specialcharacters;
  requiredentries = saininfo->countSstartype +
                    GT_DIV2(saininfo->sainseq->totallength) + 1;
  suftabentries = MAX(nonspecialentries,requiredentries);
  if (sainmode == 2U)
  {
    suftab = gt_malloc(sizeof (*suftab) * (saininfo->sainseq->totallength+1));
  } else
  {
    suftab = gt_malloc(sizeof (*suftab) * suftabentries);
  }
  gt_sain_setundefined(suftab,0,suftabentries - 1);
  gt_sain_rec_sortsuffixes(saininfo,suftab,nonspecialentries,suftabentries,
                           suftabentries,sainmode);
  gt_sain_info_delete(saininfo);
  (void) gt_sain_seq_delete(sainseq);
  gt_free(suftab);
}
