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
#include "core/unused_api.h"
#include "sfx-sain.h"
#include "stamp.h"

#define GT_SSTARLENGTH_MAX 50

struct GtSainlabels
{
  GtBitsequence *isStype;
  unsigned long countStype,
                countSstartype,
                totalSstarlength,
                totallength,
                longerthanmax,
                lendist[GT_SSTARLENGTH_MAX+1];
};

/* 
   Abstract function from encseq with the following access functions
   gt_encseq_total_length
   gt_encseq_get_encoded_char
   gt_encseq_specialcharacters
   gt_encseq_alphabetnumofchars
   gt_encseq_charcount
*/

GtSainlabels *gt_sain_labels_new(const GtEncseq *encseq)
{
  unsigned long position,
                idx,
                nextcc,
                nextSstartypepos;
  bool nextisStype = true;
  GtEncseqReader *esr;
  GtSainlabels *sainlabels;

  sainlabels = gt_malloc(sizeof *sainlabels);
  sainlabels->totalSstarlength = 0;
  sainlabels->countStype = 0;
  sainlabels->countSstartype = 0;
  sainlabels->longerthanmax = 0;
  sainlabels->totallength = gt_encseq_total_length(encseq);
  GT_INITBITTAB(sainlabels->isStype,sainlabels->totallength+1);
  GT_SETIBIT(sainlabels->isStype,sainlabels->totallength);
  nextSstartypepos = sainlabels->totallength;
  nextcc = GT_UNIQUEINT(sainlabels->totallength);
  /*printf("%lu: S\n",sainlabels->totallength);*/
  for (idx = 0; idx <= (unsigned long) GT_SSTARLENGTH_MAX; idx++)
  {
    sainlabels->lendist[idx] = 0;
  }
  esr = gt_encseq_create_reader_with_readmode(encseq,GT_READMODE_REVERSE,0);
  for (position = sainlabels->totallength-1; /* Nothing */; position--)
  {
    GtUchar cc = gt_encseq_reader_next_encoded_char(esr);
    bool currentisStype;
    unsigned long currentcc
      = ISSPECIAL(cc) ? GT_UNIQUEINT(position) : (unsigned long) cc;

    if (position < sainlabels->totallength-1 &&
        (currentcc < nextcc || (currentcc == nextcc && nextisStype)))
    {
      sainlabels->countStype++;
      currentisStype = true;
      GT_SETIBIT(sainlabels->isStype,position);
      /*printf("%lu: S\n",position);*/
    } else
    {
      currentisStype = false;
      /*printf("%lu: L\n",position);*/
    }
    if (!currentisStype && nextisStype)
    {
      unsigned long currentlen;

      sainlabels->countSstartype++;
      gt_assert(position < nextSstartypepos);
      currentlen = nextSstartypepos - position;
      sainlabels->totalSstarlength += currentlen;
      /*printf("Sstar: %lu\n",position+1);*/
      if (currentlen <= (unsigned long) GT_SSTARLENGTH_MAX)
      {
        sainlabels->lendist[currentlen]++;
      } else
      {
        sainlabels->longerthanmax++;
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
  gt_encseq_reader_delete(esr);
  return sainlabels;
}

void gt_sain_labels_delete(GtSainlabels *sainlabels)
{
  if (sainlabels != NULL)
  {
    gt_free(sainlabels->isStype);
    gt_free(sainlabels);
  }
}

void gt_sain_labels_show(const GtSainlabels *sainlabels)
{
  unsigned long idx;

  printf("S-type: %lu (%.2f)\n",sainlabels->countStype,
                      (double) sainlabels->countStype/sainlabels->totallength);
  printf("Sstar-type: %lu (%.2f)\n",sainlabels->countSstartype,
                   (double) sainlabels->countSstartype/sainlabels->totallength);
  printf("Sstar-type.length: %lu (%.2f)\n",sainlabels->totalSstarlength,
              (double) sainlabels->totalSstarlength/sainlabels->countSstartype);
  for (idx = 0; idx <= (unsigned long) GT_SSTARLENGTH_MAX; idx++)
  {
    if (sainlabels->lendist[idx] > 0)
    {
      printf("%lu %lu (%.2f)\n",idx,sainlabels->lendist[idx],
               (double) sainlabels->lendist[idx]/sainlabels->countSstartype);
    }
  }
  if (sainlabels->longerthanmax)
  {
    printf(">%d %lu (%.2f)\n",GT_SSTARLENGTH_MAX,sainlabels->longerthanmax,
                 (double) sainlabels->longerthanmax/sainlabels->countSstartype);
  }
}

static bool gt_sain_labels_isSstartype(const GtSainlabels *sainlabels,
                                       unsigned long position)
{
  gt_assert(position <= sainlabels->totallength);
  return position == sainlabels->totallength ||
         (position > 0 &&
         GT_ISIBITSET(sainlabels->isStype,position) &&
         !GT_ISIBITSET(sainlabels->isStype,position-1)) ? true : false;
}

unsigned long countDcriticalsubstrings (const GtSainlabels *sainlabels,
                                        unsigned long d)
{
  unsigned long int i = 0, j = 0;

  gt_assert(d >= 2UL);
  while (i < sainlabels->totallength)
  {
    unsigned long h;
    bool isLMS = false;

    for (h = 1UL; h <= d; h++)
    {
      if (gt_sain_labels_isSstartype(sainlabels,i+h))
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

static bool gt_sain_labels_isStype(const GtSainlabels *sainlabels,
                                   unsigned long position)
{
  return GT_ISIBITSET(sainlabels->isStype,position) ? true : false;
}

static void gt_sain_endbuckets(unsigned long *leftborder,
                               const unsigned long *bucketsize,
                               unsigned long numofchars)
{
  unsigned long charidx;

  leftborder[0] = bucketsize[0];
  for (charidx = 1UL; charidx<numofchars; charidx++)
  {
    leftborder[charidx] = leftborder[charidx-1] + bucketsize[charidx];
  }
}

static void gt_sain_startbuckets(unsigned long *leftborder,
                                 const unsigned long *bucketsize,
                                 unsigned long numofchars)
{
  unsigned long charidx;

  leftborder[0] = 0;
  for (charidx = 1UL; charidx<numofchars; charidx++)
  {
    leftborder[charidx] = leftborder[charidx-1] + bucketsize[charidx-1];
  }
}

static void insertSstarsuffixes(const GtSainlabels *sainlabels,
                                const GtEncseq *encseq,
                                unsigned long *suftab,
                                unsigned long *leftborder,
                                unsigned long regularpositions)
{
  unsigned long position;

  for (position = 0; position < sainlabels->totallength; position++)
  {
    if (gt_sain_labels_isSstartype(sainlabels,position))
    {
      unsigned long putidx;
      GtUchar cc = gt_encseq_get_encoded_char(encseq,
                                              position,
                                              GT_READMODE_FORWARD);
      gt_assert(ISNOTSPECIAL(cc) && leftborder[cc] > 0);
      putidx = --leftborder[cc];
      gt_assert(putidx + 1 < regularpositions);
      suftab[putidx+1] = position;
      /*printf("Sstar.suftab[%lu]=%lu\n",putidx+1,position);*/
    }
  }
}

static void induceLtypesuffixes(const GtSainlabels *sainlabels,
                                const GtEncseq *encseq,
                                unsigned long *suftab,
                                unsigned long *leftborder,
                                unsigned long regularpositions)
{
  unsigned long idx;

  for (idx = 0; idx < regularpositions; idx++)
  {
    unsigned long position = suftab[idx];
    if (position != ULONG_MAX &&
        position > 0 && !gt_sain_labels_isStype(sainlabels,position-1))
    {
      GtUchar cc = gt_encseq_get_encoded_char(encseq,
                                              position-1,
                                              GT_READMODE_FORWARD);
      if (ISNOTSPECIAL(cc))
      {
        unsigned long putidx = leftborder[cc]++;
        gt_assert(putidx + 1 < regularpositions);
        suftab[putidx+1] = position-1;
        /*printf("1: suftab[%lu]=%lu\n",putidx+1,position-1);*/
      }
    }
  }
}

static void induceStypesuffixes(const GtSainlabels *sainlabels,
                                const GtEncseq *encseq,
                                unsigned long *suftab,
                                unsigned long *leftborder,
                                unsigned long regularpositions)
{
  unsigned long idx;

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

        gt_assert(gt_sain_labels_isStype(sainlabels,range.start-1));
        cc = gt_encseq_get_encoded_char(encseq,range.start-1,
                                        GT_READMODE_FORWARD);
        gt_assert (ISNOTSPECIAL(cc) && leftborder[cc] > 0);
        putidx = --leftborder[cc];
        gt_assert(putidx + 1 < regularpositions);
        /*printf("suftab[%lu]=%lu in %d-bucket\n",putidx+1,range.start-1,
                 (int) cc);*/
        suftab[putidx+1] = range.start-1;
      }
    }
    gt_specialrangeiterator_delete(sri);
  }
  for (idx = regularpositions - 1; /* Nothing */; idx--)
  {
    unsigned long position = suftab[idx];
    if (position != ULONG_MAX &&
        position > 0 && gt_sain_labels_isStype(sainlabels,position-1))
    {
      GtUchar cc = gt_encseq_get_encoded_char(encseq,
                                              position-1,
                                              GT_READMODE_FORWARD);
      if (ISNOTSPECIAL(cc))
      {
        unsigned long putidx = --leftborder[cc];
        gt_assert(putidx + 1 < regularpositions);
        suftab[putidx+1] = position-1;
        /*printf("2: suftab[%lu]=%lu\n",putidx+1,position-1);*/
      }
    }
    if (idx == 0)
    {
      break;
    }
  }
}

static void sain_setundefined(unsigned long *suftab,
                              unsigned long from, unsigned long to)
{
  unsigned long idx;

  for (idx = from; idx <= to; idx++)
  {
    suftab[idx] = ULONG_MAX;
  }
}

static void moveSstar2front(const GtSainlabels *sainlabels,
                            unsigned long *suftab,
                            unsigned long regularpositions,
                            unsigned long suftabentries)
{
  unsigned long ridx, position;

  for (ridx = 0; ridx < regularpositions; ridx++)
  {
    position = suftab[ridx];
    gt_assert(position != ULONG_MAX);
    if (!gt_sain_labels_isSstartype(sainlabels,position))
    {
      break;
    }
  }
  if (ridx < sainlabels->countSstartype)
  {
    unsigned long widx;
    for (widx = ridx, ridx++; /* Nothing */; ridx++)
    {
      gt_assert (widx < ridx && ridx < regularpositions);
      position = suftab[ridx];
      gt_assert(position != ULONG_MAX);
      if (gt_sain_labels_isSstartype(sainlabels,position))
      {
        suftab[widx++] = position;
        suftab[ridx] = ULONG_MAX;
        if (widx == sainlabels->countSstartype)
        {
          break;
        }
      }
    }
  }
  gt_assert(suftabentries > 0);
  sain_setundefined(suftab,sainlabels->countSstartype, suftabentries-1);
}

static int gt_sain_compareSstarstrings(const GtEncseq *encseq,
                                       const GtSainlabels *sainlabels,
                                       unsigned long start1,
                                       unsigned long start2)
{
  bool firstcmp = true;
  gt_assert(start1 <= sainlabels->totallength &&
            start2 <= sainlabels->totallength && start1 != start2);

  while (true)
  {
    GtUchar cc;
    unsigned long cc1, cc2;

    if (start1 == sainlabels->totallength)
    {
      return -1;
    }
    if (start2 == sainlabels->totallength)
    {
      return 1;
    }
    cc = gt_encseq_get_encoded_char(encseq,start1,GT_READMODE_FORWARD);
    cc1 = ISSPECIAL(cc) ? GT_UNIQUEINT(start1) : (unsigned long) cc;
    cc = gt_encseq_get_encoded_char(encseq,start2,GT_READMODE_FORWARD);
    cc2 = ISSPECIAL(cc) ? GT_UNIQUEINT(start2) : (unsigned long) cc;
    if (cc1 < cc2)
    {
      return -1;
    }
    if (cc1 > cc2)
    {
      return 1;
    }
    if (gt_sain_labels_isStype(sainlabels,start1))
    {
      if (gt_sain_labels_isStype(sainlabels,start2))
      {
        if (firstcmp)
        {
          start1++;
          start2++;
        } else
        {
          if (gt_sain_labels_isSstartype(sainlabels,start1))
          {
            if (gt_sain_labels_isSstartype(sainlabels,start2))
            {
              /* strings are of equal length */
              return 0;
            } else
            {
              /* first is shorter than second */
              return -1;
            }
          } else
          {
            if (gt_sain_labels_isSstartype(sainlabels,start2))
            {
              /* first is longer than second */
              return 1;
            } else
            {
              start1++;
              start2++;
            }
          }
        }
      } else
      {
        /* S > L */
        return 1;
      }
    } else
    {
      if (gt_sain_labels_isStype(sainlabels,start2))
      {
        /* L < S */
        return -1;
      } else
      {
        start1++;
        start2++;
      }
    }
    firstcmp = false;
  }
}

static void assignSstarnames(const GtSainlabels *sainlabels,
                             const GtEncseq *encseq,
                             unsigned long *suftab,
                             GT_UNUSED unsigned long suftabentries)
{
  unsigned long idx, previouspos = suftab[0], currentname = 0;
  gt_assert(gt_sain_labels_isSstartype(sainlabels,previouspos));
  for (idx = 1UL; idx < sainlabels->countSstartype; idx++)
  {
    int cmp;
    unsigned long position = suftab[idx];

    gt_assert(gt_sain_labels_isSstartype(sainlabels,position));
    cmp = gt_sain_compareSstarstrings(encseq,sainlabels,previouspos,position);
    gt_assert(cmp != 1);
    if (cmp == -1)
    {
      currentname++;
    }
    /*printf("%03lu: Sstar %lu with name %lu\n",
                  sainlabels->countSstartype + GT_DIV2(position),
                  position,currentname);*/
    gt_assert(sainlabels->countSstartype + GT_DIV2(position) <
              suftabentries);
    suftab[sainlabels->countSstartype + GT_DIV2(position)] = currentname;
    previouspos = position;
  }
  printf("number of names: %lu\n",currentname);
}

static void movenames2front(const GtSainlabels *sainlabels,
                            unsigned long *suftab,
                            GT_UNUSED unsigned long suftabentries)
{
  unsigned long ridx, widx,
                maxridx = sainlabels->countSstartype +
                          GT_DIV2(sainlabels->totallength);
  for (ridx = widx = sainlabels->countSstartype; ridx <= maxridx; ridx++)
  {
    if (suftab[ridx] != ULONG_MAX)
    {
      if (widx < ridx)
      {
        gt_assert(widx < suftabentries);
        suftab[widx++] = suftab[ridx];
      } else
      {
        gt_assert(widx == ridx);
        widx++;
      }
    }
  }
  gt_assert(widx < suftabentries);
  suftab[widx++] = 0;
  gt_assert (widx == GT_MULT2(sainlabels->countSstartype));
}

void gt_sain_sortstarsuffixes(const GtEncseq *encseq)
{
  unsigned long specialcharacters, idx, d, regularpositions, numofchars,
                suftabentries, *bucketsize, *leftborder, *suftab;
  GtSainlabels *sainlabels;

  specialcharacters = gt_encseq_specialcharacters(encseq);
  numofchars = (unsigned long) gt_encseq_alphabetnumofchars(encseq);
  bucketsize = gt_malloc(sizeof (*bucketsize) * numofchars);
  for (idx = 0; idx<numofchars; idx++)
  {
    bucketsize[idx] = gt_encseq_charcount(encseq,(GtUchar) idx);
  }
  leftborder = gt_malloc(sizeof (*leftborder) * numofchars);
  gt_sain_endbuckets(leftborder,bucketsize,numofchars);
  sainlabels = gt_sain_labels_new(encseq);
  gt_sain_labels_show(sainlabels);
  for (d=2UL; d<10UL; d++)
  {
    unsigned long critical = countDcriticalsubstrings (sainlabels,d);
    printf("d=%lu,critical=%lu (%.2f)\n",d,critical,
                              (double) critical/sainlabels->totallength);
  }
  gt_assert(sainlabels->totallength >= specialcharacters);
  regularpositions = sainlabels->totallength + 1 - specialcharacters;
  suftabentries = sainlabels->countSstartype +
                  GT_DIV2(sainlabels->totallength) + 1;
  if (regularpositions > suftabentries)
  {
    suftabentries = regularpositions;
  }
  suftab = gt_malloc(sizeof (*suftab) * suftabentries);
  suftab[0] = sainlabels->totallength;
  sain_setundefined(suftab,1UL,regularpositions-1);
  insertSstarsuffixes(sainlabels, encseq, suftab, leftborder, regularpositions);
  gt_sain_startbuckets(leftborder,bucketsize,numofchars);
  induceLtypesuffixes(sainlabels, encseq, suftab, leftborder, regularpositions);
  gt_sain_endbuckets(leftborder,bucketsize,numofchars);
  induceStypesuffixes(sainlabels, encseq, suftab, leftborder, regularpositions);
  moveSstar2front(sainlabels,suftab,regularpositions, suftabentries);
  assignSstarnames(sainlabels,encseq,suftab,suftabentries);
  movenames2front(sainlabels,suftab,suftabentries);
  gt_sain_labels_delete(sainlabels);
  gt_free(leftborder);
  gt_free(bucketsize);
  gt_free(suftab);
}
