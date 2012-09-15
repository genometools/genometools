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
#include "sfx-sain.h"
#include "core/intbits.h"
#include "stamp.h"

#define GT_SSTARLENGTH_MAX 100

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

GtSainlabels *gt_sain_labels_new(const GtEncseq *encseq)
{
  unsigned long position,
                idx,
                nextcc = GT_UNIQUEINT(SEPARATOR),
                nextSstartypepos;
  bool nextisStype = true;
  GtEncseqReader *esr;
  GtSainlabels *sainlabels;

  esr = gt_encseq_create_reader_with_readmode(encseq,GT_READMODE_REVERSE,0);
  sainlabels = gt_malloc(sizeof *sainlabels);
  sainlabels->totalSstarlength = 0;
  sainlabels->countStype = 0;
  sainlabels->countSstartype = 0;
  sainlabels->longerthanmax = 0;
  sainlabels->totallength = gt_encseq_total_length(encseq);
  GT_INITBITTAB(sainlabels->isStype,sainlabels->totallength+1);
  GT_SETIBIT(sainlabels->isStype,sainlabels->totallength);
  nextSstartypepos = sainlabels->totallength;
  /*printf("%lu: S\n",sainlabels->totallength);*/
  for (idx = 0; idx<=(unsigned long) GT_SSTARLENGTH_MAX; idx++)
  {
    sainlabels->lendist[idx] = 0;
  }
  for (position = sainlabels->totallength-1; /* Nothing */; position--)
  {
    GtUchar cc = gt_encseq_reader_next_encoded_char(esr);
    bool currentisStype;
    unsigned long currentcc
      = ISSPECIAL(cc) ? GT_UNIQUEINT(cc) : (unsigned long) cc;

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
                    (double) sainlabels->countSstartype/sainlabels->countStype);
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

static bool gt_sain_labels_isStype(const GtSainlabels *sainlabels,
                                   unsigned long position)
{
  return GT_ISIBITSET(sainlabels->isStype,position) ? true : false;
}

static void gt_sain_endbuckets(unsigned long *leftborder,
                               const unsigned long *bucketsize,
                               unsigned int numofchars)
{
  unsigned int charidx;

  leftborder[0] = bucketsize[0];
  for (charidx = 1U; charidx<numofchars; charidx++)
  {
    leftborder[charidx] = leftborder[charidx-1] + bucketsize[charidx];
  }
}

static void gt_sain_startbuckets(unsigned long *leftborder,
                                 const unsigned long *bucketsize,
                                 unsigned int numofchars)
{
  unsigned int charidx;

  leftborder[0] = 0;
  for (charidx = 1U; charidx<numofchars; charidx++)
  {
    leftborder[charidx] = leftborder[charidx-1] + bucketsize[charidx-1];
  }
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
    cc1 = ISSPECIAL(cc) ? GT_UNIQUEINT(cc) : (unsigned long) cc;
    cc = gt_encseq_get_encoded_char(encseq,start2,GT_READMODE_FORWARD);
    cc2 = ISSPECIAL(cc) ? GT_UNIQUEINT(cc) : (unsigned long) cc;
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

void gt_sain_sortstarsuffixes(const GtEncseq *encseq)
{
  unsigned int charidx;
  const unsigned int numofchars = gt_encseq_alphabetnumofchars(encseq);
  const unsigned long specialcharacters = gt_encseq_specialcharacters(encseq);
  unsigned long position, idx, regularpositions, currentname = 0,
                previouspos = ULONG_MAX,
                *bucketsize, *leftborder, *suftab;
  GtSainlabels *sainlabels;

  bucketsize = gt_malloc(sizeof (*bucketsize) * numofchars);
  for (charidx = 0; charidx<numofchars; charidx++)
  {
    bucketsize[charidx] = gt_encseq_charcount(encseq,(GtUchar) charidx);
  }
  leftborder = gt_malloc(sizeof (*leftborder) * numofchars);
  gt_sain_endbuckets(leftborder,bucketsize,numofchars);
  sainlabels = gt_sain_labels_new(encseq);
  gt_assert(sainlabels->totallength >= specialcharacters);
  regularpositions = sainlabels->totallength + 1 - specialcharacters;
  suftab = gt_malloc(sizeof (*suftab) * regularpositions);
  for (idx = 0; idx < regularpositions; idx++)
  {
    suftab[idx] = ULONG_MAX;
  }
  for (position = 0; position < sainlabels->totallength; position++)
  {
    if (gt_sain_labels_isSstartype(sainlabels,position))
    {
      GtUchar cc = gt_encseq_get_encoded_char(encseq,
                                              position,
                                              GT_READMODE_FORWARD);
      if (ISNOTSPECIAL(cc))
      {
        unsigned long putidx;

        gt_assert(leftborder[cc] > 0);
        putidx = --leftborder[cc];
        gt_assert(putidx + 1 < regularpositions);
        suftab[putidx+1] = position;
        /*printf("Sstar.suftab[%lu]=%lu\n",putidx+1,position);*/
      }
    }
  }
  suftab[0] = sainlabels->totallength;
  /*printf("Sstar.suftab[0]=%lu\n",sainlabels->totallength);*/
  gt_sain_startbuckets(leftborder,bucketsize,numofchars);
  for (idx = 0; idx < regularpositions; idx++)
  {
    position = suftab[idx];
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
  gt_sain_endbuckets(leftborder,bucketsize,numofchars);
  for (idx = regularpositions - 1; /* Nothing */; idx--)
  {
    position = suftab[idx];
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
  for (idx = 0; idx < regularpositions; idx++)
  {
    position = suftab[idx];
    if (position != ULONG_MAX &&
        gt_sain_labels_isSstartype(sainlabels,position))
    {
      if (previouspos != ULONG_MAX)
      {
        int cmp = gt_sain_compareSstarstrings(encseq,
                                              sainlabels,
                                              previouspos,
                                              position);
        gt_assert(cmp != 1);
        if (cmp == -1)
        {
          currentname++;
        }
      }
      /*printf("Sstar %lu with label %lu\n",position,currentname);*/
      previouspos = position;
    }
  }
  gt_sain_labels_show(sainlabels);
  printf("number of names: %lu\n",currentname);
  gt_sain_labels_delete(sainlabels);
  gt_free(leftborder);
  gt_free(bucketsize);
  gt_free(suftab);
}
