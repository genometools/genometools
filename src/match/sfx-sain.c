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

typedef struct
{
  GtBitsequence *isStype;
  unsigned long countStype,
                countSstartype,
                totalSstarlength,
                totallength,
                longerthanmax,
                lendist[GT_SSTARLENGTH_MAX+1];
} GtSaininfo;

/*
   Abstract function from encseq with the following access functions
   gt_encseq_total_length
   gt_encseq_get_encoded_char
   gt_encseq_specialcharacters
   gt_encseq_alphabetnumofchars
   gt_encseq_charcount
*/

static GtSaininfo *gt_sain_info_new(const GtEncseq *encseq)
{
  unsigned long position,
                idx,
                nextcc,
                nextSstartypepos;
  bool nextisStype = true;
  GtEncseqReader *esr;
  GtSaininfo *saininfo;

  saininfo = gt_malloc(sizeof *saininfo);
  saininfo->totalSstarlength = 0;
  saininfo->countStype = 0;
  saininfo->countSstartype = 0;
  saininfo->longerthanmax = 0;
  saininfo->totallength = gt_encseq_total_length(encseq);
  GT_INITBITTAB(saininfo->isStype,saininfo->totallength+1);
  GT_SETIBIT(saininfo->isStype,saininfo->totallength);
  nextSstartypepos = saininfo->totallength;
  nextcc = GT_UNIQUEINT(saininfo->totallength);
  /*printf("%lu: S\n",saininfo->totallength);*/
  for (idx = 0; idx <= (unsigned long) GT_SSTARLENGTH_MAX; idx++)
  {
    saininfo->lendist[idx] = 0;
  }
  esr = gt_encseq_create_reader_with_readmode(encseq,GT_READMODE_REVERSE,0);
  for (position = saininfo->totallength-1; /* Nothing */; position--)
  {
    GtUchar cc = gt_encseq_reader_next_encoded_char(esr);
    bool currentisStype;
    unsigned long currentcc
      = ISSPECIAL(cc) ? GT_UNIQUEINT(position) : (unsigned long) cc;

    if (position < saininfo->totallength-1 &&
        (currentcc < nextcc || (currentcc == nextcc && nextisStype)))
    {
      saininfo->countStype++;
      currentisStype = true;
      GT_SETIBIT(saininfo->isStype,position);
      /*printf("%lu: S\n",position);*/
    } else
    {
      currentisStype = false;
      /*printf("%lu: L\n",position);*/
    }
    if (!currentisStype && nextisStype)
    {
      unsigned long currentlen;

      saininfo->countSstartype++;
      gt_assert(position < nextSstartypepos);
      currentlen = nextSstartypepos - position;
      saininfo->totalSstarlength += currentlen;
      /*printf("Sstar: %lu\n",position+1);*/
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
  gt_encseq_reader_delete(esr);
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

static void gt_sain_info_show(const GtSaininfo *saininfo)
{
  unsigned long idx;

  printf("S-type: %lu (%.2f)\n",saininfo->countStype,
                      (double) saininfo->countStype/saininfo->totallength);
  printf("Sstar-type: %lu (%.2f)\n",saininfo->countSstartype,
                   (double) saininfo->countSstartype/saininfo->totallength);
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
}

static bool gt_sain_info_isSstartype(const GtSaininfo *saininfo,
                                       unsigned long position)
{
  gt_assert(position <= saininfo->totallength);
  return position == saininfo->totallength ||
         (position > 0 &&
         GT_ISIBITSET(saininfo->isStype,position) &&
         !GT_ISIBITSET(saininfo->isStype,position-1)) ? true : false;
}

unsigned long countDcriticalsubstrings (const GtSaininfo *saininfo,
                                        unsigned long d)
{
  unsigned long int i = 0, j = 0;

  gt_assert(d >= 2UL);
  while (i < saininfo->totallength)
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

static bool gt_sain_info_isStype(const GtSaininfo *saininfo,
                                   unsigned long position)
{
  return GT_ISIBITSET(saininfo->isStype,position) ? true : false;
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

static void insertSstarsuffixes(const GtSaininfo *saininfo,
                                const GtEncseq *encseq,
                                unsigned long *suftab,
                                unsigned long *leftborder,
                                GT_UNUSED unsigned long regularpositions)
{
  unsigned long position;

  for (position = 0; position < saininfo->totallength; position++)
  {
    if (gt_sain_info_isSstartype(saininfo,position))
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

static void induceLtypesuffixes(const GtSaininfo *saininfo,
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
        position > 0 && !gt_sain_info_isStype(saininfo,position-1))
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

static void induceStypesuffixes(const GtSaininfo *saininfo,
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

        gt_assert(gt_sain_info_isStype(saininfo,range.start-1));
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
        position > 0 && gt_sain_info_isStype(saininfo,position-1))
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

static void moveSstar2front(const GtSaininfo *saininfo,
                            unsigned long *suftab,
                            unsigned long regularpositions,
                            unsigned long suftabentries)
{
  unsigned long ridx, position;

  for (ridx = 0; ridx < regularpositions; ridx++)
  {
    position = suftab[ridx];
    gt_assert(position != ULONG_MAX);
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
      gt_assert (widx < ridx && ridx < regularpositions);
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
  gt_assert(suftabentries > 0);
  sain_setundefined(suftab,saininfo->countSstartype, suftabentries-1);
}

static int gt_sain_compareSstarstrings(const GtEncseq *encseq,
                                       const GtSaininfo *saininfo,
                                       unsigned long start1,
                                       unsigned long start2)
{
  bool firstcmp = true;
  gt_assert(start1 <= saininfo->totallength &&
            start2 <= saininfo->totallength && start1 != start2);

  while (true)
  {
    GtUchar cc;
    unsigned long cc1, cc2;

    if (start1 == saininfo->totallength)
    {
      return -1;
    }
    if (start2 == saininfo->totallength)
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
    if (gt_sain_info_isStype(saininfo,start1))
    {
      if (gt_sain_info_isStype(saininfo,start2))
      {
        if (firstcmp)
        {
          start1++;
          start2++;
        } else
        {
          if (gt_sain_info_isSstartype(saininfo,start1))
          {
            if (gt_sain_info_isSstartype(saininfo,start2))
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
            if (gt_sain_info_isSstartype(saininfo,start2))
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
      if (gt_sain_info_isStype(saininfo,start2))
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

static void assignSstarnames(const GtSaininfo *saininfo,
                             const GtEncseq *encseq,
                             unsigned long *suftab,
                             GT_UNUSED unsigned long suftabentries)
{
  unsigned long idx, previouspos = suftab[0], currentname = 0;
  gt_assert(gt_sain_info_isSstartype(saininfo,previouspos));
  for (idx = 1UL; idx < saininfo->countSstartype; idx++)
  {
    int cmp;
    unsigned long position = suftab[idx];

    gt_assert(gt_sain_info_isSstartype(saininfo,position));
    cmp = gt_sain_compareSstarstrings(encseq,saininfo,previouspos,position);
    gt_assert(cmp != 1);
    if (cmp == -1)
    {
      currentname++;
    }
    /*printf("%03lu: Sstar %lu with name %lu\n",
                  saininfo->countSstartype + GT_DIV2(position),
                  position,currentname);*/
    gt_assert(saininfo->countSstartype + GT_DIV2(position) <
              suftabentries);
    suftab[saininfo->countSstartype + GT_DIV2(position)] = currentname;
    previouspos = position;
  }
  printf("number of names: %lu\n",currentname);
}

static void movenames2front(const GtSaininfo *saininfo,
                            unsigned long *suftab,
                            GT_UNUSED unsigned long suftabentries)
{
  unsigned long ridx, widx,
                maxridx = saininfo->countSstartype +
                          GT_DIV2(saininfo->totallength);
  for (ridx = widx = saininfo->countSstartype; ridx <= maxridx; ridx++)
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
  gt_assert (widx == GT_MULT2(saininfo->countSstartype));
}

void gt_sain_sortstarsuffixes(const GtEncseq *encseq)
{
  unsigned long specialcharacters, idx, d, regularpositions, numofchars,
                suftabentries, *bucketsize, *leftborder, *suftab;
  GtSaininfo *saininfo;

  specialcharacters = gt_encseq_specialcharacters(encseq);
  numofchars = (unsigned long) gt_encseq_alphabetnumofchars(encseq);
  bucketsize = gt_malloc(sizeof (*bucketsize) * numofchars);
  for (idx = 0; idx<numofchars; idx++)
  {
    bucketsize[idx] = gt_encseq_charcount(encseq,(GtUchar) idx);
  }
  leftborder = gt_malloc(sizeof (*leftborder) * numofchars);
  gt_sain_endbuckets(leftborder,bucketsize,numofchars);
  saininfo = gt_sain_info_new(encseq);
  gt_sain_info_show(saininfo);
  for (d=2UL; d<10UL; d++)
  {
    unsigned long critical = countDcriticalsubstrings (saininfo,d);
    printf("d=%lu,critical=%lu (%.2f)\n",d,critical,
                              (double) critical/saininfo->totallength);
  }
  gt_assert(saininfo->totallength >= specialcharacters);
  regularpositions = saininfo->totallength + 1 - specialcharacters;
  suftabentries = saininfo->countSstartype +
                  GT_DIV2(saininfo->totallength) + 1;
  if (regularpositions > suftabentries)
  {
    suftabentries = regularpositions;
  }
  suftab = gt_malloc(sizeof (*suftab) * suftabentries);
  suftab[0] = saininfo->totallength;
  sain_setundefined(suftab,1UL,regularpositions-1);
  insertSstarsuffixes(saininfo, encseq, suftab, leftborder, regularpositions);
  gt_sain_startbuckets(leftborder,bucketsize,numofchars);
  induceLtypesuffixes(saininfo, encseq, suftab, leftborder, regularpositions);
  gt_sain_endbuckets(leftborder,bucketsize,numofchars);
  induceStypesuffixes(saininfo, encseq, suftab, leftborder, regularpositions);
  moveSstar2front(saininfo,suftab,regularpositions, suftabentries);
  assignSstarnames(saininfo,encseq,suftab,suftabentries);
  movenames2front(saininfo,suftab,suftabentries);
  gt_sain_info_delete(saininfo);
  gt_free(leftborder);
  gt_free(bucketsize);
  gt_free(suftab);
}
