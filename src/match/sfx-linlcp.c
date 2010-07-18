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
#include "core/chardef.h"
#include "core/ma_api.h"

#include "core/encseq.h"
#include "compressedtab.h"
#include "sfx-linlcp.h"

unsigned long *gt_lcp13_manzini(const GtEncseq *encseq,
                                GtReadmode readmode,
                                unsigned long partwidth,
                                unsigned long totallength,
                                const Suffixsortspace *sortedsuffixes,
                                Compressedtable *inversesuftab)
{
  unsigned long pos, lcpvalue = 0, *lcptab;

  lcptab = gt_malloc(sizeof (unsigned long) * partwidth);
  lcptab[0] = 0;
  for (pos=0; pos <= totallength; pos++)
  {
    unsigned long fillpos = compressedtable_get(inversesuftab,pos);
    if (fillpos > 0 && fillpos < partwidth)
    {
      unsigned long previousstart = suffixptrget3(sortedsuffixes,fillpos-1);
      while (pos+lcpvalue < totallength &&
             previousstart+lcpvalue < totallength)
      {
        GtUchar cc1, cc2;

        cc1 = gt_encseq_get_encoded_char(encseq,pos+lcpvalue,readmode);
        cc2 = gt_encseq_get_encoded_char(encseq,previousstart+lcpvalue,
                                                readmode);
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

static unsigned long *computeocclesstab(const GtEncseq *encseq)
{
  unsigned long *occless, numofchars, idx;

  numofchars = (unsigned long)
                  gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
  occless = gt_malloc(sizeof (unsigned long) * numofchars);
  occless[0] = 0;
  for (idx = 1UL; idx < numofchars; idx++)
  {
    occless[idx] = occless[idx-1] +
                   gt_encseq_charcount(encseq,(GtUchar) (idx-1));
  }
  return occless;
}

/* for computing the ranknext-values of special positions, we only
   need the values inversesuftab[range.end] in this order,
   where range is a special range
   Now, if range.end = suffixarray[i] for some i, then
   inversesuftab[range.end] = inversesuftab[suffixarray[i]] = i.
   Thus, in case where the inversesuftab is not available,
   we obtain these values by the following function:
*/

static void setrelevantfrominversetab(Compressedtable *rightposinverse,
                                      const GtEncseq *encseq,
                                      GtReadmode readmode,
                                      const Suffixsortspace *sortedsuffixes,
                                      unsigned long partwidth)
{
  if (gt_encseq_has_specialranges(encseq))
  {
    unsigned long idx, pos;

    for (idx = 0; idx < partwidth; idx++)
    {
      pos = suffixptrget3(sortedsuffixes,idx);
      if (pos > 0)
      {
        GtUchar cc = gt_encseq_get_encoded_char(encseq,pos-1,readmode);
        if (ISSPECIAL(cc))
        {
          compressedtable_update(rightposinverse,pos,idx);
        }
      }
    }
  }
}

static unsigned long *fillrightofpartwidth(
                                         const Compressedtable *rightposinverse,
                                         const GtEncseq *encseq,
                                         GtReadmode readmode,
                                         unsigned long partwidth,
                                         unsigned long totallength)
{
  GtSpecialrangeiterator *sri;
  GtRange range;
  unsigned long realspecialranges, *rightofpartwidth = NULL;
  unsigned long countranges = 0, nextrightofpartwidth = 0;

  realspecialranges = gt_encseq_realspecialranges(encseq);
  sri = gt_specialrangeiterator_new(encseq,
                                GT_ISDIRREVERSE(readmode) ? false : true);
  while (gt_specialrangeiterator_next(sri,&range))
  {
    if (range.end < partwidth)
    {
      countranges++;
    } else
    {
      if (nextrightofpartwidth == 0)
      {
        size_t allocsize =
                     sizeof (unsigned long) * (realspecialranges - countranges);
        rightofpartwidth = gt_malloc(allocsize);
        printf("allocated %lu bytes for rightofpartwidth (%.2f)\n",
                 (unsigned long) allocsize,
                 (double) allocsize/totallength);
      }
      gt_assert(rightofpartwidth != NULL
                   && (unsigned long) nextrightofpartwidth <
                      (realspecialranges - countranges));
      rightofpartwidth[nextrightofpartwidth++]
        = compressedtable_get(rightposinverse,range.end);
    }
  }
  gt_specialrangeiterator_delete(sri);
  return rightofpartwidth;
}

static void inversesuffixarray2specialranknext(
                         const Compressedtable *rightposinverse,
                         Compressedtable *ranknext,
                         const GtEncseq *encseq,
                         GtReadmode readmode,
                         unsigned long partwidth,
                         unsigned long totallength)
{
  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;
    unsigned long specialcharacters, idx, *rightofpartwidth = NULL;
    unsigned long specialranklistindex, nextrightofpartwidth = 0;

    rightofpartwidth = fillrightofpartwidth(rightposinverse,
                                            encseq,
                                            readmode,
                                            partwidth,
                                            totallength);
    specialcharacters = gt_encseq_specialcharacters(encseq);
    specialranklistindex = partwidth;
    sri = gt_specialrangeiterator_new(encseq,
                                      GT_ISDIRREVERSE(readmode) ? false : true);
    nextrightofpartwidth = 0;
    while (gt_specialrangeiterator_next(sri,&range))
    {
      gt_assert(range.end <= totallength);
      for (idx = range.start; idx < range.end-1; idx++)
      {
        gt_assert(specialranklistindex < totallength);
        compressedtable_update(ranknext,specialranklistindex,
                               specialranklistindex + 1);
        specialranklistindex++;
      }
      gt_assert(specialranklistindex < totallength);
      if (range.end < partwidth)
      {
        compressedtable_update(ranknext,specialranklistindex,
                               compressedtable_get(rightposinverse,
                                                   range.end));
        /*
        printf("(2) set ranknext[%lu] = %lu = rightposinverse[%lu]\n",
                  (unsigned long) specialranklistindex,
                  (unsigned long) ranknext[specialranklistindex],
                  (unsigned long) range.end);
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
    gt_assert(specialranklistindex == totallength);
    gt_specialrangeiterator_delete(sri);
  }
}

static unsigned long sa2ranknext(Compressedtable *ranknext,
                                 const GtEncseq *encseq,
                                 GtReadmode readmode,
                                 unsigned long partwidth,
                                 unsigned long totallength,
                                 const Suffixsortspace *sortedsuffixes)
{
  unsigned long idx, pos, longest = 0, *occless;

  gt_assert(partwidth > 0);

  occless = computeocclesstab(encseq);
  /* now inveresuftab is not used any more, and thus the
     ranknext array (which points to ranknext can savely be stored */
  for (idx=0; idx < partwidth; idx++)
  {
    pos = suffixptrget3(sortedsuffixes,idx);
    if (pos > 0)
    {
      GtUchar cc = gt_encseq_get_encoded_char(encseq,pos-1, readmode);
      if (ISNOTSPECIAL(cc))
      {
        gt_assert(occless[cc] < (unsigned long) partwidth);
        /*
        printf("(2) cc=%u: set ranknext[%lu]=%lu\n",
                        (unsigned int) cc,
                        (unsigned long) occless[cc],
                        (unsigned long) idx);
        */
        compressedtable_update(ranknext,(unsigned long) occless[cc],idx);
        occless[cc]++;
      }
    } else
    {
      longest = idx;
    }
  }
  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;
    unsigned long specialidx, specialpos;

    sri = gt_specialrangeiterator_new(encseq,
                                  GT_ISDIRREVERSE(readmode) ? false : true);
    gt_assert(partwidth > 0); /* otherwise all lcps would be 0 */
    specialidx = partwidth;
    while (gt_specialrangeiterator_next(sri,&range))
    {
      for (specialpos = range.start; specialpos < range.end;
           specialpos++)
      {
        if (specialpos > 0)
        {
          GtUchar cc = gt_encseq_get_encoded_char(encseq,specialpos-1,
                                                         readmode);
          if (ISNOTSPECIAL(cc))
          {
            gt_assert(occless[cc] < (unsigned long) partwidth);
            compressedtable_update(ranknext,
                                   (unsigned long) occless[cc],
                                   specialidx);
            occless[cc]++;
          }
        } else
        {
          longest = partwidth;
        }
        specialidx++;
      }
    }
    if (gt_encseq_lengthofspecialsuffix(encseq))
    {
      compressedtable_update(ranknext,totallength-1,totallength);
    }
    gt_specialrangeiterator_delete(sri);
  }
  gt_free(occless);
  return longest;
}

Compressedtable *gt_lcp9_manzini(Compressedtable *spacefortab,
                                 const GtEncseq *encseq,
                                 GtReadmode readmode,
                                 unsigned long partwidth,
                                 unsigned long totallength,
                                 const Suffixsortspace *sortedsuffixes)
{
  unsigned long pos, previousstart, nextfillpos = 0, fillpos, lcpvalue = 0,
                previouscc1pos, previouscc2pos;
  Compressedtable *lcptab, *ranknext, *rightposinverse;
  GtUchar cc1, cc2;
  GtEncseqReader *esr1, *esr2;

  if (spacefortab == NULL)
  {
    rightposinverse = ranknext
                    = compressedtable_new(totallength+1,totallength);
    compressedtable_update(ranknext,totallength,totallength);
    setrelevantfrominversetab(rightposinverse,encseq,readmode,sortedsuffixes,
                              partwidth);
  } else
  {
    rightposinverse = ranknext = spacefortab;
  }
  inversesuffixarray2specialranknext(rightposinverse,ranknext,
                                     encseq,
                                     readmode,
                                     partwidth,
                                     totallength);
  fillpos = sa2ranknext(ranknext,encseq,readmode,partwidth,totallength,
                        sortedsuffixes);
  printf("longest=%lu\n",fillpos);
  lcptab = ranknext;
  /* now ranknext and lcptab point to the same memory area. After reading
     ranknext at position fillpos, the same cell is used for storing
     the determined lcp-value */
  /* exploit the fact, that pos + lcpvalue is monotone */
  esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode,0);
  cc1 = gt_encseq_reader_next_encoded_char(esr1);
  previouscc1pos = 0;
  esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  previouscc2pos = totallength;
  cc2 = 0;
  for (pos = 0; pos < totallength; pos++)
  {
    if (pos < totallength - 1)
    {
      nextfillpos = compressedtable_get(ranknext,fillpos);
    }
    if (fillpos > 0 && fillpos - 1 < partwidth)
    {
      previousstart = suffixptrget3(sortedsuffixes,fillpos-1);
      while (pos+lcpvalue < totallength &&
             previousstart+lcpvalue < totallength)
      {
        gt_assert(pos + lcpvalue >= previouscc1pos);
        while (previouscc1pos < pos + lcpvalue)
        {
          previouscc1pos++;
          cc1 = gt_encseq_reader_next_encoded_char(esr1);
        }
        if (ISSPECIAL(cc1))
        {
          break;
        }
        if (previousstart+lcpvalue < previouscc2pos ||
            previousstart+lcpvalue > previouscc2pos+1)
        {
          previouscc2pos = previousstart+lcpvalue;
          gt_encseq_reader_reinit_with_readmode(esr2, encseq, readmode,
                                                previouscc2pos);
          cc2 = gt_encseq_reader_next_encoded_char(esr2);
        } else
        {
          if (previousstart+lcpvalue == previouscc2pos+1)
          {
            previouscc2pos++;
            cc2 = gt_encseq_reader_next_encoded_char(esr2);
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
  gt_encseq_reader_delete(esr1);
  gt_encseq_reader_delete(esr2);
  return lcptab;
}
