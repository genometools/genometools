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
#include "seqpos-def.h"
#include "encseq-def.h"
#include "compressedtab.h"
#include "sfx-linlcp.h"

Seqpos *lcp13_manzini(const Encodedsequence *encseq,
                      Readmode readmode,
                      Seqpos partwidth,
                      Seqpos totallength,
                      const Seqpos *sortedsuffixes,
                      Compressedtable *inversesuftab)
{
  Seqpos pos, lcpvalue = 0, *lcptab;

  lcptab = gt_malloc(sizeof(Seqpos) * partwidth);
  lcptab[0] = 0;
  for (pos=0; pos <= totallength; pos++)
  {
    Seqpos fillpos = compressedtable_get(inversesuftab,pos);
    if (fillpos > 0 && fillpos < partwidth)
    {
      Seqpos previousstart = sortedsuffixes[fillpos-1];
      while (pos+lcpvalue < totallength &&
             previousstart+lcpvalue < totallength)
      {
        GtUchar cc1, cc2;

        cc1 = getencodedchar(encseq,pos+lcpvalue,readmode);
        cc2 = getencodedchar(encseq,previousstart+lcpvalue,readmode);
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

static unsigned long *computeocclesstab(const Encodedsequence *encseq)
{
  unsigned long *occless, numofchars, idx;

  numofchars = (unsigned long) getencseqAlphabetnumofchars(encseq);
  occless = gt_malloc(sizeof (unsigned long) * numofchars);
  occless[0] = 0;
  for (idx = 1UL; idx < numofchars; idx++)
  {
    occless[idx] = occless[idx-1] +
                   getencseqcharactercount(encseq,(GtUchar) (idx-1));
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
                                      const Encodedsequence *encseq,
                                      Readmode readmode,
                                      const Seqpos *sortedsuffixes,
                                      Seqpos partwidth)
{
  if (hasspecialranges(encseq))
  {
    Seqpos idx;

    for (idx = 0; idx < partwidth; idx++)
    {
      Seqpos pos = sortedsuffixes[idx];
      if (pos > 0)
      {
        GtUchar cc = getencodedchar(encseq,pos-1,readmode);
        if (ISSPECIAL(cc))
        {
          compressedtable_update(rightposinverse,pos,idx);
        }
      }
    }
  }
}

static Seqpos *fillrightofpartwidth(const Compressedtable *rightposinverse,
                                    const Encodedsequence *encseq,
                                    Readmode readmode,
                                    Seqpos partwidth,
                                    Seqpos totallength)
{
  Specialrangeiterator *sri;
  Sequencerange range;
  Seqpos realspecialranges, *rightofpartwidth = NULL;
  unsigned long countranges = 0, nextrightofpartwidth = 0;

  realspecialranges = getencseqrealspecialranges(encseq);
  sri = newspecialrangeiterator(encseq,ISDIRREVERSE(readmode) ? false : true);
  while (nextspecialrangeiterator(&range,sri))
  {
    if (range.rightpos < partwidth)
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
                 (double) allocsize/totallength);
      }
      gt_assert(rightofpartwidth != NULL && (Seqpos) nextrightofpartwidth <
                (realspecialranges - countranges));
      rightofpartwidth[nextrightofpartwidth++]
        = compressedtable_get(rightposinverse,range.rightpos);
    }
  }
  freespecialrangeiterator(&sri);
  return rightofpartwidth;
}

static void inversesuffixarray2specialranknext(
                         const Compressedtable *rightposinverse,
                         Compressedtable *ranknext,
                         const Encodedsequence *encseq,
                         Readmode readmode,
                         Seqpos partwidth,
                         Seqpos totallength)
{
  if (hasspecialranges(encseq))
  {
    Specialrangeiterator *sri;
    Sequencerange range;
    Seqpos specialcharacters, idx, *rightofpartwidth = NULL;
    Seqpos specialranklistindex, nextrightofpartwidth = 0;

    rightofpartwidth = fillrightofpartwidth(rightposinverse,
                                            encseq,
                                            readmode,
                                            partwidth,
                                            totallength);
    specialcharacters = getencseqspecialcharacters(encseq);
    specialranklistindex = partwidth;
    sri = newspecialrangeiterator(encseq,
                                  ISDIRREVERSE(readmode) ? false : true);
    nextrightofpartwidth = 0;
    while (nextspecialrangeiterator(&range,sri))
    {
      gt_assert(range.rightpos <= totallength);
      for (idx = range.leftpos; idx < range.rightpos-1; idx++)
      {
        gt_assert(specialranklistindex < totallength);
        compressedtable_update(ranknext,specialranklistindex,
                               specialranklistindex + 1);
        specialranklistindex++;
      }
      gt_assert(specialranklistindex < totallength);
      if (range.rightpos < partwidth)
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
    gt_assert(specialranklistindex == totallength);
    freespecialrangeiterator(&sri);
  }
}

static Seqpos sa2ranknext(Compressedtable *ranknext,
                          const Encodedsequence *encseq,
                          Readmode readmode,
                          Seqpos partwidth,
                          Seqpos totallength,
                          const Seqpos *sortedsuffixes)
{
  Seqpos idx, longest = 0;
  unsigned long *occless;

  gt_assert(partwidth > 0);

  occless = computeocclesstab(encseq);
  /* now inveresuftab is not used any more, and thus the
     ranknext array (which points to ranknext can savely be stored */
  for (idx=0; idx < partwidth; idx++)
  {
    Seqpos pos = sortedsuffixes[idx];
    if (pos > 0)
    {
      GtUchar cc = getencodedchar(encseq,pos-1, readmode);
      if (ISNOTSPECIAL(cc))
      {
        gt_assert(occless[cc] < (unsigned long) partwidth);
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
  if (hasspecialranges(encseq))
  {
    Specialrangeiterator *sri;
    Sequencerange range;
    Seqpos specialidx, specialpos;

    sri = newspecialrangeiterator(encseq,
                                  ISDIRREVERSE(readmode) ? false : true);
    gt_assert(partwidth > 0); /* otherwise all lcps would be 0 */
    specialidx = partwidth;
    while (nextspecialrangeiterator(&range,sri))
    {
      for (specialpos = range.leftpos; specialpos < range.rightpos;
           specialpos++)
      {
        if (specialpos > 0)
        {
          GtUchar cc = getencodedchar(encseq,specialpos-1,readmode);
          if (ISNOTSPECIAL(cc))
          {
            gt_assert(occless[cc] < (unsigned long) partwidth);
            compressedtable_update(ranknext,(Seqpos) occless[cc],specialidx);
            occless[cc]++;
          }
        } else
        {
          longest = partwidth;
        }
        specialidx++;
      }
    }
    if (getencseqlengthofspecialsuffix(encseq))
    {
      compressedtable_update(ranknext,totallength-1,totallength);
    }
    freespecialrangeiterator(&sri);
  }
  gt_free(occless);
  return longest;
}

Compressedtable *lcp9_manzini(Compressedtable *spacefortab,
                              const Encodedsequence *encseq,
                              Readmode readmode,
                              Seqpos partwidth,
                              Seqpos totallength,
                              const Seqpos *sortedsuffixes)
{
  Seqpos pos, previousstart, nextfillpos = 0, fillpos, lcpvalue = 0;
  Compressedtable *lcptab, *ranknext, *rightposinverse;
  Seqpos previouscc1pos, previouscc2pos;
  GtUchar cc1, cc2;
  Encodedsequencescanstate *esr1, *esr2;

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
  printf("longest=" FormatSeqpos "\n",PRINTSeqposcast(fillpos));
  lcptab = ranknext;
  /* now ranknext and lcptab point to the same memory area. After reading
     ranknext at position fillpos, the same cell is used for storing
     the determined lcp-value */
  /* exploit the fact, that pos + lcpvalue is monotone */
  esr1 = newEncodedsequencescanstate();
  initEncodedsequencescanstate(esr1,encseq,readmode,0);
  cc1 = sequentialgetencodedchar(encseq,esr1,0,readmode);
  previouscc1pos = 0;
  esr2 = newEncodedsequencescanstate();
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
      previousstart = sortedsuffixes[fillpos-1];
      while (pos+lcpvalue < totallength &&
             previousstart+lcpvalue < totallength)
      {
        gt_assert(pos + lcpvalue >= previouscc1pos);
        while (previouscc1pos < pos + lcpvalue)
        {
          previouscc1pos++;
          cc1 = sequentialgetencodedchar(encseq,
                                         esr1,
                                         previouscc1pos,
                                         readmode);
        }
        if (ISSPECIAL(cc1))
        {
          break;
        }
        if (previousstart+lcpvalue < previouscc2pos ||
            previousstart+lcpvalue > previouscc2pos+1)
        {
          previouscc2pos = previousstart+lcpvalue;
          initEncodedsequencescanstate(esr2,encseq,
                                       readmode,
                                       previouscc2pos);
          cc2 = sequentialgetencodedchar(encseq,
                                         esr2,previouscc2pos,
                                         readmode);
        } else
        {
          if (previousstart+lcpvalue == previouscc2pos+1)
          {
            previouscc2pos++;
            cc2 = sequentialgetencodedchar(encseq,
                                           esr2,
                                           previouscc2pos,
                                           readmode);
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
