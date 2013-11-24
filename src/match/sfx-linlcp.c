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
#include "core/range.h"
#include "core/mathsupport.h"
#include "core/logger.h"
#include "core/minmax.h"
#include "core/compact_ulong_store.h"
#include "esa-seqread.h"
#include "sarr-def.h"
#include "sfx-linlcp.h"

GtUword *gt_ENCSEQ_lcp13_kasai(const GtEncseq *encseq,
                               GtReadmode readmode,
                               GtUword partwidth,
                               GtUword totallength,
                               const ESASuffixptr *suftab,
                               GtCompactUlongStore *inversesuftab)
{
  GtUword pos, lcpvalue = 0, *lcptab;

  lcptab = gt_malloc(sizeof (*lcptab) * partwidth);
  lcptab[0] = 0;
  for (pos=0; pos <= totallength; pos++)
  {
    GtUword fillpos = gt_compact_ulong_store_get(inversesuftab, pos);
    if (fillpos > 0 && fillpos < partwidth)
    {
      GtUword previousstart = ESASUFFIXPTRGET(suftab, fillpos-1);
      while (pos+lcpvalue < totallength &&
             previousstart+lcpvalue < totallength)
      {
        GtUchar cc1, cc2;

        cc1 = gt_encseq_get_encoded_char(encseq, pos+lcpvalue, readmode);
        cc2 = gt_encseq_get_encoded_char(encseq, previousstart+lcpvalue,
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

unsigned int *gt_plain_lcp13_kasai(GtUword *maxlcp,
                                   const GtUchar *sequence,
                                   bool withspecial,
                                   GtUword partwidth,
                                   GtUword totallength,
                                   const unsigned int *suftab)
{
  unsigned int *lcptab, *inversesuftab, idx;
  GtUword pos, lcpvalue = 0;

  inversesuftab = gt_malloc(sizeof (*inversesuftab) * (totallength+1));
  gt_assert(totallength <= (GtUword) UINT_MAX);
  for (idx = 0; idx <= (unsigned int) totallength; idx++)
  {
    inversesuftab[suftab[idx]] = idx;
  }
  lcptab = gt_malloc(sizeof (*lcptab) * (totallength+1));
  lcptab[0] = 0;
  *maxlcp = 0;
  for (pos = 0; pos <= totallength; pos++)
  {
    GtUword fillpos = (GtUword) inversesuftab[pos];
    if (fillpos > 0 && fillpos < partwidth)
    {
      GtUword previousstart = (GtUword) suftab[fillpos-1],
              lastoffset = totallength - MAX(pos,previousstart);

      while (lcpvalue < lastoffset)
      {
        GtUchar cc1, cc2;

        cc1 = sequence[pos+lcpvalue];
        cc2 = sequence[previousstart+lcpvalue];
        if (cc1 == cc2 && (!withspecial || ISNOTSPECIAL(cc1)))
        {
          lcpvalue++;
        } else
        {
          break;
        }
      }
      gt_assert(lcpvalue <= (GtUword) UINT_MAX);
      lcptab[fillpos] = (unsigned int) lcpvalue;
      if (*maxlcp < lcpvalue)
      {
        *maxlcp = lcpvalue;
      }
    }
    if (lcpvalue > 0)
    {
      lcpvalue--;
    }
  }
  gt_free(inversesuftab);
  return lcptab;
}

unsigned int *gt_plain_lcp_phialgorithm(bool onlyplcp,
                                        GtUword *maxlcp,
                                        const GtUchar *sequence,
                                        bool withspecial,
                                        GtUword partwidth,
                                        GtUword totallength,
                                        const unsigned int *suftab)
{
  unsigned int *plcptab, *phitab, pos, suftab0, previousvalue;
  GtUword idx, lcpvalue = 0;

  phitab = gt_malloc(sizeof (*phitab) * (totallength+1));
  previousvalue = suftab[0];
  for (idx = 1UL; idx <= totallength; idx++)
  {
    unsigned int currentvalue = suftab[idx];
    phitab[currentvalue] = previousvalue;
    previousvalue = currentvalue;
  }
  plcptab = phitab; /* overlay both arrays */
  suftab0 = suftab[0];
  gt_assert(totallength <= (GtUword) UINT_MAX);
  *maxlcp = 0;
  for (pos = 0; pos < (unsigned int) totallength; pos++)
  {
    if (pos != suftab0)
    {
      const unsigned int currentphitab = phitab[pos];
      const GtUword lastoffset = totallength - MAX(pos,currentphitab);
      const GtUchar *ptr1 = sequence + pos,
                    *ptr2 = sequence + currentphitab;

      while (lcpvalue < lastoffset)
      {
        GtUchar cc1 = ptr1[lcpvalue];
        GtUchar cc2 = ptr2[lcpvalue];
        if (cc1 == cc2 && (!withspecial || ISNOTSPECIAL(cc1)))
        {
          lcpvalue++;
        } else
        {
          break;
        }
      }
      gt_assert(lcpvalue <= (GtUword) UINT_MAX);
      plcptab[pos] = (unsigned int) lcpvalue;
      if (lcpvalue > 0)
      {
        if (*maxlcp < lcpvalue)
        {
          *maxlcp = lcpvalue;
        }
        lcpvalue--;
      }
    } else
    {
      plcptab[pos] = 0;
    }
  }
  if (onlyplcp)
  {
    return plcptab;
  } else
  {
    unsigned int *lcptab = gt_malloc(sizeof (*lcptab) * (totallength+1));

    for (idx = 0; idx < partwidth; idx++)
    {
      lcptab[idx] = plcptab[suftab[idx]];
    }
    gt_free(plcptab);
    for (idx = partwidth; idx <= totallength; idx++)
    {
      lcptab[idx] = 0;
    }
    return lcptab;
  }
}

static GtUword *gt_ENCSEQ_compute_occless_tab(const GtEncseq *encseq,
                                              GtReadmode readmode)
{
  GtUword *occless;
  unsigned int charidx, numofchars;

  numofchars = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
  occless = gt_malloc(sizeof (GtUword) * numofchars);
  occless[0] = 0;
  for (charidx = 1U; charidx < numofchars; charidx++)
  {
    GtUword count;

    if (GT_ISDIRCOMPLEMENT(readmode))
    {
      count = gt_encseq_charcount(encseq,
                                  (GtUchar) GT_COMPLEMENTBASE(charidx-1));
    } else
    {
      count = gt_encseq_charcount(encseq, (GtUchar) (charidx-1));
    }
    occless[charidx] = occless[charidx-1] + count;
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

static void gt_ENCSEQ_set_relevant_from_inversetab(
                                      GtCompactUlongStore *rightposinverse,
                                      const GtEncseq *encseq,
                                      GtReadmode readmode,
                                      const ESASuffixptr *suftab,
                                      GtUword partwidth)
{
  if (gt_encseq_has_specialranges(encseq))
  {
    GtUword idx, pos;

    for (idx = 0; idx < partwidth; idx++)
    {
      pos = ESASUFFIXPTRGET(suftab, idx);
      if (pos > 0)
      {
        GtUchar cc = gt_encseq_get_encoded_char(encseq, pos-1, readmode);
        if (ISSPECIAL(cc))
        {
          gt_compact_ulong_store_update(rightposinverse, pos, idx);
        }
      }
    }
  }
}

static GtUword *gt_ENCSEQ_fill_right_of_partwidth(
                               const GtCompactUlongStore *rightposinverse,
                               const GtEncseq *encseq,
                               GtReadmode readmode,
                               GtUword partwidth,
                               GtUword totallength)
{
  GtSpecialrangeiterator *sri;
  GtRange range;
  GtUword countlargeranges, *rightofpartwidth = NULL,
                nextrightofpartwidth = 0;

  countlargeranges = gt_encseq_realspecialranges(encseq);
  sri = gt_specialrangeiterator_new(encseq,
                                    GT_ISDIRREVERSE(readmode) ? false : true);
  while (gt_specialrangeiterator_next(sri, &range))
  {
    if (GT_ISDIRREVERSE(readmode))
    {
      gt_range_reverse(totallength, &range);
    }
    if (range.end < partwidth)
    {
      gt_assert(countlargeranges > 0);
      countlargeranges--;
    } else
    {
      if (rightofpartwidth == NULL)
      {
        size_t allocsize = sizeof (*rightofpartwidth) * countlargeranges;
        rightofpartwidth = gt_malloc(allocsize);
      }
      gt_assert(nextrightofpartwidth < countlargeranges);
      rightofpartwidth[nextrightofpartwidth++]
        = gt_compact_ulong_store_get(rightposinverse, range.end);
    }
  }
  gt_specialrangeiterator_delete(sri);
  return rightofpartwidth;
}

static void gt_ENCSEQ_inversesuffixarray2specialranknext(
                         const GtCompactUlongStore *rightposinverse,
                         GtCompactUlongStore *ranknext,
                         const GtEncseq *encseq,
                         GtReadmode readmode,
                         GtUword partwidth,
                         GtUword totallength)
{
  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;
    GtUword idx, *rightofpartwidth = NULL,
                  specialranklistindex, nextrightofpartwidth = 0;

    rightofpartwidth = gt_ENCSEQ_fill_right_of_partwidth(rightposinverse,
                                                         encseq,
                                                         readmode,
                                                         partwidth,
                                                         totallength);
    specialranklistindex = partwidth;
    sri = gt_specialrangeiterator_new(encseq,
                                      GT_ISDIRREVERSE(readmode) ? false : true);
    nextrightofpartwidth = 0;
    while (gt_specialrangeiterator_next(sri, &range))
    {
      if (GT_ISDIRREVERSE(readmode))
      {
        gt_range_reverse(totallength, &range);
      }
      gt_assert(range.end <= totallength);
      for (idx = range.start; idx < range.end-1; idx++)
      {
        gt_assert(specialranklistindex < totallength);
        gt_compact_ulong_store_update(ranknext, specialranklistindex,
                                      specialranklistindex + 1);
        specialranklistindex++;
      }
      gt_assert(specialranklistindex < totallength);
      if (range.end < partwidth)
      {
        gt_compact_ulong_store_update(ranknext, specialranklistindex,
                                      gt_compact_ulong_store_get(
                                             rightposinverse, range.end));
      } else
      {
        gt_compact_ulong_store_update(ranknext, specialranklistindex,
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

static GtUword gt_ENCSEQ_sa2ranknext(GtCompactUlongStore *ranknext,
                                     const GtEncseq *encseq,
                                     GtReadmode readmode,
                                     GtUword partwidth,
                                     GtUword totallength,
                                     const ESASuffixptr *suftab)
{
  GtUword idx, pos, longest = 0, *occless;

  gt_assert(partwidth > 0);
  occless = gt_ENCSEQ_compute_occless_tab(encseq, readmode);
  /* now inveresuftab is not used any more, and thus the
     ranknext array (which points to ranknext can savely be stored */
  for (idx=0; idx < partwidth; idx++)
  {
    pos = ESASUFFIXPTRGET(suftab, idx);
    if (pos > 0)
    {
      GtUchar cc = gt_encseq_get_encoded_char(encseq, pos-1, readmode);
      if (ISNOTSPECIAL(cc))
      {
        gt_assert(occless[cc] < partwidth);
        gt_compact_ulong_store_update(ranknext, occless[cc], idx);
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
    GtUword specialidx;

    sri = gt_specialrangeiterator_new(encseq,
                                      GT_ISDIRREVERSE(readmode) ? false : true);
    gt_assert(partwidth > 0); /* otherwise all lcps would be 0 */
    specialidx = partwidth;
    while (gt_specialrangeiterator_next(sri, &range))
    {
      if (GT_ISDIRREVERSE(readmode))
      {
        gt_range_reverse(totallength, &range);
      }
      gt_assert(range.start < range.end);
      if (range.start > 0)
      {
        GtUchar cc = gt_encseq_get_encoded_char(encseq, range.start-1,
                                                readmode);
        if (ISNOTSPECIAL(cc))
        {
          gt_assert(occless[cc] < partwidth);
          gt_compact_ulong_store_update(ranknext, occless[cc], specialidx);
          occless[cc]++;
        }
      } else
      {
        longest = partwidth;
      }
      specialidx += range.end - range.start;
    }
    if ((GT_ISDIRREVERSE(readmode) &&
         gt_encseq_lengthofspecialprefix(encseq) > 0) ||
        (!GT_ISDIRREVERSE(readmode) &&
         gt_encseq_lengthofspecialsuffix(encseq) > 0))
    {
      gt_compact_ulong_store_update(ranknext, totallength-1, totallength);
    }
    gt_specialrangeiterator_delete(sri);
  }
  gt_free(occless);
  return longest;
}

GtCompactUlongStore *gt_ENCSEQ_lcp9_manzini(GtCompactUlongStore *spacefortab,
                                            const GtEncseq *encseq,
                                            GtReadmode readmode,
                                            GtUword partwidth,
                                            GtUword totallength,
                                            const ESASuffixptr *suftab)
{
  GtUword pos, previousstart, nextfillpos = 0, fillpos, lcpvalue = 0,
                previouscc1pos, previouscc2pos;
  GtCompactUlongStore *lcptab, *ranknext, *rightposinverse;
  GtUchar cc1, cc2;
  GtEncseqReader *esr1, *esr2;

  if (spacefortab == NULL)
  {
    unsigned int bitsperentry = gt_determinebitspervalue(totallength);
    rightposinverse = ranknext
                    = gt_compact_ulong_store_new(totallength+1, bitsperentry);
    gt_compact_ulong_store_update(ranknext, totallength, totallength);
    gt_ENCSEQ_set_relevant_from_inversetab(rightposinverse, encseq, readmode,
                                           suftab, partwidth);
  } else
  {
    rightposinverse = ranknext = spacefortab;
  }
  gt_ENCSEQ_inversesuffixarray2specialranknext(rightposinverse, ranknext,
                                               encseq,
                                               readmode,
                                               partwidth,
                                               totallength);
  fillpos = gt_ENCSEQ_sa2ranknext(ranknext, encseq, readmode, partwidth,
                                  totallength, suftab);
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
      nextfillpos = gt_compact_ulong_store_get(ranknext, fillpos);
    }
    if (fillpos > 0 && fillpos - 1 < partwidth)
    {
      previousstart = ESASUFFIXPTRGET(suftab, fillpos-1);
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
      gt_compact_ulong_store_update(lcptab, fillpos, lcpvalue);
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

int gt_lcptab_lightweightcheck(const char *esaindexname,
                               const GtEncseq *encseq,
                               GtReadmode readmode,
                               const ESASuffixptr *suftab,
                               GtLogger *logger,
                               GtError *err)
{
  bool haserr = false;
  Sequentialsuffixarrayreader *ssar;
  GtUword partwidth, totallength = gt_encseq_total_length(encseq),
                idx, specials = gt_encseq_specialcharacters(encseq);
  GtCompactUlongStore *lcptab = NULL;

  gt_assert(specials <= totallength);
  partwidth = totallength - specials;
  if (partwidth > 0)
  {
    lcptab = gt_ENCSEQ_lcp9_manzini(NULL,
                                    encseq,
                                    readmode,
                                    partwidth,
                                    totallength,
                                    suftab);
    gt_logger_log(logger,
                  "computed reference lcp table with manzini algorithm");
  }
  ssar = gt_newSequentialsuffixarrayreaderfromfile(esaindexname,
                                                   SARR_LCPTAB,
                                                   true,
                                                   logger,
                                                   err);
  if (ssar == NULL)
  {
    haserr = true;
  }
  for (idx = 1UL; /* Nothing */; idx++)
  {
    GtUword mlcpvalue, lcpvalue;
    int retval = gt_nextSequentiallcpvalue(&lcpvalue, ssar, err);

    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      break;
    }
    if (idx < partwidth)
    {
      mlcpvalue = gt_compact_ulong_store_get(lcptab, idx);
    } else
    {
      mlcpvalue = 0;
    }
    if (mlcpvalue != lcpvalue)
    {
      fprintf(stderr, ""GT_WU": mlcpvalue = "GT_WU" != "GT_WU" = lcpvalue\n",
                       idx, mlcpvalue, lcpvalue);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
  gt_logger_log(logger, "compare lcp-values against reference");
  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  gt_compact_ulong_store_delete(lcptab);
  return haserr ? -1 : 0;
}
