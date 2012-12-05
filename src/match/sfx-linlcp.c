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
#include "core/compact_ulong_store.h"
#include "esa-seqread.h"
#include "sarr-def.h"
#include "sfx-linlcp.h"

unsigned long *gt_lcp13_manzini(const GtEncseq *encseq,
                                GtReadmode readmode,
                                unsigned long partwidth,
                                unsigned long totallength,
                                const ESASuffixptr *suftab,
                                GtCompactUlongStore *inversesuftab)
{
  unsigned long pos, lcpvalue = 0, *lcptab;

  lcptab = gt_malloc(sizeof (unsigned long) * partwidth);
  lcptab[0] = 0;
  for (pos=0; pos <= totallength; pos++)
  {
    unsigned long fillpos = gt_compact_ulong_store_get(inversesuftab,pos);
    if (fillpos > 0 && fillpos < partwidth)
    {
      unsigned long previousstart = ESASUFFIXPTRGET(suftab,fillpos-1);
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

static unsigned long *computeocclesstab(const GtEncseq *encseq,
                                        GtReadmode readmode)
{
  unsigned long *occless;
  unsigned int charidx, numofchars;

  numofchars = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
  occless = gt_malloc(sizeof (unsigned long) * numofchars);
  occless[0] = 0;
  for (charidx = 1U; charidx < numofchars; charidx++)
  {
    unsigned long count;

    if (GT_ISDIRCOMPLEMENT(readmode))
    {
      count = gt_encseq_charcount(encseq,
                                  (GtUchar) GT_COMPLEMENTBASE(charidx-1));
    } else
    {
      count = gt_encseq_charcount(encseq,(GtUchar) (charidx-1));
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

static void setrelevantfrominversetab(GtCompactUlongStore *rightposinverse,
                                      const GtEncseq *encseq,
                                      GtReadmode readmode,
                                      const ESASuffixptr *suftab,
                                      unsigned long partwidth)
{
  if (gt_encseq_has_specialranges(encseq))
  {
    unsigned long idx, pos;

    for (idx = 0; idx < partwidth; idx++)
    {
      pos = ESASUFFIXPTRGET(suftab,idx);
      if (pos > 0)
      {
        GtUchar cc = gt_encseq_get_encoded_char(encseq,pos-1,readmode);
        if (ISSPECIAL(cc))
        {
          gt_compact_ulong_store_update(rightposinverse,pos,idx);
        }
      }
    }
  }
}

static unsigned long *fillrightofpartwidth(
                                     const GtCompactUlongStore *rightposinverse,
                                     const GtEncseq *encseq,
                                     GtReadmode readmode,
                                     unsigned long partwidth,
                                     unsigned long totallength)
{
  GtSpecialrangeiterator *sri;
  GtRange range;
  unsigned long countlargeranges, *rightofpartwidth = NULL,
                nextrightofpartwidth = 0;

  countlargeranges = gt_encseq_realspecialranges(encseq);
  sri = gt_specialrangeiterator_new(encseq,
                                    GT_ISDIRREVERSE(readmode) ? false : true);
  while (gt_specialrangeiterator_next(sri,&range))
  {
    if (GT_ISDIRREVERSE(readmode))
    {
      gt_range_reverse(totallength,&range);
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
        /*printf("allocated %lu bytes for rightofpartwidth (%.2f)\n",
            (unsigned long) allocsize, (double) allocsize/totallength);*/
      }
      gt_assert(nextrightofpartwidth < countlargeranges);
      rightofpartwidth[nextrightofpartwidth++]
        = gt_compact_ulong_store_get(rightposinverse,range.end);
    }
  }
  gt_specialrangeiterator_delete(sri);
  return rightofpartwidth;
}

static void inversesuffixarray2specialranknext(
                         const GtCompactUlongStore *rightposinverse,
                         GtCompactUlongStore *ranknext,
                         const GtEncseq *encseq,
                         GtReadmode readmode,
                         unsigned long partwidth,
                         unsigned long totallength)
{
  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;
    unsigned long idx, *rightofpartwidth = NULL,
                  specialranklistindex, nextrightofpartwidth = 0;

    rightofpartwidth = fillrightofpartwidth(rightposinverse,
                                            encseq,
                                            readmode,
                                            partwidth,
                                            totallength);
    specialranklistindex = partwidth;
    sri = gt_specialrangeiterator_new(encseq,
                                      GT_ISDIRREVERSE(readmode) ? false : true);
    nextrightofpartwidth = 0;
    while (gt_specialrangeiterator_next(sri,&range))
    {
      if (GT_ISDIRREVERSE(readmode))
      {
        gt_range_reverse(totallength,&range);
      }
      gt_assert(range.end <= totallength);
      for (idx = range.start; idx < range.end-1; idx++)
      {
        gt_assert(specialranklistindex < totallength);
        gt_compact_ulong_store_update(ranknext,specialranklistindex,
                                      specialranklistindex + 1);
        specialranklistindex++;
      }
      gt_assert(specialranklistindex < totallength);
      if (range.end < partwidth)
      {
        gt_compact_ulong_store_update(ranknext,specialranklistindex,
                                      gt_compact_ulong_store_get(
                                             rightposinverse,range.end));
      } else
      {
        gt_compact_ulong_store_update(ranknext,specialranklistindex,
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

static unsigned long sa2ranknext(GtCompactUlongStore *ranknext,
                                 const GtEncseq *encseq,
                                 GtReadmode readmode,
                                 unsigned long partwidth,
                                 unsigned long totallength,
                                 const ESASuffixptr *suftab)
{
  unsigned long idx, pos, longest = 0, *occless;

  gt_assert(partwidth > 0);

  occless = computeocclesstab(encseq,readmode);
  /* now inveresuftab is not used any more, and thus the
     ranknext array (which points to ranknext can savely be stored */
  for (idx=0; idx < partwidth; idx++)
  {
    pos = ESASUFFIXPTRGET(suftab,idx);
    if (pos > 0)
    {
      GtUchar cc = gt_encseq_get_encoded_char(encseq,pos-1, readmode);
      if (ISNOTSPECIAL(cc))
      {
        gt_assert(occless[cc] < partwidth);
        gt_compact_ulong_store_update(ranknext,occless[cc],idx);
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
    unsigned long specialidx;

    sri = gt_specialrangeiterator_new(encseq,
                                      GT_ISDIRREVERSE(readmode) ? false : true);
    gt_assert(partwidth > 0); /* otherwise all lcps would be 0 */
    specialidx = partwidth;
    while (gt_specialrangeiterator_next(sri,&range))
    {
      if (GT_ISDIRREVERSE(readmode))
      {
        gt_range_reverse(totallength,&range);
      }
      gt_assert(range.start < range.end);
      if (range.start > 0)
      {
        GtUchar cc = gt_encseq_get_encoded_char(encseq,range.start-1,readmode);
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
      gt_compact_ulong_store_update(ranknext,totallength-1,totallength);
    }
    gt_specialrangeiterator_delete(sri);
  }
  gt_free(occless);
  return longest;
}

GtCompactUlongStore *gt_lcp9_manzini(GtCompactUlongStore *spacefortab,
                                 const GtEncseq *encseq,
                                 GtReadmode readmode,
                                 unsigned long partwidth,
                                 unsigned long totallength,
                                 const ESASuffixptr *suftab)
{
  unsigned long pos, previousstart, nextfillpos = 0, fillpos, lcpvalue = 0,
                previouscc1pos, previouscc2pos;
  GtCompactUlongStore *lcptab, *ranknext, *rightposinverse;
  GtUchar cc1, cc2;
  GtEncseqReader *esr1, *esr2;

  if (spacefortab == NULL)
  {
    unsigned int bitsperentry = gt_determinebitspervalue(totallength);
    rightposinverse = ranknext
                    = gt_compact_ulong_store_new(totallength+1,bitsperentry);
    gt_compact_ulong_store_update(ranknext,totallength,totallength);
    setrelevantfrominversetab(rightposinverse,encseq,readmode,suftab,
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
                        suftab);
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
      nextfillpos = gt_compact_ulong_store_get(ranknext,fillpos);
    }
    if (fillpos > 0 && fillpos - 1 < partwidth)
    {
      previousstart = ESASUFFIXPTRGET(suftab,fillpos-1);
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
      gt_compact_ulong_store_update(lcptab,fillpos,lcpvalue);
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

static unsigned long gt_check_for_range_occurrence(const ESASuffixptr *suftab,
                                                   unsigned long suffix,
                                                   unsigned long start,
                                                   unsigned long end)
{
  unsigned long idx;

  for (idx = start; idx <= end; idx++)
  {
    unsigned long position = ESASUFFIXPTRGET(suftab,idx);

    if (suffix == position)
    {
      return idx;
    }
  }
  return ULONG_MAX;
}

typedef struct
{
  unsigned long start, end;
  GtUchar firstchar;
} GtRangewithchar;

/* The following funktion implements the linear time algorithm of
   @INPROCEEDINGS{BUR:KAER:2003,
   author = {Burkhardt, S. and K{\"a}rkk{\"a}inen, J.},
   title = {{Fast Lightweight Suffix Array Construction and Checking}},
   booktitle = {{Proceedings of the 14th Annual Symposium on Combinatorial
                 Pattern Matching (CPM)}},
   year = {2003},
   editor = {{Baeza-Yates, R. and Ch{\'a}vez, E. and Crochemore, M.}},
   volume = {2676},
   series = {LNCS},
   pages = {200-210},
   publisher = {Springer-Verlag}
   }
   to check the following suffix-order condition of the sorted suffix array:
   For all characters c, if SA[i,j] contains the suffixes starting
   with charcter c, then SA[i]+1, SA[i+1]+1, \ldots, SA[j]+1 occur
   in SA in this order (but not consecutively in general).
   The running time of the algorithm is independent of the alphabet size.
   The main problem is that it requires random access to the sequence
   which slows it down.
*/

static void gt_suftab_bk_suffixorder(const GtEncseq *encseq,
                                     GtReadmode readmode,
                                     unsigned long totallength,
                                     unsigned int numofchars,
                                     const ESASuffixptr *suftab,
                                     const GtRangewithchar *rangestore,
                                     unsigned int numofranges)
{
  unsigned int rangeidx;
  unsigned long idx, *nexttab = gt_calloc((size_t) numofchars,sizeof(*nexttab));

  for (rangeidx = 0; rangeidx < numofranges; rangeidx++)
  {
    nexttab[rangestore[rangeidx].firstchar] = rangestore[rangeidx].start;
  }
  for (idx = 0; idx < totallength; idx++)
  {
    unsigned long position = ESASUFFIXPTRGET(suftab,idx);

    if (position > 0)
    {
      GtUchar cc = gt_encseq_get_encoded_char(encseq,position - 1,readmode);
      if (ISNOTSPECIAL(cc))
      {
        unsigned long checkpos;

        checkpos = ESASUFFIXPTRGET(suftab,nexttab[(int) cc]) + 1;
        if (checkpos != position)
        {
          fprintf(stderr,"idx=%lu,checkpos=%lu,position=%lu\n",
                          idx,checkpos,position);
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
        nexttab[(int) cc]++;
      }
    }
  }
  gt_free(nexttab);
}

/* The following function checks the suffix-order condition described
   above using an O(\sigma n) algorithm where \sigma is the alphabet size
   and n is the length of the suffix array. The algorithm does not
   access the sequence and performs \sigma linear scans of the suffix array.
   This makes it faster than the previous method. It is probably possible
   to perform the check in one linear scan.
*/

static void gt_suftab_sk_suffixorder(unsigned long totallength,
                                     unsigned int numofchars,
                                     const ESASuffixptr *suftab,
                                     const GtRangewithchar *rangestore,
                                     unsigned int numofranges)
{
  unsigned int rangeidx;
  unsigned long numofcomparisons = 0;
  double ratio;

  for (rangeidx = 0; rangeidx < numofranges; rangeidx++)
  {
    unsigned long idx, start = 0;

    for (idx = rangestore[rangeidx].start; idx <= rangestore[rangeidx].end;
         idx++)
    {
      unsigned long position = ESASUFFIXPTRGET(suftab,idx);

      if (position + 1 <= totallength)
      {
        unsigned long found = gt_check_for_range_occurrence(suftab,
                                                            position + 1,
                                                            start,
                                                            totallength);
        if (found == ULONG_MAX)
        {
          fprintf(stderr,"Cannot find position+1=%lu in range [%lu,%lu]\n",
                            position+1,start,totallength);
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
        numofcomparisons += found - start + 1;
        start = found + 1;
      }
    }
  }
  ratio = (double) numofcomparisons/totallength;
  if (gt_double_compare(ratio,(double) numofchars) > 0)
  {
    fprintf(stderr,"gt_double_compare(%.2f,%u) = %d > 0 not exected\n",
                    ratio,numofchars,
                    gt_double_compare(ratio,(double) numofchars));
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

void gt_suftab_lightweightcheck(const GtEncseq *encseq,
                               GtReadmode readmode,
                               unsigned long totallength,
                               const ESASuffixptr *suftab,
                               GtLogger *logger)
{
  unsigned long idx, countbitsset = 0, previouspos = 0,
                firstspecial = totallength, rangestart = 0;
  unsigned int numofchars, charidx, rangeidx = 0, numofranges;
  GtBitsequence *startposoccurs;
  GtUchar previouscc = 0;
  GtRangewithchar *rangestore;
  const bool skcheck = true;

  GT_INITBITTAB(startposoccurs,totallength+1);
  numofchars = gt_encseq_alphabetnumofchars(encseq);
  rangestore = gt_malloc(sizeof(*rangestore) * numofchars);
  for (idx = 0; idx < totallength; idx++)
  {
    unsigned long position = ESASUFFIXPTRGET(suftab,idx);
    GtUchar cc;

    if (GT_ISIBITSET(startposoccurs,position))
    {
      fprintf(stderr,"ERROR: suffix with startpos %lu already occurs\n",
              ESASUFFIXPTRGET(suftab,idx));
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    GT_SETIBIT(startposoccurs,position);
    countbitsset++;
    cc = gt_encseq_get_encoded_char(encseq,position,readmode);
    if (idx > 0)
    {
      if (ISSPECIAL(cc))
      {
        if (firstspecial == totallength)
        {
          firstspecial = idx;
          gt_assert(rangeidx < numofchars);
          rangestore[rangeidx].start = rangestart;
          rangestore[rangeidx].end = idx-1;
          rangestore[rangeidx++].firstchar = previouscc;
        }
        if (ISSPECIAL(previouscc))
        {
          if (previouspos > position)
          {
            fprintf(stderr,"incorrect order: %lu = %lu=SPECIAL > SPECIAL=%lu "
                           " = %lu\n",
                      idx-1,position,previouspos,idx);
            exit(GT_EXIT_PROGRAMMING_ERROR);
          }
        }
      } else
      {
        if (ISSPECIAL(previouscc))
        {
          fprintf(stderr,"incorrect order: %lu=%lu=SPECIAL > %u=%lu=%lu\n",
                    idx-1,position,(unsigned int) cc,previouspos,idx);
          exit(GT_EXIT_PROGRAMMING_ERROR);
        } else
        {
          if (previouscc > cc)
          {
            fprintf(stderr,"incorrect order: %lu = %lu=%u > %u=%lu=%lu\n",
                      idx-1,position,(unsigned int) previouscc,
                      (unsigned int) cc,previouspos,idx);
            exit(GT_EXIT_PROGRAMMING_ERROR);
          } else
          {
            if (previouscc < cc)
            {
              gt_assert(rangeidx < numofchars);
              rangestore[rangeidx].start = rangestart;
              rangestore[rangeidx].end = idx-1;
              rangestore[rangeidx++].firstchar = previouscc;
              rangestart = idx;
            }
          }
        }
      }
    } else
    {
      if (ISSPECIAL(cc))
      {
        firstspecial = 0;
      }
    }
    previouscc = cc;
    previouspos = position;
  }
  if (countbitsset != totallength)
  {
    fprintf(stderr,"ERROR: only %lu of %lu suffixes occur\n",countbitsset,
                                totallength);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_free(startposoccurs);
  if (firstspecial == totallength)
  {
    gt_assert(firstspecial > 0 && rangeidx < numofchars);
    rangestore[rangeidx].start = rangestart;
    rangestore[rangeidx].end = firstspecial-1;
    rangestore[rangeidx++].firstchar = previouscc;
  }
  numofranges = rangeidx;
  for (charidx = 0, rangeidx = 0; charidx < numofchars; charidx++)
  {
    unsigned long count;

    if (GT_ISDIRCOMPLEMENT(readmode))
    {
      count = gt_encseq_charcount(encseq,(GtUchar) GT_COMPLEMENTBASE(charidx));
    } else
    {
      count = gt_encseq_charcount(encseq,(GtUchar) charidx);
    }
    if (count != 0)
    {
      gt_assert(rangestore[rangeidx].firstchar == (GtUchar) charidx);
      gt_assert(rangestore[rangeidx].end - rangestore[rangeidx].start + 1
                == count);
      rangeidx++;
    }
  }
  gt_logger_log(logger,"suftab-check, first phase done");
  if (skcheck)
  {
    gt_suftab_sk_suffixorder(totallength,
                             numofchars,
                             suftab,
                             rangestore,
                             numofranges);
    gt_logger_log(logger,"suftab-check, second phase (sk-method) done");
  } else
  {
    gt_suftab_bk_suffixorder(encseq,
                             readmode,
                             totallength,
                             numofchars,
                             suftab,
                             rangestore,
                             numofranges);
    gt_logger_log(logger,"suftab-check, second phase (bk-method) done");
  }
  gt_free(rangestore);
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
  unsigned long partwidth, totallength = gt_encseq_total_length(encseq),
                idx, specials = gt_encseq_specialcharacters(encseq);
  GtCompactUlongStore *lcptab = NULL;

  gt_assert(specials <= totallength);
  partwidth = totallength - specials;
  if (partwidth > 0)
  {
    lcptab = gt_lcp9_manzini(NULL,
                           encseq,
                           readmode,
                           partwidth,
                           totallength,
                           suftab);
    gt_logger_log(logger,"computed reference lcp table with manzini algorithm");
  }
  ssar = gt_newSequentialsuffixarrayreaderfromfile(esaindexname,
                                                   SARR_LCPTAB,
                                                   SEQ_scan,
                                                   logger,
                                                   err);
  if (ssar == NULL)
  {
    haserr = true;
  }
  for (idx = 1UL; /* Nothing */; idx++)
  {
    unsigned long mlcpvalue, lcpvalue;
    int retval = gt_nextSequentiallcpvalue(&lcpvalue,ssar,err);

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
      mlcpvalue = gt_compact_ulong_store_get(lcptab,idx);
    } else
    {
      mlcpvalue = 0;
    }
    if (mlcpvalue != lcpvalue)
    {
      fprintf(stderr,"%lu: mlcpvalue = %lu != %lu = lcpvalue\n",
                       idx,mlcpvalue,lcpvalue);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
  gt_logger_log(logger,"compare lcp-values against reference");
  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  gt_compact_ulong_store_delete(lcptab);
  return haserr ? -1 : 0;
}
