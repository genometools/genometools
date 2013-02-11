/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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
#include <stdio.h>
#include "core/assert_api.h"
#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/minmax.h"
#include "core/arraydef.h"
#include "core/unused_api.h"
#include "core/types_api.h"
#include "core/encseq.h"
#include "bcktab.h"
#include "kmer2string.h"
#include "sfx-radixsort.h"
#include "sfx-bltrie.h"
#include "sfx-copysort.h"
#include "sfx-bentsedg.h"
#include "sfx-lcpvalues.h"
#include "sfx-suftaborder.h"
#include "sfx-suffixgetset.h"
#include "sfx-shortreadsort.h"

#define ACCESSCHARRAND(POS)    gt_encseq_get_encoded_char(bsr->encseq,\
                                                          POS,bsr->readmode)
#define ACCESSCHARSEQ(ESR)     gt_encseq_reader_next_encoded_char(ESR)
#define ISNOTEND(POS)          ((POS) < bsr->totallength &&\
                                ISNOTSPECIAL(ACCESSCHARRAND(POS)))

#define DEREFSTOPPOSSEQ(VAR,POS,STOPPOS,ESR)\
        (((POS) < (STOPPOS) && ISNOTSPECIAL(VAR = ACCESSCHARSEQ(ESR))) ?\
        ((unsigned long) VAR) : GT_UNIQUEINT(POS))

#define DEREFSEQ(VAR,POS,ESR) DEREFSTOPPOSSEQ(VAR,POS,bsr->totallength,ESR)

#define BS_SWAPARRAY(TMP,SUBBUCKETLEFT,IDX1,IDX2)\
        if ((IDX1) != (IDX2))\
        {\
          TMP = gt_suffixsortspace_get(bsr->sssp,SUBBUCKETLEFT,IDX1);\
          gt_suffixsortspace_set(bsr->sssp,SUBBUCKETLEFT,IDX1,\
                                 gt_suffixsortspace_get(bsr->sssp,\
                                                        SUBBUCKETLEFT,\
                                                        IDX2));\
          gt_suffixsortspace_set(bsr->sssp,SUBBUCKETLEFT,IDX2,TMP);\
        }

#define STACKTOP\
        bsr->mkvauxstack.spaceGtMKVstack[bsr->mkvauxstack.nextfreeGtMKVstack]

#define GT_BSR_UPDATELCP(MINVAL,MAXVAL,LCP)\
        if ((MINVAL) > (LCP))\
        {\
          MINVAL = LCP;\
        }\
        GT_UPDATE_MAX(MAXVAL,LCP)

#define CMPCHARBYCHARPTR2INT(VAR,SUBBUCKETLEFT,TMPVAR,IDX)\
        VAR = (((cptr = gt_suffixsortspace_get(bsr->sssp,SUBBUCKETLEFT,IDX)+\
                        depth)\
                < bsr->totallength &&\
                ISNOTSPECIAL(TMPVAR = ACCESSCHARRAND(cptr)))\
                    ? ((unsigned long) TMPVAR) : GT_UNIQUEINT(cptr))

typedef GtEndofTwobitencoding GtSfxcmp;

#define PTR2INTSTOREPOS(TMPVAR,SUBBUCKETLEFT,IDX,POSASSIGNMENT)\
        {\
          unsigned long pos\
            = gt_suffixsortspace_get(bsr->sssp,SUBBUCKETLEFT,IDX);\
          POSASSIGNMENT;\
          TMPVAR.referstartpos = pos;\
          if (pos + depth < bsr->totallength)\
          {\
            pos += depth;\
            (void) gt_encseq_extract2bitencwithtwobitencodingstoppos(&(TMPVAR),\
                                                                 bsr->esr1,\
                                                                 bsr->encseq,\
                                                                 bsr->readmode,\
                                                                 pos);\
          } else\
          {\
            TMPVAR.tbe = 0;\
            TMPVAR.unitsnotspecial = 0;\
          }\
        }

#define PTR2INT(TMPVAR,SUBBUCKETLEFT,IDX)\
        PTR2INTSTOREPOS(TMPVAR,SUBBUCKETLEFT,IDX,/* Nothing */)

#define Sfxdocompare(COMMONUNITS,X,Y)\
        ret##X##Y = gt_encseq_compare_pairof_twobitencodings(bsr->fwd,\
                                                             bsr->complement,\
                                                             COMMONUNITS,&X,&Y)

#define GtSfxcmpEQUAL(X,Y)      (ret##X##Y == 0)
#define GtSfxcmpSMALLER(X,Y)    (ret##X##Y < 0)
#define GtSfxcmpGREATER(X,Y)    (ret##X##Y > 0)

typedef struct
{
  unsigned long subbucketleft,
                width,
                depth;
} GtMKVstack;

typedef struct
{
  GtEndofTwobitencoding etbe;
  unsigned long suftaboffset;
} GtMedianinfo;

typedef GtMedianinfo GtMedianElem;

typedef struct
{
  unsigned long suffix;
  unsigned char lcpwithpivot;
  char cmpresult;
} GtCountingsortinfo;

GT_DECLAREARRAYSTRUCT(GtMKVstack);

typedef struct
{
  const GtEncseq *encseq;
  GtEncseqReader *esr1, /* XXX be careful with threads */
                 *esr2;
  GtReadmode readmode;
  bool fwd, complement;
  unsigned long totallength;
  GtArrayGtMKVstack mkvauxstack; /* XXX be careful with threads */
  GtLcpvalues *tableoflcpvalues;
  GtMedianinfo *medianinfospace;
  GtCountingsortinfo *countingsortinfo;
  const Sfxstrategy *sfxstrategy;
  unsigned int sortmaxdepth,
               prefixlength;
  GtBlindtrie *blindtrie;
  unsigned long leftlcpdist[GT_UNITSIN2BITENC],
                rightlcpdist[GT_UNITSIN2BITENC];
  GtSuffixsortspace *sssp;
  GtProcessunsortedsuffixrange processunsortedsuffixrange;
  void *processunsortedsuffixrangeinfo;
  bool *equalwithprevious;
  size_t sizeofworkspace;
  unsigned long countinsertionsort,
                counttqsort,
                countshortreadsort,
                countradixsort,
                countcountingsort,
                countbltriesort,
                srs_maxremain; /* only relevant for short read sort */
  unsigned long radixsortminwidth,
                radixsortmaxwidth,
                shortreadsort_maxwidth;
  GtShortreadsortworkinfo *srsw;
  const GtTwobitencoding *twobitencoding;
  GtRadixsortstringinfo *rsi;
} GtBentsedgresources;

#ifdef WITHCHECKSTARTPOINTER
static unsigned int checkstartpointorder(const unsigned long *left,
                                         const unsigned long *right)
{
  const unsigned long *ptr;
  bool ascending;

  gt_assert(left < right);
  gt_assert(*left != *(left+1));
  ascending = (*left < *(left+1)) ? true : false;
  for (ptr = left+1; ptr < right; ptr++)
  {
    gt_assert(*ptr != *(ptr+1));
    if (*ptr < *(ptr+1))
    {
      if (!ascending)
      {
        return 0;
      }
    } else
    {
      if (*ptr > *(ptr+1))
      {
        if (ascending)
        {
          return 0;
        }
      }
    }
  }
  return ascending ? 1U : 2U;
}
#endif

static unsigned long medianof3cmpcharbychar(const GtBentsedgresources *bsr,
                                            unsigned long subbucketleft,
                                            unsigned long depth,
                                            unsigned long a,
                                            unsigned long b,
                                            unsigned long c)
{
  unsigned long vala, valb, valc, cptr;
  GtUchar tmpavar, tmpbvar;

  CMPCHARBYCHARPTR2INT(vala,subbucketleft,tmpavar,a);
  CMPCHARBYCHARPTR2INT(valb,subbucketleft,tmpbvar,b);
  if (vala == valb)
  {
    return a;
  }
  CMPCHARBYCHARPTR2INT(valc,subbucketleft,tmpavar,c);
  if (vala == valc || valb == valc)
  {
    return c;
  }
  return vala < valb ?
        (valb < valc ? b : (vala < valc ? c : a))
      : (valb > valc ? b : (vala < valc ? a : c));
}

static unsigned long medianof3(const GtBentsedgresources *bsr,
                               unsigned long subbucketleft,
                               unsigned long depth,
                               unsigned long a,
                               unsigned long b,
                               unsigned long c)
{
  GtSfxcmp vala, valb, valc;
  GtCommonunits commonunits;
  int retvalavalb, retvalavalc, retvalbvalc;

  PTR2INT(vala,subbucketleft,a);
  PTR2INT(valb,subbucketleft,b);
  Sfxdocompare(&commonunits,vala,valb);
  if (GtSfxcmpEQUAL(vala,valb))
  {
    return a;
  }
  PTR2INT(valc,subbucketleft,c);
  Sfxdocompare(&commonunits,vala,valc);
  if (GtSfxcmpEQUAL(vala,valc))
  {
    return c;
  }
  Sfxdocompare(&commonunits,valb,valc);
  if (GtSfxcmpEQUAL(valb,valc))
  {
    return c;
  }
  return GtSfxcmpSMALLER(vala,valb) ?
        (GtSfxcmpSMALLER(valb,valc) ? b : (GtSfxcmpSMALLER(vala,valc) ? c : a))
      : (GtSfxcmpGREATER(valb,valc) ? b : (GtSfxcmpSMALLER(vala,valc) ? a : c));
}

static void bs_insertionsort(GtBentsedgresources *bsr,
                             unsigned long subbucketleft,
                             unsigned long width,
                             unsigned long offset)
{
  unsigned long sval1, sval2, pm, pl, startpos1, startpos2, lcplen = 0;
  int retval;
  GtCommonunits commonunits;

  bsr->countinsertionsort++;
  for (pm = 1UL; pm < width; pm++)
  {
    for (pl = pm; pl > 0; pl--)
    {
      sval1 = gt_suffixsortspace_get(bsr->sssp,subbucketleft,pl-1);
      sval2 = gt_suffixsortspace_get(bsr->sssp,subbucketleft,pl);
      if (bsr->sfxstrategy->cmpcharbychar)
      {
        startpos1 = sval1 + offset;
        if (startpos1 < bsr->totallength)
        {
          gt_encseq_reader_reinit_with_readmode(bsr->esr1, bsr->encseq,
                                                bsr->readmode, startpos1);
        }
        startpos2 = sval2 + offset;
        if (startpos2 < bsr->totallength)
        {
          gt_encseq_reader_reinit_with_readmode(bsr->esr2, bsr->encseq,
                                                bsr->readmode, startpos2);
        }
        for (;;)
        {
          unsigned long ccs, cct;
          GtUchar tmp1, tmp2;

          ccs = DEREFSEQ(tmp1,startpos1,bsr->esr1);
          cct = DEREFSEQ(tmp2,startpos2,bsr->esr2);
          if (ccs != cct)
          {
            lcplen = startpos2 - sval2;
            retval = (ccs < cct) ? -1 : 1;
            break;
          }
          startpos1++;
          startpos2++;
        }
      } else
      {
#ifdef SKDEBUG
        printf("%s[%lu,%lu] at offset %lu\n",__func__,sval1,sval2,offset);
        gt_encseq_showatstartpos(stdout,
                                 bsr->fwd,
                                 bsr->complement,
                                 bsr->encseq,
                                 sval1);
        gt_encseq_showatstartpos(stdout,
                                 bsr->fwd,
                                 bsr->complement,
                                 bsr->encseq,
                                 sval2);
#endif
        retval = gt_encseq_compare_viatwobitencoding(&commonunits,
                                                     bsr->encseq,
                                                     bsr->encseq,
                                                     bsr->readmode,
                                                     bsr->esr1,
                                                     bsr->esr2,
                                                     sval1,
                                                     sval2,
                                                     offset,
                                                     0);
        lcplen = commonunits.finaldepth;
      }
      gt_assert(retval != 0);
      if (bsr->tableoflcpvalues != NULL)
      {
        if (pl < pm && retval > 0)
        {
          gt_lcptab_update(bsr->tableoflcpvalues,subbucketleft,pl+1,
                           gt_lcptab_getvalue(bsr->tableoflcpvalues,
                                              subbucketleft,pl));
        }
        gt_lcptab_update(bsr->tableoflcpvalues,subbucketleft,pl,lcplen);
      }
      if (retval < 0)
      {
        break;
      }
      gt_suffixsortspace_set(bsr->sssp,subbucketleft,pl-1,sval2);
      gt_suffixsortspace_set(bsr->sssp,subbucketleft,pl,sval1);
    }
  }
}

static void bs_insertionsortmaxdepth(GtBentsedgresources *bsr,
                                     unsigned long subbucketleft,
                                     unsigned long width,
                                     unsigned long offset,
                                     unsigned long sortmaxdepth)
{
  unsigned long sval1, sval2, pm, pl, startpos1, startpos2,
                lcplen = 0, idx = 0;
  int retval;
  bool tempb;
  GtCommonunits commonunits;

  bsr->countinsertionsort++;
  for (pm = 1UL; pm < width; pm++)
  {
    for (pl = pm; pl > 0; pl--)
    {
      sval1 = gt_suffixsortspace_get(bsr->sssp,subbucketleft,pl-1);
      sval2 = gt_suffixsortspace_get(bsr->sssp,subbucketleft,pl);
      if (bsr->sfxstrategy->cmpcharbychar)
      {
        unsigned long endpos1, endpos2;

        endpos1 = sval1+sortmaxdepth;
        if (endpos1 > bsr->totallength)
        {
          endpos1 = bsr->totallength;
        }
        endpos2 = sval2+sortmaxdepth;
        if (endpos2 > bsr->totallength)
        {
          endpos2 = bsr->totallength;
        }
        startpos1 = sval1+offset;
        if (startpos1 < bsr->totallength)
        {
          gt_encseq_reader_reinit_with_readmode(bsr->esr1, bsr->encseq,
                                                bsr->readmode, startpos1);
        }
        startpos2 = sval2+offset;
        if (startpos2 < bsr->totallength)
        {
          gt_encseq_reader_reinit_with_readmode(bsr->esr2, bsr->encseq,
                                                bsr->readmode, startpos2);
        }
        for (;;)
        {
          unsigned long ccs, cct;
          GtUchar tmp1, tmp2;

          ccs = DEREFSTOPPOSSEQ(tmp1,startpos1,endpos1,bsr->esr1);
          cct = DEREFSTOPPOSSEQ(tmp2,startpos2,endpos2,bsr->esr2);
          lcplen = startpos2 - sval2;
          if (lcplen == sortmaxdepth)
          {
            retval = 0;
            break;
          }
          gt_assert(lcplen < sortmaxdepth);
          if (ccs != cct)
          {
            retval = (ccs < cct) ? -1 : 1;
            break;
          }
          startpos1++;
          startpos2++;
        }
      } else
      {
        gt_assert(offset < sortmaxdepth);
        retval = gt_encseq_compare_viatwobitencoding(&commonunits,
                                                     bsr->encseq,
                                                     bsr->encseq,
                                                     bsr->readmode,
                                                     bsr->esr1,
                                                     bsr->esr2,
                                                     sval1,
                                                     sval2,
                                                     offset,
                                                     sortmaxdepth);
        lcplen = commonunits.finaldepth;
        gt_assert(lcplen <= sortmaxdepth);
        if (lcplen == sortmaxdepth)
        {
          gt_assert(retval == 0);
        }
      }
#ifdef SKDEBUG
      printf("cmp %lu and %lu: retval = %d, lcplen = %lu\n",
             sval1, sval2, retval, (unsigned long) lcplen);
#endif
      if (bsr->tableoflcpvalues != NULL)
      {
        if (pl < pm && retval > 0)
        {
          gt_lcptab_update(bsr->tableoflcpvalues,subbucketleft,pl+1,
                           gt_lcptab_getvalue(bsr->tableoflcpvalues,
                                              subbucketleft,pl));
        }
        gt_lcptab_update(bsr->tableoflcpvalues,subbucketleft,pl,lcplen);
      }
      if (retval < 0)
      {
        break;
      }
      idx = pl;
      if (retval == 0)
      {
        gt_assert(idx > 0);
        bsr->equalwithprevious[idx] = true;
        break;
      }
      gt_suffixsortspace_set(bsr->sssp,subbucketleft,pl-1,sval2);
      gt_suffixsortspace_set(bsr->sssp,subbucketleft,pl,sval1);
      tempb = bsr->equalwithprevious[idx-1];
      bsr->equalwithprevious[idx-1] = bsr->equalwithprevious[idx];
      bsr->equalwithprevious[idx] = tempb;
    }
  }
  if (idx > 0)
  {
    unsigned long equalsrangewidth = 0;
    unsigned long bucketleftidx
     = gt_suffixsortspace_bucketleftidx_get(bsr->sssp);
#ifdef SKDEBUG
    printf("ordered suffix %lu\n",gt_suffixsortspace_get(bsr->sssp,
                                                         subbucketleft,0));
#endif
    for (idx = 1UL; idx < width; idx++)
    {
#ifdef SKDEBUG
      printf("ordered suffix %lu, equalwithprevious=%s\n",
              gt_suffixsortspace_get(bsr->sssp,subbucketleft,idx),
              bsr->equalwithprevious[idx] ? "true" : "false");
#endif
      if (bsr->equalwithprevious[idx])
      {
        bsr->equalwithprevious[idx] = false;
        equalsrangewidth++;
      } else
      {
        if (equalsrangewidth > 0)
        {
#ifdef SKDEBUG
          printf("process interval of width %lu\n",
                 equalsrangewidth + 1);
#endif
          if (bsr->processunsortedsuffixrange != NULL)
          {
            bsr->processunsortedsuffixrange(bsr->processunsortedsuffixrangeinfo,
                                            bucketleftidx + subbucketleft
                                                          + idx - 1
                                                          - equalsrangewidth,
                                            equalsrangewidth + 1, sortmaxdepth);
          }
          equalsrangewidth = 0;
        }
      }
    }
    if (equalsrangewidth > 0)
    {
#ifdef SKDEBUG
      printf("process interval of width %lu\n",
             equalsrangewidth + 1);
#endif
      if (bsr->processunsortedsuffixrange != NULL)
      {
        bsr->processunsortedsuffixrange(bsr->processunsortedsuffixrangeinfo,
                                        bucketleftidx + subbucketleft + width
                                                      - 1 - equalsrangewidth,
                                        equalsrangewidth + 1, sortmaxdepth);
      }
    }
  }
}

#define DOMEDIANCOMPARE(A,B)\
        gt_encseq_compare_pairof_twobitencodings(fwd,complement,&commonunits,\
                                                 &((A)->etbe),&((B)->etbe))

#define GT_MedianElemGREATER(A,B)  (DOMEDIANCOMPARE(A,B) > 0)

#define GT_MedianElemSWAP(A,B)     {\
                                     register GtMedianElem tmp = *(A);\
                                                           *(A) = *(B);\
                                                           *(B) = tmp;\
                                   }

/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */

static GtMedianElem *quickmedian (bool fwd,bool complement,
                                  GtMedianElem *arr,unsigned long width)
{
  GtMedianElem *low, *high, *median, *middle, *ll, *hh;
  GtCommonunits commonunits;

  gt_assert(width > 0);
  low = arr;
  high = arr + width - 1;
  median = low + GT_DIV2(width);
  for (;;)
  {
    if (high <= low)                   /* One element only */
    {
      return median;
    }
    if (high == low + 1)
    {                                  /* Two elements only */
      if (GT_MedianElemGREATER(low,high))
      {
        GT_MedianElemSWAP (low, high);
      }
      return median;
    }

    /* Find median of low, middle and high items; swap into position low */
    middle = low + GT_DIV2(high - low + 1);
    if (GT_MedianElemGREATER(middle,high))
    {
      GT_MedianElemSWAP (middle, high);
    }
    if (GT_MedianElemGREATER(low,high))
    {
      GT_MedianElemSWAP (low, high);
    }
    if (GT_MedianElemGREATER(middle,low))
    {
      GT_MedianElemSWAP (middle, low);
    }
    /* Swap low item (now in position middle) into position (low+1) */
    GT_MedianElemSWAP (middle, low + 1);

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;)
    {
      do
      {
        ll++;
      } while (GT_MedianElemGREATER(low,ll));
      do
      {
        hh--;
      } while  (GT_MedianElemGREATER(hh,low));
      if (hh < ll)
      {
        break;
      }
      GT_MedianElemSWAP (ll, hh);
    }

    /* Swap middle item (in position low) back into correct position */
    GT_MedianElemSWAP (low, hh);

    /* Re-set active partition */
    if (hh <= median)
    {
      low = ll;
    }
    if (hh >= median)
    {
      high = hh - 1;
    }
  }
}

#ifdef WITHcheckmedian

static void checkmedian(bool fwd,
                        bool complement,
                        const GtMedianinfo *median,
                        const GtMedianinfo *space,
                        unsigned long width)
{
  unsigned long sum1, sum2, idx, smaller = 0, larger = 0, equal = 0, equalpart;
  unsigned int commonunits;
  int cmp;

  for (idx = 0; idx < width; idx++)
  {
    cmp = DOMEDIANCOMPARE(space + idx,median);
    if (cmp > 0)
    {
     larger++;
    } else
    {
      if (cmp < 0)
      {
        smaller++;
      } else
      {
        equal++;
      }
    }
  }
  if (smaller == larger)
  {
    return;
  }
  for (equalpart = 0; equalpart < equal; equalpart++)
  {
    sum1 = smaller + equalpart;
    sum2 = larger + (equal-1) - equalpart;
    if (sum1 < sum2)
    {
      if (sum1 + 1 == sum2)
      {
        return;
      }
    } else
    {
      if (sum1 > sum2)
      {
        if (sum1 == sum2 + 1)
        {
          return;
        }
      } else
      {
        return;
      }
    }
  }
  fprintf(stderr,"problem with equal=%lu,smaller=%lu,larger=%lu\n",
                  equal,smaller,larger);
  exit(GT_EXIT_PROGRAMMING_ERROR);
}
#endif

static unsigned long realmedian(const GtBentsedgresources *bsr,
                                unsigned long subbucketleft,
                                unsigned long width,
                                unsigned long depth)
{
  GtMedianinfo *medianptr;
  unsigned long idx;

  for (idx = 0; idx < width; idx++)
  {
    bsr->medianinfospace[idx].suftaboffset = idx;
    PTR2INT(bsr->medianinfospace[idx].etbe,subbucketleft,idx);
  }
  medianptr = quickmedian(bsr->fwd,bsr->complement,bsr->medianinfospace,width);
/*
  checkmedian(bsr->fwd,bsr->complement,medianptr,medianinfospace,width);
*/
  gt_assert(medianptr != NULL);
  return medianptr->suftaboffset;
}

#define MINMEDIANOF9WIDTH 31UL

static unsigned long cmpcharbychardelivermedian(const GtBentsedgresources *bsr,
                                                unsigned long subbucketleft,
                                                unsigned long width,
                                                unsigned long depth)
{
  unsigned long pl = 0,
                pm = GT_DIV2(width),
                pr = width - 1;

  if (width >= MINMEDIANOF9WIDTH)
  { /* On big arrays, pseudomedian of 9 */
    unsigned long offset, doubleoffset;
    offset = GT_DIV8(width);
    doubleoffset = GT_MULT2(offset);
    pl = medianof3cmpcharbychar(bsr,subbucketleft,depth,pl,pl+offset,
                                pl+doubleoffset);
    pm = medianof3cmpcharbychar(bsr,subbucketleft,depth,pm-offset,
                                pm,pm+offset);
    pr = medianof3cmpcharbychar(bsr,subbucketleft,depth,
                                pr-doubleoffset,pr-offset,
                                pr);
  }
  return medianof3cmpcharbychar(bsr,subbucketleft,depth,pl,pm,pr);
}

static unsigned long blockcmpdelivermedian(const GtBentsedgresources *bsr,
                                           unsigned long subbucketleft,
                                           unsigned long width,
                                           unsigned long depth,
                                           unsigned long maxwidthrealmedian)
{
  unsigned long pl = 0,
                pm = GT_DIV2(width),
                pr = width - 1;

  if (width >= MINMEDIANOF9WIDTH)
  {
    if (width > maxwidthrealmedian)
    { /* On big arrays, pseudomedian of 9 */
      unsigned long offset, doubleoffset;
      offset = GT_DIV8(width);
      doubleoffset = GT_MULT2(offset);
      pl = medianof3(bsr,subbucketleft,depth,pl,pl+offset,
                     pl+doubleoffset);
      pm = medianof3(bsr,subbucketleft,depth,pm-offset,pm,pm+offset);
      pr = medianof3(bsr,subbucketleft,depth,pr-doubleoffset,
                     pr-offset,pr);
      pm = medianof3(bsr,subbucketleft,depth,pl,pm,pr);
    } else /* width <= maxwidthrealmedian */
    {
      pm = realmedian(bsr, subbucketleft,width, depth);
    }
  } else
  {
    pm = medianof3(bsr,subbucketleft,depth,pl,pm,pr);
  }
  return pm;
}

/*
static void showcountingsortinfo(const GtCountingsortinfo *countingsortinfo,
                              unsigned long idx)
{
  printf("countingsortinfo[%lu]=(%lu,",idx,
          (unsigned long) countingsortinfo[idx].suffix);
  printf("%lu,",(unsigned long) countingsortinfo[idx].lcpwithpivot);
  printf("%d)\n",countingsortinfo[idx].cmpresult);
}
*/

#undef CHECKFORWHOLELEAFS
#ifdef CHECKFORWHOLELEAFS
static bool gt_containswholeleaf(const GtBentsedgresources *bsr,
                                 unsigned long subbucketleft,
                                 unsigned long width)
{
  unsigned long position, idx;

  for (idx = 0;  idx < width; idx++)
  {
    position = gt_suffixsortspace_get(bsr->sssp,subbucketleft,idx);
    if (position == 0 || gt_encseq_position_is_separator(bsr->encseq,
                                                         position - 1,
                                                         bsr->readmode))
    {
      return true;
    }
  }
  return false;
}
static unsigned long saved_intervals = 0, saved_width = 0;
#endif

static bool multistrategysort(GtBentsedgresources *bsr,
                              unsigned long subbucketleft,
                              unsigned long width,
                              unsigned long depth,
                              unsigned long sortmaxdepth)
{
  gt_assert(width > 1UL);
  if (width <= bsr->sfxstrategy->maxinsertionsort)
  {
    if (sortmaxdepth == 0)
    {
      bs_insertionsort(bsr,subbucketleft,width,depth);
    } else
    {
      bs_insertionsortmaxdepth(bsr,subbucketleft,width,depth,sortmaxdepth);
    }
    return true;
  }
  if (width <= bsr->sfxstrategy->maxbltriesort)
  {
    gt_blindtrie_suffixsort(bsr->blindtrie,
                            subbucketleft,
                            bsr->tableoflcpvalues,
                            width,
                            depth,
                            sortmaxdepth,
                            bsr->processunsortedsuffixrangeinfo,
                            bsr->processunsortedsuffixrange);
    bsr->countbltriesort++;
    return true;
  }
  return false;
}

static void subsort_bentleysedgewick(GtBentsedgresources *bsr,
                                     unsigned long subbucketleft,
                                     unsigned long width,
                                     unsigned long depth)
{
  gt_assert(width > 1UL);
  if (width > 1UL)
  {
#ifdef SKDEBUG
    if (depth > 0)
    {
      gt_checkifprefixesareidentical(__FILE__,
                                     __LINE__,
                                     bsr->encseq,
                                     bsr->readmode,
                                     bsr->sssp,
                                     subbucketleft,
                                     width,
                                     depth);
    }
#endif
#ifdef CHECKFORWHOLELEAFS
    if (bsr->sfxstrategy->spmopt > 0 &&
        !gt_containswholeleaf(bsr,subbucketleft,width))
    {
      saved_intervals++;
      saved_width += width;
      return;
    }
#endif
    if (bsr->rsi != NULL && width >= bsr->radixsortminwidth
                         && width <= bsr->radixsortmaxwidth)
    {
      gt_sfx_radixsort_str(bsr->rsi,
                           depth,
                           bsr->sortmaxdepth,
                           subbucketleft,
                           width,
                           bsr->sssp,
                           bsr->tableoflcpvalues);
      bsr->countradixsort++;
      return;
    }
    if (bsr->srsw != NULL &&
        !bsr->sfxstrategy->noshortreadsort &&
        width <= bsr->shortreadsort_maxwidth)
    {
      gt_shortreadsort_sssp_sort(bsr->srsw,
                                 bsr->encseq,
                                 bsr->srs_maxremain,
                                 bsr->readmode,
                                 bsr->esr1,
                                 bsr->sssp,
                                 subbucketleft,
                                 width,
                                 depth,
                                 (unsigned long) bsr->sortmaxdepth);
      if (bsr->sortmaxdepth > 0)
      {
        gt_shortreadsort_sssp_add_unsorted(
                              bsr->srsw,
                              gt_suffixsortspace_bucketleftidx_get(bsr->sssp),
                              subbucketleft,
                              width,
                              (unsigned long) bsr->sortmaxdepth,
                              bsr->processunsortedsuffixrange,
                              bsr->processunsortedsuffixrangeinfo);
      }
      bsr->countshortreadsort++;
      return;
    }
    if (bsr->sortmaxdepth > 0 && depth >= (unsigned long) bsr->sortmaxdepth)
    {
      if (bsr->processunsortedsuffixrange != NULL)
      {
        bsr->processunsortedsuffixrange(
                         bsr->processunsortedsuffixrangeinfo,
                         gt_suffixsortspace_bucketleftidx_get(bsr->sssp) +
                         subbucketleft,
                         width,depth);
      }
      return;
    }
    if (multistrategysort(bsr,subbucketleft,width,depth,
                          (unsigned long) bsr->sortmaxdepth))
    {
      return;
    }
    /* push */
    GT_CHECKARRAYSPACE(&bsr->mkvauxstack,GtMKVstack,1024);
    STACKTOP.subbucketleft = subbucketleft;
    STACKTOP.width = width;
    STACKTOP.depth = depth;
    bsr->mkvauxstack.nextfreeGtMKVstack++;
  }
}

static void sarrcountingsort(GtBentsedgresources *bsr,
                             unsigned long subbucketleft,
                             unsigned long width,
                             const GtSfxcmp *pivotcmpbits,
                             unsigned long pivotidx,
                             unsigned long depth)
{
  int cmp;
  unsigned int maxsmallerwithlcp = 0, maxlargerwithlcp = 0;
  GtCommonunits commonunits;
  GtEndofTwobitencoding etbecurrent;
  unsigned long idx, smaller = 0, larger = 0,
                insertindex, end, equaloffset, currentwidth;
  GtCountingsortinfo *csiptr;
  /* const bool cmpcharbychar = false; */

  bsr->countcountingsort++;
  for (idx = 0; idx < width; idx++)
  {
    if (idx != pivotidx)
    {
      PTR2INTSTOREPOS(etbecurrent,subbucketleft,idx,
                      bsr->countingsortinfo[idx].suffix = pos);
      cmp = gt_encseq_compare_pairof_twobitencodings(bsr->fwd,
                                                     bsr->complement,
                                                     &commonunits,
                                                     &etbecurrent,
                                                     pivotcmpbits);
      gt_assert(commonunits.common <= (unsigned int) GT_UNITSIN2BITENC);
      bsr->countingsortinfo[idx].lcpwithpivot = commonunits.common;
      if (cmp > 0)
      {
        gt_assert(commonunits.common < (unsigned int) GT_UNITSIN2BITENC);
        bsr->rightlcpdist[commonunits.common]++;
        if (maxlargerwithlcp < commonunits.common)
        {
          maxlargerwithlcp = commonunits.common;
        }
        bsr->countingsortinfo[idx].cmpresult = (char) 1;
        larger++;
      } else
      {
        if (cmp < 0)
        {
          gt_assert(commonunits.common < (unsigned int) GT_UNITSIN2BITENC);
          bsr->leftlcpdist[commonunits.common]++;
          if (maxsmallerwithlcp < commonunits.common)
          {
            maxsmallerwithlcp = commonunits.common;
          }
          bsr->countingsortinfo[idx].cmpresult = (char) -1;
          smaller++;
        } else
        {
          gt_assert(commonunits.common == (unsigned int) GT_UNITSIN2BITENC);
          bsr->countingsortinfo[idx].cmpresult = 0;
        }
      }
    } else
    {
      bsr->countingsortinfo[idx].suffix
        = gt_suffixsortspace_get(bsr->sssp,subbucketleft,idx);
      bsr->countingsortinfo[idx].lcpwithpivot = (unsigned char)
                                                GT_UNITSIN2BITENC;
      bsr->countingsortinfo[idx].cmpresult = (char) 0;
    }
  }
  for (idx = 1UL; idx <= (unsigned long) maxsmallerwithlcp; idx++)
  {
    bsr->leftlcpdist[idx] += bsr->leftlcpdist[idx-1];
  }
  for (idx = 1UL; idx <= (unsigned long) maxlargerwithlcp; idx++)
  {
    bsr->rightlcpdist[idx] += bsr->rightlcpdist[idx-1];
  }
  equaloffset = width - larger;
  for (csiptr = bsr->countingsortinfo + width -1;
       csiptr >= bsr->countingsortinfo;
       csiptr--)
  {
    switch (csiptr->cmpresult)
    {
      case -1:
        insertindex = --(bsr->leftlcpdist[csiptr->lcpwithpivot]);
        gt_suffixsortspace_set(bsr->sssp,subbucketleft,insertindex,
                               csiptr->suffix);
        break;
      case 0:
        gt_suffixsortspace_set(bsr->sssp,subbucketleft,--equaloffset,
                               csiptr->suffix);
        break;
      case 1:
        insertindex = --(bsr->rightlcpdist[csiptr->lcpwithpivot]);
        gt_suffixsortspace_set(bsr->sssp,subbucketleft,
                               width - 1 - insertindex,
                               csiptr->suffix);
        break;
    }
  }
  for (idx = 0; idx <= (unsigned long) maxsmallerwithlcp; idx++)
  {
    if (idx < (unsigned long) maxsmallerwithlcp)
    {
      end = bsr->leftlcpdist[idx+1];
    } else
    {
      end = smaller;
    }
    if (bsr->leftlcpdist[idx] + 1 < end) /* at least two elements */
    {
      currentwidth = end - bsr->leftlcpdist[idx];
      subsort_bentleysedgewick(bsr,
                               subbucketleft + bsr->leftlcpdist[idx],
                               currentwidth,
                               depth + idx);
    }
    if (bsr->tableoflcpvalues != NULL && bsr->leftlcpdist[idx] < end)
    { /* at least one element */
      gt_lcptab_update(bsr->tableoflcpvalues,subbucketleft,end,depth + idx);
    }
    bsr->leftlcpdist[idx] = 0;
  }
  if (width - smaller - larger > 1UL)
  {
    currentwidth = width - smaller - larger;
    subsort_bentleysedgewick(bsr,
                             subbucketleft + smaller,
                             currentwidth,
                             depth + GT_UNITSIN2BITENC);
  }
  for (idx = 0; idx <= (unsigned long) maxlargerwithlcp; idx++)
  {
    if (idx < (unsigned long) maxlargerwithlcp)
    {
      end = bsr->rightlcpdist[idx+1];
    } else
    {
      end = larger;
    }
    if (bsr->rightlcpdist[idx] + 1 < end) /* at least two elements */
    {
      currentwidth = end - bsr->rightlcpdist[idx];
      subsort_bentleysedgewick(bsr,
                               /* for largest idx use first bucket */
                               subbucketleft + width - end,
                               currentwidth,
                               depth + idx);
    }
    if (bsr->tableoflcpvalues != NULL && bsr->rightlcpdist[idx] < end)
    { /* at least one element */
      gt_assert(width >= end);
      gt_lcptab_update(bsr->tableoflcpvalues,subbucketleft,width - end,
                    depth + idx);
    }
    bsr->rightlcpdist[idx] = 0;
  }
}

static inline void vectorswap(GtSuffixsortspace *sssp,
                              unsigned long subbucketleft1,
                              unsigned long subbucketleft2,
                              unsigned long width)
{
  unsigned long idx, tmp;

  for (idx = 0; idx < width; idx++)
  {
    tmp = gt_suffixsortspace_get(sssp,subbucketleft1,idx);
    gt_suffixsortspace_set(sssp,subbucketleft1,idx,
                           gt_suffixsortspace_get(sssp,subbucketleft2,idx));
    gt_suffixsortspace_set(sssp,subbucketleft2,idx,tmp);
  }
}

static void gt_sort_bentleysedgewick(GtBentsedgresources *bsr,
                                     unsigned long width,
                                     unsigned long depth)
{
  bsr->mkvauxstack.nextfreeGtMKVstack = 0;
/*#ifndef NDEBUG
  gt_suffixsortspace_checkorder(bsr->sssp,0,width);
#endif*/
  subsort_bentleysedgewick(bsr, 0, width, depth);
  while (bsr->mkvauxstack.nextfreeGtMKVstack > 0)
  {
    unsigned long leftplusw, pa, pb, pc, pd, pm, bucketright,
                  cptr, temp, pivotcmpcharbychar = 0, valcmpcharbychar,
                  wtmp, subbucketleft;
    unsigned int smallermaxlcp, greatermaxlcp, smallerminlcp, greaterminlcp;
    GtSfxcmp pivotcmpbits, val;
    int retvalpivotcmpbits;
    GtUchar tmpvar;
    GtCommonunits commonunits;
    const int commonunitsequal = bsr->sfxstrategy->cmpcharbychar
                                 ? 1
                                 : GT_UNITSIN2BITENC;

    /* pop */
    bsr->mkvauxstack.nextfreeGtMKVstack--;
    subbucketleft = STACKTOP.subbucketleft;
    width = STACKTOP.width;
    depth = STACKTOP.depth;
    bucketright = width - 1;

    if (bsr->sfxstrategy->cmpcharbychar)
    {
      pm = cmpcharbychardelivermedian(bsr, subbucketleft, width, depth);
      BS_SWAPARRAY(temp, subbucketleft, 0, pm);
      CMPCHARBYCHARPTR2INT(pivotcmpcharbychar, subbucketleft,tmpvar,0);
    } else
    {
      pm = blockcmpdelivermedian(bsr,
                                 subbucketleft,
                                 width,
                                 depth,
                                 bsr->sfxstrategy->maxwidthrealmedian);
      if (width <= bsr->sfxstrategy->maxcountingsort &&
          width >= MINMEDIANOF9WIDTH)
      {
        PTR2INT(pivotcmpbits,subbucketleft,pm);
        gt_assert(width >= bsr->sfxstrategy->maxbltriesort);
        sarrcountingsort(bsr,
                         subbucketleft,
                         width,
                         &pivotcmpbits,
                         pm,
                         depth);
        /* new values for subbucketleft, bucketright, depth */
        continue;
      }
      BS_SWAPARRAY(temp, subbucketleft, 0, pm);
      PTR2INT(pivotcmpbits,subbucketleft,0UL);
    }
    bsr->counttqsort++;
    /* now pivot element is at index subbucketleft */
    /* all elements to be compared are between pb and pc */
    /* pa is the position at which the next element smaller than the
       pivot element is inserted at */
    /* pd is the position at which the next element greater than the
       pivot element is inserted at */
    pa = pb = 1UL;
    pc = pd = bucketright;
    if (bsr->sfxstrategy->cmpcharbychar)
    {
      smallerminlcp = greaterminlcp = smallermaxlcp = greatermaxlcp = 0;
      for (;;)
      {
        while (pb <= pc)
        {
          CMPCHARBYCHARPTR2INT(valcmpcharbychar,subbucketleft,tmpvar,pb);
          if (valcmpcharbychar > pivotcmpcharbychar)
          {
            break;
          }
          if (valcmpcharbychar == pivotcmpcharbychar)
          {
            BS_SWAPARRAY(temp, subbucketleft, pa, pb);
            pa++;
          }
          pb++;
        }
        while (pb <= pc)
        {
          CMPCHARBYCHARPTR2INT(valcmpcharbychar,subbucketleft,tmpvar,pc);
          if (valcmpcharbychar < pivotcmpcharbychar)
          { /* stop for elements < pivot */
            break;
          }
          if (valcmpcharbychar == pivotcmpcharbychar)
          {
            /* exchange equal element and element at index pd */
            BS_SWAPARRAY(temp, subbucketleft, pc, pd);
            pd--;
          }
          pc--;
        }
        if (pb > pc)
        { /* no elements to compare to pivot */
          break;
        }
        BS_SWAPARRAY(temp, subbucketleft, pb, pc);
        pb++;
        pc--;
      }
    } else
    {
      smallermaxlcp = greatermaxlcp = 0;
      smallerminlcp = greaterminlcp = (unsigned int) GT_UNITSIN2BITENC;
      for (;;)
      {
        /* look for elements identical or smaller than pivot from left */
        while (pb <= pc)
        {
          PTR2INT(val,subbucketleft,pb);
          Sfxdocompare(&commonunits,val,pivotcmpbits);
          if (GtSfxcmpGREATER(val,pivotcmpbits))
          { /* stop for elements val > pivot */
            GT_BSR_UPDATELCP(greaterminlcp,greatermaxlcp,commonunits.common);
            break;
          }
          if (GtSfxcmpEQUAL(val,pivotcmpbits))
          {
            /* exchange equal element and element at index pa */
            BS_SWAPARRAY(temp, subbucketleft, pa, pb);
            pa++;
          } else /* smaller */
          {
            GT_BSR_UPDATELCP(smallerminlcp,smallermaxlcp,commonunits.common);
          }
          pb++;
        }
        /* look for elements identical or greater than pivot from right */
        while (pb <= pc)
        {
          PTR2INT(val,subbucketleft,pc);
          Sfxdocompare(&commonunits,val,pivotcmpbits);
          if (GtSfxcmpSMALLER(val,pivotcmpbits))
          { /* stop for elements val < pivot */
            GT_BSR_UPDATELCP(smallerminlcp,smallermaxlcp,commonunits.common);
            break;
          }
          if (GtSfxcmpEQUAL(val,pivotcmpbits))
          {
            /* exchange equal element and element at index pa */
            BS_SWAPARRAY(temp, subbucketleft, pc, pd);
            pd--;
          } else /* greater */
          {
            GT_BSR_UPDATELCP(greaterminlcp,greatermaxlcp,commonunits.common);
          }
          pc--;
        }
        if (pb > pc)
        { /* interval is empty */
          break;
        }
        BS_SWAPARRAY(temp, subbucketleft, pb, pc);
        pb++;
        pc--;
      }
    }
    gt_assert(pb >= pa);
    wtmp = MIN(pa,pb-pa);
    /* move w elements at the left to the middle */
    gt_assert(pb >= wtmp && wtmp <= pb - wtmp);
    vectorswap(bsr->sssp, subbucketleft, subbucketleft+pb-wtmp, wtmp);
    gt_assert(pd >= pc);
    gt_assert(bucketright >= pd);
    wtmp = MIN(pd-pc, bucketright-pd);
    /* move w elements at the right to the middle */
    gt_assert(bucketright + 1 >= wtmp && pb + wtmp <= bucketright+1-wtmp);
    vectorswap(bsr->sssp, subbucketleft+pb, subbucketleft+bucketright+1-wtmp,
               wtmp);

    /* all elements equal to the pivot are now in the middle namely in the
       range [subbucketleft + (pb-pa) and bucketright - (pd-pc)] */
    /* hence we have to sort the elements in the intervals
       [subbucketleft..subbucketleft+(pb-pa)-1] and
       [bucketright-(pd-pc)+1..bucketright] */

    gt_assert(pb >= pa);
    if ((wtmp = pb-pa) > 0)
    {
      leftplusw = wtmp;
      if (bsr->tableoflcpvalues != NULL)
      {
        /*
          left part has suffix with lcp up to length smallermaxlcp w.r.t.
          to the pivot. This lcp belongs to a suffix on the left
          which is at a minimum distance to the pivot and thus to an
          element in the final part of the left side.
        */
        gt_lcptab_update(bsr->tableoflcpvalues,subbucketleft,leftplusw,
                         depth + smallermaxlcp);
      }
      if (wtmp > 1UL)
      {
        subsort_bentleysedgewick(bsr,
                                 subbucketleft,
                                 wtmp,
                                 depth + smallerminlcp);
      }
    } else
    {
      leftplusw = 0;
    }

    cptr = gt_suffixsortspace_get(bsr->sssp,subbucketleft,leftplusw) + depth;
    if (ISNOTEND(cptr) && (wtmp = bucketright-(pd-pb)-leftplusw) > 1UL)
    {
      subsort_bentleysedgewick(bsr,
                               subbucketleft + leftplusw,
                               wtmp,
                               depth+commonunitsequal);
    }

    gt_assert(pd >= pc);
    if ((wtmp = (unsigned long) (pd-pc)) > 0)
    {
      if (bsr->tableoflcpvalues != NULL)
      {
        /*
          right part has suffix with lcp up to length largermaxlcp w.r.t.
          to the pivot. This lcp belongs to a suffix on the right
          which is at a minimum distance to the pivot and thus to an
          element in the first part of the right side.
        */
        gt_lcptab_update(bsr->tableoflcpvalues,
                         subbucketleft,bucketright - wtmp + 1,
                         depth + greatermaxlcp);
      }
      if (wtmp > 1UL)
      {
        subsort_bentleysedgewick(bsr,
                                 subbucketleft + bucketright - wtmp + 1,
                                 wtmp,
                                 depth + greaterminlcp);
      }
    }
  }
}

static void bentsedgresources_init(GtBentsedgresources *bsr,
                                   GtSuffixsortspace *suffixsortspace,
                                   const GtEncseq *encseq,
                                   GtReadmode readmode,
                                   unsigned int prefixlength,
                                   unsigned int sortmaxdepth,
                                   const Sfxstrategy *sfxstrategy,
                                   bool withlcps)
{
  unsigned long idx;

  bsr->readmode = readmode;
  bsr->totallength = gt_encseq_total_length(encseq);
  bsr->sfxstrategy = sfxstrategy;
  bsr->sssp = suffixsortspace;
  bsr->rsi = NULL;
  gt_suffixsortspace_bucketleftidx_set(bsr->sssp,0);
  bsr->encseq = encseq;
  bsr->fwd = GT_ISDIRREVERSE(bsr->readmode) ? false : true;
  bsr->complement = GT_ISDIRCOMPLEMENT(bsr->readmode) ? true : false;
  bsr->tableoflcpvalues = NULL;
  bsr->prefixlength = prefixlength;
  bsr->esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  bsr->esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  bsr->radixsortminwidth = gt_radixsort_str_minwidth();
  bsr->radixsortmaxwidth = 0;
  if (gt_encseq_accesstype_get(bsr->encseq) == GT_ACCESS_TYPE_EQUALLENGTH)
  {
    bsr->twobitencoding = gt_encseq_twobitencoding_export(bsr->encseq);
  } else
  {
    bsr->twobitencoding = NULL;
  }
  if (sortmaxdepth > 0)
  {
    if (sortmaxdepth > prefixlength)
    {
      bsr->srs_maxremain = sortmaxdepth - (unsigned long) prefixlength;
    } else
    {
      bsr->srs_maxremain = 0;
    }
  } else
  {
    if (gt_encseq_lengthoflongestnonspecial(encseq) >
        (unsigned long) prefixlength)
    {
      bsr->srs_maxremain = gt_encseq_lengthoflongestnonspecial(encseq) -
                           prefixlength;
    } else
    {
      bsr->srs_maxremain = 0;
    }
  }
  bsr->sizeofworkspace = gt_size_of_sort_workspace (sfxstrategy);
  bsr->shortreadsort_maxwidth = gt_shortreadsort_maxwidth(false,
                                                          bsr->srs_maxremain,
                                                          bsr->sizeofworkspace);
  GT_INITARRAY(&bsr->mkvauxstack,GtMKVstack);
  bsr->countingsortinfo = NULL;
  bsr->medianinfospace = NULL;
  bsr->blindtrie = NULL;
  bsr->equalwithprevious = NULL;
  bsr->srsw = NULL;
  if (!sfxstrategy->cmpcharbychar)
  {
    const bool withmediumsizelcps = (sortmaxdepth > 0 && !withlcps)
                                    ? true : false;
    if (sortmaxdepth == 0 ||
        gt_encseq_has_twobitencoding_stoppos_support(encseq))
    {
      bsr->srsw = gt_shortreadsort_new(0,bsr->srs_maxremain,readmode,false,
                                       withmediumsizelcps);
    }
    for (idx = 0; idx < (unsigned long) GT_UNITSIN2BITENC; idx++)
    {
      bsr->leftlcpdist[idx] = bsr->rightlcpdist[idx] = 0;
    }
    bsr->countingsortinfo = gt_malloc(sizeof (*bsr->countingsortinfo) *
                                      sfxstrategy->maxcountingsort);
    if (sfxstrategy->maxwidthrealmedian >= MINMEDIANOF9WIDTH)
    {
      bsr->medianinfospace = gt_malloc(sizeof (*bsr->medianinfospace) *
                                       sfxstrategy->maxwidthrealmedian);
    } else
    {
      bsr->medianinfospace = NULL;
    }
  }
  bsr->blindtrie = gt_blindtrie_new(bsr->sssp,
                                    sfxstrategy->maxbltriesort,
                                    0, /* the nodenumberincrement */
                                    encseq,
                                    sfxstrategy->cmpcharbychar,
                                    bsr->esr1,
                                    bsr->esr2,
                                    readmode);
  if (sortmaxdepth > 0 && bsr->sfxstrategy->maxinsertionsort >= 2UL)
  {
    bsr->equalwithprevious = gt_malloc(sizeof (*bsr->equalwithprevious) *
                                       bsr->sfxstrategy->maxinsertionsort);
    for (idx=0; idx < bsr->sfxstrategy->maxinsertionsort; idx++)
    {
      bsr->equalwithprevious[idx] = false;
    }
  }
  bsr->sortmaxdepth = sortmaxdepth;
  bsr->processunsortedsuffixrangeinfo = NULL;
  bsr->processunsortedsuffixrange = NULL;
  bsr->countinsertionsort = 0;
  bsr->counttqsort = 0;
  bsr->countcountingsort = 0;
  bsr->countbltriesort = 0;
  bsr->countshortreadsort = 0;
  bsr->countradixsort = 0;
}

size_t gt_size_of_sort_workspace (const Sfxstrategy *sfxstrategy)
{
  size_t sumsize = 0;

  if (!sfxstrategy->cmpcharbychar)
  {
    sumsize += sizeof (GtCountingsortinfo) * sfxstrategy->maxcountingsort;
    if (sfxstrategy->maxwidthrealmedian >= MINMEDIANOF9WIDTH)
    {
      sumsize += sizeof (GtMedianinfo) * sfxstrategy->maxwidthrealmedian;
    }
  }
  sumsize += gt_blindtrie_size(sfxstrategy->maxbltriesort);
  return sumsize;
}

static void bentsedgresources_delete(GtBentsedgresources *bsr, GtLogger *logger)
{
  gt_free(bsr->countingsortinfo);
  bsr->countingsortinfo = NULL;
  gt_free(bsr->medianinfospace);
  bsr->medianinfospace = NULL;
  gt_shortreadsort_delete(bsr->srsw);
  bsr->srsw = NULL;
  gt_blindtrie_delete(bsr->blindtrie);
  gt_encseq_reader_delete(bsr->esr1);
  gt_encseq_reader_delete(bsr->esr2);
  gt_free(bsr->equalwithprevious);
  GT_FREEARRAY(&bsr->mkvauxstack,GtMKVstack);
  gt_radixsort_str_delete(bsr->rsi);
  gt_logger_log(logger,"countinsertionsort=%lu",bsr->countinsertionsort);
  gt_logger_log(logger,"countbltriesort=%lu",bsr->countbltriesort);
  gt_logger_log(logger,"countcountingsort=%lu",bsr->countcountingsort);
  gt_logger_log(logger,"countshortreadsort=%lu",bsr->countshortreadsort);
  gt_logger_log(logger,"countradixsort=%lu",bsr->countradixsort);
  gt_logger_log(logger,"counttqsort=%lu",bsr->counttqsort);
}

/*
  The following function is called in sfx-suffixer.c and sorts all buckets by
  different suffix comparison methods without the help of other sorting
  information. GtSuffixsortspace contains the sortspace which is accessed
  by some negative offset.
*/

void gt_sortallbuckets(GtSuffixsortspace *suffixsortspace,
                       unsigned long numberofsuffixes,
                       GtBucketspec2 *bucketspec2,
                       const GtEncseq *encseq,
                       GtReadmode readmode,
                       GtCodetype mincode,
                       GtCodetype maxcode,
                       GtBcktab *bcktab,
                       unsigned int numofchars,
                       unsigned int prefixlength,
                       GtOutlcpinfo *outlcpinfo,
                       unsigned int sortmaxdepth,
                       const Sfxstrategy *sfxstrategy,
                       GtProcessunsortedsuffixrange processunsortedsuffixrange,
                       void *processunsortedsuffixrangeinfo,
                       unsigned long long *bucketiterstep,
                       GtLogger *logger)
{
  GtCodetype code;
  unsigned int rightchar = (unsigned int) (mincode % numofchars);
  GtBucketspecification bucketspec;
  GtBentsedgresources bsr;

  bentsedgresources_init(&bsr,
                         suffixsortspace,
                         encseq,
                         readmode,
                         prefixlength,
                         sortmaxdepth,
                         sfxstrategy,
                         outlcpinfo != NULL ? true : false);
  gt_bcktab_determinemaxsize(bcktab, mincode, maxcode, numberofsuffixes);
  if (bsr.sfxstrategy->withradixsort &&
      gt_encseq_accesstype_get(bsr.encseq) == GT_ACCESS_TYPE_EQUALLENGTH &&
      bsr.readmode == GT_READMODE_FORWARD)
  {
    bsr.rsi = gt_radixsort_str_new(bsr.twobitencoding,
                                   gt_encseq_is_mirrored(encseq)
                                     ? GT_DIV2(bsr.totallength - 1)
                                     : bsr.totallength,
                                   1 + gt_encseq_equallength(bsr.encseq),
                                   gt_bcktab_nonspecialsmaxsize(bcktab));
    bsr.radixsortmaxwidth = gt_radixsort_str_maxwidth(bsr.rsi);
  }
  if (outlcpinfo != NULL)
  {
    bsr.tableoflcpvalues = gt_Outlcpinfo_resizereservoir(outlcpinfo,bcktab);
    if (bsr.srsw != NULL)
    {
      gt_shortreadsort_assigntableoflcpvalues(bsr.srsw,bsr.tableoflcpvalues);
    }
  }
  bsr.processunsortedsuffixrangeinfo = processunsortedsuffixrangeinfo;
  bsr.processunsortedsuffixrange = processunsortedsuffixrange;
  for (code = mincode; code <= maxcode; code++)
  {
    if (bucketspec2 != NULL)
    {
      if (gt_copysort_checkhardwork(bucketspec2,code))
      {
        rightchar = (unsigned int) (code % numofchars);
      } else
      {
        continue;
      }
    }
    (*bucketiterstep)++;
    rightchar = gt_bcktab_calcboundsparts(&bucketspec,
                                          bcktab,
                                          code,
                                          maxcode,
                                          numberofsuffixes,
                                          rightchar);
    gt_Outlcpinfo_prebucket(outlcpinfo,code,bucketspec.left);
    if (bucketspec.nonspecialsinbucket > 0)
    {
      if (bucketspec.nonspecialsinbucket > 1UL)
      {
        gt_suffixsortspace_bucketleftidx_set(bsr.sssp,bucketspec.left);
        gt_sort_bentleysedgewick(&bsr,bucketspec.nonspecialsinbucket,
                                 (unsigned long) prefixlength);
        gt_suffixsortspace_bucketleftidx_set(bsr.sssp,0);
      }
      gt_Outlcpinfo_nonspecialsbucket(outlcpinfo,
                                      prefixlength,
                                      bsr.sssp,
                                      bsr.tableoflcpvalues,
                                      &bucketspec,
                                      code);
    }
    gt_Outlcpinfo_postbucket(outlcpinfo,
                             prefixlength,
                             bsr.sssp,
                             bcktab,
                             &bucketspec,
                             code);
  }
#ifdef CHECKFORWHOLELEAFS
  printf("saved_intervals=%lu,saved_width=%lu (%.2f)\n",
          saved_intervals,saved_width,100.0 *
                                      (double) saved_width/numberofsuffixes);
  saved_intervals = 0;
  saved_width = 0;
#endif
  bentsedgresources_delete(&bsr, logger);
}

void gt_sortallsuffixesfromstart(GtSuffixsortspace *suffixsortspace,
                                 unsigned long numberofsuffixes,
                                 const GtEncseq *encseq,
                                 GtReadmode readmode,
                                 GtOutlcpinfo *outlcpinfo,
                                 unsigned int sortmaxdepth,
                                 const Sfxstrategy *sfxstrategy,
                                 GtProcessunsortedsuffixrange
                                   processunsortedsuffixrange,
                                 void *processunsortedsuffixrangeinfo,
                                 GtLogger *logger)
{
  GtBentsedgresources bsr;

  if (numberofsuffixes > 1UL)
  {
    bentsedgresources_init(&bsr,
                           suffixsortspace,
                           encseq,
                           readmode,
                           0,
                           sortmaxdepth,
                           sfxstrategy,
                           outlcpinfo != NULL ? true : false);
    if (outlcpinfo != NULL)
    {
      bsr.tableoflcpvalues = gt_Outlcpinfo_resizereservoir(outlcpinfo,NULL);
      if (bsr.srsw != NULL)
      {
        gt_shortreadsort_assigntableoflcpvalues(bsr.srsw,bsr.tableoflcpvalues);
      }
    }
    bsr.processunsortedsuffixrangeinfo = processunsortedsuffixrangeinfo;
    bsr.processunsortedsuffixrange = processunsortedsuffixrange;
    gt_suffixsortspace_bucketleftidx_set(bsr.sssp,0);
    gt_sort_bentleysedgewick(&bsr,numberofsuffixes,0);
    gt_suffixsortspace_bucketleftidx_set(bsr.sssp,0);
    bentsedgresources_delete(&bsr, logger);
  }
}
