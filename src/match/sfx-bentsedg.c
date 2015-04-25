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
#ifdef GT_THREADS_ENABLED
#include "core/thread_api.h"
#endif

#define ACCESSCHARRAND(POS)    gt_encseq_get_encoded_char(bsr->encseq,\
                                                          POS,bsr->readmode)
#define ACCESSCHARSEQ(ESR)     gt_encseq_reader_next_encoded_char(ESR)
#define ISNOTEND(POS)          ((POS) < bsr->totallength &&\
                                ISNOTSPECIAL(ACCESSCHARRAND(POS)))

#define DEREFSTOPPOSSEQ(VAR,POS,STOPPOS,ESR)\
        (((POS) < (STOPPOS) && ISNOTSPECIAL(VAR = ACCESSCHARSEQ(ESR))) ?\
        ((GtUword) VAR) : GT_UNIQUEINT(POS))

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
                    ? ((GtUword) TMPVAR) : GT_UNIQUEINT(cptr))

typedef GtEndofTwobitencoding GtSfxcmp;

#define PTR2INTSTOREPOS(TMPVAR,SUBBUCKETLEFT,IDX,POSASSIGNMENT)\
        {\
          GtUword pos\
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

#define GtSfxcmpEQUAL(X,Y)      ret##X##Y == 0
#define GtSfxcmpSMALLER(X,Y)    ret##X##Y < 0
#define GtSfxcmpGREATER(X,Y)    ret##X##Y > 0

typedef struct
{
  GtUword subbucketleft,
                width,
                depth;
} GtMKVstack;

typedef struct
{
  GtEndofTwobitencoding etbe;
  GtUword suftaboffset;
} GtMedianinfo;

typedef GtMedianinfo GtMedianElem;

typedef struct
{
  GtUword suffix;
  unsigned char lcpwithpivot;
  signed char cmpresult;
} GtCountingsortinfo;

GT_DECLAREARRAYSTRUCT(GtMKVstack);

typedef struct
{
  const GtEncseq *encseq;
  GtEncseqReader *esr1,
                 *esr2;
  GtReadmode readmode;
  bool fwd, complement;
  GtArrayGtMKVstack mkvauxstack;
  GtLcpvalues *tableoflcpvalues;
  const Sfxstrategy *sfxstrategy;
  GtSuffixsortspace *sssp;
  GtProcessunsortedsuffixrange processunsortedsuffixrange;
  void *processunsortedsuffixrangeinfo;
  unsigned int sortmaxdepth,
               prefixlength;
  GtUword countinsertionsort,
          counttqsort,
          countshortreadsort,
          countradixsort,
          countcountingsort,
          countbltriesort,
          srs_maxremain, /* only relevant for short read sort */
          radixsortminwidth,
          radixsortmaxwidth,
          totallength,
          shortreadsort_maxwidth,
          leftlcpdist[GT_UNITSIN2BITENC],
          rightlcpdist[GT_UNITSIN2BITENC];
  bool *equalwithprevious;
  GtBlindtrie *blindtrie;
  GtMedianinfo *medianinfospace;
  GtCountingsortinfo *countingsortinfo;
  GtShortreadsortworkinfo *srsw;
  const GtTwobitencoding *twobitencoding;
  GtRadixsortstringinfo *rsi;
} GtBentsedgresources;

#ifdef WITHCHECKSTARTPOINTER
static unsigned int checkstartpointorder(const GtUword *left,
                                         const GtUword *right)
{
  const GtUword *ptr;
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

static GtUword medianof3cmpcharbychar(const GtBentsedgresources *bsr,
                                            GtUword subbucketleft,
                                            GtUword depth,
                                            GtUword a,
                                            GtUword b,
                                            GtUword c)
{
  GtUword vala, valb, valc, cptr;
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

static GtUword medianof3(const GtBentsedgresources *bsr,
                               GtUword subbucketleft,
                               GtUword depth,
                               GtUword a,
                               GtUword b,
                               GtUword c)
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
                             GtUword subbucketleft,
                             GtUword width,
                             GtUword offset)
{
  GtUword sval1, sval2, pm, pl, startpos1, startpos2, lcplen = 0;
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
          GtUword ccs, cct;
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
        printf("%s["GT_WU","GT_WU"] at offset"
               GT_WU"\n",__func__,sval1,sval2,offset);
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
                                     GtUword subbucketleft,
                                     GtUword width,
                                     GtUword offset,
                                     GtUword sortmaxdepth)
{
  GtUword sval1, sval2, pm, pl, startpos1, startpos2,
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
        GtUword endpos1, endpos2;

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
          GtUword ccs, cct;
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
      printf("cmp "GT_WU" and "GT_WU": retval = %d, lcplen = "GT_WU"\n",
             sval1, sval2, retval, (GtUword) lcplen);
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
    GtUword equalsrangewidth = 0,
            bucketleftidx = gt_suffixsortspace_bucketleftidx_get(bsr->sssp);

#ifdef SKDEBUG
    printf("ordered suffix "GT_WU"\n",gt_suffixsortspace_get(bsr->sssp,
                                                         subbucketleft,0));
#endif
    for (idx = 1UL; idx < width; idx++)
    {
#ifdef SKDEBUG
      printf("ordered suffix "GT_WU", equalwithprevious=%s\n",
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
          printf("process interval of width "GT_WU"\n",
                 equalsrangewidth + 1);
#endif
          if (bsr->processunsortedsuffixrange != NULL)
          {
            bsr->processunsortedsuffixrange(bsr->processunsortedsuffixrangeinfo,
                                            bsr->sssp,
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
      printf("process interval of width "GT_WU"\n",
             equalsrangewidth + 1);
#endif
      if (bsr->processunsortedsuffixrange != NULL)
      {
        bsr->processunsortedsuffixrange(bsr->processunsortedsuffixrangeinfo,
                                        bsr->sssp,
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
                                  GtMedianElem *arr,GtUword width)
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
                        GtUword width)
{
  GtUword sum1, sum2, idx, smaller = 0, larger = 0, equal = 0, equalpart;
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
  fprintf(stderr,"problem with equal="GT_WU",smaller="GT_WU",larger="GT_WU"\n",
                  equal,smaller,larger);
  exit(GT_EXIT_PROGRAMMING_ERROR);
}
#endif

static GtUword realmedian(const GtBentsedgresources *bsr,
                                GtUword subbucketleft,
                                GtUword width,
                                GtUword depth)
{
  GtMedianinfo *medianptr;
  GtUword idx;

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

static GtUword cmpcharbychardelivermedian(const GtBentsedgresources *bsr,
                                                GtUword subbucketleft,
                                                GtUword width,
                                                GtUword depth)
{
  GtUword pl = 0,
                pm = GT_DIV2(width),
                pr = width - 1;

  if (width >= MINMEDIANOF9WIDTH)
  { /* On big arrays, pseudomedian of 9 */
    GtUword offset, doubleoffset;
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

static GtUword blockcmpdelivermedian(const GtBentsedgresources *bsr,
                                           GtUword subbucketleft,
                                           GtUword width,
                                           GtUword depth,
                                           GtUword maxwidthrealmedian)
{
  GtUword pl = 0,
                pm = GT_DIV2(width),
                pr = width - 1;

  if (width >= MINMEDIANOF9WIDTH)
  {
    if (width > maxwidthrealmedian)
    { /* On big arrays, pseudomedian of 9 */
      GtUword offset, doubleoffset;
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
                              GtUword idx)
{
  printf("countingsortinfo["GT_WU"]=("GT_WU",",idx,
          (GtUword) countingsortinfo[idx].suffix);
  printf(""GT_WU",",(GtUword) countingsortinfo[idx].lcpwithpivot);
  printf("%d)\n",countingsortinfo[idx].cmpresult);
}
*/

static bool multistrategysort(GtBentsedgresources *bsr,
                              GtUword subbucketleft,
                              GtUword width,
                              GtUword depth,
                              GtUword sortmaxdepth)
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
                            bsr->processunsortedsuffixrange,
                            bsr->processunsortedsuffixrangeinfo);
    bsr->countbltriesort++;
    return true;
  }
  return false;
}

static void subsort_bentleysedgewick(GtBentsedgresources *bsr,
                                     GtUword subbucketleft,
                                     GtUword width,
                                     GtUword depth)
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
                                 (GtUword) bsr->sortmaxdepth);
      if (bsr->sortmaxdepth > 0)
      {
        gt_shortreadsort_sssp_add_unsorted(
                              bsr->srsw,
                              gt_suffixsortspace_bucketleftidx_get(bsr->sssp),
                              subbucketleft,
                              width,
                              (GtUword) bsr->sortmaxdepth,
                              bsr->sssp,
                              bsr->processunsortedsuffixrange,
                              bsr->processunsortedsuffixrangeinfo);
      }
      bsr->countshortreadsort++;
      return;
    }
    if (bsr->sortmaxdepth > 0 && depth >= (GtUword) bsr->sortmaxdepth)
    {
      if (bsr->processunsortedsuffixrange != NULL)
      {
        bsr->processunsortedsuffixrange(
                         bsr->processunsortedsuffixrangeinfo,
                         bsr->sssp,
                         gt_suffixsortspace_bucketleftidx_get(bsr->sssp) +
                         subbucketleft,
                         width,depth);
      }
      return;
    }
    if (multistrategysort(bsr,subbucketleft,width,depth,
                          (GtUword) bsr->sortmaxdepth))
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
                             GtUword subbucketleft,
                             GtUword width,
                             const GtSfxcmp *pivotcmpbits,
                             GtUword pivotidx,
                             GtUword depth)
{
  int cmp;
  unsigned int maxsmallerwithlcp = 0, maxlargerwithlcp = 0;
  GtCommonunits commonunits;
  GtEndofTwobitencoding etbecurrent;
  GtUword idx, smaller = 0, larger = 0,
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
  for (idx = 1UL; idx <= (GtUword) maxsmallerwithlcp; idx++)
  {
    bsr->leftlcpdist[idx] += bsr->leftlcpdist[idx-1];
  }
  for (idx = 1UL; idx <= (GtUword) maxlargerwithlcp; idx++)
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
  for (idx = 0; idx <= (GtUword) maxsmallerwithlcp; idx++)
  {
    if (idx < (GtUword) maxsmallerwithlcp)
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
  for (idx = 0; idx <= (GtUword) maxlargerwithlcp; idx++)
  {
    if (idx < (GtUword) maxlargerwithlcp)
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
                              GtUword subbucketleft1,
                              GtUword subbucketleft2,
                              GtUword width)
{
  GtUword idx, tmp;

  for (idx = 0; idx < width; idx++)
  {
    tmp = gt_suffixsortspace_get(sssp,subbucketleft1,idx);
    gt_suffixsortspace_set(sssp,subbucketleft1,idx,
                           gt_suffixsortspace_get(sssp,subbucketleft2,idx));
    gt_suffixsortspace_set(sssp,subbucketleft2,idx,tmp);
  }
}

static void gt_sort_bentleysedgewick(GtBentsedgresources *bsr,
                                     GtUword bucketleftidx,
                                     GtUword width,
                                     GtUword depth)
{
  bsr->mkvauxstack.nextfreeGtMKVstack = 0;
  gt_suffixsortspace_bucketrange_set(bsr->sssp,bucketleftidx,width);
  subsort_bentleysedgewick(bsr, 0, width, depth);
  while (bsr->mkvauxstack.nextfreeGtMKVstack > 0)
  {
    GtUword leftplusw, pa, pb, pc, pd, pm, bucketright,
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
    if ((wtmp = (GtUword) (pd-pc)) > 0)
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
  gt_suffixsortspace_bucketrange_reset(bsr->sssp);
}

static GtBentsedgresources *bentsedgresources_new(
                                   GtSuffixsortspace *suffixsortspace,
                                   const GtEncseq *encseq,
                                   GtReadmode readmode,
                                   unsigned int prefixlength,
                                   const GtBcktab *bcktab,
                                   unsigned int sortmaxdepth,
                                   const Sfxstrategy *sfxstrategy,
                                   bool withlcps)
{
  GtBentsedgresources *bsr = (GtBentsedgresources *) gt_malloc(sizeof *bsr);

  bsr->encseq = encseq;
  bsr->readmode = readmode;
  bsr->totallength = gt_encseq_total_length(encseq);
  bsr->sfxstrategy = sfxstrategy;
  bsr->sssp = suffixsortspace;
  bsr->rsi = NULL;
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
      bsr->srs_maxremain = sortmaxdepth - (GtUword) prefixlength;
    } else
    {
      bsr->srs_maxremain = 0;
    }
  } else
  {
    GtUword longestnsp = gt_encseq_lengthoflongestnonspecial(encseq);
    if (longestnsp > (GtUword) prefixlength)
    {
      bsr->srs_maxremain = longestnsp - prefixlength;
    } else
    {
      bsr->srs_maxremain = 0;
    }
  }
  bsr->shortreadsort_maxwidth
    = gt_shortreadsort_maxwidth(false,
                                bsr->srs_maxremain,
                                gt_size_of_sort_workspace (sfxstrategy));
  GT_INITARRAY(&bsr->mkvauxstack,GtMKVstack);
  bsr->countingsortinfo = NULL;
  bsr->medianinfospace = NULL;
  bsr->blindtrie = NULL;
  bsr->equalwithprevious = NULL;
  bsr->srsw = NULL;
  if (!sfxstrategy->cmpcharbychar)
  {
    int idx;
    const bool withmediumsizelcps = (sortmaxdepth > 0 && !withlcps)
                                    ? true : false;
    if (sortmaxdepth == 0 ||
        gt_encseq_has_twobitencoding_stoppos_support(encseq))
    {
      bsr->srsw = gt_shortreadsort_new(0,bsr->srs_maxremain,readmode,false,
                                       withmediumsizelcps);
    }
    for (idx = 0; idx < GT_UNITSIN2BITENC; idx++)
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
    GtUword idx;

    bsr->equalwithprevious = gt_malloc(sizeof (*bsr->equalwithprevious) *
                                       bsr->sfxstrategy->maxinsertionsort);
    for (idx=0; idx < bsr->sfxstrategy->maxinsertionsort; idx++)
    {
      bsr->equalwithprevious[idx] = false;
    }
  }
  bsr->sortmaxdepth = sortmaxdepth;
  bsr->processunsortedsuffixrange = NULL;
  bsr->processunsortedsuffixrangeinfo = NULL;
  bsr->countinsertionsort = 0;
  bsr->counttqsort = 0;
  bsr->countcountingsort = 0;
  bsr->countbltriesort = 0;
  bsr->countshortreadsort = 0;
  bsr->countradixsort = 0;
  if (bcktab != NULL &&
      sfxstrategy->withradixsort &&
      gt_encseq_accesstype_get(encseq) == GT_ACCESS_TYPE_EQUALLENGTH &&
      readmode == GT_READMODE_FORWARD)
  {
    bsr->rsi = gt_radixsort_str_new(bsr->twobitencoding,
                                    gt_encseq_is_mirrored(encseq)
                                      ? GT_DIV2(bsr->totallength - 1)
                                      : bsr->totallength,
                                    1 + gt_encseq_equallength(encseq),
                                    gt_bcktab_nonspecialsmaxsize(bcktab));
    bsr->radixsortmaxwidth = gt_radixsort_str_maxwidth(bsr->rsi);
  }
  return bsr;
}

static void gt_bentsedgresources_addlcpinfo(GtBentsedgresources *bsr,
                                            GtOutlcpinfo *outlcpinfo,
                                            const GtBcktab *bcktab)
{
  if (outlcpinfo != NULL)
  {
    bsr->tableoflcpvalues = gt_Outlcpinfo_resizereservoir(outlcpinfo,bcktab);
    if (bsr->srsw != NULL)
    {
      gt_shortreadsort_assigntableoflcpvalues(bsr->srsw,bsr->tableoflcpvalues);
    }
  }
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
  gt_logger_log(logger,"countinsertionsort="GT_WU"",bsr->countinsertionsort);
  gt_logger_log(logger,"countbltriesort="GT_WU"",bsr->countbltriesort);
  gt_logger_log(logger,"countcountingsort="GT_WU"",bsr->countcountingsort);
  gt_logger_log(logger,"countshortreadsort="GT_WU"",bsr->countshortreadsort);
  gt_logger_log(logger,"countradixsort="GT_WU"",bsr->countradixsort);
  gt_logger_log(logger,"counttqsort="GT_WU"",bsr->counttqsort);
  gt_free(bsr);
}

/*
  The following function is called in sfx-suffixer.c and sorts all buckets by
  different suffix comparison methods without the help of other sorting
  information. GtSuffixsortspace contains the sortspace which is accessed
  by some negative offset.
*/

void gt_sortallbuckets(GtSuffixsortspace *suffixsortspace,
                       GtUword numberofsuffixes,
                       GtBucketspec2 *bucketspec2,
                       const GtEncseq *encseq,
                       GtReadmode readmode,
                       GtCodetype mincode,
                       GtCodetype maxcode,
                       const GtBcktab *bcktab,
                       unsigned int numofchars,
                       unsigned int prefixlength,
                       GtOutlcpinfo *outlcpinfo,
                       unsigned int sortmaxdepth,
                       const Sfxstrategy *sfxstrategy,
                       GtProcessunsortedsuffixrange processunsortedsuffixrange,
                       void *processunsortedsuffixrangeinfo,
                       GtUint64 *bucketiterstep,
                       GtLogger *logger)
{
  GtCodetype code;
  unsigned int rightchar = (unsigned int) (mincode % numofchars);
  GtBucketspecification bucketspec;
  GtBentsedgresources *bsr = bentsedgresources_new(suffixsortspace,
                                                   encseq,
                                                   readmode,
                                                   prefixlength,
                                                   bcktab,
                                                   sortmaxdepth,
                                                   sfxstrategy,
                                                   outlcpinfo != NULL ? true
                                                                      : false);
  gt_bentsedgresources_addlcpinfo(bsr,outlcpinfo,bcktab);
  bsr->processunsortedsuffixrangeinfo = processunsortedsuffixrangeinfo;
  bsr->processunsortedsuffixrange = processunsortedsuffixrange;
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
        gt_sort_bentleysedgewick(bsr,bucketspec.left,
                                 bucketspec.nonspecialsinbucket,
                                 (GtUword) prefixlength);
      }
      gt_Outlcpinfo_nonspecialsbucket(outlcpinfo,
                                      prefixlength,
                                      bsr->sssp,
                                      bsr->tableoflcpvalues,
                                      &bucketspec,
                                      code);
    }
    gt_Outlcpinfo_postbucket(outlcpinfo,
                             prefixlength,
                             bsr->sssp,
                             bcktab,
                             &bucketspec,
                             code);
  }
  bentsedgresources_delete(bsr, logger);
}

void gt_sortallsuffixesfromstart(GtSuffixsortspace *suffixsortspace,
                                 GtUword numberofsuffixes,
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
  if (numberofsuffixes > 1UL)
  {
    GtBentsedgresources *bsr = bentsedgresources_new(suffixsortspace,
                                                     encseq,
                                                     readmode,
                                                     0,
                                                     NULL,
                                                     sortmaxdepth,
                                                     sfxstrategy,
                                                     outlcpinfo != NULL
                                                       ? true : false);
    gt_bentsedgresources_addlcpinfo(bsr,outlcpinfo,NULL);
    bsr->processunsortedsuffixrange = processunsortedsuffixrange;
    bsr->processunsortedsuffixrangeinfo = processunsortedsuffixrangeinfo;
    gt_sort_bentleysedgewick(bsr,0,numberofsuffixes,0);
    bentsedgresources_delete(bsr, logger);
  }
}

#ifdef GT_THREADS_ENABLED
#ifdef GT_THREADS_PARTITION
typedef struct
{
  unsigned int numofchars,
               prefixlength;
  const GtBcktab *bcktab;

  GtCodetype mincode,
             maxcode;
  GtUword totalwidth;
  GtBentsedgresources *bsr;
  unsigned int thread_num;
  GtThread *thread;
} GtBentsedg_partition_thread_info;

static void *gt_bentsedg_partition_thread_caller(void *data)
{
  GtCodetype code;
  GtBentsedg_partition_thread_info *thinfo
    = (GtBentsedg_partition_thread_info *) data;
  unsigned int rightchar;

  rightchar = (unsigned int) (thinfo->mincode % thinfo->numofchars);
  for (code = thinfo->mincode; code <= thinfo->maxcode; code++)
  {
    GtBucketspecification bucketspec;
    rightchar = gt_bcktab_calcboundsparts(&bucketspec,
                                          thinfo->bcktab,
                                          code,
                                          thinfo->maxcode,
                                          thinfo->totalwidth,
                                          rightchar);
    if (bucketspec.nonspecialsinbucket > 1UL)
    {
      gt_sort_bentleysedgewick(thinfo->bsr,bucketspec.left,
                               bucketspec.nonspecialsinbucket,
                               (GtUword) thinfo->prefixlength);
    }
  }
  return NULL;
}

void gt_threaded_partition_sortallbuckets(GtSuffixsortspace *suffixsortspace,
                       const GtSuftabparts *partition_for_threads,
                       const GtEncseq *encseq,
                       GtReadmode readmode,
                       const GtBcktab *bcktab,
                       unsigned int numofchars,
                       unsigned int prefixlength,
                       unsigned int sortmaxdepth,
                       const Sfxstrategy *sfxstrategy,
                       GtProcessunsortedsuffixrange processunsortedsuffixrange,
                       void *processunsortedsuffixrangeinfo,
                       GtLogger *logger)
{
  unsigned int tp, thread_parts;
  bool haserr = false;
  GtBentsedg_partition_thread_info *th_tab;
  GtSuffixsortspace **sssp_tab;

  gt_assert(partition_for_threads != NULL);
  thread_parts = gt_suftabparts_numofparts(partition_for_threads);
  gt_assert(thread_parts > 1U);
  th_tab = gt_malloc(sizeof *th_tab * thread_parts);
  sssp_tab = gt_malloc(sizeof *sssp_tab * thread_parts);
  for (tp = 0; !haserr && tp < thread_parts; tp++)
  {
    th_tab[tp].thread_num = tp;
    th_tab[tp].numofchars = numofchars;
    th_tab[tp].prefixlength = prefixlength;
    th_tab[tp].bcktab = bcktab;
    th_tab[tp].mincode = gt_suftabparts_minindex(tp,partition_for_threads);
    th_tab[tp].maxcode = gt_suftabparts_maxindex(tp,partition_for_threads);
    th_tab[tp].totalwidth = gt_suftabparts_sumofwidth(tp,partition_for_threads);
    if (tp == 0)
    {
      sssp_tab[tp] = suffixsortspace;
    } else
    {
      sssp_tab[tp] = gt_suffixsortspace_clone(suffixsortspace,tp,logger);
    }
    th_tab[tp].bsr = bentsedgresources_new(sssp_tab[tp],
                                           encseq,
                                           readmode,
                                           prefixlength,
                                           bcktab,
                                           sortmaxdepth,
                                           sfxstrategy,
                                           false);
    th_tab[tp].bsr->processunsortedsuffixrange
      = processunsortedsuffixrange;
    th_tab[tp].bsr->processunsortedsuffixrangeinfo
      = processunsortedsuffixrangeinfo;
    th_tab[tp].thread
      = gt_thread_new (gt_bentsedg_partition_thread_caller,th_tab + tp,NULL);
    if (th_tab[tp].thread == NULL)
    {
      haserr = true;
    }
  }
  for (tp = 0; tp < thread_parts; tp++)
  {
    if (!haserr)
    {
      gt_thread_join(th_tab[tp].thread);
    }
    gt_thread_delete(th_tab[tp].thread);
  }
  for (tp = 0; tp < thread_parts; tp++)
  {
    bentsedgresources_delete(th_tab[tp].bsr, logger);
  }
  gt_suffixsortspace_delete_cloned(sssp_tab,thread_parts);
  gt_free(sssp_tab);
  gt_free(th_tab);
  gt_assert(!haserr);
}
#else

typedef struct
{
  GtUword totalwidth, bucketnumber;
  const GtBcktab *bcktab;
  GtCodetype code, mincode, maxcode;
  unsigned int rightchar;
  GtMutex *mutex;
} GtBentsedgIterator;

static GtBentsedgIterator *gt_BentsedgIterator_new(GtCodetype mincode,
                                                   GtCodetype maxcode,
                                                   GtUword totalwidth,
                                                   unsigned int numofchars,
                                                   const GtBcktab *bcktab)
{
  GtBentsedgIterator *bentsedg_iterator = gt_malloc(sizeof *bentsedg_iterator);

  bentsedg_iterator->code = mincode;
  bentsedg_iterator->mincode = mincode;
  bentsedg_iterator->maxcode = maxcode;
  bentsedg_iterator->totalwidth = totalwidth;
  bentsedg_iterator->rightchar = (unsigned int) (mincode % numofchars);
  bentsedg_iterator->bcktab = bcktab;
  bentsedg_iterator->bucketnumber = 0;
  bentsedg_iterator->mutex = gt_mutex_new();
  return bentsedg_iterator;
}

static bool gt_BentsedgIterator_next(GtBucketspecification *bucketspec,
                                     GtBentsedgIterator *bs_it)
{
  if (bs_it->code <= bs_it->maxcode)
  {
    bs_it->rightchar = gt_bcktab_calcboundsparts(bucketspec,
                                                 bs_it->bcktab,
                                                 bs_it->code,
                                                 bs_it->maxcode,
                                                 bs_it->totalwidth,
                                                 bs_it->rightchar);
    bs_it->code++;
    return true;
  } else
  {
    return false;
  }
}
static void gt_BentsedgIterator_delete(GtBentsedgIterator *bs_it)
{
  if (bs_it != NULL)
  {
    gt_assert(bs_it->bucketnumber == bs_it->maxcode - bs_it->mincode + 1);
    gt_mutex_delete(bs_it->mutex);
    gt_free(bs_it);
  }
}

typedef struct
{
  GtUword bucketnumber;
} GtPriorecord;

typedef struct
{
  GtUword capacity, numofelements;
  GtPriorecord *elements;
} GtBentsedgePrioqueue;

static GtBentsedgePrioqueue *gt_bs_priority_queue_new(GtUword numofelements)
{
  GtBentsedgePrioqueue *pq = gt_malloc(sizeof *pq);

  pq->elements = gt_malloc(sizeof (*pq->elements) * numofelements);
  pq->capacity = numofelements;
  pq->numofelements = 0;
  return pq;
}

static bool gt_bs_priority_queue_is_empty(const GtBentsedgePrioqueue *pq)
{
  gt_assert(pq != NULL);
  return pq->numofelements == 0 ? true : false;
}

static void gt_bs_priority_queue_add(GtBentsedgePrioqueue *pq,
                                     GtPriorecord priorecord)
{
  GtPriorecord *ptr;

  gt_assert(pq != NULL);
  if (pq->numofelements >= pq->capacity)
  {
    pq->capacity *= 2;
    pq->elements = gt_realloc(pq->elements,sizeof *pq->elements * pq->capacity);
  }
  /* store elements in reverse order, i.e.\ with the minimum element
     at the last index; shift elements to the right until an element larger
     or equal than the key is found. */
  for (ptr = pq->elements + pq->numofelements;
       ptr > pq->elements && (ptr-1)->bucketnumber < priorecord.bucketnumber;
       ptr--)
  {
    *ptr = *(ptr-1);
  }
  pq->numofelements++;
  gt_assert(ptr >= pq->elements && ptr < pq->elements + pq->numofelements);
  *ptr = priorecord;
}

static void gt_bs_priority_queue_delete_min(GtBentsedgePrioqueue *pq)
{
  gt_assert(pq != NULL && !gt_bs_priority_queue_is_empty(pq));
  pq->numofelements--;
}

const GtUword gt_bs_priority_queue_min_bucketnumber(
                                              const GtBentsedgePrioqueue *pq)
{
  gt_assert(pq != NULL && !gt_bs_priority_queue_is_empty(pq));
  return pq->elements[pq->numofelements-1].bucketnumber;
}

static void gt_bs_priority_queue_delete(GtBentsedgePrioqueue *pq)
{
  if (pq != NULL)
  {
    gt_free(pq->elements);
    gt_free(pq);
  }
}

typedef struct
{
  GtBentsedgePrioqueue *queue;
  GtUword nextrequest;
  GtMutex *mutex;
} GtBentsedgSynchronizer;

static GtBentsedgSynchronizer *gt_bendsedgSynchronizer_new(void)
{
  GtBentsedgSynchronizer *bs_sync = gt_malloc(sizeof *bs_sync);

  bs_sync->queue = gt_bs_priority_queue_new(gt_jobs);
  bs_sync->nextrequest = 0;
  bs_sync->mutex = gt_mutex_new();
  return bs_sync;
}

static void gt_bendsedgSynchronizer_process(GtBentsedgSynchronizer *bs_sync,
                                            GtUword bucketnumber)
{
  gt_assert(bs_sync != NULL);
  if (bs_sync->nextrequest == bucketnumber)
  {
    bs_sync->nextrequest++;
    while (!gt_bs_priority_queue_is_empty(bs_sync->queue) &&
           gt_bs_priority_queue_min_bucketnumber(bs_sync->queue) ==
           bs_sync->nextrequest)
    {
      gt_bs_priority_queue_delete_min(bs_sync->queue);
      bs_sync->nextrequest++;
    }
  } else
  {
    GtPriorecord priorecord;
    gt_assert(bs_sync->nextrequest < bucketnumber);
    priorecord.bucketnumber = bucketnumber;
    gt_bs_priority_queue_add(bs_sync->queue,priorecord);
  }
}

static void gt_bendsedgSynchronizer_delete(GtBentsedgSynchronizer *bs_sync)
{
  if (bs_sync != NULL)
  {
    gt_assert(gt_bs_priority_queue_is_empty(bs_sync->queue));
    gt_mutex_delete(bs_sync->mutex);
    gt_bs_priority_queue_delete(bs_sync->queue);
    gt_free(bs_sync);
  }
}

typedef struct
{
  GtBentsedgresources *bsr;
  unsigned int prefixlength, thread_num;
  GtBentsedgIterator *bs_it; /* shared, _next-function needs a mutex */
  GtBentsedgSynchronizer *bs_sync; /* shared _process-function needs a mutex */
  GtThread *thread;
} GtBentsedg_stream_thread_info;

static void *gt_bentsedg_stream_thread_caller(void *data)
{
  GtBentsedg_stream_thread_info *thinfo
    = (GtBentsedg_stream_thread_info *) data;

  while (true)
  {
    GtBucketspecification bucketspec;
    GtUword bucketnumber;

    gt_mutex_lock(thinfo->bs_it->mutex);
    if (!gt_BentsedgIterator_next(&bucketspec,thinfo->bs_it))
    {
      gt_mutex_unlock(thinfo->bs_it->mutex);
      break;
    }
    bucketnumber = thinfo->bs_it->bucketnumber++;
    gt_mutex_unlock(thinfo->bs_it->mutex);
    if (bucketspec.nonspecialsinbucket > 1UL)
    {
      gt_sort_bentleysedgewick(thinfo->bsr,bucketspec.left,
                               bucketspec.nonspecialsinbucket,
                               (GtUword) thinfo->prefixlength);
    }
    gt_mutex_lock(thinfo->bs_sync->mutex);
    gt_bendsedgSynchronizer_process(thinfo->bs_sync,bucketnumber);
    gt_mutex_unlock(thinfo->bs_sync->mutex);
  }
  return NULL;
}

void gt_threaded_stream_sortallbuckets(GtSuffixsortspace *suffixsortspace,
                       const GtEncseq *encseq,
                       GtReadmode readmode,
                       const GtBcktab *bcktab,
                       GtCodetype mincode,
                       GtCodetype maxcode,
                       GtUword sumofwidth,
                       unsigned int numofchars,
                       unsigned int prefixlength,
                       unsigned int sortmaxdepth,
                       const Sfxstrategy *sfxstrategy,
                       GtProcessunsortedsuffixrange processunsortedsuffixrange,
                       void *processunsortedsuffixrangeinfo,
                       GtLogger *logger)
{
  GtBentsedgIterator *bs_it;
  GtBentsedgSynchronizer *bs_sync;
  unsigned int tp;
  bool haserr = false;
  GtBentsedg_stream_thread_info *th_tab;
  GtSuffixsortspace **sssp_tab;

  gt_assert(gt_jobs > 1U);
  th_tab = gt_malloc(sizeof *th_tab * gt_jobs);
  sssp_tab = gt_malloc(sizeof *sssp_tab * gt_jobs);
  bs_it = gt_BentsedgIterator_new(mincode,maxcode,sumofwidth,numofchars,bcktab);
  bs_sync = gt_bendsedgSynchronizer_new();
  for (tp = 0; !haserr && tp < gt_jobs; tp++)
  {
    th_tab[tp].thread_num = tp;
    th_tab[tp].prefixlength = prefixlength;
    if (tp == 0)
    {
      sssp_tab[tp] = suffixsortspace;
    } else
    {
      sssp_tab[tp] = gt_suffixsortspace_clone(suffixsortspace,tp,logger);
    }
    th_tab[tp].bsr = bentsedgresources_new(sssp_tab[tp],
                                           encseq,
                                           readmode,
                                           prefixlength,
                                           bcktab,
                                           sortmaxdepth,
                                           sfxstrategy,
                                           false);
    th_tab[tp].bsr->processunsortedsuffixrange
      = processunsortedsuffixrange;
    th_tab[tp].bsr->processunsortedsuffixrangeinfo
      = processunsortedsuffixrangeinfo;
    th_tab[tp].bs_it = bs_it;
    th_tab[tp].bs_sync = bs_sync;
    th_tab[tp].thread
      = gt_thread_new (gt_bentsedg_stream_thread_caller,th_tab + tp,NULL);
    if (th_tab[tp].thread == NULL)
    {
      haserr = true;
    }
  }
  for (tp = 0; tp < gt_jobs; tp++)
  {
    if (!haserr)
    {
      gt_thread_join(th_tab[tp].thread);
    }
    gt_thread_delete(th_tab[tp].thread);
  }
  for (tp = 0; tp < gt_jobs; tp++)
  {
    bentsedgresources_delete(th_tab[tp].bsr, logger);
  }
  gt_suffixsortspace_delete_cloned(sssp_tab,gt_jobs);
  gt_BentsedgIterator_delete(bs_it);
  gt_bendsedgSynchronizer_delete(bs_sync);
  gt_free(sssp_tab);
  gt_free(th_tab);
  gt_assert(!haserr);
}
#endif
#endif
