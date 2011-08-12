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
#include "core/stack-inlined.h"
#include "bcktab.h"
#include "kmer2string.h"
#include "sfx-bltrie.h"
#include "sfx-copysort.h"
#include "sfx-bentsedg.h"
#include "sfx-lcpvalues.h"
#include "sfx-suftaborder.h"
#include "sfx-suffixgetset.h"

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

#define UPDATEMAXLCP(MAXVAL,LCP)\
        if ((MAXVAL) < (LCP))\
        {\
          MAXVAL = LCP;\
        }

#define UPDATELCP(MINVAL,MAXVAL,LCP)\
        if ((MINVAL) > (LCP))\
        {\
          MINVAL = LCP;\
        }\
        UPDATEMAXLCP(MAXVAL,LCP)

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
            TMPVAR.position = pos;\
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

#ifdef _LP64
#define GT_NUMOFTBEVALUEFOR100 3
#else
#define GT_NUMOFTBEVALUEFOR100 6
#endif

typedef struct
{
  unsigned long suffix;
  GtTwobitencoding tbe[GT_NUMOFTBEVALUEFOR100];
  unsigned int unitsnotspecial;
} GtShortreadsort;

GT_DECLAREARRAYSTRUCT(GtMKVstack);

typedef struct
{
  const GtEncseq *encseq;
  GtEncseqReader *esr1, /* XXX be carefull with threads */
                 *esr2;
  GtReadmode readmode;
  bool fwd, complement;
  unsigned long totallength;
  GtArrayGtMKVstack mkvauxstack; /* XXX be carefull with treads */
  GtLcpvalues *tableoflcpvalues;
  GtMedianinfo *medianinfospace;
  GtCountingsortinfo *countingsortinfo;
  GtShortreadsort *shortreadsortinfo;
  uint16_t *shortreadsortrefs;
  const Sfxstrategy *sfxstrategy;
  unsigned int sortmaxdepth,
               prefixlength;
  Blindtrie *blindtrie;
  unsigned long leftlcpdist[GT_UNITSIN2BITENC],
                rightlcpdist[GT_UNITSIN2BITENC];
  GtSuffixsortspace *sssp;
  GtProcessunsortedsuffixrange processunsortedsuffixrange;
  void *processunsortedsuffixrangeinfo;
  bool *equalwithprevious;
  unsigned long countinsertionsort,
                counttqsort,
                countshortreadsort,
                countcountingsort,
                countbltriesort;
  unsigned long lcplen;
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
  unsigned long sval1, sval2, pm, pl, startpos1, startpos2, temp,
                lcpindex, lcplen = 0;
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
        printf("gt_encseq_compare_viatwobitencoding[%lu,%lu] "
               "at offset %lu\n",
                sval1,sval2,offset);
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
        lcpindex = subbucketleft+pl;
        if (pl < pm && retval > 0)
        {
          lcptab_update(bsr->tableoflcpvalues,lcpindex+1,
                        lcpsubtab_getvalue(bsr->tableoflcpvalues,lcpindex));
        }
        lcptab_update(bsr->tableoflcpvalues,lcpindex,lcplen);
      }
      if (retval < 0)
      {
        break;
      }
      BS_SWAPARRAY(temp,subbucketleft,pl,pl-1);
    }
  }
}

static void bs_insertionsortmaxdepth(GtBentsedgresources *bsr,
                                     unsigned long subbucketleft,
                                     unsigned long width,
                                     unsigned long offset,
                                     unsigned long maxdepth)
{
  unsigned long sval1, sval2, pm, pl, startpos1, startpos2, temp,
                lcpindex, lcplen = 0, idx = 0;
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

        endpos1 = sval1+maxdepth;
        if (endpos1 > bsr->totallength)
        {
          endpos1 = bsr->totallength;
        }
        endpos2 = sval2+maxdepth;
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
          if (lcplen == maxdepth)
          {
            retval = 0;
            break;
          }
          gt_assert(lcplen < maxdepth);
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
        gt_assert(offset < maxdepth);
        retval = gt_encseq_compare_viatwobitencoding(&commonunits,
                                                     bsr->encseq,
                                                     bsr->readmode,
                                                     bsr->esr1,
                                                     bsr->esr2,
                                                     sval1,
                                                     sval2,
                                                     offset,
                                                     maxdepth);
        lcplen = commonunits.finaldepth;
        gt_assert(lcplen <= maxdepth);
        if (lcplen == maxdepth)
        {
          gt_assert(retval == 0);
        }
      }
#ifdef SKDEBUG
      printf("cmp %lu and %lu: retval = %d, lcplen = %lu\n",
             sval1, sval2, retval, (unsigned long) lcplen);
#endif
      if (bsr->tableoflcpvalues != NULL && retval != 0)
      {
        lcpindex = subbucketleft + pl;
        if (pl < pm && retval > 0)
        {
          lcptab_update(bsr->tableoflcpvalues,lcpindex+1,
                        lcpsubtab_getvalue(bsr->tableoflcpvalues,lcpindex));
        }
        lcptab_update(bsr->tableoflcpvalues,lcpindex,lcplen);
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
      BS_SWAPARRAY(temp,subbucketleft,pl,pl-1);
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
          bsr->processunsortedsuffixrange(
                              bsr->processunsortedsuffixrangeinfo,
                              bucketleftidx + subbucketleft + idx - 1
                                            - equalsrangewidth,
                              equalsrangewidth + 1, maxdepth);
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
      bsr->processunsortedsuffixrange(
                           bsr->processunsortedsuffixrangeinfo,
                           bucketleftidx + subbucketleft + width - 1
                                         - equalsrangewidth,
                           equalsrangewidth + 1, maxdepth);
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
                              unsigned long maxdepth)
{
  gt_assert(width > 1UL);
  if (width <= bsr->sfxstrategy->maxinsertionsort)
  {
    if (maxdepth == 0)
    {
      bs_insertionsort(bsr,subbucketleft,width,depth);
    } else
    {
      bs_insertionsortmaxdepth(bsr,subbucketleft,width,depth,maxdepth);
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
                            maxdepth,
                            bsr->processunsortedsuffixrangeinfo,
                            bsr->processunsortedsuffixrange);
    bsr->countbltriesort++;
    return true;
  }
  return false;
}

static int compareshortreadsortinfo(const GtShortreadsort *aq,
                                    const GtShortreadsort *bq,
                                    GtBentsedgresources *bsr)
{
  int idx, retval;
  unsigned int maxprefix;
  GtCommonunits commonunits;

  for (idx=0, maxprefix = (unsigned int) GT_UNITSIN2BITENC;
       idx<GT_NUMOFTBEVALUEFOR100;
       idx++, maxprefix+=(unsigned int) GT_UNITSIN2BITENC)
  {
    if (aq->unitsnotspecial >= maxprefix &&
        bq->unitsnotspecial >= maxprefix)
    {
      if (aq->tbe[idx] != bq->tbe[idx])
      {
        retval = gt_encseq_compare_pairof_different_twobitencodings(
                                                     bsr->fwd,
                                                     bsr->complement,
                                                     &commonunits,
                                                     aq->tbe[idx],
                                                     bq->tbe[idx]);
        bsr->lcplen = (unsigned long)
                      (maxprefix - GT_UNITSIN2BITENC + commonunits.common);
        return retval;
      }
    } else
    {
      GtEndofTwobitencoding tbe_a, tbe_b;

      tbe_a.tbe = aq->tbe[idx];
      tbe_b.tbe = bq->tbe[idx];
      tbe_a.position = aq->suffix;
      tbe_b.position = bq->suffix;
      tbe_a.unitsnotspecial
        = aq->unitsnotspecial >= maxprefix
           ? maxprefix
           : aq->unitsnotspecial + GT_UNITSIN2BITENC - maxprefix;
      tbe_b.unitsnotspecial
        = bq->unitsnotspecial >= maxprefix
           ? maxprefix
           : bq->unitsnotspecial + GT_UNITSIN2BITENC - maxprefix;
      retval = gt_encseq_compare_pairof_twobitencodings(bsr->fwd,
                                                        bsr->complement,
                                                        &commonunits,
                                                        &tbe_a,
                                                        &tbe_b);
      bsr->lcplen = (unsigned long) (maxprefix - GT_UNITSIN2BITENC +
                                     commonunits.common);
      return retval;
    }
  }
  gt_assert(false);
  return 0;
}

#ifdef  QSORTNAME
#undef  QSORTNAME
#endif

#define QSORTNAME(NAME) shortread_##NAME

#define shortread_ARRAY_GET(ARR,IDX)\
        (unsigned long) data->shortreadsortrefs[IDX]

#define shortread_ARRAY_SET(ARR,IDX,VALUE)\
        data->shortreadsortrefs[IDX] = (uint16_t) VALUE

typedef GtBentsedgresources * QSORTNAME(Datatype);

static int QSORTNAME(qsortcmparr) (unsigned long a,
                                   unsigned long b,
                                   const QSORTNAME(Datatype) data)
{
  return compareshortreadsortinfo(
                      data->shortreadsortinfo + QSORTNAME(ARRAY_GET)(NULL,a),
                      data->shortreadsortinfo + QSORTNAME(ARRAY_GET)(NULL,b),
                      data);
}

typedef unsigned long QSORTNAME(Sorttype);

/*
 * Qsort routine from Bentley & McIlroy's ``Engineering a Sort Function''.
 */

#ifndef GT_QSORT_ARR_SWAP
#define GT_QSORT_ARR_SWAP(ARR,A,B)\
        if ((A) != (B))\
        {\
          tmp = QSORTNAME(ARRAY_GET)(ARR,A);\
          QSORTNAME(ARRAY_SET)(ARR,A,QSORTNAME(ARRAY_GET)(ARR,B));\
          QSORTNAME(ARRAY_SET)(ARR,B,tmp);\
        }
#endif

#ifndef GT_QSORT_ARR_VECSWAP
#define GT_QSORT_ARR_VECSWAP(ARR,A,B,N)\
        aidx = A;\
        bidx = B;\
        while ((N)-- > 0)\
        {\
          tmp = QSORTNAME(ARRAY_GET)(ARR,aidx);\
          QSORTNAME(ARRAY_SET)(ARR,aidx,QSORTNAME(ARRAY_GET)(ARR,bidx));\
          QSORTNAME(ARRAY_SET)(ARR,bidx,tmp);\
          aidx++;\
          bidx++;\
        }
#endif

static inline unsigned long QSORTNAME(gt_inlined_qsort_arr_r_med3)
                     (unsigned long a, unsigned long b, unsigned long c,
                      QSORTNAME(Datatype) data)
{
  return QSORTNAME(qsortcmparr) (a, b, data) < 0
           ? (QSORTNAME(qsortcmparr) (b, c, data) < 0
                ? b
                : (QSORTNAME(qsortcmparr) (a, c, data) < 0
                     ? c : a))
           : (QSORTNAME(qsortcmparr) (b, c, data) > 0
                ? b
                : (QSORTNAME(qsortcmparr) (a, c, data) < 0
                     ? a
                     : c));
}

#ifndef GT_STACK_INTERVALARRAYTOBESORTED_DEFINED
typedef struct
{
  unsigned long startindex,
                len;
} Intervalarrtobesorted;

GT_STACK_DECLARESTRUCT(Intervalarrtobesorted,32UL);
#define GT_STACK_INTERVALARRAYTOBESORTED_DEFINED
#endif

static void QSORTNAME(gt_inlinedarr_qsort_r) (
                                   unsigned long insertionsortthreshold,
                                   bool handlenotswapped,
                                   unsigned long len,
                                   QSORTNAME(Datatype) data,
                                   unsigned long depth,
                                   unsigned long subbucketleft)
{
  unsigned long tmp, pa, pb, pc, pd, pl, pm, pn, aidx, bidx, s,
                smallermaxlcp, greatermaxlcp;
  int r;
  bool swapped;
  GtStackIntervalarrtobesorted intervalstack;
  Intervalarrtobesorted current;

  GT_STACK_INIT(&intervalstack,32UL);
  current.startindex = 0;
  current.len = len;
  GT_STACK_PUSH(&intervalstack,current);
  if (insertionsortthreshold <= 2UL)
  {
    insertionsortthreshold = 6UL;
  }
  while (!GT_STACK_ISEMPTY(&intervalstack))
  {
    swapped = false;
    current = GT_STACK_POP(&intervalstack);
    if (current.len <= insertionsortthreshold)
    {
      for (pm = current.startindex + 1;
           pm < current.startindex + current.len; pm++)
      {
        for (pl = pm; pl > current.startindex; pl--)
        {
          r = QSORTNAME(qsortcmparr) (pl - 1, pl, data);
          if (data->tableoflcpvalues != NULL)
          {
            unsigned long lcpindex = subbucketleft + pl;
            if (pl < pm && r > 0)
            {
              lcptab_update(data->tableoflcpvalues,lcpindex+1,
                            lcpsubtab_getvalue(data->tableoflcpvalues,
                                               lcpindex));
            }
            lcptab_update(data->tableoflcpvalues,lcpindex,depth+data->lcplen);
          }
          if (r <= 0)
          {
            break;
          }
          GT_QSORT_ARR_SWAP (arr, pl, pl - 1);
        }
      }
      continue;
    }
    pm = current.startindex + GT_DIV2 (current.len);
    if (current.len > 7UL)
    {
      pl = current.startindex;
      pn = current.startindex + current.len - 1;
      if (current.len > 40UL)
      {
        s = GT_DIV8 (current.len);
        pl = QSORTNAME(gt_inlined_qsort_arr_r_med3) (pl, pl + s,
                                                     pl + GT_MULT2 (s), data);
        gt_assert(pm >= s);
        pm = QSORTNAME(gt_inlined_qsort_arr_r_med3) (pm - s, pm,
                                                     pm + s, data);
        gt_assert(pn >= GT_MULT2(s));
        pn = QSORTNAME(gt_inlined_qsort_arr_r_med3) (pn - GT_MULT2 (s),
                                                     pn - s, pn, data);
      }
      pm = QSORTNAME(gt_inlined_qsort_arr_r_med3) (pl, pm, pn, data);
    }
    GT_QSORT_ARR_SWAP (arr, current.startindex, pm);
    pa = pb = current.startindex + 1;
    pc = pd = current.startindex + current.len - 1;
    smallermaxlcp = greatermaxlcp = 0;
    while (1)
    {
      while (pb <= pc)
      {
        r = QSORTNAME(qsortcmparr) (pb, current.startindex, data);
        if (r > 0)
        {
          UPDATEMAXLCP(greatermaxlcp,data->lcplen);
          break;
        }
        if (r == 0)
        {
          swapped = true;
          GT_QSORT_ARR_SWAP (arr, pa, pb);
          pa++;
        } else
        {
          UPDATEMAXLCP(smallermaxlcp,data->lcplen);
        }
        pb++;
      }
      while (pb <= pc)
      {
        r = QSORTNAME(qsortcmparr) (pc, current.startindex, data);
        if (r < 0)
        {
          UPDATEMAXLCP(smallermaxlcp,data->lcplen);
          break;
        }
        if (r == 0)
        {
          swapped = true;
          GT_QSORT_ARR_SWAP (arr, pc, pd);
          gt_assert(pd > 0);
          pd--;
        } else
        {
          UPDATEMAXLCP(greatermaxlcp,data->lcplen);
        }
        gt_assert(pc > 0);
        pc--;
      }
      if (pb > pc)
      {
        break;
      }
      GT_QSORT_ARR_SWAP (arr, pb, pc);
      swapped = true;
      pb++;
      gt_assert(pc > 0);
      pc--;
    }
    /* The following switch is not explained in the above mentioned
       paper and therefore we ignore it. */
    if (handlenotswapped && !swapped)
    {                                  /* Switch to insertion sort */
      gt_assert(current.len <= 40UL);
      for (pm = current.startindex + 1;
           pm < current.startindex + current.len; pm++)
      {
        for (pl = pm; pl > current.startindex; pl--)
        {
          r = QSORTNAME(qsortcmparr) (pl - 1, pl, data);
          if (r <= 0)
          {
            break;
          }
          GT_QSORT_ARR_SWAP (arr, pl, pl - 1);
        }
      }
      continue;
    }
    pn = current.startindex + current.len;
    gt_assert(pa >= current.startindex);
    gt_assert(pb >= pa);
    s = MIN ((unsigned long) (pa - current.startindex),
             (unsigned long) (pb - pa));
    gt_assert(pb >= s);
    GT_QSORT_ARR_VECSWAP (arr, current.startindex, pb - s, s);
    gt_assert(pd >= pc);
    gt_assert(pn > pd);
    s = MIN ((unsigned long) (pd - pc), (unsigned long) (pn - pd - 1));
    gt_assert(pn > s);
    GT_QSORT_ARR_VECSWAP (arr, pb, pn - s, s);
    gt_assert(pb >= pa);
    if ((s = (unsigned long) (pb - pa)) > 0)
    {
      if (data->tableoflcpvalues != NULL)
      {
        /*
          left part has suffix with lcp up to length smallermaxlcp w.r.t.
          to the pivot. This lcp belongs to a suffix on the left
          which is at a minimum distance to the pivot and thus to an
          element in the final part of the left side.
        */
        lcptab_update(data->tableoflcpvalues,subbucketleft+current.startindex+s,
                      depth+smallermaxlcp);
      }
      if (s > 1UL)
      {
        current.len = s;
        GT_STACK_PUSH(&intervalstack,current);
      }
    }
    gt_assert(pd >= pc);
    if ((s = (unsigned long) (pd - pc)) > 0)
    {
      if (data->tableoflcpvalues != NULL)
      {
        /*
          right part has suffix with lcp up to length largermaxlcp w.r.t.
          to the pivot. This lcp belongs to a suffix on the right
          which is at a minimum distance to the pivot and thus to an
          element in the first part of the right side.
        */
        gt_assert(pn >= s);
        lcptab_update(data->tableoflcpvalues, subbucketleft + pn - s,
                      depth + greatermaxlcp);
      }
      if (s > 1UL)
      {
        gt_assert(pn >= s);
        current.startindex = pn - s;
        current.len = s;
        GT_STACK_PUSH(&intervalstack,current);
      }
    }
  }
  GT_STACK_DELETE(&intervalstack);
}

static void sarrshortreadsort(GtBentsedgresources *bsr,
                              unsigned long subbucketleft,
                              unsigned long width,
                              unsigned long depth)
{
  unsigned long idx, pos;
  GtSuffixsortspace_exportptr *exportptr;

  gt_assert(width <= (unsigned long) bsr->sfxstrategy->maxshortreadsort);
  exportptr = gt_suffixsortspace_exportptr(subbucketleft, bsr->sssp);
  if (exportptr->ulongtabsectionptr != NULL)
  {
    for (idx = 0; idx < width; idx++)
    {
      pos = exportptr->ulongtabsectionptr[idx];
      bsr->shortreadsortinfo[idx].suffix = pos;
      bsr->shortreadsortinfo[idx].unitsnotspecial
        = gt_encseq_extract2bitencvector(bsr->shortreadsortinfo[idx].tbe,
                                         GT_NUMOFTBEVALUEFOR100,
                                         bsr->encseq,
                                         bsr->esr1,
                                         bsr->readmode,
                                         pos+depth);
    }
  } else
  {
    for (idx = 0; idx < width; idx++)
    {
      pos = (unsigned long) exportptr->uinttabsectionptr[idx];
      bsr->shortreadsortinfo[idx].suffix = pos;
      bsr->shortreadsortinfo[idx].unitsnotspecial
        = gt_encseq_extract2bitencvector(bsr->shortreadsortinfo[idx].tbe,
                                         GT_NUMOFTBEVALUEFOR100,
                                         bsr->encseq,
                                         bsr->esr1,
                                         bsr->readmode,
                                         pos+depth);
    }
  }
  QSORTNAME(gt_inlinedarr_qsort_r) (6UL, false, width,bsr,depth,subbucketleft);
  if (exportptr->ulongtabsectionptr != NULL)
  {
    for (idx = 0; idx < width; idx++)
    {
      exportptr->ulongtabsectionptr[idx]
        = bsr->shortreadsortinfo[bsr->shortreadsortrefs[idx]].suffix;
      bsr->shortreadsortrefs[idx] = (uint16_t) idx;
      if (exportptr->ulongtabsectionptr[idx] == 0)
      {
        gt_suffixsortspace_updatelongest(bsr->sssp,idx);
      }
    }
  } else
  {
    for (idx = 0; idx < width; idx++)
    {
      exportptr->uinttabsectionptr[idx]
        = (uint32_t) bsr->shortreadsortinfo[bsr->shortreadsortrefs[idx]].suffix;
      bsr->shortreadsortrefs[idx] = (uint16_t) idx;
      if (exportptr->uinttabsectionptr[idx] == 0)
      {
        gt_suffixsortspace_updatelongest(bsr->sssp,idx);
      }
    }
  }
  gt_suffixsortspace_export_done(bsr->sssp);
  bsr->countshortreadsort++;
}

static bool allowforshortreadsort(const GtBentsedgresources *bsr)
{
  return (gt_encseq_accesstype_get(bsr->encseq) == GT_ACCESS_TYPE_EQUALLENGTH &&
          gt_encseq_equallength(bsr->encseq) == 100UL &&
          bsr->prefixlength >= 4U /* as GT_NUMOFTBEVALUEFOR100 *
                                        GT_UNITSIN2BITENC + 4 >= 100 */
          && bsr->fwd && !bsr->complement)
          /* currently !fwd or complement does not work = > exclude for now */
          ? true : false;
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
    /* generalize the following for cases with small maximal length */
    if (allowforshortreadsort(bsr) &&
        width <= (unsigned long) bsr->sfxstrategy->maxshortreadsort)
    {
      sarrshortreadsort(bsr,subbucketleft,width,depth);
      return;
    }
    if (bsr->sortmaxdepth > 0 && depth >= (unsigned long) bsr->sortmaxdepth)
    {
      bsr->processunsortedsuffixrange(bsr->processunsortedsuffixrangeinfo,
                         gt_suffixsortspace_bucketleftidx_get(bsr->sssp) +
                         subbucketleft,
                         width,depth);
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
      lcptab_update(bsr->tableoflcpvalues,subbucketleft+end,depth + idx);
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
      lcptab_update(bsr->tableoflcpvalues,subbucketleft + width - end,
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
      PTR2INT(pivotcmpbits,subbucketleft,0);
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
            UPDATELCP(greaterminlcp,greatermaxlcp,commonunits.common);
            break;
          }
          if (GtSfxcmpEQUAL(val,pivotcmpbits))
          {
            /* exchange equal element and element at index pa */
            BS_SWAPARRAY(temp, subbucketleft, pa, pb);
            pa++;
          } else /* smaller */
          {
            UPDATELCP(smallerminlcp,smallermaxlcp,commonunits.common);
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
            UPDATELCP(smallerminlcp,smallermaxlcp,commonunits.common);
            break;
          }
          if (GtSfxcmpEQUAL(val,pivotcmpbits))
          {
            /* exchange equal element and element at index pa */
            BS_SWAPARRAY(temp, subbucketleft, pc, pd);
            pd--;
          } else /* greater */
          {
            UPDATELCP(greaterminlcp,greatermaxlcp,commonunits.common);
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
    vectorswap(bsr->sssp, subbucketleft, subbucketleft+pb-wtmp, wtmp);
    gt_assert(pd >= pc);
    gt_assert(bucketright >= pd);
    wtmp = MIN(pd-pc, bucketright-pd);
    /* move w elements at the right to the middle */
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
        lcptab_update(bsr->tableoflcpvalues,subbucketleft + leftplusw,
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
        lcptab_update(bsr->tableoflcpvalues,
                      subbucketleft + bucketright - wtmp + 1,
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

static void initBentsedgresources(GtBentsedgresources *bsr,
                                  GtSuffixsortspace *suffixsortspace,
                                  const GtEncseq *encseq,
                                  GtReadmode readmode,
                                  unsigned int prefixlength,
                                  unsigned int sortmaxdepth,
                                  const Sfxstrategy *sfxstrategy)
{
  unsigned long idx;

  bsr->readmode = readmode;
  bsr->totallength = gt_encseq_total_length(encseq);
  bsr->sfxstrategy = sfxstrategy;
  bsr->sssp = suffixsortspace;
  gt_suffixsortspace_bucketleftidx_set(bsr->sssp,0);
  bsr->encseq = encseq;
  bsr->fwd = GT_ISDIRREVERSE(bsr->readmode) ? false : true;
  bsr->complement = GT_ISDIRCOMPLEMENT(bsr->readmode) ? true : false;
  bsr->tableoflcpvalues = NULL;
  bsr->prefixlength = prefixlength;
  for (idx = 0; idx < (unsigned long) GT_UNITSIN2BITENC; idx++)
  {
    bsr->leftlcpdist[idx] = bsr->rightlcpdist[idx] = 0;
  }
  bsr->esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  bsr->esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  GT_INITARRAY(&bsr->mkvauxstack,GtMKVstack);
  if (sfxstrategy->cmpcharbychar)
  {
    bsr->countingsortinfo = NULL;
    bsr->shortreadsortinfo = NULL;
    bsr->shortreadsortrefs = NULL;
    bsr->medianinfospace = NULL;
  } else
  {
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
    if (allowforshortreadsort(bsr))
    {
      bsr->shortreadsortinfo = gt_malloc(sizeof (*bsr->shortreadsortinfo) *
                                         (sfxstrategy->maxshortreadsort+1));
      bsr->shortreadsortrefs = gt_malloc(sizeof (*bsr->shortreadsortrefs) *
                                          (sfxstrategy->maxshortreadsort+1));
      for (idx = 0; idx <= sfxstrategy->maxshortreadsort; idx++)
      {
        bsr->shortreadsortrefs[idx] = (uint16_t) idx;
      }
    } else
    {
      bsr->shortreadsortinfo = NULL;
      bsr->shortreadsortrefs = NULL;
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
  bsr->processunsortedsuffixrangeinfo = NULL;
  bsr->processunsortedsuffixrange = NULL;
  bsr->sortmaxdepth = sortmaxdepth;
  if (sortmaxdepth > 0 && bsr->sfxstrategy->maxinsertionsort >= 2UL)
  {
    bsr->equalwithprevious = gt_malloc(sizeof (*bsr->equalwithprevious) *
                                       bsr->sfxstrategy->maxinsertionsort);
    for (idx=0; idx < bsr->sfxstrategy->maxinsertionsort; idx++)
    {
      bsr->equalwithprevious[idx] = false;
    }
  } else
  {
    bsr->equalwithprevious = NULL;
  }
  bsr->countinsertionsort = 0;
  bsr->counttqsort = 0;
  bsr->countcountingsort = 0;
  bsr->countbltriesort = 0;
  bsr->countshortreadsort = 0;
}

static void bentsedgresources_delete(GtBentsedgresources *bsr, GtLogger *logger)
{
  gt_free(bsr->countingsortinfo);
  bsr->countingsortinfo = NULL;
  gt_free(bsr->medianinfospace);
  bsr->medianinfospace = NULL;
  gt_free(bsr->shortreadsortinfo);
  bsr->shortreadsortinfo = NULL;
  gt_free(bsr->shortreadsortrefs);
  bsr->shortreadsortrefs = NULL;
  gt_blindtrie_delete(bsr->blindtrie);
  gt_encseq_reader_delete(bsr->esr1);
  gt_encseq_reader_delete(bsr->esr2);
  gt_free(bsr->equalwithprevious);
  GT_FREEARRAY(&bsr->mkvauxstack,GtMKVstack);
  gt_logger_log(logger,"countinsertionsort=%lu",bsr->countinsertionsort);
  gt_logger_log(logger,"countbltriesort=%lu",bsr->countbltriesort);
  gt_logger_log(logger,"countcountingsort=%lu",bsr->countcountingsort);
  gt_logger_log(logger,"countshortreadsort=%lu",bsr->countshortreadsort);
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
                       Outlcpinfo *outlcpinfo,
                       unsigned int sortmaxdepth,
                       const Sfxstrategy *sfxstrategy,
                       GtProcessunsortedsuffixrange
                         processunsortedsuffixrange,
                       void *processunsortedsuffixrangeinfo,
                       unsigned long long *bucketiterstep,
                       GtLogger *logger)
{
  GtCodetype code;
  unsigned int rightchar = (unsigned int) (mincode % numofchars);
  GtBucketspecification bucketspec;
  GtBentsedgresources bsr;

  initBentsedgresources(&bsr,
                        suffixsortspace,
                        encseq,
                        readmode,
                        prefixlength,
                        sortmaxdepth,
                        sfxstrategy);
  gt_bcktab_determinemaxsize(bcktab, mincode, maxcode, numberofsuffixes);
  if (outlcpinfo != NULL)
  {
    bsr.tableoflcpvalues = gt_Outlcpinfo_resizereservoir(outlcpinfo,bcktab);
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
                                 Outlcpinfo *outlcpinfo,
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
    initBentsedgresources(&bsr,
                          suffixsortspace,
                          encseq,
                          readmode,
                          0,
                          sortmaxdepth,
                          sfxstrategy);
    if (outlcpinfo != NULL)
    {
      bsr.tableoflcpvalues = gt_Outlcpinfo_resizereservoir(outlcpinfo,NULL);
    }
    bsr.processunsortedsuffixrangeinfo = processunsortedsuffixrangeinfo;
    bsr.processunsortedsuffixrange = processunsortedsuffixrange;
    gt_suffixsortspace_bucketleftidx_set(bsr.sssp,0);
    gt_sort_bentleysedgewick(&bsr,numberofsuffixes,0);
    gt_suffixsortspace_bucketleftidx_set(bsr.sssp,0);
    bentsedgresources_delete(&bsr, logger);
  }
}
