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

#include "libgtcore/chardef.h"
#include "libgtcore/symboldef.h"
#include "libgtcore/unused.h"
#include "libgtcore/arraydef.h"
#include "sarr-def.h"
#include "seqpos-def.h"
#include "esa-splititv.h"
#include "spacedef.h"
#include "esa-limdfs.h"

#define UNDEFINDEX      (patternlength+1)

#define SKDEBUG

typedef struct
{
  unsigned long Pv,    /* the plus-vector for Myers Algorithm */
                Mv;    /* the minus-vector for Myers Algorithm */
  unsigned long maxleqk;  /* \(\max\{i\in[0,m]\mid D(i)\leq k\}\) where
                          \(m\) is the length of the pattern, \(k\) is the
                          distance threshold, and \(D\) is
                          the current distance column */
#ifdef SKDEBUG
  unsigned long scorevalue;    /* the score for the given depth */
#endif
} Myerscolumn;

#ifdef SKDEBUG
static void verifycolumnvalues(unsigned long patternlength,
                               unsigned long maxdistance,
                               const Myerscolumn *col,
                               unsigned long startscore)
{
  unsigned long idx, score = startscore, minscore, mask;
  unsigned long maxleqkindex;

  if (score <= maxdistance)
  {
    maxleqkindex = 0;
    minscore = score;
  } else
  {
    maxleqkindex = UNDEFINDEX;
    minscore = 0;
  }
  for (idx=1UL, mask = 1UL; idx <= patternlength; idx++, mask <<= 1)
  {
    if (col->Pv & mask)
    {
      score++;
    } else
    {
      if (col->Mv & mask)
      {
        score--;
      }
    }
    if (score <= maxdistance)
    {
      maxleqkindex = idx;
      minscore = score;
    }
  }
  if (maxleqkindex != col->maxleqk)
  {
    fprintf(stderr,"correct maxleqkindex = %lu != %lu = col->maxleqk\n",
                   maxleqkindex,
                   col->maxleqk);
    exit(EXIT_FAILURE);
  }
  if (maxleqkindex != UNDEFINDEX)
  {
    if (minscore != col->scorevalue)
    {
      fprintf(stderr,"correct score = %lu != %lu = col->score\n",
                   minscore,
                   col->scorevalue);
      exit(EXIT_FAILURE);
    }
  }
}
#endif

static void initeqsvector(unsigned long *eqsvector,
                          unsigned long eqslen,
                          const Uchar *u,
                          unsigned long ulen)
{
  unsigned long *vptr, shiftmask;
  const Uchar *uptr;

  for (vptr = eqsvector; vptr < eqsvector + eqslen; vptr++)
  {
    *vptr = 0;
  }
  for (uptr = u, shiftmask = 1UL;
       uptr < u + ulen && shiftmask != 0;
       uptr++, shiftmask <<= 1)
  {
    assert (*uptr != (Uchar) SEPARATOR);
    if (*uptr != (Uchar) WILDCARD)
    {
      eqsvector[(unsigned long) *uptr] |= shiftmask;
    }
  }
}

static void nextEDcolumn(const unsigned long *eqsvector,
                         unsigned long patternlength,
                         unsigned long maxdistance,
                         Myerscolumn *outcol,
                         Uchar currentchar,
                         const Myerscolumn *incol)
{
  unsigned long Eq, Xv, Xh, Ph, Mh, /* as in Myers Paper */
                backmask;           /* only one bit is on */
  unsigned long idx;                /* a counter */
  unsigned long score;              /* current score */

  Eq = eqsvector[(unsigned long) currentchar];
  Xv = Eq | incol->Mv;
  Xh = (((Eq & incol->Pv) + incol->Pv) ^ incol->Pv) | Eq;

  Ph = incol->Mv | ~ (Xh | incol->Pv);
  Mh = incol->Pv & Xh;

  Ph = (Ph << 1) | 1UL;
  outcol->Pv = (Mh << 1) | ~ (Xv | Ph);
  outcol->Mv = Ph & Xv;
  /* printf("incol->maxleqk %ld\n",(Showsint) incol->maxleqk); */
#ifdef SKDEBUG
  if ((unsigned long) incol->maxleqk == patternlength)
  {
    fprintf(stderr,"incol->maxleqk = %lu = patternlength not allowed\n",
            patternlength);
    exit(EXIT_FAILURE);
  }
  if (incol->maxleqk == UNDEFINDEX)
  {
    fprintf(stderr,"incol->maxleqk = UNDEFINDEX not allowed\n");
    exit(EXIT_FAILURE);
  }
#endif
  backmask = 1UL << incol->maxleqk;
  if (Eq & backmask || Mh & backmask)
  {
    outcol->maxleqk = incol->maxleqk + 1UL;
#ifdef SKDEBUG
    outcol->scorevalue = incol->scorevalue;
#endif
  } else
  {
    if (Ph & backmask)
    {
      score = maxdistance+1;
      outcol->maxleqk = UNDEFINDEX;
      if (incol->maxleqk > 0)
      {
        for (idx = incol->maxleqk - 1, backmask >>= 1;
             /* Nothing */;
             backmask >>= 1)
        {
          if (outcol->Pv & backmask)
          {
            score--;
            if (score <= maxdistance)
            {
              outcol->maxleqk = idx;
#ifdef SKDEBUG
              outcol->scorevalue = score;
#endif
              break;
            }
          } else
          {
            if (outcol->Mv & backmask)
            {
              score++;
            }
          }
          if (idx > 0)
          {
            idx--;
          } else
          {
            break;
          }
        }
      }
    } else
    {
      outcol->maxleqk = incol->maxleqk;
#ifdef SKDEBUG
      outcol->scorevalue = incol->scorevalue;
#endif
    }
  }
}

typedef struct
{
  Lcpinterval lcpitv;
  Myerscolumn column;
} Lcpintervalwithinfo;

DECLAREARRAYSTRUCT(Lcpintervalwithinfo);

struct Limdfsresources
{
  unsigned long *eqsvector;
  Rightboundwithchar *rbwc;
  ArrayLcpintervalwithinfo stack;
  Uchar alphasize;
};

Limdfsresources *newLimdfsresources(unsigned int mapsize)
{
  Limdfsresources *limdfsresources;

  ALLOCASSIGNSPACE(limdfsresources,NULL,Limdfsresources,1);
  ALLOCASSIGNSPACE(limdfsresources->eqsvector,NULL,unsigned long,mapsize-1);
  ALLOCASSIGNSPACE(limdfsresources->rbwc,NULL,Rightboundwithchar,mapsize-1);
  INITARRAY(&limdfsresources->stack,Lcpintervalwithinfo);
  assert(mapsize-1 <= UCHAR_MAX);
  limdfsresources->alphasize = (Uchar) (mapsize-1);
  return limdfsresources;
}

void freeLimdfsresources(Limdfsresources **ptrlimdfsresources)
{
  Limdfsresources *limdfsresources = *ptrlimdfsresources;

  FREESPACE(limdfsresources->eqsvector);
  FREESPACE(limdfsresources->rbwc);
  FREEARRAY(&limdfsresources->stack,Lcpintervalwithinfo);
  FREESPACE(*ptrlimdfsresources);
}

static bool possiblypush(Limdfsresources *limdfsresources,
                         unsigned long patternlength,
                         unsigned long maxdistance,
                         Lcpintervalwithinfo *stackptr,
                         Seqpos offset,
                         Seqpos lbound,
                         Seqpos rbound,
                         Uchar inchar,
                         const Myerscolumn *incol)
{
  stackptr->lcpitv.offset = offset;
  stackptr->lcpitv.left = lbound;
  stackptr->lcpitv.right = rbound;
  nextEDcolumn(limdfsresources->eqsvector,
               patternlength,
               maxdistance,
               &stackptr->column,
               inchar,
               incol);
#ifdef SKDEBUG
  verifycolumnvalues(patternlength,
                     maxdistance,
                     &stackptr->column,
                     (unsigned long) offset);
#endif
  if (stackptr->column.maxleqk == UNDEFINDEX)
  {
    return true;
  }
  if (stackptr->column.maxleqk == patternlength)
  {
    printf("success " FormatSeqpos " " FormatSeqpos "\n",
            PRINTSeqposcast(stackptr->lcpitv.left),
            PRINTSeqposcast(stackptr->lcpitv.right));
    return true;
  }
  return false;
}

void esalimiteddfs(Limdfsresources *limdfsresources,
                   const Encodedsequence *encseq,
                   const Seqpos *suftab,
                   const Uchar *pattern,
                   unsigned long patternlength,
                   unsigned long maxdistance)
{
  Lcpintervalwithinfo *stackptr;
  unsigned long idx, rboundscount;
  Uchar extendchar;
  Seqpos lbound, rbound, offset, totallength = getencseqtotallength(encseq);
  Myerscolumn previouscolumn;
  bool remstack;

  limdfsresources->stack.nextfreeLcpintervalwithinfo = 0;
  assert(maxdistance < patternlength);
  previouscolumn.Pv = ~0UL;
  previouscolumn.Mv = 0UL;
  previouscolumn.maxleqk = maxdistance;
  initeqsvector(limdfsresources->eqsvector,
                (unsigned long) limdfsresources->alphasize,
                pattern,patternlength);
#ifdef SKDEBUG
  previouscolumn.scorevalue = maxdistance;
#endif
  rboundscount = lcpintervalsplitwithoutspecial(limdfsresources->rbwc,
                                                limdfsresources->alphasize,
                                                encseq,
                                                suftab,
                                                0,
                                                0,
                                                totallength);
  for (idx=0; idx < rboundscount; idx++)
  {
    lbound = limdfsresources->rbwc[idx].bound;
    rbound = limdfsresources->rbwc[idx+1].bound-1;
    assert(lbound <= rbound);
    if (lbound < rbound)
    {
      assert(limdfsresources->stack.spaceLcpintervalwithinfo != NULL);
      GETNEXTFREEINARRAY(stackptr,&limdfsresources->stack,
                         Lcpintervalwithinfo,128);
      remstack = possiblypush(limdfsresources,
                              patternlength,
                              maxdistance,
                              stackptr,
                              (Seqpos) 1,
                              lbound,
                              rbound,
                              limdfsresources->rbwc[idx].inchar,
                              &previouscolumn);
      if (remstack)
      {
        limdfsresources->stack.nextfreeLcpintervalwithinfo--;
      }
    } else
    {
      printf("singleton " FormatSeqpos "\n",PRINTSeqposcast(lbound));
    }
  }
  while (limdfsresources->stack.nextfreeLcpintervalwithinfo > 0)
  {
    assert(limdfsresources->stack.spaceLcpintervalwithinfo != NULL);
    stackptr = limdfsresources->stack.spaceLcpintervalwithinfo +
               limdfsresources->stack.nextfreeLcpintervalwithinfo - 1;
    extendchar = lcpintervalextendlcp(encseq,
                                      suftab,
                                      &stackptr->lcpitv,
                                      limdfsresources->alphasize);
    previouscolumn = stackptr->column;
    if (extendchar < limdfsresources->alphasize)
    {
      nextEDcolumn(limdfsresources->eqsvector,
                   patternlength,
                   maxdistance,
                   &stackptr->column,
                   extendchar,
                   &previouscolumn);
      stackptr->lcpitv.offset++;
#ifdef SKDEBUG
      verifycolumnvalues(patternlength,
                         maxdistance,
                         &stackptr->column,
                         (unsigned long) stackptr->lcpitv.offset);
#endif
    } else
    {
      rboundscount = lcpintervalsplitwithoutspecial(limdfsresources->rbwc,
                                                    limdfsresources->alphasize,
                                                    encseq,
                                                    suftab,
                                                    stackptr->lcpitv.offset,
                                                    stackptr->lcpitv.left,
                                                    stackptr->lcpitv.right);
      for (idx=0; idx < rboundscount; idx++)
      {
        lbound = limdfsresources->rbwc[idx].bound;
        rbound = limdfsresources->rbwc[idx+1].bound-1;
        offset = stackptr->lcpitv.offset + 1;
        assert(lbound <= rbound);
        if (lbound < rbound)
        {
          GETNEXTFREEINARRAY(stackptr,&limdfsresources->stack,
                             Lcpintervalwithinfo,128);
          remstack = possiblypush(limdfsresources,
                                  patternlength,
                                  maxdistance,
                                  stackptr,
                                  offset,
                                  lbound,
                                  rbound,
                                  limdfsresources->rbwc[idx].inchar,
                                  &previouscolumn);
          if (remstack)
          {
            limdfsresources->stack.nextfreeLcpintervalwithinfo--;
          }
        } else
        {
          printf("singleton " FormatSeqpos "\n",PRINTSeqposcast(lbound));
        }
      }
    }
  }
}
