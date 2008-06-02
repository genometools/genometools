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

#undef SKDEBUG

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

static void showmaxleqvalue(FILE *fp,unsigned long maxleqk,
                            unsigned long patternlength)
{
  if (maxleqk == UNDEFINDEX)
  {
    fprintf(fp,"undefined");
  } else
  {
    fprintf(fp,"%lu",maxleqk);
  }
}

static void verifycolumnvalues(unsigned long patternlength,
                               unsigned long maxdistance,
                               const Myerscolumn *col,
                               unsigned long startscore)
{
  unsigned long idx, score, minscore, mask, bfmaxleqk;

  if (startscore <= maxdistance)
  {
    bfmaxleqk = 0;
    minscore = startscore;
  } else
  {
    bfmaxleqk = UNDEFINDEX;
    minscore = 0;
  }
  score = startscore;
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
      bfmaxleqk = idx;
      minscore = score;
    }
  }
  if (bfmaxleqk != col->maxleqk)
  {
    fprintf(stderr,"correct maxleqk = ");
    showmaxleqvalue(stderr,bfmaxleqk,patternlength);
    fprintf(stderr," != ");
    showmaxleqvalue(stderr,col->maxleqk,patternlength);
    fprintf(stderr," = col->maxleqk\n");
    exit(EXIT_FAILURE);
  }
  if (bfmaxleqk != UNDEFINDEX && minscore != col->scorevalue)
  {
    fprintf(stderr,"correct score = %lu != %lu = col->score\n",
                 minscore,
                 col->scorevalue);
    exit(EXIT_FAILURE);
  }
}

static void showcolumn(const Myerscolumn *col,unsigned long score,
                       unsigned long patternlength)
{
  if (col->maxleqk == UNDEFINDEX)
  {
    printf("[]");
  } else
  {
    unsigned long idx, backmask;

    printf("[%lu",score);
    for (idx=1UL, backmask = 1UL; idx<=col->maxleqk; idx++, backmask <<= 1)
    {
      if (col->Pv & backmask)
      {
        score++;
      } else
      {
        if (col->Mv & backmask)
        {
          score--;
        }
      }
      printf(",%lu",score);
    }
    printf("] with maxleqk=%lu",col->maxleqk);
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

  assert(incol->maxleqk != UNDEFINDEX && incol->maxleqk != patternlength);
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
  if (incol->maxleqk == patternlength)
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

static void inplacenextEDcolumn(const unsigned long *eqsvector,
                                unsigned long patternlength,
                                unsigned long maxdistance,
                                Myerscolumn *col,
                                Uchar currentchar)
{
  unsigned long Eq, Xv, Xh, Ph, Mh, /* as in Myers Paper */
                backmask;           /* only one bit is on */
  unsigned long idx;                /* a counter */
  unsigned long score;              /* current score */

#ifdef SKDEBUG
  if (col->maxleqk == patternlength)
  {
    fprintf(stderr,"col->maxleqk = %lu = patternlength not allowed\n",
            patternlength);
    exit(EXIT_FAILURE);
  }
  if (col->maxleqk == UNDEFINDEX)
  {
    fprintf(stderr,"col->maxleqk = UNDEFINDEX not allowed\n");
    exit(EXIT_FAILURE);
  }
#endif
  Eq = eqsvector[(unsigned long) currentchar];
  Xv = Eq | col->Mv;
  Xh = (((Eq & col->Pv) + col->Pv) ^ col->Pv) | Eq;

  Ph = col->Mv | ~ (Xh | col->Pv);
  Mh = col->Pv & Xh;

  Ph = (Ph << 1) | 1UL;
  col->Pv = (Mh << 1) | ~ (Xv | Ph);
  col->Mv = Ph & Xv;
  backmask = 1UL << col->maxleqk;
  if (Eq & backmask || Mh & backmask)
  {
    col->maxleqk++;
  } else
  {
    if (Ph & backmask)
    {
      unsigned long tmpmaxleqk = UNDEFINDEX;
      score = maxdistance+1;
      if (col->maxleqk > 0)
      {
        for (idx = col->maxleqk - 1, backmask >>= 1;
             /* Nothing */;
             backmask >>= 1)
        {
          if (col->Pv & backmask)
          {
            score--;
            if (score <= maxdistance)
            {
              tmpmaxleqk = idx;
#ifdef SKDEBUG
              col->scorevalue = score;
#endif
              break;
            }
          } else
          {
            if (col->Mv & backmask)
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
      col->maxleqk = tmpmaxleqk;
    }
  }
}

static bool iternextEDcolumn(unsigned long *matchpref,
                             const unsigned long *eqsvector,
                             unsigned long patternlength,
                             unsigned long maxdistance,
                             const Myerscolumn *col,
                             Seqpos startpos,
                             Readmode readmode,
                             const Encodedsequence *encseq,
                             Seqpos totallength)
{
  Seqpos pos;
  Uchar cc;
  Myerscolumn currentcol = *col;

  for (pos = startpos; pos < totallength; pos++)
  {
    cc = getencodedchar(encseq,pos,readmode);
    if (cc == (Uchar) SEPARATOR)
    {
      break;
    }
    inplacenextEDcolumn(eqsvector,
                        patternlength,
                        maxdistance,
                        &currentcol,
                        cc);
    if (currentcol.maxleqk == UNDEFINDEX)
    {
      break;
    }
    if (currentcol.maxleqk == patternlength)
    {
      *matchpref = (unsigned long) (pos - startpos + 1);
      return true;
    }
  }
  return false;
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
  const Seqpos *suftab;
  Uchar alphasize;
};

Limdfsresources *newLimdfsresources(unsigned int mapsize,const Seqpos *suftab)
{
  Limdfsresources *limdfsresources;

  ALLOCASSIGNSPACE(limdfsresources,NULL,Limdfsresources,1);
  ALLOCASSIGNSPACE(limdfsresources->eqsvector,NULL,unsigned long,mapsize-1);
  ALLOCASSIGNSPACE(limdfsresources->rbwc,NULL,Rightboundwithchar,mapsize);
  INITARRAY(&limdfsresources->stack,Lcpintervalwithinfo);
  assert(mapsize-1 <= UCHAR_MAX);
  limdfsresources->alphasize = (Uchar) (mapsize-1);
  limdfsresources->suftab = suftab;
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

static void showmatch(unsigned long patternlength,
                      unsigned long maxdistance,
                      Seqpos startpos, Seqpos offset, unsigned long matchpref)
{
  Seqpos matchlen = offset + matchpref;

  printf("match " FormatSeqpos " " FormatSeqpos "+%lu\n",
         PRINTSeqposcast(startpos), PRINTSeqposcast(offset), matchpref);
  assert((Seqpos) (patternlength - maxdistance) <= matchlen &&
         matchlen <= (Seqpos) (patternlength + maxdistance));
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
#ifdef SKDEBUG
  printf("nextEDcol(");
  assert(offset > 0);
  showcolumn(incol,(unsigned long) (offset-1),patternlength);
#endif
  nextEDcolumn(limdfsresources->eqsvector,
               patternlength,
               maxdistance,
               &stackptr->column,
               inchar,
               incol);
#ifdef SKDEBUG
  printf(",%u)=",(unsigned int) inchar);
  showcolumn(&stackptr->column,(unsigned long) offset,patternlength);
  printf("\n");
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
    Seqpos idx;

    for (idx = lbound; idx <= rbound; idx++)
    {
      showmatch(patternlength,maxdistance,limdfsresources->suftab[idx],
                offset,0);
    }
    return true;
  }
  return false;
}

void esalimiteddfs(Limdfsresources *limdfsresources,
                   const Encodedsequence *encseq,
                   const Alphabet *alpha,
                   Readmode readmode,
                   const Uchar *pattern,
                   unsigned long patternlength,
                   unsigned long maxdistance)
{
  Lcpintervalwithinfo *stackptr;
  unsigned long matchpref = 0, idx, rboundscount;
  Uchar extendchar;
  Seqpos lbound, rbound, bound, offset,
         totallength = getencseqtotallength(encseq);
  Myerscolumn previouscolumn;
  bool remstack;

  printf("# patternlength=%lu\n",patternlength);
  printf("# maxdistance=%lu\n",maxdistance);
  printf("# tag=");
  showsymbolstringgeneric(stdout,alpha,pattern,patternlength);
  printf("\n");
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
                                                limdfsresources->alphasize+1,
                                                encseq,
                                                limdfsresources->suftab,
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
      assert(lbound == rbound);
      assert (previouscolumn.maxleqk != UNDEFINDEX &&
              previouscolumn.maxleqk != patternlength);
      /*
      printf("(2) singleton " FormatSeqpos "\n",PRINTSeqposcast(lbound));
      */
      if (iternextEDcolumn(&matchpref,
                           limdfsresources->eqsvector,
                           patternlength,
                           maxdistance,
                           &previouscolumn,
                           limdfsresources->suftab[lbound],
                           readmode,
                           encseq,
                           totallength))
      {
        showmatch(patternlength,maxdistance,
                  limdfsresources->suftab[lbound],
                  (Seqpos) 1,matchpref);
      }
    }
  }
  while (limdfsresources->stack.nextfreeLcpintervalwithinfo > 0)
  {
    assert(limdfsresources->stack.spaceLcpintervalwithinfo != NULL);
    stackptr = limdfsresources->stack.spaceLcpintervalwithinfo +
               limdfsresources->stack.nextfreeLcpintervalwithinfo - 1;
#ifdef SKDEBUG
    printf("top=(offset=%lu,%lu,%lu) with ",
                (unsigned long) stackptr->lcpitv.offset,
                (unsigned long) stackptr->lcpitv.left,
                (unsigned long) stackptr->lcpitv.right);
    showcolumn(&stackptr->column,
               (unsigned long) stackptr->lcpitv.offset,patternlength);
    printf("\n");
#endif
    extendchar = lcpintervalextendlcp(encseq,
                                      limdfsresources->suftab,
                                      &stackptr->lcpitv,
                                      limdfsresources->alphasize);
    previouscolumn = stackptr->column;
    if (extendchar < limdfsresources->alphasize)
    {
#ifdef SKDEBUG
      printf("nextEDcol(");
      showcolumn(&previouscolumn,
                 (unsigned long) stackptr->lcpitv.offset,patternlength);
      printf(",%u)=",(unsigned int) extendchar);
#endif
      nextEDcolumn(limdfsresources->eqsvector,
                   patternlength,
                   maxdistance,
                   &stackptr->column,
                   extendchar,
                   &previouscolumn);
      stackptr->lcpitv.offset++;
      if (stackptr->column.maxleqk == UNDEFINDEX)
      {
        assert(limdfsresources->stack.nextfreeLcpintervalwithinfo > 0);
        limdfsresources->stack.nextfreeLcpintervalwithinfo--;
      }
      if (stackptr->column.maxleqk == patternlength)
      {
        for (bound = stackptr->lcpitv.left; bound <= stackptr->lcpitv.right;
             bound++)
        {
          showmatch(patternlength,maxdistance,
                    limdfsresources->suftab[bound],stackptr->lcpitv.offset,0);
        }
        assert(limdfsresources->stack.nextfreeLcpintervalwithinfo > 0);
        limdfsresources->stack.nextfreeLcpintervalwithinfo--;
      }
#ifdef SKDEBUG
      showcolumn(&stackptr->column,
                 (unsigned long) stackptr->lcpitv.offset,patternlength);
      printf("\n");
      verifycolumnvalues(patternlength,
                         maxdistance,
                         &stackptr->column,
                         (unsigned long) stackptr->lcpitv.offset);
#endif
    } else
    {
      rboundscount = lcpintervalsplitwithoutspecial(
                              limdfsresources->rbwc,
                              limdfsresources->alphasize+1,
                              encseq,
                              limdfsresources->suftab,
                              stackptr->lcpitv.offset,
                              stackptr->lcpitv.left,
                              stackptr->lcpitv.right);
      offset = stackptr->lcpitv.offset;
      assert(limdfsresources->stack.nextfreeLcpintervalwithinfo > 0);
      limdfsresources->stack.nextfreeLcpintervalwithinfo--;
      for (idx=0; idx < rboundscount; idx++)
      {
        lbound = limdfsresources->rbwc[idx].bound;
        rbound = limdfsresources->rbwc[idx+1].bound-1;
        assert(lbound <= rbound);
        if (lbound < rbound)
        {
          GETNEXTFREEINARRAY(stackptr,&limdfsresources->stack,
                             Lcpintervalwithinfo,128);
          remstack = possiblypush(limdfsresources,
                                  patternlength,
                                  maxdistance,
                                  stackptr,
                                  offset+1,
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
          assert(lbound == rbound);
          /*
          printf("(1) singleton " FormatSeqpos "\n",PRINTSeqposcast(lbound));
          */
          if (iternextEDcolumn(&matchpref,
                               limdfsresources->eqsvector,
                               patternlength,
                               maxdistance,
                               &previouscolumn,
                               limdfsresources->suftab[lbound] + offset,
                               readmode,
                               encseq,
                               totallength))
          {
            showmatch(patternlength,maxdistance,
                      limdfsresources->suftab[lbound],offset,matchpref);
          }
        }
      }
    }
  }
}
