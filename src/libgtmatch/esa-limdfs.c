/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "spacedef.h"
#include "divmodmul.h"
#include "esa-splititv.h"
#include "esa-limdfs.h"
#include "eis-voiditf.h"
#include "defined-types.h"
#include "stamp.h"

#define UNDEFINDEX      (patternlength+1)

#undef SKDEBUG

typedef struct
{
  unsigned long Pv,       /* the plus-vector for Myers Algorithm */
                Mv,       /* the minus-vector for Myers Algorithm */
                maxleqk;  /* \(\max\{i\in[0,m]\mid D(i)\leq k\}\) where
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

static void showMyerscolumn(const Myerscolumn *col,unsigned long score,
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

#ifdef SKDEBUG

static void showinterval(bool withesa,const Lcpinterval *itv)
{
  Seqpos width;

  width = withesa ? (itv->right - itv->left + 1) : itv->right - itv->left;
  printf("(%lu,width=%lu)",(unsigned long) itv->offset,(unsigned long) width);
  printf("(%lu,%lu)",(unsigned long) itv->left,
                     (unsigned long) (withesa ? itv->right : itv->right-1));
}

#define SHOWSTACKTOP(STACKPTR)\
        printf("top=");\
        showinterval(limdfsresources->withesa,&STACKPTR->lcpitv);\
        showMyerscolumn(&STACKPTR->column,\
               (unsigned long) STACKPTR->lcpitv.offset,patternlength);\
        printf("\n")
#else
#define SHOWSTACKTOP(STACKPTR) /* Nothing */
#endif

static void nextEDcolumn(const unsigned long *eqsvector,
                         unsigned long patternlength,
                         unsigned long maxdistance,
                         Myerscolumn *outcol,
                         UNUSED unsigned long startscore,
                         Uchar currentchar,
                         const Myerscolumn *incol)
{
  unsigned long Eq = 0, Xv, Xh, Ph, Mh, /* as in Myers Paper */
                backmask,               /* only one bit is on */
                idx,                    /* a counter */
                score;                  /* current score */

  assert(incol->maxleqk != UNDEFINDEX && incol->maxleqk != patternlength);
  assert(currentchar != (Uchar) SEPARATOR);
  if (currentchar != (Uchar) WILDCARD)
  {
    Eq = eqsvector[(unsigned long) currentchar];
  }
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
#ifdef SKDEBUG
  verifycolumnvalues(patternlength,
                     maxdistance,
                     outcol,
                     startscore+1);
#endif
}

static void inplacenextEDcolumn(const unsigned long *eqsvector,
                                unsigned long patternlength,
                                unsigned long maxdistance,
                                Myerscolumn *col,
                                Uchar currentchar)
{
  unsigned long Eq = 0, Xv, Xh, Ph, Mh, /* as in Myers Paper */
                backmask,           /* only one bit is on */
                idx,                /* a counter */
                score;              /* current score */

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
  if (currentchar != (Uchar) WILDCARD)
  {
    Eq = eqsvector[(unsigned long) currentchar];
  }
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

typedef struct
{
  Lcpinterval lcpitv;
  Myerscolumn column;
} Lcpintervalwithinfo;

DECLAREARRAYSTRUCT(Lcpintervalwithinfo);

struct Limdfsresources
{
  unsigned long *eqsvector;
  ArrayBoundswithchar bwci;
  ArrayLcpintervalwithinfo stack;
  Uchar alphasize;
  Seqpos totallength;
  void (*processmatch)(void *,bool,Seqpos,Seqpos,Seqpos);
  void *processmatchinfo;
  const void *genericindex;
  bool withesa, nospecials;
  Seqpos *rangeOccs;
};

Limdfsresources *newLimdfsresources(const void *genericindex,
                                    bool withesa,
                                    bool nospecials,
                                    unsigned int mapsize,
                                    Seqpos totallength,
                                    void (*processmatch)(void *,bool,Seqpos,
                                                         Seqpos,Seqpos),
                                    void *processmatchinfo)
{
  Limdfsresources *limdfsresources;

  ALLOCASSIGNSPACE(limdfsresources,NULL,Limdfsresources,1);
  ALLOCASSIGNSPACE(limdfsresources->eqsvector,NULL,unsigned long,mapsize-1);
  ALLOCASSIGNSPACE(limdfsresources->bwci.spaceBoundswithchar,NULL,
                   Boundswithchar,mapsize);
  limdfsresources->bwci.nextfreeBoundswithchar = 0;
  limdfsresources->bwci.allocatedBoundswithchar
    = (unsigned long) mapsize;
  INITARRAY(&limdfsresources->stack,Lcpintervalwithinfo);
  assert(mapsize-1 <= UCHAR_MAX);
  limdfsresources->alphasize = (Uchar) (mapsize-1);
  limdfsresources->processmatch = processmatch;
  limdfsresources->processmatchinfo = processmatchinfo;
  limdfsresources->genericindex = genericindex;
  limdfsresources->totallength = totallength;
  limdfsresources->withesa = withesa;
  limdfsresources->nospecials = nospecials;
  if (withesa)
  {
    limdfsresources->rangeOccs = NULL;
  } else
  {
    ALLOCASSIGNSPACE(limdfsresources->rangeOccs,NULL,Seqpos,
                     MULT2(limdfsresources->alphasize));
  }
  return limdfsresources;
}

const void *getgenericindexfromresource(Limdfsresources *limdfsresources)
{
  return limdfsresources->genericindex;
}

static void initMyerscolumn(Myerscolumn *column,unsigned long maxdistance)
{
  column->Pv = ~0UL;
  column->Mv = 0UL;
  column->maxleqk = maxdistance;
#ifdef SKDEBUG
  column->scorevalue = maxdistance;
#endif
}

static void initlcpinfostack(ArrayLcpintervalwithinfo *stack,
                             Seqpos left,
                             Seqpos right,
                             unsigned long maxdistance)
{
  Lcpintervalwithinfo *stackptr;

  stack->nextfreeLcpintervalwithinfo = 0;
  GETNEXTFREEINARRAY(stackptr,stack,Lcpintervalwithinfo,128);
  stackptr->lcpitv.offset = 0;
  stackptr->lcpitv.left = left;
  stackptr->lcpitv.right = right;
  initMyerscolumn(&stackptr->column,maxdistance);
}

void freeLimdfsresources(Limdfsresources **ptrlimdfsresources)
{
  Limdfsresources *limdfsresources = *ptrlimdfsresources;

  FREESPACE(limdfsresources->eqsvector);
  FREEARRAY(&limdfsresources->bwci,Boundswithchar);
  FREEARRAY(&limdfsresources->stack,Lcpintervalwithinfo);
  FREESPACE(limdfsresources->rangeOccs);
  FREESPACE(*ptrlimdfsresources);
}

/* enumerate the suffixes in an LCP-interval */

static void esa_overinterval(Limdfsresources *limdfsresources,
                             const Lcpinterval *itv)
{
  Seqpos idx;
  const Suffixarray *suffixarray
    = (const Suffixarray *) limdfsresources->genericindex;

  for (idx = itv->left; idx <= itv->right; idx++)
  {
    limdfsresources->processmatch(limdfsresources->processmatchinfo,
                                  true,
                                  limdfsresources->totallength,
                                  suffixarray->suftab[idx],
                                  itv->offset);
  }
}

static void pck_overinterval(Limdfsresources *limdfsresources,
                             const Lcpinterval *itv)
{
  Bwtseqpositioniterator *bspi;
  Seqpos dbstartpos;

  bspi = newBwtseqpositioniterator (limdfsresources->genericindex,
                                    itv->left,itv->right);
  assert(itv->left < itv->right);
  while (nextBwtseqpositioniterator(&dbstartpos,bspi))
  {
    limdfsresources->processmatch(limdfsresources->processmatchinfo,
                                  false,
                                  limdfsresources->totallength,
                                  dbstartpos + itv->offset,
                                  itv->offset);
  }
  freeBwtseqpositioniterator(&bspi);
}

Definedunsignedlong esa_findshortestmatch(const Encodedsequence *encseq,
                                          bool nospecials,
                                          unsigned long alphasize,
                                          const Uchar *pattern,
                                          unsigned long patternlength,
                                          unsigned long maxdistance,
                                          Seqpos startpos)
{
  Seqpos pos, totallength = getencseqtotallength(encseq);
  unsigned long *eqsvector;
  Myerscolumn currentcol;
  Uchar cc;
  Definedunsignedlong result;

  ALLOCASSIGNSPACE(eqsvector,NULL,unsigned long,alphasize);
  initeqsvector(eqsvector,alphasize,pattern,patternlength);
  initMyerscolumn(&currentcol,maxdistance);
  for (pos = startpos; /* Nothing */; pos++)
  {
    assert(pos - startpos <= (Seqpos) (patternlength + maxdistance));
    cc = getencodedchar(encseq,pos,Forwardmode);
    if (nospecials && cc == (Uchar) WILDCARD)
    {
      FREESPACE(eqsvector);
      result.defined = false;
      result.valueunsignedlong = 0;
      return result;
    }
    inplacenextEDcolumn(eqsvector,
                        patternlength,
                        maxdistance,
                        &currentcol,
                        cc);
    assert (currentcol.maxleqk != UNDEFINDEX);
    if (currentcol.maxleqk == patternlength || pos == totallength-1)
    {
      break;
    }
  }
  FREESPACE(eqsvector);
  result.defined = true;
  result.valueunsignedlong = (unsigned long) (pos - startpos + 1);
  return result;
}

/* iterate myers algorithm over a sequence context */

static void esa_overcontext(Limdfsresources *limdfsresources,
                            unsigned long patternlength,
                            unsigned long maxdistance,
                            const Myerscolumn *col,
                            Seqpos left,
                            Seqpos offset)
{
  Seqpos pos, startpos;
  Uchar cc;
  Myerscolumn currentcol = *col;
  const Suffixarray *suffixarray
    = (const Suffixarray *) limdfsresources->genericindex;

  startpos = suffixarray->suftab[left];
#ifdef SKDEBUG
  printf("retrieve context of startpos=%lu\n",(unsigned long) startpos);
#endif
  for (pos = startpos + offset - 1; pos < limdfsresources->totallength; pos++)
  {
    cc = getencodedchar(suffixarray->encseq,pos,suffixarray->readmode);
    if (cc != (Uchar) SEPARATOR &&
        (!limdfsresources->nospecials || cc != (Uchar) WILDCARD))
    {
#ifdef SKDEBUG
      printf("cc=%u\n",cc);
#endif
      inplacenextEDcolumn(limdfsresources->eqsvector,
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
        limdfsresources->processmatch(limdfsresources->processmatchinfo,
                                      true,
                                      limdfsresources->totallength,
                                      startpos,
                                      pos - startpos + 1);
        break;
      }
    } else
    {
      break;
    }
  }
}

static void pck_overcontext(Limdfsresources *limdfsresources,
                            unsigned long patternlength,
                            unsigned long maxdistance,
                            const Myerscolumn *col,
                            Uchar inchar,
                            Seqpos left,
                            Seqpos offset)
{
  Uchar cc;
  unsigned long matchlength;
  Myerscolumn currentcol = *col;
  bool processinchar = true;
  Bwtseqcontextiterator *bsci
    = newBwtseqcontextiterator(limdfsresources->genericindex,left);

#ifdef SKDEBUG
  printf("retrieve context of left=%lu\n",(unsigned long) left);
#endif
  for (matchlength = 0; /* nothing */; matchlength++)
  {
    if (processinchar)
    {
      cc = inchar;
      processinchar = false;
    } else
    {
      cc = nextBwtseqcontextiterator(bsci);
    }
    if (cc != (Uchar) SEPARATOR &&
        (!limdfsresources->nospecials || cc != (Uchar) WILDCARD))
    {
#ifdef SKDEBUG
      printf("cc=%u\n",cc);
#endif
      inplacenextEDcolumn(limdfsresources->eqsvector,
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
        Seqpos startpos = bwtseqfirstmatch(limdfsresources->genericindex,left);

        limdfsresources->processmatch(limdfsresources->processmatchinfo,
                                      false,
                                      limdfsresources->totallength,
                                      startpos + offset,
                                      offset + matchlength);
        break;
      }
    } else
    {
      break;
    }
  }
  freeBwtseqcontextiterator(&bsci);
}

static bool pushandpossiblypop(Limdfsresources *limdfsresources,
                               unsigned long patternlength,
                               unsigned long maxdistance,
                               Lcpintervalwithinfo *stackptr,
                               const Lcpinterval *child,
                               Uchar inchar,
                               const Myerscolumn *incol)
{
  stackptr->lcpitv = *child;

#ifdef SKDEBUG
  printf("(2) nextEDcol(");
  showMyerscolumn(incol,(unsigned long) (child->offset-1),patternlength);
  printf(",%u)=",(unsigned int) inchar);
#endif
  nextEDcolumn(limdfsresources->eqsvector,
               patternlength,
               maxdistance,
               &stackptr->column,
               (unsigned long) (child->offset-1),
               inchar,
               incol);
#ifdef SKDEBUG
  showMyerscolumn(&stackptr->column,(unsigned long) child->offset,
                  patternlength);
  printf("\n");
#endif
  if (stackptr->column.maxleqk == UNDEFINDEX)
  {
    return true;
  }
  if (stackptr->column.maxleqk == patternlength)
  {
    if (limdfsresources->withesa)
    {
      esa_overinterval(limdfsresources,child);
    } else
    {
      pck_overinterval(limdfsresources,child);
    }
    return true;
  }
  return false;
}

static void processchildinterval(Limdfsresources *limdfsresources,
                                 unsigned long patternlength,
                                 unsigned long maxdistance,
                                 const Lcpinterval *child,
                                 Uchar inchar,
                                 const Myerscolumn *previouscolumn)
{
  if (child->left + 1 < child->right || (limdfsresources->withesa &&
                                        child->left + 1 == child->right))
  {
    Lcpintervalwithinfo *stackptr;

    GETNEXTFREEINARRAY(stackptr,&limdfsresources->stack,
                       Lcpintervalwithinfo,128);
    if (pushandpossiblypop(limdfsresources,
                           patternlength,
                           maxdistance,
                           stackptr,
                           child,
                           inchar,
                           previouscolumn))
    {
      limdfsresources->stack.nextfreeLcpintervalwithinfo--;
    }
  } else
  {
    if (limdfsresources->withesa)
    {
      esa_overcontext(limdfsresources,
                      patternlength,
                      maxdistance,
                      previouscolumn,
                      child->left,
                      child->offset);
    } else
    {
      pck_overcontext(limdfsresources,
                      patternlength,
                      maxdistance,
                      previouscolumn,
                      inchar,
                      child->left,
                      child->offset);
    }
  }
}

static void esa_splitandprocess(Limdfsresources *limdfsresources,
                                unsigned long patternlength,
                                unsigned long maxdistance,
                                const Lcpinterval *parent,
                                const Myerscolumn *previouscolumn)
{
  Seqpos firstnonspecial;
  Uchar extendchar;
  unsigned long idx;
  const Suffixarray *suffixarray
    = (const Suffixarray *) limdfsresources->genericindex;

  extendchar = lcpintervalextendlcp(suffixarray->encseq,
                                    suffixarray->readmode,
                                    suffixarray->suftab,
                                    limdfsresources->totallength,
                                    parent,
                                    limdfsresources->alphasize);
  if (extendchar < limdfsresources->alphasize)
  {
    limdfsresources->bwci.spaceBoundswithchar[0].lbound = parent->left;
    limdfsresources->bwci.spaceBoundswithchar[0].rbound = parent->right;
    limdfsresources->bwci.spaceBoundswithchar[0].inchar = extendchar;
    limdfsresources->bwci.nextfreeBoundswithchar = 1UL;
  } else
  {
    limdfsresources->bwci.nextfreeBoundswithchar = 0;
    lcpintervalsplitwithoutspecial(&limdfsresources->bwci,
                                   suffixarray->encseq,
                                   suffixarray->readmode,
                                   limdfsresources->totallength,
                                   suffixarray->suftab,
                                   parent);
  }
  firstnonspecial = parent->left;
  for (idx = 0; idx < limdfsresources->bwci.nextfreeBoundswithchar; idx++)
  {
    Lcpinterval child;
    Uchar inchar = limdfsresources->bwci.spaceBoundswithchar[idx].inchar;

    child.offset = parent->offset+1;
    child.left = limdfsresources->bwci.spaceBoundswithchar[idx].lbound;
    child.right = limdfsresources->bwci.spaceBoundswithchar[idx].rbound;
#ifdef SKDEBUG
    printf("%u-child of ",(unsigned int) inchar);
    showinterval(limdfsresources->withesa,parent);
    printf(" is ");
    showinterval(limdfsresources->withesa,&child);
    printf("\n");
#endif
    processchildinterval(limdfsresources,
                         patternlength,
                         maxdistance,
                         &child,
                         inchar,
                         previouscolumn);
    firstnonspecial = child.right+1;
  }
  if (!limdfsresources->nospecials)
  {
    Seqpos bound;

    for (bound=firstnonspecial; bound <= parent->right; bound++)
    {
      esa_overcontext(limdfsresources,
                      patternlength,
                      maxdistance,
                      previouscolumn,
                      bound,
                      parent->offset+1);
    }
  }
}

static void pck_splitandprocess(Limdfsresources *limdfsresources,
                                unsigned long patternlength,
                                unsigned long maxdistance,
                                const Lcpinterval *parent,
                                const Myerscolumn *previouscolumn)
{
  unsigned long idx;

  bwtrangesplitwithoutspecial(&limdfsresources->bwci,
                              limdfsresources->rangeOccs,
                              (unsigned long) limdfsresources->alphasize,
                              limdfsresources->genericindex,
                              parent);
  for (idx = 0; idx < limdfsresources->bwci.nextfreeBoundswithchar; idx++)
  {
    Uchar inchar;
    Lcpinterval child;

    inchar = limdfsresources->bwci.spaceBoundswithchar[idx].inchar;
    child.offset = parent->offset+1;
    child.left = limdfsresources->bwci.spaceBoundswithchar[idx].lbound;
    child.right = limdfsresources->bwci.spaceBoundswithchar[idx].rbound;
#ifdef SKDEBUG
    printf("%u-child of ",(unsigned int) inchar);
    showinterval(limdfsresources->withesa,parent);
    printf(" is ");
    showinterval(limdfsresources->withesa,&child);
    printf("\n");
#endif
    processchildinterval(limdfsresources,
                         patternlength,
                         maxdistance,
                         &child,
                         inchar,
                         previouscolumn);
  }
}

void esalimiteddfs(Limdfsresources *limdfsresources,
                   const Uchar *pattern,
                   unsigned long patternlength,
                   unsigned long maxdistance)
{
  Lcpintervalwithinfo *stackptr;
  Lcpinterval parent;
  Myerscolumn previouscolumn;

  assert(maxdistance < patternlength);
  initeqsvector(limdfsresources->eqsvector,
                (unsigned long) limdfsresources->alphasize,
                pattern,patternlength);
  initlcpinfostack(&limdfsresources->stack,
                   0,
                   limdfsresources->withesa
                     ? limdfsresources->totallength
                     : limdfsresources->totallength+1,
                   maxdistance);
  while (limdfsresources->stack.nextfreeLcpintervalwithinfo > 0)
  {
    assert(limdfsresources->stack.spaceLcpintervalwithinfo != NULL);
    stackptr = limdfsresources->stack.spaceLcpintervalwithinfo +
               limdfsresources->stack.nextfreeLcpintervalwithinfo - 1;
    SHOWSTACKTOP(stackptr);
    previouscolumn = stackptr->column;
    parent = stackptr->lcpitv; /* save, since stackptr is changed in split */
    /* split interval */
    assert(limdfsresources->stack.nextfreeLcpintervalwithinfo > 0);
    limdfsresources->stack.nextfreeLcpintervalwithinfo--;
    if (limdfsresources->withesa)
    {
      esa_splitandprocess(limdfsresources,
                          patternlength,
                          maxdistance,
                          &parent,
                          &previouscolumn);
    } else
    {
      pck_splitandprocess(limdfsresources,
                          patternlength,
                          maxdistance,
                          &parent,
                          &previouscolumn);
    }
  }
}

/*
  Specification steve: only output sequences which occur at most
  t times, where t is a parameter given by the user.
  Output all matches involving a prefix of the pattern and the current
  path with up to k error (k=2 for tagsize around 25).
  Output exakt matching statistics for each suffix of the pattern
*/
