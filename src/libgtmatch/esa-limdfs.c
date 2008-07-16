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

#define UNDEFMAXLEQK      (mti->patternlength+1)
#define SUCCESSMAXLEQK    mti->patternlength

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

typedef struct
{
  unsigned long patternlength,
                maxdistance,
                *eqsvector;
} Matchtaskinfo;

static size_t deliversizeofdfsstate(void)
{
  return sizeof (Myerscolumn);
}

typedef unsigned long Aliasdfsstate;

#define DECLAREPTRDFSSTATE(V)\
        Aliasdfsstate * V

#ifdef SKDEBUG

static void showmaxleqvalue(FILE *fp,unsigned long maxleqk,
                            const Matchtaskinfo *mti)
{
  if (maxleqk == UNDEFMAXLEQK)
  {
    fprintf(fp,"undefined");
  } else
  {
    fprintf(fp,"%lu",maxleqk);
  }
}

static void apm_showLimdfsstate(const DECLAREPTRDFSSTATE(aliascol),
                                unsigned long score,
                                const void *dfsconstinfo)
{
  const Matchtaskinfo *mti = (const Matchtaskinfo *) dfsconstinfo;
  const Myerscolumn *col = (const Myerscolumn *) aliascol;

  if (col->maxleqk == UNDEFMAXLEQK)
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

static void verifycolumnvalues(const Matchtaskinfo *mti,
                               const Myerscolumn *col,
                               unsigned long startscore)
{
  unsigned long idx, score, minscore, mask, bfmaxleqk;

  if (startscore <= mti->maxdistance)
  {
    bfmaxleqk = 0;
    minscore = startscore;
  } else
  {
    bfmaxleqk = UNDEFMAXLEQK;
    minscore = 0;
  }
  score = startscore;
  for (idx=1UL, mask = 1UL; idx <= mti->patternlength; idx++, mask <<= 1)
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
    if (score <= mti->maxdistance)
    {
      bfmaxleqk = idx;
      minscore = score;
    }
  }
  if (bfmaxleqk != col->maxleqk)
  {
    fprintf(stderr,"correct maxleqk = ");
    showmaxleqvalue(stderr,bfmaxleqk,mti);
    fprintf(stderr," != ");
    showmaxleqvalue(stderr,col->maxleqk,mti);
    fprintf(stderr," = col->maxleqk\n");
    exit(EXIT_FAILURE);
  }
  if (bfmaxleqk != UNDEFMAXLEQK && minscore != col->scorevalue)
  {
    fprintf(stderr,"correct score = %lu != %lu = col->score\n",
                 minscore,
                 col->scorevalue);
    exit(EXIT_FAILURE);
  }
}

#endif

static void apm_initdfsconstinfo(void *dfsconstinfo,
                                 unsigned int alphasize,
                                 const Uchar *pattern,
                                 unsigned long patternlength,
                                 unsigned long maxdistance)
{
  unsigned long *eptr, shiftmask;
  const Uchar *pptr;
  Matchtaskinfo *mti = (Matchtaskinfo *) dfsconstinfo;

  for (eptr = mti->eqsvector; eptr < mti->eqsvector + alphasize; eptr++)
  {
    *eptr = 0;
  }
  for (pptr = pattern, shiftmask = 1UL;
       pptr < pattern + patternlength && shiftmask != 0;
       pptr++, shiftmask <<= 1)
  {
    assert (*pptr != (Uchar) SEPARATOR);
    if (*pptr != (Uchar) WILDCARD)
    {
      mti->eqsvector[(unsigned long) *pptr] |= shiftmask;
    }
  }
  mti->patternlength = patternlength;
  mti->maxdistance = maxdistance;
}

static void *apm_allocatedfsconstinfo(unsigned int alphasize)
{
  Matchtaskinfo *mti = ma_malloc(sizeof(Matchtaskinfo));
  mti->eqsvector = ma_malloc(sizeof(*mti->eqsvector) * alphasize);
  return mti;
}

static void apm_freedfsconstinfo(void **dfsconstinfo)
{
  Matchtaskinfo *mti = (Matchtaskinfo *) *dfsconstinfo;

  ma_free(mti->eqsvector);
  ma_free(mti);
  *dfsconstinfo = NULL;
}

static unsigned long apm_limdfsnextstep(const DECLAREPTRDFSSTATE(aliascolumn),
                                        const void *dfsconstinfo)
{
  const Matchtaskinfo *mti = (Matchtaskinfo *) dfsconstinfo;
  const Myerscolumn *limdfsstate = (const Myerscolumn *) aliascolumn;

  if (limdfsstate->maxleqk == UNDEFMAXLEQK)
  {
    return 0;
  }
  if (limdfsstate->maxleqk == SUCCESSMAXLEQK)
  {
    return mti->patternlength+1;
  }
  return 1UL;
}

/* implement types Limdfsstate and all function appearing until here */

static void apm_nextDfsstate(const void *dfsconstinfo,
                             DECLAREPTRDFSSTATE(aliasoutcol),
                             UNUSED unsigned long startscore,
                             Uchar currentchar,
                             const DECLAREPTRDFSSTATE(aliasincol))
{
  unsigned long Eq = 0, Xv, Xh, Ph, Mh, /* as in Myers Paper */
                backmask,               /* only one bit is on */
                idx,                    /* a counter */
                score;                  /* current score */
  const Matchtaskinfo *mti = (const Matchtaskinfo *) dfsconstinfo;
  Myerscolumn *outcol = (Myerscolumn *) aliasoutcol;
  const Myerscolumn *incol = (const Myerscolumn *) aliasincol;

  assert(incol->maxleqk != UNDEFMAXLEQK && incol->maxleqk != SUCCESSMAXLEQK);
  assert(currentchar != (Uchar) SEPARATOR);
  if (currentchar != (Uchar) WILDCARD)
  {
    Eq = mti->eqsvector[(unsigned long) currentchar];
  }
  Xv = Eq | incol->Mv;
  Xh = (((Eq & incol->Pv) + incol->Pv) ^ incol->Pv) | Eq;

  Ph = incol->Mv | ~ (Xh | incol->Pv);
  Mh = incol->Pv & Xh;

  Ph = (Ph << 1) | 1UL;
  outcol->Pv = (Mh << 1) | ~ (Xv | Ph);
  outcol->Mv = Ph & Xv;
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
      score = mti->maxdistance+1;
      outcol->maxleqk = UNDEFMAXLEQK;
      if (incol->maxleqk > 0)
      {
        for (idx = incol->maxleqk - 1, backmask >>= 1;
             /* Nothing */;
             backmask >>= 1)
        {
          if (outcol->Pv & backmask)
          {
            score--;
            if (score <= mti->maxdistance)
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
  verifycolumnvalues(mti,outcol,startscore+1);
#endif
}

static void apm_inplacenextDfsstate(const void *dfsconstinfo,
                                    DECLAREPTRDFSSTATE(limdfsstate),
                                    Uchar currentchar)
{
  unsigned long Eq = 0, Xv, Xh, Ph, Mh, /* as in Myers Paper */
                backmask,           /* only one bit is on */
                idx,                /* a counter */
                score;              /* current score */
  const Matchtaskinfo *mti = (const Matchtaskinfo *) dfsconstinfo;
  Myerscolumn *col = (Myerscolumn *) limdfsstate;

  assert(col->maxleqk != UNDEFMAXLEQK && col->maxleqk != SUCCESSMAXLEQK);
  if (currentchar != (Uchar) WILDCARD)
  {
    Eq = mti->eqsvector[(unsigned long) currentchar];
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
      unsigned long tmpmaxleqk = UNDEFMAXLEQK;
      score = mti->maxdistance+1;
      if (col->maxleqk > 0)
      {
        for (idx = col->maxleqk - 1, backmask >>= 1;
             /* Nothing */;
             backmask >>= 1)
        {
          if (col->Pv & backmask)
          {
            score--;
            if (score <= mti->maxdistance)
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

static void apm_copyDfsstate(DECLAREPTRDFSSTATE(aliasdest),
                             const DECLAREPTRDFSSTATE(aliassrc))
{
  memcpy(aliasdest,aliassrc,sizeof (Myerscolumn));
}

static void apm_initlimdfsstate(DECLAREPTRDFSSTATE(aliascolumn),
                                const void *dfsconstinfo)
{
  Myerscolumn *column = (Myerscolumn *) aliascolumn;

  column->Pv = ~0UL;
  column->Mv = 0UL;
  column->maxleqk = ((Matchtaskinfo *) dfsconstinfo)->maxdistance;
#ifdef SKDEBUG
  column->scorevalue = ((Matchtaskinfo *) dfsconstinfo)->maxdistance;
#endif
}

#define DECLAREDFSSTATE(V)\
        Aliasdfsstate V[4]

typedef struct
{
  Lcpinterval lcpitv;
  DECLAREDFSSTATE(aliasstate);
} Lcpintervalwithinfo;

DECLAREARRAYSTRUCT(Lcpintervalwithinfo);

struct Limdfsresources
{
  void *dfsconstinfo;
  ArrayBoundswithchar bwci;
  ArrayLcpintervalwithinfo stack;
  Uchar alphasize;
  Seqpos totallength;
  void (*processmatch)(void *,bool,Seqpos,Seqpos,Seqpos,unsigned long);
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
                                                         Seqpos,Seqpos,
                                                         unsigned long),
                                    void *processmatchinfo)
{
  Limdfsresources *limdfsresources;

  ALLOCASSIGNSPACE(limdfsresources,NULL,Limdfsresources,1);
  ALLOCASSIGNSPACE(limdfsresources->bwci.spaceBoundswithchar,NULL,
                   Boundswithchar,mapsize);
  limdfsresources->bwci.nextfreeBoundswithchar = 0;
  limdfsresources->bwci.allocatedBoundswithchar = (unsigned long) mapsize;
  INITARRAY(&limdfsresources->stack,Lcpintervalwithinfo);
  assert(mapsize-1 <= MAXALPHABETCHARACTER);
  limdfsresources->alphasize = (Uchar) (mapsize-1);
  limdfsresources->processmatch = processmatch;
  limdfsresources->processmatchinfo = processmatchinfo;
  limdfsresources->genericindex = genericindex;
  limdfsresources->totallength = totallength;
  limdfsresources->withesa = withesa;
  limdfsresources->nospecials = nospecials;
  /* Application specific */
  limdfsresources->dfsconstinfo
    = apm_allocatedfsconstinfo((unsigned int) limdfsresources->alphasize);
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

static void initlcpinfostack(ArrayLcpintervalwithinfo *stack,
                             Seqpos left,
                             Seqpos right,
                             const void *dfsconstinfo)
{
  Lcpintervalwithinfo *stackptr;

  stack->nextfreeLcpintervalwithinfo = 0;
  GETNEXTFREEINARRAY(stackptr,stack,Lcpintervalwithinfo,128);
  stackptr->lcpitv.offset = 0;
  stackptr->lcpitv.left = left;
  stackptr->lcpitv.right = right;
  apm_initlimdfsstate(stackptr->aliasstate,dfsconstinfo);
}

void freeLimdfsresources(Limdfsresources **ptrlimdfsresources)
{
  Limdfsresources *limdfsresources = *ptrlimdfsresources;

  apm_freedfsconstinfo(&limdfsresources->dfsconstinfo);
  FREEARRAY(&limdfsresources->bwci,Boundswithchar);
  FREEARRAY(&limdfsresources->stack,Lcpintervalwithinfo);
  FREESPACE(limdfsresources->rangeOccs);
  FREESPACE(*ptrlimdfsresources);
}

/* enumerate the suffixes in an LCP-interval */

static void esa_overinterval(Limdfsresources *limdfsresources,
                             const Lcpinterval *itv,
                             unsigned long pprefixlen)
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
                                  itv->offset,
                                  pprefixlen);
  }
}

static void pck_overinterval(Limdfsresources *limdfsresources,
                             const Lcpinterval *itv,
                             unsigned long pprefixlen)
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
                                  itv->offset,
                                  pprefixlen);
  }
  freeBwtseqpositioniterator(&bspi);
}

Definedunsignedlong esa_findshortestmatch(const Encodedsequence *encseq,
                                          bool nospecials,
                                          unsigned int alphasize,
                                          const Uchar *pattern,
                                          unsigned long patternlength,
                                          unsigned long maxdistance,
                                          Seqpos startpos)
{
  Seqpos pos, totallength = getencseqtotallength(encseq);
  void *dfsconstinfo;
  Matchtaskinfo *mti;
  Myerscolumn currentcol;
  Uchar cc;
  Definedunsignedlong result;

  dfsconstinfo = apm_allocatedfsconstinfo(alphasize);
  mti = (Matchtaskinfo *) dfsconstinfo;
  apm_initdfsconstinfo(dfsconstinfo,alphasize,pattern,patternlength,
                       maxdistance);
  apm_initlimdfsstate((Aliasdfsstate *) &currentcol,dfsconstinfo);
  for (pos = startpos; /* Nothing */; pos++)
  {
    assert(pos - startpos <= (Seqpos) (patternlength + maxdistance));
    cc = getencodedchar(encseq,pos,Forwardmode);
    if (nospecials && cc == (Uchar) WILDCARD)
    {
      apm_freedfsconstinfo(&dfsconstinfo);
      result.defined = false;
      result.valueunsignedlong = 0;
      return result;
    }
    apm_inplacenextDfsstate(dfsconstinfo,(Aliasdfsstate *) &currentcol,cc);
    assert (currentcol.maxleqk != UNDEFMAXLEQK);
    if (currentcol.maxleqk == SUCCESSMAXLEQK || pos == totallength-1)
    {
      break;
    }
  }
  apm_freedfsconstinfo(&dfsconstinfo);
  result.defined = true;
  result.valueunsignedlong = (unsigned long) (pos - startpos + 1);
  return result;
}

/* iterate myers algorithm over a sequence context */

static void esa_overcontext(Limdfsresources *limdfsresources,
                            const DECLAREPTRDFSSTATE(dfsstate),
                            Seqpos left,
                            Seqpos offset)
{
  Seqpos pos, startpos;
  Uchar cc;
  DECLAREDFSSTATE(currentdfsstate);
  unsigned long pprefixlen;
  const Suffixarray *suffixarray
    = (const Suffixarray *) limdfsresources->genericindex;

  apm_copyDfsstate(currentdfsstate,dfsstate);
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
      printf("cc=%u\n",(unsigned int) cc);
#endif
      apm_inplacenextDfsstate(limdfsresources->dfsconstinfo,
                              currentdfsstate,cc);
      pprefixlen = apm_limdfsnextstep(currentdfsstate,
                                      limdfsresources->dfsconstinfo);
      if (pprefixlen == 0) /* failure */
      {
        break;
      }
      if (pprefixlen > 1UL)
      {
        limdfsresources->processmatch(limdfsresources->processmatchinfo,
                                      true,
                                      limdfsresources->totallength,
                                      startpos,
                                      pos - startpos + 1,
                                      pprefixlen-1);
        break;
      }
    } else
    {
      break;
    }
  }
}

static void pck_overcontext(Limdfsresources *limdfsresources,
                            const DECLAREPTRDFSSTATE(dfsstate),
                            Seqpos left,
                            Seqpos offset,
                            Uchar inchar)
{
  Uchar cc;
  unsigned long contextlength, pprefixlen;
  DECLAREDFSSTATE(currentdfsstate);
  bool processinchar = true;
  Bwtseqcontextiterator *bsci
    = newBwtseqcontextiterator(limdfsresources->genericindex,left);

  apm_copyDfsstate(currentdfsstate,dfsstate);
#ifdef SKDEBUG
  printf("retrieve context for left = %lu\n",(unsigned long) left);
#endif
  for (contextlength = 0; /* nothing */; contextlength++)
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
      printf("cc=%u\n",(unsigned int) cc);
#endif
      apm_inplacenextDfsstate(limdfsresources->dfsconstinfo,
                              currentdfsstate,cc);
      pprefixlen = apm_limdfsnextstep(currentdfsstate,
                                      limdfsresources->dfsconstinfo);
      if (pprefixlen == 0)
      {
        break;
      }
      if (pprefixlen > 1UL) /* check for success */
      {
        Seqpos startpos = bwtseqfirstmatch(limdfsresources->genericindex,left);

        limdfsresources->processmatch(limdfsresources->processmatchinfo,
                                      false,
                                      limdfsresources->totallength,
                                      startpos + offset,
                                      offset + contextlength,
                                      pprefixlen-1);
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
                               Lcpintervalwithinfo *stackptr,
                               const Lcpinterval *child,
                               Uchar inchar,
                               const DECLAREPTRDFSSTATE(indfsstate))
{
  unsigned long pprefixlen;
  stackptr->lcpitv = *child;

#ifdef SKDEBUG
  printf("(2) apm_nextDfsstate(");
  apm_showLimdfsstate(indfsstate,(unsigned long) (child->offset-1),
                      limdfsresources->dfsconstinfo);
  printf(",%u)=",(unsigned int) inchar);
#endif
  apm_nextDfsstate(limdfsresources->dfsconstinfo,
               stackptr->aliasstate,
               (unsigned long) (child->offset-1),
               inchar,
               indfsstate);
#ifdef SKDEBUG
  apm_showLimdfsstate(stackptr->aliasstate,(unsigned long) child->offset,
                      limdfsresources->dfsconstinfo);
  printf("\n");
#endif
  pprefixlen = apm_limdfsnextstep(stackptr->aliasstate,
                                  limdfsresources->dfsconstinfo);
  if (pprefixlen == 0)
  {
    return true;
  }
  if (pprefixlen > 1UL)
  {
    (limdfsresources->withesa ? esa_overinterval : pck_overinterval)
      (limdfsresources,child,pprefixlen-1);
    return true;
  }
  return false;
}

static void processchildinterval(Limdfsresources *limdfsresources,
                                 const Lcpinterval *child,
                                 Uchar inchar,
                                 const DECLAREPTRDFSSTATE(previousdfsstate))
{
  if (child->left + 1 < child->right || (limdfsresources->withesa &&
                                        child->left + 1 == child->right))
  {
    Lcpintervalwithinfo *stackptr;

    GETNEXTFREEINARRAY(stackptr,&limdfsresources->stack,
                       Lcpintervalwithinfo,128);
    if (pushandpossiblypop(limdfsresources,
                           stackptr,
                           child,
                           inchar,
                           previousdfsstate))
    {
      limdfsresources->stack.nextfreeLcpintervalwithinfo--;
    }
  } else
  {
    if (limdfsresources->withesa)
    {
      esa_overcontext(limdfsresources,
                      previousdfsstate,
                      child->left,
                      child->offset);
    } else
    {
      pck_overcontext(limdfsresources,
                      previousdfsstate,
                      child->left,
                      child->offset,
                      inchar);
    }
  }
}

#ifdef SKDEBUG

static void showLCPinterval(bool withesa,const Lcpinterval *itv)
{
  Seqpos width;

  width = withesa ? (itv->right - itv->left + 1) : itv->right - itv->left;
  printf("(%lu,width=%lu)",(unsigned long) itv->offset,(unsigned long) width);
  printf("(%lu,%lu)",(unsigned long) itv->left,
                     (unsigned long) (withesa ? itv->right : itv->right-1));
}

#endif

static void esa_splitandprocess(Limdfsresources *limdfsresources,
                                const Lcpintervalwithinfo *parentwithinfo)
{
  Seqpos firstnonspecial;
  Uchar extendchar;
  unsigned long idx;
  const Suffixarray *suffixarray
    = (const Suffixarray *) limdfsresources->genericindex;
  const Lcpinterval *parent = &parentwithinfo->lcpitv;

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
    showLCPinterval(limdfsresources->withesa,parent);
    printf(" is ");
    showLCPinterval(limdfsresources->withesa,&child);
    printf("\n");
#endif
    processchildinterval(limdfsresources,
                         &child,
                         inchar,
                         parentwithinfo->aliasstate);
    firstnonspecial = child.right+1;
  }
  if (!limdfsresources->nospecials)
  {
    Seqpos bound;

    for (bound=firstnonspecial; bound <= parent->right; bound++)
    {
      esa_overcontext(limdfsresources,
                      parentwithinfo->aliasstate,
                      bound,
                      parent->offset+1);
    }
  }
}

static void pck_splitandprocess(Limdfsresources *limdfsresources,
                                const Lcpintervalwithinfo *parentwithinfo)
{
  unsigned long idx;
  Seqpos sumwidth = 0;
  const Lcpinterval *parent = &parentwithinfo->lcpitv;

  bwtrangesplitwithoutspecial(&limdfsresources->bwci,
                              limdfsresources->rangeOccs,
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
    sumwidth += child.right - child.left;
#ifdef SKDEBUG
    printf("%u-child of ",(unsigned int) inchar);
    showLCPinterval(limdfsresources->withesa,parent);
    printf(" is ");
    showLCPinterval(limdfsresources->withesa,&child);
    printf("\n");
#endif
    processchildinterval(limdfsresources,
                         &child,
                         inchar,
                         parentwithinfo->aliasstate);
  }
  if (!limdfsresources->nospecials)
  {
    Seqpos bound;
    for (bound = parent->left + sumwidth; bound < parent->right; bound++)
    {
      Uchar cc = bwtseqgetsymbol(bound,limdfsresources->genericindex);

      if (cc != (Uchar) SEPARATOR)
      {
        pck_overcontext(limdfsresources,
                        parentwithinfo->aliasstate,
                        bound,
                        parent->offset+1,
                        cc);
      }
    }
  }
}

#ifdef SKDEBUG
#define SHOWSTACKTOP(STACKPTR)\
        printf("top=");\
        showLCPinterval(limdfsresources->withesa,&STACKPTR->lcpitv);\
        apm_showLimdfsstate(STACKPTR->aliasstate,\
                            (unsigned long) STACKPTR->lcpitv.offset,\
                            limdfsresources->dfsconstinfo);\
        printf("\n")
#else
#define SHOWSTACKTOP(STACKPTR) /* Nothing */
#endif

void indexbasedapproxpatternmatching(Limdfsresources *limdfsresources,
                                     const Uchar *pattern,
                                     unsigned long patternlength,
                                     unsigned long maxdistance)
{
  Lcpintervalwithinfo *stackptr, parentwithinfo;

  printf("deliversizeofdfsstate()=%lu\n",
          (unsigned long) deliversizeofdfsstate());
  printf("sizeof (parentwithinfo.aliasstate)=%lu\n",
          (unsigned long) sizeof (parentwithinfo.aliasstate));
  assert(deliversizeofdfsstate() <= sizeof (parentwithinfo.aliasstate));
  assert(maxdistance < patternlength);
  apm_initdfsconstinfo(limdfsresources->dfsconstinfo,
                       (unsigned int) limdfsresources->alphasize,
                       pattern,patternlength,maxdistance);
  initlcpinfostack(&limdfsresources->stack,
                   0,
                   limdfsresources->withesa
                     ? limdfsresources->totallength
                     : limdfsresources->totallength+1,
                   limdfsresources->dfsconstinfo);
  while (limdfsresources->stack.nextfreeLcpintervalwithinfo > 0)
  {
    assert(limdfsresources->stack.spaceLcpintervalwithinfo != NULL);
    stackptr = limdfsresources->stack.spaceLcpintervalwithinfo +
               limdfsresources->stack.nextfreeLcpintervalwithinfo - 1;
    SHOWSTACKTOP(stackptr);
    parentwithinfo = *stackptr; /* make a copy */
    assert(limdfsresources->stack.nextfreeLcpintervalwithinfo > 0);
    limdfsresources->stack.nextfreeLcpintervalwithinfo--;
    (limdfsresources->withesa ? esa_splitandprocess : pck_splitandprocess)
        (limdfsresources,&parentwithinfo);
  }
}

/*
  Specification steve: only output sequences which occur at most
  t times, where t is a parameter given by the user.
  Output all matches involving a prefix of the pattern and the current
  path with up to k error (k=2 for tagsize around 25).
  Output exakt matching statistics for each suffix of the pattern
*/
