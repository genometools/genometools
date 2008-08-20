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

#include <string.h>
#include "libgtcore/chardef.h"
#include "libgtcore/symboldef.h"
#include "libgtcore/unused.h"
#include "libgtcore/arraydef.h"
#include "sarr-def.h"
#include "seqpos-def.h"
#include "spacedef.h"
#include "divmodmul.h"
#include "esa-splititv.h"
#include "eis-voiditf.h"
#include "defined-types.h"
#include "esa-mmsearch-def.h"
#include "absdfstrans-imp.h"
#include "idx-limdfs.h"
#include "esa-minunique.pr"

#ifdef SKDEBUG
#define DECLAREDFSSTATE(V)\
        Aliasdfsstate V[4]
#else
#define DECLAREDFSSTATE(V)\
        Aliasdfsstate V[3]
#endif

typedef struct
{
  Seqpos offset,
         leftbound,
         rightbound;
  Codetype code;
  Uchar inchar;
} Indexbounds;

typedef struct
{
  Indexbounds lcpitv;
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
  Processmatch processmatch;
  void *processmatchinfo;
  Processresult processresult;
  void *patterninfo;
  const void *genericindex;
  bool withesa, nowildcards;
  unsigned long maxintervalwidth;
  DECLAREDFSSTATE(currentdfsstate);
  Seqpos *rangeOccs;
  unsigned int maxdepth;
  const Matchbound **mbtab;
  const Encodedsequence *encseq;
  ArraySeqpos mstatspos;
  Uchar *currentpathspace;
  unsigned long maxpathlength;
};

Limdfsresources *newLimdfsresources(const void *genericindex,
                                    const Matchbound **mbtab,
                                    unsigned int maxdepth,
                                    const Encodedsequence *encseq,
                                    bool withesa,
                                    bool nowildcards,
                                    unsigned long maxintervalwidth,
                                    unsigned int mapsize,
                                    Seqpos totallength,
                                    unsigned long maxpathlength,
                                    Processmatch processmatch,
                                    void *processmatchinfo,
                                    Processresult processresult,
                                    void *patterninfo,
                                    const AbstractDfstransformer *adfst)
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
  limdfsresources->processresult = processresult;
  limdfsresources->patterninfo = patterninfo;
  limdfsresources->genericindex = genericindex;
  limdfsresources->totallength = totallength;
  limdfsresources->withesa = withesa;
  limdfsresources->nowildcards = nowildcards;
  limdfsresources->encseq = encseq;
  limdfsresources->maxdepth = maxdepth;
  limdfsresources->mbtab = mbtab;
  limdfsresources->maxintervalwidth = maxintervalwidth;
  limdfsresources->maxpathlength = maxpathlength;
  if (maxpathlength > 0)
  {
    ALLOCASSIGNSPACE(limdfsresources->currentpathspace,NULL,Uchar,
                     maxpathlength);
  } else
  {
    limdfsresources->currentpathspace = NULL;
  }
  /* Application specific */
  limdfsresources->dfsconstinfo
    = adfst->allocatedfsconstinfo((unsigned int) limdfsresources->alphasize);
  if (withesa)
  {
    limdfsresources->rangeOccs = NULL;
  } else
  {
    ALLOCASSIGNSPACE(limdfsresources->rangeOccs,NULL,Seqpos,
                     MULT2(limdfsresources->alphasize));
  }
  INITARRAY(&limdfsresources->mstatspos,Seqpos);
  if (maxintervalwidth > 0)
  {
    ALLOCASSIGNSPACE(limdfsresources->mstatspos.spaceSeqpos,NULL,Seqpos,
                     maxintervalwidth);
    limdfsresources->mstatspos.allocatedSeqpos = maxintervalwidth;
  }
  return limdfsresources;
}

static void initlcpinfostack(ArrayLcpintervalwithinfo *stack,
                             Seqpos leftbound,
                             Seqpos rightbound,
                             void *dfsconstinfo,
                             const AbstractDfstransformer *adfst)
{
  Lcpintervalwithinfo *stackptr;

  stack->nextfreeLcpintervalwithinfo = 0;
  GETNEXTFREEINARRAY(stackptr,stack,Lcpintervalwithinfo,128);
  stackptr->lcpitv.offset = 0;
  stackptr->lcpitv.leftbound = leftbound;
  stackptr->lcpitv.rightbound = rightbound;
  stackptr->lcpitv.code = 0;
  adfst->initLimdfsstate(stackptr->aliasstate,dfsconstinfo);
}

void freeLimdfsresources(Limdfsresources **ptrlimdfsresources,
                         const AbstractDfstransformer *adfst)
{
  Limdfsresources *limdfsresources = *ptrlimdfsresources;

  adfst->freedfsconstinfo(&limdfsresources->dfsconstinfo);
  FREEARRAY(&limdfsresources->bwci,Boundswithchar);
  FREEARRAY(&limdfsresources->stack,Lcpintervalwithinfo);
  FREESPACE(limdfsresources->rangeOccs);
  FREESPACE(limdfsresources->currentpathspace);
  FREEARRAY(&limdfsresources->mstatspos,Seqpos);
  FREESPACE(*ptrlimdfsresources);
}

/* enumerate the suffixes in an LCP-interval */

static void gen_esa_overinterval(const void *voidsuffixarray,
                                 Processmatch processmatch,
                                 void *processmatchinfo,
                                 bool rcmatch,
                                 const Indexbounds *itv,
                                 unsigned long pprefixlen,
                                 UNUSED Seqpos totallength,
                                 const Uchar *dbsubstring)
{
  const Suffixarray *suffixarray = (const Suffixarray *) voidsuffixarray;
  Seqpos idx;

  for (idx = itv->leftbound; idx <= itv->rightbound; idx++)
  {
    processmatch(processmatchinfo,
                 rcmatch,
                 suffixarray->suftab[idx],
                 itv->offset,
                 dbsubstring,
                 pprefixlen);
  }
}

static void esa_overinterval(Limdfsresources *limdfsresources,
                             bool rcmatch,
                             const Indexbounds *itv,
                             unsigned long pprefixlen)
{
  gen_esa_overinterval((const Suffixarray *) limdfsresources->genericindex,
                       limdfsresources->processmatch,
                       limdfsresources->processmatchinfo,
                       rcmatch,
                       itv,
                       pprefixlen,
                       limdfsresources->totallength,
                       limdfsresources->currentpathspace);
}

static void gen_pck_overinterval(const void *voidbwtseq,
                                 Processmatch processmatch,
                                 void *processmatchinfo,
                                 bool rcmatch,
                                 const Indexbounds *itv,
                                 unsigned long pprefixlen,
                                 Seqpos totallength,
                                 const Uchar *dbsubstring)
{
  Bwtseqpositioniterator *bspi;
  Seqpos dbstartpos;

  assert(itv->leftbound < itv->rightbound);
  bspi = newBwtseqpositioniterator (voidbwtseq,itv->leftbound,itv->rightbound);
  while (nextBwtseqpositioniterator(&dbstartpos,bspi))
  {
    assert(totallength >= (dbstartpos + itv->offset));
    processmatch(processmatchinfo,
                 rcmatch,
                 totallength - (dbstartpos + itv->offset),
                 itv->offset,
                 dbsubstring,
                 pprefixlen);
  }
  freeBwtseqpositioniterator(&bspi);
}

static void pck_overinterval(Limdfsresources *limdfsresources,
                             bool rcmatch,
                             const Indexbounds *itv,
                             unsigned long pprefixlen)
{
  gen_pck_overinterval(limdfsresources->genericindex,
                       limdfsresources->processmatch,
                       limdfsresources->processmatchinfo,
                       rcmatch,
                       itv,
                       pprefixlen,
                       limdfsresources->totallength,
                       limdfsresources->currentpathspace);
}

static void storemstatsposition(void *processinfo,
                                UNUSED bool rcmatch,
                                Seqpos dbstartpos,
                                UNUSED Seqpos dblen,
                                UNUSED const Uchar *dbsubstring,
                                UNUSED unsigned long pprefixlen)
{
  ArraySeqpos *mstatspos = (ArraySeqpos *) processinfo;

  STOREINARRAY(mstatspos,Seqpos,32,dbstartpos);
}

static int comparepositions(const void *a, const void *b)
{
  if (*((Seqpos *) a) < *((Seqpos *) b))
  {
    return -1;
  }
  return 1;
}

ArraySeqpos *fromitv2sortedmatchpositions(Limdfsresources *limdfsresources,
                                          bool rcmatch,
                                          Seqpos leftbound,
                                          Seqpos rightbound,
                                          unsigned long offset)
{
  Indexbounds itv;

  limdfsresources->mstatspos.nextfreeSeqpos = 0;
  itv.leftbound = leftbound;
  itv.rightbound = rightbound;
  itv.offset = (Seqpos) offset;
  (limdfsresources->withesa ? gen_esa_overinterval : gen_pck_overinterval)
    (limdfsresources->genericindex,
     storemstatsposition,
     &limdfsresources->mstatspos,
     rcmatch,
     &itv,
     offset,
     limdfsresources->totallength,
     limdfsresources->currentpathspace);
  qsort(limdfsresources->mstatspos.spaceSeqpos,
        (size_t) limdfsresources->mstatspos.nextfreeSeqpos,
        sizeof (Seqpos), comparepositions);
  return &limdfsresources->mstatspos;
}

/* iterate myers algorithm over a sequence context */

static void esa_overcontext(Limdfsresources *limdfsresources,
                            bool rcmatch,
                            const DECLAREPTRDFSSTATE(dfsstate),
                            Seqpos leftbound,
                            Seqpos offset,
                            const AbstractDfstransformer *adfst)
{
  Seqpos pos, startpos;
  Uchar cc;
  unsigned long pprefixlen;
  const Suffixarray *suffixarray
    = (const Suffixarray *) limdfsresources->genericindex;

  memcpy(limdfsresources->currentdfsstate,dfsstate,adfst->sizeofdfsstate);
  startpos = suffixarray->suftab[leftbound];
#ifdef SKDEBUG
  printf("retrieve context of startpos=%lu\n",(unsigned long) startpos);
#endif
  for (pos = startpos + offset - 1; pos < limdfsresources->totallength; pos++)
  {
    cc = getencodedchar(suffixarray->encseq,pos,suffixarray->readmode);
    if (cc != (Uchar) SEPARATOR &&
        (!limdfsresources->nowildcards || cc != (Uchar) WILDCARD))
    {
#ifdef SKDEBUG
      printf("cc=%u\n",(unsigned int) cc);
#endif
      adfst->inplacenextDfsstate(limdfsresources->dfsconstinfo,
                                 limdfsresources->currentdfsstate,
                                 (unsigned long) (pos - startpos + 1),
                                 cc);
      pprefixlen = adfst->limdfsnextstep(limdfsresources->currentdfsstate,
                                         leftbound,
                                         leftbound,
                                         (Seqpos) 1,
                                         (unsigned long) (pos - startpos + 1),
                                         limdfsresources->dfsconstinfo);
      if (pprefixlen == 0) /* failure */
      {
        break;
      }
      if (pprefixlen > 1UL)
      {
        limdfsresources->processmatch(limdfsresources->processmatchinfo,
                                      rcmatch,
                                      startpos,
                                      pos - startpos + 1,
                                      limdfsresources->currentpathspace,
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
                            bool rcmatch,
                            const DECLAREPTRDFSSTATE(dfsstate),
                            Seqpos leftbound,
                            Seqpos offset,
                            Uchar inchar,
                            const AbstractDfstransformer *adfst)
{
  Uchar cc;
  unsigned long contextlength, pprefixlen;
  bool processinchar = true;
  Seqpos bound = leftbound;
  Bwtseqcontextiterator *bsci;

  bsci = newBwtseqcontextiterator(limdfsresources->genericindex,bound);
  memcpy(limdfsresources->currentdfsstate,dfsstate,adfst->sizeofdfsstate);
#ifdef SKDEBUG
  printf("retrieve context for bound = %lu\n",(unsigned long) bound);
#endif
  for (contextlength = 0; /* nothing */; contextlength++)
  {
    if (processinchar)
    {
      cc = inchar;
      processinchar = false;
    } else
    {
      cc = nextBwtseqcontextiterator(&bound,bsci);
    }
    if (cc != (Uchar) SEPARATOR &&
        (!limdfsresources->nowildcards || cc != (Uchar) WILDCARD))
    {
#ifdef SKDEBUG
      printf("cc=%u\n",(unsigned int) cc);
#endif
      assert(offset - 1 + contextlength
             < (Seqpos) limdfsresources->maxpathlength);
      limdfsresources->currentpathspace[offset-1+contextlength] = cc;
      adfst->inplacenextDfsstate(limdfsresources->dfsconstinfo,
                                 limdfsresources->currentdfsstate,
                                 (unsigned long) (offset + contextlength),
                                 cc);
      pprefixlen = adfst->limdfsnextstep(limdfsresources->currentdfsstate,
                                         bound,
                                         bound+1,
                                         (Seqpos) 1,
                                         (unsigned long) (offset+contextlength),
                                         limdfsresources->dfsconstinfo);
      if (pprefixlen == 0)
      {
        break;
      }
      if (pprefixlen > 1UL) /* check for success */
      {
        Seqpos startpos = bwtseqfirstmatch(limdfsresources->genericindex,
                                           leftbound);
        limdfsresources->processmatch(limdfsresources->processmatchinfo,
                                      rcmatch,
                                      limdfsresources->totallength -
                                           (startpos + offset),
                                      offset + contextlength,
                                      limdfsresources->currentpathspace,
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
                               bool rcmatch,
                               Lcpintervalwithinfo *stackptr,
                               const Indexbounds *child,
                               Uchar inchar,
                               const DECLAREPTRDFSSTATE(indfsstate),
                               const AbstractDfstransformer *adfst)
{
  unsigned long pprefixlen;
  Seqpos width;
  stackptr->lcpitv = *child;

#ifdef SKDEBUG
  printf("(2) nextDfsstate(");
  adfst->showLimdfsstate(indfsstate,(unsigned long) (child->offset-1),
                         limdfsresources->dfsconstinfo);
  printf(",%u)=",(unsigned int) inchar);
#endif
  adfst->nextDfsstate(limdfsresources->dfsconstinfo,
                      stackptr->aliasstate,
                      (unsigned long) child->offset,
                      inchar,
                      indfsstate);
#ifdef SKDEBUG
  adfst->showLimdfsstate(stackptr->aliasstate,(unsigned long) child->offset,
                         limdfsresources->dfsconstinfo);
  printf("\n");
#endif
  if (limdfsresources->withesa)
  {
    width = child->rightbound - child->leftbound + 1;
  } else
  {
    width = child->rightbound - child->leftbound;
  }
  pprefixlen = adfst->limdfsnextstep(stackptr->aliasstate,
                                     child->leftbound,
                                     child->rightbound,
                                     width,
                                     (unsigned long) child->offset,
                                     limdfsresources->dfsconstinfo);
  if (pprefixlen == 0)
  {
    return true;
  }
  if (pprefixlen > 1UL)
  {
    (limdfsresources->withesa ? esa_overinterval : pck_overinterval)
      (limdfsresources,rcmatch,child,pprefixlen-1);
    return true;
  }
  return false;
}

static void processchildinterval(Limdfsresources *limdfsresources,
                                 bool rcmatch,
                                 const Indexbounds *child,
                                 Uchar inchar,
                                 const DECLAREPTRDFSSTATE(previousdfsstate),
                                 const AbstractDfstransformer *adfst)
{
  if (child->leftbound + 1 < child->rightbound ||
      (limdfsresources->withesa && child->leftbound + 1 == child->rightbound))
  {
    Lcpintervalwithinfo *stackptr;

    GETNEXTFREEINARRAY(stackptr,&limdfsresources->stack,
                       Lcpintervalwithinfo,128);
    if (pushandpossiblypop(limdfsresources,
                           rcmatch,
                           stackptr,
                           child,
                           inchar,
                           previousdfsstate,
                           adfst))
    {
      limdfsresources->stack.nextfreeLcpintervalwithinfo--;
    }
  } else
  {
    if (limdfsresources->withesa)
    {
      esa_overcontext(limdfsresources,
                      rcmatch,
                      previousdfsstate,
                      child->leftbound,
                      child->offset,
                      adfst);
    } else
    {
      pck_overcontext(limdfsresources,
                      rcmatch,
                      previousdfsstate,
                      child->leftbound,
                      child->offset,
                      inchar,
                      adfst);
    }
  }
}

#ifdef SKDEBUG

static void showLCPinterval(bool withesa,const Indexbounds *itv)
{
  Seqpos width;

  width = withesa ? (itv->rightbound - itv->leftbound + 1)
                  : itv->rightbound - itv->leftbound;
  printf("(%lu,width=%lu)",(unsigned long) itv->offset,(unsigned long) width);
  printf("(%lu,%lu)",(unsigned long) itv->leftbound,
                     (unsigned long)
                     (withesa ? itv->rightbound : itv->rightbound-1));
}

#endif

static void esa_splitandprocess(Limdfsresources *limdfsresources,
                                bool rcmatch,
                                const Lcpintervalwithinfo *parentwithinfo,
                                const AbstractDfstransformer *adfst)
{
  Seqpos firstnonspecial;
  Uchar extendchar;
  unsigned long idx;
  const Suffixarray *suffixarray
    = (const Suffixarray *) limdfsresources->genericindex;
  const Indexbounds *parent = &parentwithinfo->lcpitv;

  extendchar = lcpintervalextendlcp(suffixarray->encseq,
                                    suffixarray->readmode,
                                    suffixarray->suftab,
                                    limdfsresources->totallength,
                                    limdfsresources->alphasize,
                                    parent->offset,
                                    parent->leftbound,
                                    parent->rightbound);
  if (extendchar < limdfsresources->alphasize)
  {
    limdfsresources->bwci.spaceBoundswithchar[0].lbound = parent->leftbound;
    limdfsresources->bwci.spaceBoundswithchar[0].rbound = parent->rightbound;
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
                                   parent->offset,
                                   parent->leftbound,
                                   parent->rightbound);
  }
  firstnonspecial = parent->leftbound;
  for (idx = 0; idx < limdfsresources->bwci.nextfreeBoundswithchar; idx++)
  {
    Indexbounds child;
    Uchar inchar = limdfsresources->bwci.spaceBoundswithchar[idx].inchar;

    child.offset = parent->offset+1;
    child.leftbound = limdfsresources->bwci.spaceBoundswithchar[idx].lbound;
    child.rightbound = limdfsresources->bwci.spaceBoundswithchar[idx].rbound;
    child.code = 0; /* not used, but we better defined it */
    child.inchar = inchar;
#ifdef SKDEBUG
    printf("%u-child of ",(unsigned int) inchar);
    showLCPinterval(limdfsresources->withesa,parent);
    printf(" is ");
    showLCPinterval(limdfsresources->withesa,&child);
    printf("\n");
#endif
    processchildinterval(limdfsresources,
                         rcmatch,
                         &child,
                         inchar,
                         parentwithinfo->aliasstate,
                         adfst);
    firstnonspecial = child.rightbound+1;
  }
  if (!limdfsresources->nowildcards)
  {
    Seqpos bound;

    for (bound=firstnonspecial; bound <= parent->rightbound; bound++)
    {
      esa_overcontext(limdfsresources,
                      rcmatch,
                      parentwithinfo->aliasstate,
                      bound,
                      parent->offset+1,
                      adfst);
    }
  }
}

static void smalldepthbwtrangesplitwithoutspecial(ArrayBoundswithchar *bwci,
                                                  const Matchbound **mbtab,
                                                  Uchar alphasize,
                                                  Codetype parentcode,
                                                  unsigned long childdepth)
{
  Codetype childcode;
  const Matchbound *mbptr;

  assert(childdepth > 0);
  bwci->nextfreeBoundswithchar = 0;
  childcode = parentcode * alphasize;
  for (mbptr = mbtab[childdepth] + childcode;
       mbptr < mbtab[childdepth] + childcode + alphasize;
       mbptr++)
  {
    if (mbptr->lowerbound < mbptr->upperbound)
    {
      bwci->spaceBoundswithchar[bwci->nextfreeBoundswithchar].inchar
        = (Uchar) (mbptr - (mbtab[childdepth] + childcode));
      bwci->spaceBoundswithchar[bwci->nextfreeBoundswithchar].lbound
        = mbptr->lowerbound;
      bwci->spaceBoundswithchar[bwci->nextfreeBoundswithchar++].rbound
        = mbptr->upperbound;
    }
  }
}

static void pck_splitandprocess(Limdfsresources *limdfsresources,
                                bool rcmatch,
                                const Lcpintervalwithinfo *parentwithinfo,
                                const AbstractDfstransformer *adfst)
{
  unsigned long idx;
  Seqpos sumwidth = 0;
  const Indexbounds *parent = &parentwithinfo->lcpitv;
  Codetype startcode;

  if (parent->offset < (Seqpos) limdfsresources->maxdepth)
  {
    smalldepthbwtrangesplitwithoutspecial(&limdfsresources->bwci,
                                          limdfsresources->mbtab,
                                          limdfsresources->alphasize,
                                          parent->code,
                                          (unsigned long) (parent->offset + 1));
    startcode = parent->code * limdfsresources->alphasize;
  } else
  {
    bwtrangesplitwithoutspecial(&limdfsresources->bwci,
                                limdfsresources->rangeOccs,
                                limdfsresources->genericindex,
                                parent->leftbound,
                                parent->rightbound);
    startcode = 0;
  }
  for (idx = 0; idx < limdfsresources->bwci.nextfreeBoundswithchar; idx++)
  {
    Uchar inchar;
    Indexbounds child;

    inchar = limdfsresources->bwci.spaceBoundswithchar[idx].inchar;
    child.offset = parent->offset+1;
    child.leftbound = limdfsresources->bwci.spaceBoundswithchar[idx].lbound;
    child.rightbound = limdfsresources->bwci.spaceBoundswithchar[idx].rbound;
    assert(inchar < limdfsresources->alphasize);
    child.code = startcode + inchar;
    child.inchar = inchar;
    sumwidth += child.rightbound - child.leftbound;
#ifdef SKDEBUG
    printf("%u-child of ",(unsigned int) inchar);
    showLCPinterval(limdfsresources->withesa,parent);
    printf(" is ");
    showLCPinterval(limdfsresources->withesa,&child);
    printf("\n");
#endif
    processchildinterval(limdfsresources,
                         rcmatch,
                         &child,
                         inchar,
                         parentwithinfo->aliasstate,
                         adfst);
  }
  if (!limdfsresources->nowildcards)
  {
    Seqpos bound;
    for (bound = parent->leftbound + sumwidth;
         bound < parent->rightbound; bound++)
    {
      Uchar cc = bwtseqgetsymbol(bound,limdfsresources->genericindex);

      if (cc != (Uchar) SEPARATOR)
      {
        pck_overcontext(limdfsresources,
                        rcmatch,
                        parentwithinfo->aliasstate,
                        bound,
                        parent->offset+1,
                        cc,
                        adfst);
      }
    }
  }
}

#ifdef SKDEBUG
#define SHOWSTACKTOP(STACKPTR)\
        printf("top=");\
        showLCPinterval(limdfsresources->withesa,&STACKPTR->lcpitv);\
        adfst->showLimdfsstate(STACKPTR->aliasstate,\
                               (unsigned long) STACKPTR->lcpitv.offset,\
                               limdfsresources->dfsconstinfo);\
        printf("\n")
#else
#define SHOWSTACKTOP(STACKPTR) /* Nothing */
#endif

void indexbasedapproxpatternmatching(Limdfsresources *limdfsresources,
                                     bool rcmatch,
                                     const Uchar *pattern,
                                     unsigned long patternlength,
                                     unsigned long maxdistance,
                                     unsigned long maxintervalwidth,
                                     bool skpp,
                                     const AbstractDfstransformer *adfst)
{
  Lcpintervalwithinfo *stackptr, parentwithinfo;

  /*
  printf("sizeofdfsstate()=%lu\n",(unsigned long) adfst->sizeofdfsstate);
  printf("sizeof (parentwithinfo.aliasstate)=%lu\n",
          (unsigned long) sizeof (parentwithinfo.aliasstate));
  */
  assert(adfst->sizeofdfsstate <= sizeof (parentwithinfo.aliasstate));
  assert(maxdistance < patternlength);
  adfst->initdfsconstinfo(limdfsresources->dfsconstinfo,
                          (unsigned int) limdfsresources->alphasize,
                          pattern,patternlength,maxdistance,
                          maxintervalwidth,
                          skpp,
                          false);
  initlcpinfostack(&limdfsresources->stack,
                   0,
                   limdfsresources->withesa
                     ? limdfsresources->totallength
                     : limdfsresources->totallength+1,
                   limdfsresources->dfsconstinfo,
                   adfst);
  while (limdfsresources->stack.nextfreeLcpintervalwithinfo > 0)
  {
    assert(limdfsresources->stack.spaceLcpintervalwithinfo != NULL);
    stackptr = limdfsresources->stack.spaceLcpintervalwithinfo +
               limdfsresources->stack.nextfreeLcpintervalwithinfo - 1;
    SHOWSTACKTOP(stackptr);
    parentwithinfo = *stackptr; /* make a copy */
    if (parentwithinfo.lcpitv.offset > 0)
    {
      assert(parentwithinfo.lcpitv.offset-1 <
             (Seqpos) limdfsresources->maxpathlength);
      limdfsresources->currentpathspace[parentwithinfo.lcpitv.offset-1]
                     = parentwithinfo.lcpitv.inchar;
    }
    assert(limdfsresources->stack.nextfreeLcpintervalwithinfo > 0);
    limdfsresources->stack.nextfreeLcpintervalwithinfo--;
    (limdfsresources->withesa ? esa_splitandprocess : pck_splitandprocess)
        (limdfsresources,rcmatch,&parentwithinfo,adfst);
  }
  if (adfst->extractdfsconstinfo != NULL)
  {
    adfst->extractdfsconstinfo(limdfsresources->processresult,
                               limdfsresources,
                               limdfsresources->patterninfo,
                               limdfsresources->dfsconstinfo);
  }
}

unsigned long genericmstats(const Limdfsresources *limdfsresources,
                            const Uchar *qstart,
                            const Uchar *qend)
{
  return (limdfsresources->withesa ?
             suffixarraymstats : voidpackedindexmstatsforward)
                     (limdfsresources->genericindex,
                      0,
                      0,
                      limdfsresources->totallength,
                      NULL,
                      qstart,
                      qend);
}

static void esa_exactpatternmatching(const void *genericindex,
                                     bool rcmatch,
                                     const Uchar *pattern,
                                     unsigned long patternlength,
                                     const Uchar *dbsubstring,
                                     Processmatch processmatch,
                                     void *processmatchinfo)
{
  const Suffixarray *suffixarray = (const Suffixarray *) genericindex;
  MMsearchiterator *mmsi;
  Seqpos dbstartpos, totallength = getencseqtotallength(suffixarray->encseq);

  mmsi = newmmsearchiterator(suffixarray->encseq,
                             suffixarray->suftab,
                             0,  /* leftbound */
                             totallength, /* rightbound */
                             0, /* offset */
                             suffixarray->readmode,
                             pattern,
                             patternlength);
  while (nextmmsearchiterator(&dbstartpos,mmsi))
  {
    processmatch(processmatchinfo,rcmatch,dbstartpos,
                 (Seqpos) patternlength,dbsubstring,patternlength);
  }
  freemmsearchiterator(&mmsi);
}

void indexbasedexactpatternmatching(const Limdfsresources *limdfsresources,
                                    bool rcmatch,
                                    const Uchar *pattern,
                                    unsigned long patternlength,
                                    Processmatch processmatch,
                                    void *processmatchinfo)
{
  if (limdfsresources->withesa)
  {
    esa_exactpatternmatching(limdfsresources->genericindex,
                             rcmatch,
                             pattern,
                             patternlength,
                             limdfsresources->currentpathspace,
                             processmatch,
                             processmatchinfo);
  } else
  {
    pck_exactpatternmatching(limdfsresources->genericindex,
                             rcmatch,
                             pattern,
                             patternlength,
                             limdfsresources->totallength,
                             limdfsresources->currentpathspace,
                             processmatch,
                             processmatchinfo);
  }
}

Seqpos bound2startpos(const Limdfsresources *limdfsresources,
                      Seqpos bound,unsigned long matchlength)
{
  if (limdfsresources->withesa)
  {
    return ((const Suffixarray *) limdfsresources->genericindex)->suftab[bound];
  }
  return voidpackedfindfirstmatchconvert(limdfsresources->genericindex,
                                         bound,matchlength);
}

Uchar limdfsgetencodedchar(const Limdfsresources *limdfsresources,
                           Seqpos pos,
                           Readmode readmode)
{
  assert(limdfsresources->encseq != NULL);

  return getencodedchar(limdfsresources->encseq,pos,readmode);
}

Seqpos getlastbound(const Limdfsresources *limdfsresources,Seqpos rightbound)
{
  if (limdfsresources->withesa)
  {
    return rightbound;
  }
  return rightbound - 1;
}

bool intervalwidthleq(const Limdfsresources *limdfsresources,
                      Seqpos leftbound,Seqpos rightbound)
{
  Seqpos width;

  if (limdfsresources->withesa)
  {
    if (leftbound > rightbound)
    {
      width = 0;
    } else
    {
      width = rightbound - leftbound + 1;
    }
  } else
  {
    if (leftbound >= rightbound)
    {
      width = 0;
    } else
    {
      width = rightbound - leftbound;
    }
  }
  if (width > 0 && width <= (Seqpos) limdfsresources->maxintervalwidth)
  {
    return true;
  }
  return false;
}

/*
  Specification steve: only output sequences which occur at most
  t times, where t is a parameter given by the user.
  Output all matches involving a prefix of the pattern and the current
  path with up to k error (k=2 for tagsize around 25).
  Output exakt matching statistics for each suffix of the pattern
*/
