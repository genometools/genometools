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
#include "idx-limdfs.h"

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
  DECLAREDFSSTATE(currentdfsstate);
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
                                    void *processmatchinfo,
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
  limdfsresources->genericindex = genericindex;
  limdfsresources->totallength = totallength;
  limdfsresources->withesa = withesa;
  limdfsresources->nospecials = nospecials;
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
  return limdfsresources;
}

const void *getgenericindexfromresource(Limdfsresources *limdfsresources)
{
  return limdfsresources->genericindex;
}

static void initlcpinfostack(ArrayLcpintervalwithinfo *stack,
                             Seqpos left,
                             Seqpos right,
                             const void *dfsconstinfo,
                             const AbstractDfstransformer *adfst)
{
  Lcpintervalwithinfo *stackptr;

  stack->nextfreeLcpintervalwithinfo = 0;
  GETNEXTFREEINARRAY(stackptr,stack,Lcpintervalwithinfo,128);
  stackptr->lcpitv.offset = 0;
  stackptr->lcpitv.left = left;
  stackptr->lcpitv.right = right;
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

/* iterate myers algorithm over a sequence context */

static void esa_overcontext(Limdfsresources *limdfsresources,
                            const DECLAREPTRDFSSTATE(dfsstate),
                            Seqpos left,
                            Seqpos offset,
                            const AbstractDfstransformer *adfst)
{
  Seqpos pos, startpos;
  Uchar cc;
  unsigned long pprefixlen;
  const Suffixarray *suffixarray
    = (const Suffixarray *) limdfsresources->genericindex;

  memcpy(limdfsresources->currentdfsstate,dfsstate,adfst->sizeofdfsstate);
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
      adfst->inplacenextDfsstate(limdfsresources->dfsconstinfo,
                                 limdfsresources->currentdfsstate,cc);
      pprefixlen = adfst->limdfsnextstep(limdfsresources->currentdfsstate,
                                         (Seqpos) 1,
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
                            Uchar inchar,
                            const AbstractDfstransformer *adfst)
{
  Uchar cc;
  unsigned long contextlength, pprefixlen;
  bool processinchar = true;
  Bwtseqcontextiterator *bsci
    = newBwtseqcontextiterator(limdfsresources->genericindex,left);

  memcpy(limdfsresources->currentdfsstate,dfsstate,adfst->sizeofdfsstate);
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
      adfst->inplacenextDfsstate(limdfsresources->dfsconstinfo,
                                 limdfsresources->currentdfsstate,cc);
      pprefixlen = adfst->limdfsnextstep(limdfsresources->currentdfsstate,
                                         (Seqpos) 1,
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
                      (unsigned long) (child->offset-1),
                      inchar,
                      indfsstate);
#ifdef SKDEBUG
  adfst->showLimdfsstate(stackptr->aliasstate,(unsigned long) child->offset,
                         limdfsresources->dfsconstinfo);
  printf("\n");
#endif
  if (limdfsresources->withesa)
  {
    width = child->right - child->left + 1;
  } else
  {
    width = child->right - child->left;
  }
  pprefixlen = adfst->limdfsnextstep(stackptr->aliasstate,
                                     width,
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
                                 const DECLAREPTRDFSSTATE(previousdfsstate),
                                 const AbstractDfstransformer *adfst)
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
                      previousdfsstate,
                      child->left,
                      child->offset,
                      adfst);
    } else
    {
      pck_overcontext(limdfsresources,
                      previousdfsstate,
                      child->left,
                      child->offset,
                      inchar,
                      adfst);
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
                                const Lcpintervalwithinfo *parentwithinfo,
                                const AbstractDfstransformer *adfst)
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
                         parentwithinfo->aliasstate,
                         adfst);
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
                      parent->offset+1,
                      adfst);
    }
  }
}

static void pck_splitandprocess(Limdfsresources *limdfsresources,
                                const Lcpintervalwithinfo *parentwithinfo,
                                const AbstractDfstransformer *adfst)
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
                         parentwithinfo->aliasstate,
                         adfst);
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
                                     const Uchar *pattern,
                                     unsigned long patternlength,
                                     unsigned long maxdistance,
                                     Seqpos maxintervalwidth,
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
                          maxintervalwidth);
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
    assert(limdfsresources->stack.nextfreeLcpintervalwithinfo > 0);
    limdfsresources->stack.nextfreeLcpintervalwithinfo--;
    (limdfsresources->withesa ? esa_splitandprocess : pck_splitandprocess)
        (limdfsresources,&parentwithinfo,adfst);
  }
}

/*
  Specification steve: only output sequences which occur at most
  t times, where t is a parameter given by the user.
  Output all matches involving a prefix of the pattern and the current
  path with up to k error (k=2 for tagsize around 25).
  Output exakt matching statistics for each suffix of the pattern
*/
