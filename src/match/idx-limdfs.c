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
#include "core/chardef.h"
#include "core/symboldef.h"
#include "core/unused_api.h"
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
#include "stamp.h"
#include "esa-map.h"
#include "esa-minunique.pr"

#define DECLAREDFSSTATE(V)\
        Aliasdfsstate V[5]

struct Genericindex
{
  Suffixarray *suffixarray;
  Seqpos totallength;
  void *packedindex;
  bool withesa;
  const Mbtab **mbtab;      /* only relevant for packedindex */
  unsigned int maxdepth;    /* maximaldepth of boundaries */
};

void genericindex_delete(Genericindex *genericindex)
{
  if (genericindex == NULL)
  {
    return;
  }
  freesuffixarray(genericindex->suffixarray);
  gt_free(genericindex->suffixarray);
  if (genericindex->packedindex != NULL)
  {
    deletevoidBWTSeq(genericindex->packedindex);
  }
  gt_free(genericindex);
}

const Encodedsequence *genericindex_getencseq(const Genericindex *genericindex)
{
  return genericindex->suffixarray->encseq;
}

Genericindex *genericindex_new(const GtStr *indexname,
                               bool withesa,
                               bool withencseq,
                               int userdefinedmaxdepth,
                               Verboseinfo *verboseinfo,
                               GtError *err)
{
  unsigned int demand = 0;
  bool haserr = false;
  Genericindex *genericindex;

  genericindex = gt_malloc(sizeof(*genericindex));
  if (withesa)
  {
    demand |= SARR_SUFTAB;
  }
  if (withencseq)
  {
    demand |= SARR_ESQTAB;
  }
  genericindex->withesa = withesa;
  genericindex->suffixarray = gt_malloc(sizeof(*genericindex->suffixarray));
  if (mapsuffixarray(genericindex->suffixarray,
                     demand,
                     indexname,
                     verboseinfo,
                     err) != 0)
  {
    haserr = true;
    genericindex->totallength = 0;
  } else
  {
    genericindex->totallength = getencseqtotallength(genericindex->suffixarray
                                                                 ->encseq);
  }
  if (!haserr)
  {
    if (withesa && genericindex->suffixarray->readmode != Forwardmode)
    {
      gt_error_set(err,"using option -esa you can only process index "
                       "in forward mode");
      haserr = true;
    } else
    {
      if (!withesa && genericindex->suffixarray->readmode != Reversemode)
      {
        gt_error_set(err,"with option -pck you can only process index "
                         "in reverse mode");
        haserr = true;
      }
    }
  }
  genericindex->packedindex = NULL;
  genericindex->mbtab = NULL;
  genericindex->maxdepth = 0;
  if (!haserr && !withesa)
  {
    genericindex->packedindex = loadvoidBWTSeqForSA(indexname,
                                                    genericindex->suffixarray,
                                                    genericindex->totallength,
                                                    true, err);
    if (genericindex->packedindex == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr && !withesa)
  {
    genericindex->mbtab = bwtseq2mbtab(genericindex->packedindex);
    genericindex->maxdepth = bwtseq2maxdepth(genericindex->packedindex);
    if (userdefinedmaxdepth >= 0 &&
        genericindex->maxdepth > (unsigned int) userdefinedmaxdepth)
    {
      genericindex->maxdepth = (unsigned int) userdefinedmaxdepth;
    }
  }
  if (haserr)
  {
    genericindex_delete(genericindex);
    return NULL;
  }
  return genericindex;
}

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
  Lcpintervalwithinfo parentwithinfo;
  Uchar alphasize;
  Processmatch processmatch;
  void *processmatchinfo;
  Processresult processresult;
  void *patterninfo;
  const Genericindex *genericindex;
  bool nowildcards;
  unsigned long maxintervalwidth;
  DECLAREDFSSTATE(currentdfsstate);
  Seqpos *rangeOccs;
  const Encodedsequence *encseq;
  ArraySeqpos mstatspos;
  Uchar *currentpathspace;
  unsigned long maxpathlength;
  Seqpos numberofmatches;
};

Limdfsresources *newLimdfsresources(const Genericindex *genericindex,
                                    bool nowildcards,
                                    unsigned long maxintervalwidth,
                                    unsigned long maxpathlength,
                                    Processmatch processmatch,
                                    void *processmatchinfo,
                                    Processresult processresult,
                                    void *patterninfo,
                                    const AbstractDfstransformer *adfst)
{
  Limdfsresources *limdfsresources;
  unsigned int numofchars;
  const Encodedsequence *encseq;

  encseq = genericindex->suffixarray->encseq;
  numofchars = getencseqAlphabetnumofchars(encseq);
  ALLOCASSIGNSPACE(limdfsresources,NULL,Limdfsresources,1);
  ALLOCASSIGNSPACE(limdfsresources->bwci.spaceBoundswithchar,NULL,
                   Boundswithchar,numofchars+1);
  limdfsresources->bwci.nextfreeBoundswithchar = 0;
  limdfsresources->bwci.allocatedBoundswithchar
    = (unsigned long) (numofchars+1);
  INITARRAY(&limdfsresources->stack,Lcpintervalwithinfo);
  gt_assert(numofchars <= MAXALPHABETCHARACTER);
  limdfsresources->alphasize = (Uchar) numofchars;
  limdfsresources->processmatch = processmatch;
  limdfsresources->processmatchinfo = processmatchinfo;
  limdfsresources->processresult = processresult;
  limdfsresources->patterninfo = patterninfo;
  limdfsresources->genericindex = genericindex;
  limdfsresources->nowildcards = nowildcards;
  limdfsresources->encseq = encseq;
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
  if (genericindex->withesa)
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
  if (adfst->initLimdfsstackelem != NULL)
  {
    adfst->initLimdfsstackelem(limdfsresources->currentdfsstate);
    adfst->initLimdfsstackelem(limdfsresources->parentwithinfo.aliasstate);
  }
  return limdfsresources;
}

static Lcpintervalwithinfo *allocateStackspace(ArrayLcpintervalwithinfo *stack,
                                               const AbstractDfstransformer
                                                      *adfst)
{
  if (stack->nextfreeLcpintervalwithinfo >= stack->allocatedLcpintervalwithinfo)
  {
    unsigned long idx;
    const unsigned long addelems = 128UL;

    ALLOCASSIGNSPACE(stack->spaceLcpintervalwithinfo,
                     stack->spaceLcpintervalwithinfo,
                     Lcpintervalwithinfo,
                     stack->allocatedLcpintervalwithinfo + addelems);
    if (adfst->initLimdfsstackelem != NULL)
    {
      for (idx = stack->allocatedLcpintervalwithinfo;
           idx < stack->allocatedLcpintervalwithinfo + addelems; idx++)
      {
        adfst->initLimdfsstackelem(stack->spaceLcpintervalwithinfo[idx].
                                   aliasstate);
      }
    }
    stack->allocatedLcpintervalwithinfo += addelems;
  }
  return stack->spaceLcpintervalwithinfo + stack->nextfreeLcpintervalwithinfo++;
}

static void initlcpinfostack(ArrayLcpintervalwithinfo *stack,
                             Seqpos leftbound,
                             Seqpos rightbound,
                             void *dfsconstinfo,
                             const AbstractDfstransformer *adfst)
{
  Lcpintervalwithinfo *stackptr;

  stack->nextfreeLcpintervalwithinfo = 0;
  stackptr = allocateStackspace(stack,adfst);
  stackptr->lcpitv.offset = 0;
  stackptr->lcpitv.leftbound = leftbound;
  stackptr->lcpitv.rightbound = rightbound;
  stackptr->lcpitv.code = 0;
  if (adfst->initLimdfsstate != NULL)
  {
    adfst->initLimdfsstate(stackptr->aliasstate,dfsconstinfo);
  }
}

void freeLimdfsresources(Limdfsresources **ptrlimdfsresources,
                         const AbstractDfstransformer *adfst)
{
  Limdfsresources *limdfsresources = *ptrlimdfsresources;

  adfst->freedfsconstinfo(&limdfsresources->dfsconstinfo);
  FREEARRAY(&limdfsresources->bwci,Boundswithchar);
  if (adfst->freeLimdfsstackelem != NULL)
  {
    unsigned long idx;

    for (idx = 0; idx < limdfsresources->stack.allocatedLcpintervalwithinfo;
         idx++)
    {
      adfst->freeLimdfsstackelem(limdfsresources->stack.
                                 spaceLcpintervalwithinfo[idx].aliasstate);
    }
    adfst->freeLimdfsstackelem(limdfsresources->currentdfsstate);
    adfst->freeLimdfsstackelem(limdfsresources->parentwithinfo.aliasstate);
  }
  FREEARRAY(&limdfsresources->stack,Lcpintervalwithinfo);
  FREESPACE(limdfsresources->rangeOccs);
  FREESPACE(limdfsresources->currentpathspace);
  FREEARRAY(&limdfsresources->mstatspos,Seqpos);
  FREESPACE(*ptrlimdfsresources);
}

/* enumerate the suffixes in an LCP-interval */

static void gen_esa_overinterval(const Genericindex *genericindex,
                                 Processmatch processmatch,
                                 void *processmatchinfo,
                                 const Indexbounds *itv,
                                 unsigned long pprefixlen,
                                 unsigned long distance,
                                 GT_UNUSED Seqpos totallength,
                                 const Uchar *dbsubstring)
{
  Seqpos idx;

  for (idx = itv->leftbound; idx <= itv->rightbound; idx++)
  {
    processmatch(processmatchinfo,
                 genericindex->suffixarray->suftab[idx],
                 itv->offset,
                 dbsubstring,
                 pprefixlen,
                 distance);
  }
}

static void esa_overinterval(Limdfsresources *limdfsresources,
                             const Indexbounds *itv,
                             unsigned long pprefixlen,
                             unsigned long distance)
{
  gen_esa_overinterval(limdfsresources->genericindex,
                       limdfsresources->processmatch,
                       limdfsresources->processmatchinfo,
                       itv,
                       pprefixlen,
                       distance,
                       limdfsresources->genericindex->totallength,
                       limdfsresources->currentpathspace);
  limdfsresources->numberofmatches += (itv->rightbound - itv->leftbound + 1);
}

static void gen_pck_overinterval(const Genericindex *genericindex,
                                 Processmatch processmatch,
                                 void *processmatchinfo,
                                 const Indexbounds *itv,
                                 unsigned long pprefixlen,
                                 unsigned long distance,
                                 Seqpos totallength,
                                 const Uchar *dbsubstring)
{
  Bwtseqpositioniterator *bspi;
  Seqpos dbstartpos;

  gt_assert(itv->leftbound < itv->rightbound);
  bspi = newBwtseqpositioniterator (genericindex->packedindex,
                                    itv->leftbound,itv->rightbound);
  while (nextBwtseqpositioniterator(&dbstartpos,bspi))
  {
    gt_assert(totallength >= (dbstartpos + itv->offset));
    processmatch(processmatchinfo,
                 totallength - (dbstartpos + itv->offset),
                 itv->offset,
                 dbsubstring,
                 pprefixlen,
                 distance);
  }
  freeBwtseqpositioniterator(&bspi);
}

static void pck_overinterval(Limdfsresources *limdfsresources,
                             const Indexbounds *itv,
                             unsigned long pprefixlen,
                             unsigned long distance)
{
  gen_pck_overinterval(limdfsresources->genericindex,
                       limdfsresources->processmatch,
                       limdfsresources->processmatchinfo,
                       itv,
                       pprefixlen,
                       distance,
                       limdfsresources->genericindex->totallength,
                       limdfsresources->currentpathspace);
  limdfsresources->numberofmatches += (itv->rightbound - itv->leftbound);
}

static void storemstatsposition(void *processinfo,
                                Seqpos dbstartpos,
                                GT_UNUSED Seqpos dblen,
                                GT_UNUSED const Uchar *dbsubstring,
                                GT_UNUSED unsigned long pprefixlen,
                                GT_UNUSED unsigned long distance)
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
                                          Seqpos leftbound,
                                          Seqpos rightbound,
                                          unsigned long offset)
{
  Indexbounds itv;

  limdfsresources->mstatspos.nextfreeSeqpos = 0;
  itv.leftbound = leftbound;
  itv.rightbound = rightbound;
  itv.offset = (Seqpos) offset;
  (limdfsresources->genericindex->withesa
       ? gen_esa_overinterval
       : gen_pck_overinterval)
    (limdfsresources->genericindex,
     storemstatsposition,
     &limdfsresources->mstatspos,
     &itv,
     offset,
     0,
     limdfsresources->genericindex->totallength,
     limdfsresources->currentpathspace);
  qsort(limdfsresources->mstatspos.spaceSeqpos,
        (size_t) limdfsresources->mstatspos.nextfreeSeqpos,
        sizeof (Seqpos), comparepositions);
  if (limdfsresources->genericindex->withesa)
  {
    limdfsresources->numberofmatches += (rightbound - leftbound + 1);
  } else
  {
    limdfsresources->numberofmatches += (rightbound - leftbound);
  }
  return &limdfsresources->mstatspos;
}

/* iterate myers algorithm over a sequence context */

static void esa_overcontext(Limdfsresources *limdfsresources,
                            Seqpos leftbound,
                            Seqpos offset,
                            const AbstractDfstransformer *adfst)
{
  Seqpos pos, startpos;
  Uchar cc;
  Limdfsresult limdfsresult;

  if (adfst->copyLimdfsstate == NULL)
  {
    memcpy(limdfsresources->currentdfsstate,
           limdfsresources->parentwithinfo.aliasstate,
           adfst->sizeofdfsstate);
  } else
  {
    adfst->copyLimdfsstate(limdfsresources->currentdfsstate,
                           limdfsresources->parentwithinfo.aliasstate,
                           limdfsresources->dfsconstinfo,
                           true);
  }
  startpos = limdfsresources->genericindex->suffixarray->suftab[leftbound];
#ifdef SKDEBUG
  printf("retrieve context of startpos=%lu\n",(unsigned long) startpos);
#endif
  for (pos = startpos + offset - 1;
       pos < limdfsresources->genericindex->totallength; pos++)
  {
    cc = getencodedchar(limdfsresources->genericindex->suffixarray->encseq,
                        pos,
                        limdfsresources->genericindex->suffixarray->readmode);
    if (cc != (Uchar) SEPARATOR &&
        (!limdfsresources->nowildcards || cc != (Uchar) WILDCARD))
    {
#ifdef SKDEBUG
      printf("cc=%u\n",(unsigned int) cc);
#endif
      adfst->inplacenextLimdfsstate(limdfsresources->dfsconstinfo,
                                    limdfsresources->currentdfsstate,
                                    (unsigned long) (pos - startpos + 1),
                                    cc);
      adfst->fullmatchLimdfsstate(&limdfsresult,
                                  limdfsresources->currentdfsstate,
                                  leftbound,
                                  leftbound,
                                  (Seqpos) 1,
                                  (unsigned long) (pos-startpos+1),
                                  limdfsresources->dfsconstinfo);
      if (limdfsresult.status == Limdfsstop)
      {
        break;
      }
      if (limdfsresult.status == Limdfssuccess)
      {
        limdfsresources->processmatch(limdfsresources->processmatchinfo,
                                      startpos,
                                      pos - startpos + 1,
                                      limdfsresources->currentpathspace,
                                      limdfsresult.pprefixlen,
                                      limdfsresult.distance);
        limdfsresources->numberofmatches++;
        break;
      }
    } else
    {
      break; /* failure */
    }
  }
}

static void addpathchar(Limdfsresources *limdfsresources,unsigned long idx,
                        Uchar cc)
{
  gt_assert(idx < limdfsresources->maxpathlength);
  limdfsresources->currentpathspace[idx] = cc;
}

static void pck_overcontext(Limdfsresources *limdfsresources,
                            Seqpos leftbound,
                            Seqpos offset,
                            Uchar inchar,
                            const AbstractDfstransformer *adfst)
{
  Uchar cc;
  unsigned long contextlength;
  bool processinchar = true;
  Seqpos bound = leftbound;
  Bwtseqcontextiterator *bsci;
  Limdfsresult limdfsresult;

  bsci = newBwtseqcontextiterator(limdfsresources->genericindex->packedindex,
                                  bound);
  if (adfst->copyLimdfsstate == NULL)
  {
    memcpy(limdfsresources->currentdfsstate,
           limdfsresources->parentwithinfo.aliasstate,adfst->sizeofdfsstate);
  } else
  {
    adfst->copyLimdfsstate(limdfsresources->currentdfsstate,
                           limdfsresources->parentwithinfo.aliasstate,
                           limdfsresources->dfsconstinfo,
                           true);
  }
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
      if (limdfsresources->maxpathlength > 0)
      {
        addpathchar(limdfsresources,
                    (unsigned long) (offset - 1 + contextlength),cc);
      }
      adfst->inplacenextLimdfsstate(limdfsresources->dfsconstinfo,
                                    limdfsresources->currentdfsstate,
                                    (unsigned long) (offset + contextlength),
                                    cc);
      adfst->fullmatchLimdfsstate(&limdfsresult,
                                  limdfsresources->currentdfsstate,
                                  bound,
                                  bound+1,
                                  (Seqpos) 1,
                                  (unsigned long) (offset+contextlength),
                                  limdfsresources->dfsconstinfo);
      if (limdfsresult.status == Limdfsstop)
      {
        break;
      }
      if (limdfsresult.status == Limdfssuccess)
      {
        Seqpos startpos = bwtseqfirstmatch(limdfsresources->genericindex->
                                                            packedindex,
                                           leftbound);
        limdfsresources->processmatch(limdfsresources->processmatchinfo,
                                      limdfsresources->genericindex->
                                                       totallength -
                                           (startpos + offset),
                                      offset + contextlength,
                                      limdfsresources->currentpathspace,
                                      limdfsresult.pprefixlen,
                                      limdfsresult.distance);
        limdfsresources->numberofmatches++;
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
                               const Indexbounds *child,
                               const AbstractDfstransformer *adfst)
{
  Limdfsresult limdfsresult;
  Seqpos width;

  stackptr->lcpitv = *child;
#ifdef SKDEBUG
  printf("(2) nextLimdfsstate(");
  adfst->showLimdfsstate(limdfsresources->parentwithinfo.aliasstate,
                         (unsigned long) (child->offset-1),
                         limdfsresources->dfsconstinfo);
  printf(",%u)=",(unsigned int) child->inchar);
#endif
  adfst->nextLimdfsstate(limdfsresources->dfsconstinfo,
                         stackptr->aliasstate,
                         (unsigned long) child->offset,
                         child->inchar,
                         limdfsresources->parentwithinfo.aliasstate);
#ifdef SKDEBUG
  adfst->showLimdfsstate(stackptr->aliasstate,
                         (unsigned long) child->offset,
                         limdfsresources->dfsconstinfo);
  printf("\n");
#endif
  if (limdfsresources->genericindex->withesa)
  {
    width = child->rightbound - child->leftbound + 1;
  } else
  {
    width = child->rightbound - child->leftbound;
  }
  adfst->fullmatchLimdfsstate(&limdfsresult,
                              stackptr->aliasstate,
                              child->leftbound,
                              child->rightbound,
                              width,
                              (unsigned long) child->offset,
                              limdfsresources->dfsconstinfo);
  if (limdfsresult.status == Limdfsstop)
  {
    return true; /* stop traversal without match */
  }
  if (limdfsresult.status == Limdfssuccess)
  {
    (limdfsresources->genericindex->withesa
         ? esa_overinterval
         : pck_overinterval)
      (limdfsresources,child,limdfsresult.pprefixlen,limdfsresult.distance);
    /* success with match of length pprefixlen - 1 */
    return true;
  }
  return false; /* continue with depth first traversal */
}

static void processchildinterval(Limdfsresources *limdfsresources,
                                 const Indexbounds *child,
                                 const AbstractDfstransformer *adfst)
{
  if (child->leftbound + 1 < child->rightbound ||
      (limdfsresources->genericindex->withesa &&
       child->leftbound + 1 == child->rightbound))
  {
    Lcpintervalwithinfo *stackptr;

    stackptr = allocateStackspace(&limdfsresources->stack,adfst);
    if (pushandpossiblypop(limdfsresources, stackptr, child, adfst))
    {
      limdfsresources->stack.nextfreeLcpintervalwithinfo--;
    }
  } else
  {
    if (limdfsresources->genericindex->withesa)
    {
      esa_overcontext(limdfsresources,
                      child->leftbound,
                      child->offset,
                      adfst);
    } else
    {
      pck_overcontext(limdfsresources,
                      child->leftbound,
                      child->offset,
                      child->inchar,
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
                                const AbstractDfstransformer *adfst)
{
  Seqpos firstnonspecial;
  Uchar extendchar;
  unsigned long idx;
  const Indexbounds *parent = &limdfsresources->parentwithinfo.lcpitv;

  extendchar = lcpintervalextendlcp(
                       limdfsresources->genericindex->suffixarray->encseq,
                       limdfsresources->genericindex->suffixarray->readmode,
                       limdfsresources->genericindex->suffixarray->suftab,
                       limdfsresources->genericindex->totallength,
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
    lcpintervalsplitwithoutspecial(
                     &limdfsresources->bwci,
                     limdfsresources->genericindex->suffixarray->encseq,
                     limdfsresources->genericindex->suffixarray->readmode,
                     limdfsresources->genericindex->totallength,
                     limdfsresources->genericindex->suffixarray->suftab,
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
    child.code = 0; /* not used, but we better define it */
    child.inchar = inchar;
#ifdef SKDEBUG
    printf("%u-child of ",(unsigned int) inchar);
    showLCPinterval(limdfsresources->genericindex->withesa,parent);
    printf(" is ");
    showLCPinterval(limdfsresources->genericindex->withesa,&child);
    printf("\n");
#endif
    processchildinterval(limdfsresources, &child, adfst);
    firstnonspecial = child.rightbound+1;
  }
  if (!limdfsresources->nowildcards)
  {
    Seqpos bound;

    for (bound=firstnonspecial; bound <= parent->rightbound; bound++)
    {
      esa_overcontext(limdfsresources,
                      bound,
                      parent->offset+1,
                      adfst);
    }
  }
}

static void smalldepthbwtrangesplitwithoutspecial(ArrayBoundswithchar *bwci,
                                                  const Mbtab **mbtab,
                                                  Uchar alphasize,
                                                  Codetype parentcode,
                                                  unsigned long childdepth)
{
  Codetype childcode;
  const Mbtab *mbptr;

  gt_assert(childdepth > 0);
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

unsigned long exactmatchuptomaxdepth(const Mbtab *mbptr,
                                     const Mbtab **mbtab,
                                     Uchar alphasize,
                                     unsigned int maxdepth,
                                     const Uchar *sequence,
                                     unsigned long seqlen)
{
  Codetype code = 0;
  unsigned int depth = 0;
  Uchar cc;

  gt_assert(seqlen > 0);
  for (depth = 0, code = 0;  depth <= maxdepth;
       depth++, code = code * alphasize + cc)
  {
    if (seqlen <= (unsigned long) depth)
    {
      return seqlen;
    }
    cc = sequence[depth];
    gt_assert(ISNOTSPECIAL(cc));
    mbptr = mbtab[depth] + code + cc;
    if (mbptr->lowerbound >= mbptr->upperbound)
    {
      break;
    }
  }
  return (unsigned long) depth;
}

static void pck_splitandprocess(Limdfsresources *limdfsresources,
                                const AbstractDfstransformer *adfst)
{
  unsigned long idx;
  Seqpos sumwidth = 0;
  const Indexbounds *parent = &limdfsresources->parentwithinfo.lcpitv;
  Codetype startcode;

  if (parent->offset < (Seqpos) limdfsresources->genericindex->maxdepth)
  {
    smalldepthbwtrangesplitwithoutspecial(&limdfsresources->bwci,
                                          limdfsresources->genericindex->mbtab,
                                          limdfsresources->alphasize,
                                          parent->code,
                                          (unsigned long) (parent->offset + 1));
    startcode = parent->code * limdfsresources->alphasize;
  } else
  {
    bwtrangesplitwithoutspecial(&limdfsresources->bwci,
                                limdfsresources->rangeOccs,
                                limdfsresources->genericindex->packedindex,
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
    gt_assert(inchar < limdfsresources->alphasize);
    child.code = startcode + inchar;
    child.inchar = inchar;
    if (limdfsresources->maxpathlength > 0)
    {
      addpathchar(limdfsresources,(unsigned long) parent->offset,inchar);
    }
    sumwidth += child.rightbound - child.leftbound;
#ifdef SKDEBUG
    printf("%u-child of ",(unsigned int) inchar);
    showLCPinterval(limdfsresources->genericindex->withesa,parent);
    printf(" is ");
    showLCPinterval(limdfsresources->genericindex->withesa,&child);
    printf("\n");
#endif
    processchildinterval(limdfsresources, &child, adfst);
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
        showLCPinterval(limdfsresources->genericindex->withesa,\
                        &STACKPTR->lcpitv);\
        adfst->showLimdfsstate(STACKPTR->aliasstate,\
                               (unsigned long) STACKPTR->lcpitv.offset,\
                               limdfsresources->dfsconstinfo);\
        printf("\n")
#else
#define SHOWSTACKTOP(STACKPTR) /* Nothing */
#endif

static void runlimdfs(Limdfsresources *limdfsresources,
                      const AbstractDfstransformer *adfst)
{
  Lcpintervalwithinfo *stackptr;

  gt_assert(adfst->sizeofdfsstate <=
            sizeof (limdfsresources->parentwithinfo.aliasstate));
  limdfsresources->numberofmatches = 0;
  initlcpinfostack(&limdfsresources->stack,
                   0,
                   limdfsresources->genericindex->withesa
                     ? limdfsresources->genericindex->totallength
                     : limdfsresources->genericindex->totallength+1,
                   limdfsresources->dfsconstinfo,
                   adfst);
  while (limdfsresources->stack.nextfreeLcpintervalwithinfo > 0)
  {
    gt_assert(limdfsresources->stack.spaceLcpintervalwithinfo != NULL);
    printf("nextfreeLcpintervalwithinfo=%lu\n",
               limdfsresources->stack.nextfreeLcpintervalwithinfo);
    stackptr = limdfsresources->stack.spaceLcpintervalwithinfo +
               limdfsresources->stack.nextfreeLcpintervalwithinfo - 1;
    SHOWSTACKTOP(stackptr);
    /* no make a copy of the top most stack element to be used as source */
    if (adfst->copyLimdfsstate == NULL)
    {
      limdfsresources->parentwithinfo = *stackptr; /* make a copy */
    } else
    {
      limdfsresources->parentwithinfo.lcpitv = stackptr->lcpitv;
      adfst->copyLimdfsstate(limdfsresources->parentwithinfo.aliasstate,
                             stackptr->aliasstate,
                             limdfsresources->dfsconstinfo,
                             true);
    }
    if (limdfsresources->parentwithinfo.lcpitv.offset > 0 &&
        limdfsresources->maxpathlength > 0)
    {
      addpathchar(limdfsresources,
                  (unsigned long)
                   (limdfsresources->parentwithinfo.lcpitv.offset-1),
                  limdfsresources->parentwithinfo.lcpitv.inchar);
    }
    gt_assert(limdfsresources->stack.nextfreeLcpintervalwithinfo > 0);
    limdfsresources->stack.nextfreeLcpintervalwithinfo--;
    (limdfsresources->genericindex->withesa
        ? esa_splitandprocess
        : pck_splitandprocess) (limdfsresources,adfst);
  }
  if (adfst->extractdfsconstinfo != NULL)
  {
    adfst->extractdfsconstinfo(limdfsresources->processresult,
                               limdfsresources,
                               limdfsresources->patterninfo,
                               limdfsresources->dfsconstinfo);
  }
}

bool indexbasedapproxpatternmatching(Limdfsresources *limdfsresources,
                                     const Uchar *pattern,
                                     unsigned long patternlength,
                                     unsigned long maxdistance,
                                     unsigned long maxintervalwidth,
                                     bool skpp,
                                     const AbstractDfstransformer *adfst)
{
  adfst->initdfsconstinfo(limdfsresources->dfsconstinfo,
                          (unsigned int) limdfsresources->alphasize,
                          pattern,
                          patternlength,
                          maxdistance,
                          maxintervalwidth,
                          skpp);
  runlimdfs(limdfsresources,adfst);
  return (limdfsresources->numberofmatches > 0) ? true : false;
}

void indexbasedmstats(Limdfsresources *limdfsresources,
                      const Uchar *pattern,
                      unsigned long patternlength,
                      const AbstractDfstransformer *adfst)
{
  adfst->initdfsconstinfo(limdfsresources->dfsconstinfo,
                          (unsigned int) limdfsresources->alphasize,
                          pattern,
                          patternlength);
  runlimdfs(limdfsresources,adfst);
}

void indexbasedspacedseeds(Limdfsresources *limdfsresources,
                           const Uchar *pattern,
                           Bitstring seedbitvector,
                           unsigned long seedweight,
                           const AbstractDfstransformer *adfst)
{
  adfst->initdfsconstinfo(limdfsresources->dfsconstinfo,
                          (unsigned int) limdfsresources->alphasize,
                          pattern,
                          seedbitvector,
                          seedweight);
  runlimdfs(limdfsresources,adfst);
}

void indexbasedlocali(Limdfsresources *limdfsresources,
                      long matchscore,
                      long mismatchscore,
                      long gapstart,
                      long gapextend,
                      unsigned long threshold,
                      const Uchar *query,
                      unsigned long querylength,
                      const AbstractDfstransformer *adfst)
{
  adfst->initdfsconstinfo(limdfsresources->dfsconstinfo,
                          (unsigned int) limdfsresources->alphasize,
                          matchscore,
                          mismatchscore,
                          gapstart,
                          gapextend,
                          query,
                          querylength,
                          threshold);
  runlimdfs(limdfsresources,adfst);
}

unsigned long genericmstats(const Limdfsresources *limdfsresources,
                            const Uchar *qstart,
                            const Uchar *qend)
{
  if (limdfsresources->genericindex->withesa)
  {
    return suffixarraymstats (limdfsresources->genericindex->suffixarray,
                              0,
                              0,
                              limdfsresources->genericindex->totallength,
                              NULL,
                              qstart,
                              qend);
  }
  return voidpackedindexmstatsforward(limdfsresources->genericindex->
                                                       packedindex,
                                      0,
                                      0,
                                      limdfsresources->genericindex->
                                                       totallength,
                                      NULL,
                                      qstart,
                                      qend);
}

static bool esa_exactpatternmatching(const Suffixarray *suffixarray,
                                     const Uchar *pattern,
                                     unsigned long patternlength,
                                     const Uchar *dbsubstring,
                                     Processmatch processmatch,
                                     void *processmatchinfo)
{
  MMsearchiterator *mmsi;
  Seqpos dbstartpos, totallength = getencseqtotallength(suffixarray->encseq);
  bool nomatches;

  mmsi = newmmsearchiterator(suffixarray->encseq,
                             suffixarray->suftab,
                             0,  /* leftbound */
                             totallength, /* rightbound */
                             0, /* offset */
                             suffixarray->readmode,
                             pattern,
                             patternlength);

  nomatches = isemptymmsearchiterator(mmsi);
  while (nextmmsearchiterator(&dbstartpos,mmsi))
  {
    processmatch(processmatchinfo,dbstartpos,
                 (Seqpos) patternlength,dbsubstring,patternlength,0);
  }
  freemmsearchiterator(&mmsi);
  return nomatches ? false : true;
}

bool indexbasedexactpatternmatching(const Limdfsresources *limdfsresources,
                                    const Uchar *pattern,
                                    unsigned long patternlength)
{
  if (limdfsresources->genericindex->withesa)
  {
    return esa_exactpatternmatching(limdfsresources->genericindex->suffixarray,
                                    pattern,
                                    patternlength,
                                    limdfsresources->currentpathspace,
                                    limdfsresources->processmatch,
                                    limdfsresources->processmatchinfo);
  } else
  {
    return pck_exactpatternmatching(limdfsresources->genericindex->packedindex,
                                    pattern,
                                    patternlength,
                                    limdfsresources->genericindex->totallength,
                                    limdfsresources->currentpathspace,
                                    limdfsresources->processmatch,
                                    limdfsresources->processmatchinfo);
  }
}

Seqpos bound2startpos(const Limdfsresources *limdfsresources,
                      Seqpos bound,unsigned long matchlength)
{
  if (limdfsresources->genericindex->withesa)
  {
    return limdfsresources->genericindex->suffixarray->suftab[bound];
  }
  return voidpackedfindfirstmatchconvert(limdfsresources->genericindex,
                                         bound,matchlength);
}

Uchar limdfsgetencodedchar(const Limdfsresources *limdfsresources,
                           Seqpos pos,
                           Readmode readmode)
{
  gt_assert(limdfsresources->encseq != NULL);

  return getencodedchar(limdfsresources->encseq,pos,readmode);
}

Seqpos getlastbound(const Limdfsresources *limdfsresources,Seqpos rightbound)
{
  if (limdfsresources->genericindex->withesa)
  {
    return rightbound;
  }
  return rightbound - 1;
}

bool intervalwidthleq(const Limdfsresources *limdfsresources,
                      Seqpos leftbound,Seqpos rightbound)
{
  Seqpos width;

  if (limdfsresources->genericindex->withesa)
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
