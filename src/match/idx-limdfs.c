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
#include "core/divmodmul.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "core/defined-types.h"
#include "core/ma_api.h"
#include "sarr-def.h"
#include "esa-splititv.h"
#include "eis-voiditf.h"
#include "esa-mmsearch.h"
#include "absdfstrans-imp.h"
#include "idx-limdfs.h"
#include "esa-map.h"
#include "idxlocalidp.h"
#include "esa-minunique.h"

#define DECLAREDFSSTATE(V)\
        Aliasdfsstate V[5]

struct Genericindex
{
  Suffixarray *suffixarray;
  unsigned long totallength;
  FMindex *packedindex;
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
  gt_freesuffixarray(genericindex->suffixarray);
  gt_free(genericindex->suffixarray);
  if (genericindex->packedindex != NULL)
  {
    gt_deletevoidBWTSeq(genericindex->packedindex);
  }
  gt_free(genericindex);
}

const GtEncseq *genericindex_getencseq(const Genericindex *genericindex)
{
  gt_assert(genericindex->suffixarray->encseq != NULL);
  return genericindex->suffixarray->encseq;
}

const Suffixarray *genericindex_getsuffixarray(const Genericindex
                                                *genericindex)
{
  gt_assert(genericindex->suffixarray != NULL);
  return genericindex->suffixarray;
}

Genericindex *genericindex_new(const char *indexname,
                               bool withesa,
                               bool withencseq,
                               bool withdestab,
                               bool withssptab,
                               int userdefinedmaxdepth,
                               GtLogger *logger,
                               GtError *err)
{
  unsigned int demand = 0;
  bool haserr = false;
  Genericindex *genericindex;

  genericindex = gt_malloc(sizeof (*genericindex));
  if (withesa)
  {
    demand |= SARR_SUFTAB;
  }
  if (withencseq)
  {
    demand |= SARR_ESQTAB;
  }
  if (withdestab)
  {
    demand |= SARR_DESTAB;
  }
  if (withssptab)
  {
    demand |= SARR_SSPTAB;
  }
  genericindex->withesa = withesa;
  genericindex->suffixarray = gt_malloc(sizeof (*genericindex->suffixarray));
  if (gt_mapsuffixarray(genericindex->suffixarray,
                        demand,
                        indexname,
                        logger,
                        err) != 0)
  {
    haserr = true;
    genericindex->totallength = 0;
  } else
  {
    genericindex->totallength
      = gt_encseq_total_length(genericindex->suffixarray->encseq);
  }
  if (!haserr)
  {
    if (withesa && genericindex->suffixarray->readmode != GT_READMODE_FORWARD)
    {
      gt_error_set(err,"using option -esa you can only process index "
                       "in forward mode");
      haserr = true;
    } else
    {
      if (!withesa
            && genericindex->suffixarray->readmode != GT_READMODE_REVERSE)
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
    genericindex->packedindex = gt_loadvoidBWTSeqForSA(indexname,true, err);
    if (genericindex->packedindex == NULL)
    {
      gt_assert(gt_error_is_set(err));
      haserr = true;
    }
  }
  if (!haserr && !withesa)
  {
    genericindex->mbtab = gt_bwtseq2mbtab(genericindex->packedindex);
    genericindex->maxdepth = gt_bwtseq2maxdepth(genericindex->packedindex);
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
  unsigned long offset,
                leftbound,
                rightbound;
  GtCodetype code;
  GtUchar inchar;
} Indexbounds;

typedef struct
{
  Indexbounds lcpitv;
  DECLAREDFSSTATE(aliasstate);
  /* the following two components are only required if keepexpandedonstack
     is true */
  bool keeponstack;
  unsigned long previousstackelem;
} Lcpintervalwithinfo;

GT_DECLAREARRAYSTRUCT(Lcpintervalwithinfo);

struct Limdfsresources
{
  Limdfsconstinfo *dfsconstinfo;
  GtArrayBoundswithchar bwci;
  GtArrayLcpintervalwithinfo stack;
  Lcpintervalwithinfo copyofparent;
  DECLAREDFSSTATE(copyofcopyofparentstate);
  unsigned long parentindex;
  bool keepexpandedonstack;
  GtUchar alphasize;
  void *patterninfo;
  const Genericindex *genericindex;
  bool nowildcards;
  unsigned long maxintervalwidth;
  unsigned long *rangeOccs;
  const GtEncseq *encseq;
  GtArrayGtUlong mstatspos;
  GtUchar *currentpathspace;
  unsigned long allocatedpathspace;
  unsigned long numberofmatches;
  ProcessIdxMatch processmatch;
  void *processmatchinfo;
  Processresult processresult;
};

Limdfsresources *gt_newLimdfsresources(const Genericindex *genericindex,
                                    bool nowildcards,
                                    unsigned long maxintervalwidth,
                                    unsigned long maxpathlength,
                                    bool keepexpandedonstack,
                                    ProcessIdxMatch processmatch,
                                    void *processmatchinfo,
                                    Processresult processresult,
                                    void *patterninfo,
                                    const AbstractDfstransformer *adfst)
{
  Limdfsresources *limdfsresources;
  unsigned int numofchars;
  const GtEncseq *encseq;

  encseq = genericindex->suffixarray->encseq;
  numofchars = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
  limdfsresources = gt_malloc(sizeof *limdfsresources);
  limdfsresources->bwci.spaceBoundswithchar
    = gt_malloc(sizeof *limdfsresources->bwci.spaceBoundswithchar
                * (numofchars+1));
  limdfsresources->bwci.nextfreeBoundswithchar = 0;
  limdfsresources->bwci.allocatedBoundswithchar
    = (unsigned long) (numofchars+1);
  GT_INITARRAY(&limdfsresources->stack,Lcpintervalwithinfo);
  gt_assert(numofchars <= GT_MAXALPHABETCHARACTER);
  limdfsresources->alphasize = (GtUchar) numofchars;
  limdfsresources->processmatch = processmatch;
  limdfsresources->processmatchinfo = processmatchinfo;
  limdfsresources->processresult = processresult;
  limdfsresources->patterninfo = patterninfo;
  limdfsresources->genericindex = genericindex;
  limdfsresources->nowildcards = nowildcards;
  limdfsresources->encseq = encseq;
  limdfsresources->maxintervalwidth = maxintervalwidth;
  limdfsresources->keepexpandedonstack = keepexpandedonstack;
  if (maxpathlength > 0)
  {
    limdfsresources->currentpathspace
      = gt_malloc(sizeof *limdfsresources->currentpathspace
                  * maxpathlength);
    limdfsresources->allocatedpathspace = maxpathlength;
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
    limdfsresources->rangeOccs
      = gt_malloc(sizeof *limdfsresources->rangeOccs
                  * GT_MULT2(limdfsresources->alphasize));
  }
  GT_INITARRAY(&limdfsresources->mstatspos,GtUlong);
  if (maxintervalwidth > 0)
  {
    limdfsresources->mstatspos.spaceGtUlong
      = gt_malloc(sizeof *limdfsresources->mstatspos.spaceGtUlong
                  * maxintervalwidth);
    limdfsresources->mstatspos.allocatedGtUlong = maxintervalwidth;
  }
  if (adfst->initLimdfsstackelem != NULL &&
      !limdfsresources->keepexpandedonstack)
  {
    adfst->initLimdfsstackelem(limdfsresources->copyofcopyofparentstate);
    adfst->initLimdfsstackelem(limdfsresources->copyofparent.aliasstate);
  }
  return limdfsresources;
}

static void tracethestackelems(GtIdxMatch *match,
                               Limdfsresources *limdfsresources,
                               unsigned long pprefixlen,
                               const Lcpintervalwithinfo *runptr)
{
  unsigned long previous = 0;

  gt_reinitLocaliTracebackstate(limdfsresources->dfsconstinfo,
                             runptr->lcpitv.offset,
                             pprefixlen);
  do
  {
    if (previous > 0)
    {
      gt_assert(previous - 1 == runptr->lcpitv.offset);
    }
    previous = runptr->lcpitv.offset;
    gt_assert(previous > 0);
    gt_assert(runptr->previousstackelem <
              limdfsresources->stack.nextfreeLcpintervalwithinfo);
    gt_processelemLocaliTracebackstate(limdfsresources->dfsconstinfo,
                                    runptr->lcpitv.inchar,runptr->aliasstate);
    runptr = limdfsresources->stack.spaceLcpintervalwithinfo +
             runptr->previousstackelem;
  } while (runptr->lcpitv.offset > 0);
  match->alignment
    = gt_completealignmentfromLocaliTracebackstate(&match->querylen,
                                                limdfsresources->dfsconstinfo);
  gt_assert(pprefixlen >= match->querylen);
  match->querystartpos = pprefixlen - match->querylen;
}

static Lcpintervalwithinfo *allocateStackspace(Limdfsresources *limdfsresources,
                                               const AbstractDfstransformer
                                                      *adfst)
{
  GtArrayLcpintervalwithinfo *stack = &limdfsresources->stack;

  if (stack->nextfreeLcpintervalwithinfo >= stack->allocatedLcpintervalwithinfo)
  {
    unsigned long idx;
    const unsigned long addelems = 32UL;

    stack->spaceLcpintervalwithinfo
      = gt_realloc(stack->spaceLcpintervalwithinfo,
                   sizeof *stack->spaceLcpintervalwithinfo
                   * (stack->allocatedLcpintervalwithinfo + addelems));
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

static void initlcpinfostack(Limdfsresources *limdfsresources,
                             const AbstractDfstransformer *adfst)
{
  Lcpintervalwithinfo *stackptr;

  limdfsresources->stack.nextfreeLcpintervalwithinfo = 0;
  stackptr = allocateStackspace(limdfsresources,adfst);
  stackptr->lcpitv.offset = 0;
  stackptr->lcpitv.leftbound = 0;
  stackptr->lcpitv.rightbound = limdfsresources->genericindex->withesa
                                ? limdfsresources->genericindex->totallength
                                : limdfsresources->genericindex->totallength+1;
  stackptr->lcpitv.code = 0;
  if (limdfsresources->keepexpandedonstack)
  {
    stackptr->previousstackelem = 0; /* Just to let it be defined */
    stackptr->keeponstack = true;
  }
  if (adfst->initrootLimdfsstate != NULL)
  {
    adfst->initrootLimdfsstate(stackptr->aliasstate,
                               limdfsresources->dfsconstinfo);
  }
}

void gt_freeLimdfsresources(Limdfsresources **ptrlimdfsresources,
                         const AbstractDfstransformer *adfst)
{
  Limdfsresources *limdfsresources = *ptrlimdfsresources;

  adfst->freedfsconstinfo(&limdfsresources->dfsconstinfo);
  GT_FREEARRAY(&limdfsresources->bwci,Boundswithchar);
  if (adfst->freeLimdfsstackelem != NULL)
  {
    unsigned long idx;

    for (idx = 0; idx < limdfsresources->stack.allocatedLcpintervalwithinfo;
         idx++)
    {
      adfst->freeLimdfsstackelem(limdfsresources->stack.
                                 spaceLcpintervalwithinfo[idx].aliasstate);
    }
    if (!limdfsresources->keepexpandedonstack)
    {
      adfst->freeLimdfsstackelem(limdfsresources->copyofcopyofparentstate);
      adfst->freeLimdfsstackelem(limdfsresources->copyofparent.aliasstate);
    }
  }
  GT_FREEARRAY(&limdfsresources->stack,Lcpintervalwithinfo);
  gt_free(limdfsresources->rangeOccs);
  gt_free(limdfsresources->currentpathspace);
  GT_FREEARRAY(&limdfsresources->mstatspos,GtUlong);
  gt_free(*ptrlimdfsresources);
}

/* enumerate the suffixes in an LCP-interval */

static void gen_esa_overinterval(const Genericindex *genericindex,
                                 ProcessIdxMatch processmatch,
                                 void *processmatchinfo,
                                 const Indexbounds *itv,
                                 GT_UNUSED unsigned long totallength,
                                 GtIdxMatch *match)
{
  unsigned long idx;

  for (idx = itv->leftbound; idx <= itv->rightbound; idx++)
  {
    match->dbstartpos = ESASUFFIXPTRGET(genericindex->suffixarray->suftab,idx);
    /* call processmatch */
    processmatch(processmatchinfo,match);
  }
}

static void esa_overinterval(Limdfsresources *limdfsresources,
                             const Indexbounds *itv,
                             GtIdxMatch *match)
{
  gen_esa_overinterval(limdfsresources->genericindex,
                       limdfsresources->processmatch,
                       limdfsresources->processmatchinfo,
                       itv,
                       limdfsresources->genericindex->totallength,
                       match);
  limdfsresources->numberofmatches += (itv->rightbound - itv->leftbound + 1);
}

static void gen_pck_overinterval(const Genericindex *genericindex,
                                 ProcessIdxMatch processmatch,
                                 void *processmatchinfo,
                                 const Indexbounds *itv,
                                 unsigned long totallength,
                                 GtIdxMatch *match)
{
  Bwtseqpositioniterator *bspi;
  unsigned long dbstartpos;

  gt_assert(itv->leftbound < itv->rightbound);
  bspi = gt_Bwtseqpositioniterator_new(genericindex->packedindex,
                                       itv->leftbound,itv->rightbound);
  while (gt_Bwtseqpositioniterator_next(&dbstartpos,bspi))
  {
    gt_assert(totallength >= (dbstartpos + itv->offset));
    /* call processmatch */
    match->dbstartpos = totallength - (dbstartpos + itv->offset);
    processmatch(processmatchinfo,match);
  }
  gt_Bwtseqpositioniterator_delete(bspi);
}

static void pck_overinterval(Limdfsresources *limdfsresources,
                             const Indexbounds *itv,
                             GtIdxMatch *match)
{
  gen_pck_overinterval(limdfsresources->genericindex,
                       limdfsresources->processmatch,
                       limdfsresources->processmatchinfo,
                       itv,
                       limdfsresources->genericindex->totallength,
                       match);
  limdfsresources->numberofmatches += (itv->rightbound - itv->leftbound);
}

static void storemstatsposition(void *processinfo,const GtIdxMatch *match)
{
  GtArrayGtUlong *mstatspos = (GtArrayGtUlong *) processinfo;

  GT_STOREINARRAY(mstatspos,GtUlong,32,match->dbstartpos);
}

static int comparepositions(const void *a, const void *b)
{
  if (*((unsigned long *) a) < *((unsigned long *) b))
  {
    return -1;
  }
  return 1;
}

GtArrayGtUlong *gt_fromitv2sortedmatchpositions(
                                             Limdfsresources *limdfsresources,
                                             unsigned long leftbound,
                                             unsigned long rightbound,
                                             unsigned long offset)
{
  Indexbounds itv;
  GtIdxMatch match;

  limdfsresources->mstatspos.nextfreeGtUlong = 0;
  itv.leftbound = leftbound;
  itv.rightbound = rightbound;
  itv.offset = (unsigned long) offset;
  match.dbabsolute = true;
  match.dblen = itv.offset;
  match.dbsubstring = limdfsresources->currentpathspace;
  match.querystartpos = 0;
  match.querylen = offset;
  match.distance = 0;
  match.alignment = NULL;
  (limdfsresources->genericindex->withesa
       ? gen_esa_overinterval
       : gen_pck_overinterval)
    (limdfsresources->genericindex,
     storemstatsposition,
     &limdfsresources->mstatspos,
     &itv,
     limdfsresources->genericindex->totallength,
     &match);
  qsort(limdfsresources->mstatspos.spaceGtUlong,
        (size_t) limdfsresources->mstatspos.nextfreeGtUlong,
        sizeof (unsigned long), comparepositions);
  if (limdfsresources->genericindex->withesa)
  {
    limdfsresources->numberofmatches += (rightbound - leftbound + 1);
  } else
  {
    limdfsresources->numberofmatches += (rightbound - leftbound);
  }
  return &limdfsresources->mstatspos;
}

static void initparentcopy(Limdfsresources *limdfsresources,
                           const AbstractDfstransformer *adfst)
{
  if (!limdfsresources->keepexpandedonstack)
  {
    if (adfst->copyLimdfsstate == NULL)
    {
      memcpy(limdfsresources->copyofcopyofparentstate,
             limdfsresources->copyofparent.aliasstate,
             adfst->sizeofdfsstate);
    } else
    {
      adfst->copyLimdfsstate(limdfsresources->copyofcopyofparentstate,
                             limdfsresources->copyofparent.aliasstate,
                             limdfsresources->dfsconstinfo);
    }
  }
}

static Lcpintervalwithinfo *expandsingleton(Limdfsresources *limdfsresources,
                                            unsigned long *resetvalue,
                                            bool notfirst,
                                            GtUchar cc,
                                            unsigned long currentdepth,
                                            const AbstractDfstransformer *adfst)
{
  if (limdfsresources->keepexpandedonstack)
  {
    Lcpintervalwithinfo *instate;
    Lcpintervalwithinfo *outstate;

    outstate = allocateStackspace(limdfsresources,adfst);
    outstate->keeponstack = true;
    outstate->lcpitv.offset = (unsigned long) currentdepth;
    outstate->lcpitv.inchar = cc;
    if (notfirst)
    {
      instate = outstate-1;
      gt_assert(limdfsresources->stack.nextfreeLcpintervalwithinfo >= 2UL);
      outstate->previousstackelem
        = limdfsresources->stack.nextfreeLcpintervalwithinfo - 2;
    } else
    {
      gt_assert(limdfsresources->parentindex <
                limdfsresources->stack.nextfreeLcpintervalwithinfo-1);
      *resetvalue = limdfsresources->stack.nextfreeLcpintervalwithinfo-1;
      instate = limdfsresources->stack.spaceLcpintervalwithinfo +
                limdfsresources->parentindex;
      outstate->previousstackelem = limdfsresources->parentindex;
    }
    gt_assert(instate < outstate);
    adfst->nextLimdfsstate(limdfsresources->dfsconstinfo,
                           outstate->aliasstate,
                           currentdepth,
                           cc,
                           instate->aliasstate);
    return outstate;
  }
  adfst->inplacenextLimdfsstate(limdfsresources->dfsconstinfo,
                                limdfsresources->copyofcopyofparentstate,
                                currentdepth,
                                cc);
  return NULL;
}

static void addpathchar(Limdfsresources *limdfsresources,unsigned long idx,
                        GtUchar cc)
{
  if (limdfsresources->currentpathspace != NULL)
  {
    gt_assert(idx < limdfsresources->allocatedpathspace);
    limdfsresources->currentpathspace[idx] = cc;
  }
}

/* iterate transformation algorithm over a sequence context */

static void esa_overcontext(Limdfsresources *limdfsresources,
                            const Indexbounds *child,
                            const AbstractDfstransformer *adfst)
{
  unsigned long pos, startpos;
  unsigned long resetvalue = 0;
  GtUchar cc;
  Limdfsresult limdfsresult;
  GtIdxMatch match;

  initparentcopy(limdfsresources,adfst);
  startpos = ESASUFFIXPTRGET(limdfsresources->genericindex->suffixarray->suftab,
                             child->leftbound);
#ifdef SKDEBUG
  printf("retrieve context of startpos=%lu\n",(unsigned long) startpos);
#endif
  for (pos = startpos + child->offset - 1;
       pos < limdfsresources->genericindex->totallength; pos++)
  {
    cc = gt_encseq_get_encoded_char(
                          limdfsresources->genericindex->suffixarray->encseq,
                          pos,
                          limdfsresources->genericindex->suffixarray->readmode);
    if (cc != (GtUchar) SEPARATOR &&
        (!limdfsresources->nowildcards || cc != (GtUchar) WILDCARD))
    {
#ifdef SKDEBUG
      printf("cc=%u\n",(unsigned int) cc);
#endif
      Lcpintervalwithinfo *outstate;

      addpathchar(limdfsresources,(unsigned long) (pos - startpos),cc);
      outstate = expandsingleton(limdfsresources,
                                 &resetvalue,
                                 (pos > startpos + child->offset - 1)
                                    ? true : false,
                                 cc,
                                 (unsigned long) (pos - startpos + 1),
                                 adfst);
      adfst->fullmatchLimdfsstate(&limdfsresult,
                                  limdfsresources->keepexpandedonstack
                                    ?  outstate->aliasstate
                                    :  limdfsresources->copyofcopyofparentstate,
                                  child->leftbound,
                                  child->leftbound,
                                  (unsigned long) 1,
                                  (unsigned long) (pos-startpos+1),
                                  limdfsresources->dfsconstinfo);
      if (limdfsresult.status == Limdfsstop)
      {
        break;
      }
      if (limdfsresult.status == Limdfssuccess)
      {
        match.dbabsolute = true;
        match.dbstartpos = startpos;
        match.dblen = pos - startpos + 1;
        match.dbsubstring = limdfsresources->currentpathspace;
        match.querylen = limdfsresult.pprefixlen;
        match.distance = limdfsresult.distance;
        if (limdfsresources->keepexpandedonstack)
        {
          tracethestackelems(&match,limdfsresources,limdfsresult.pprefixlen,
                             outstate);
        } else
        {
          match.querystartpos = 0;
          match.alignment = NULL;
        }
        /* call processmatch */
        limdfsresources->processmatch(limdfsresources->processmatchinfo,&match);
        limdfsresources->numberofmatches++;
        break;
      }
    } else
    {
      break; /* failure */
    }
  }
  if (limdfsresources->keepexpandedonstack)
  {
    gt_assert(resetvalue > 0);
    limdfsresources->stack.nextfreeLcpintervalwithinfo = resetvalue;
  }
}

/*
   the following function iterates to determine the smallest length
   at which the expanded pattern matches
*/

static void pck_overcontext(Limdfsresources *limdfsresources,
                            const Indexbounds *child,
                            const AbstractDfstransformer *adfst)
{
  GtUchar cc;
  unsigned long contextlength, resetvalue = 0;
  bool processinchar = true;
  unsigned long bound;
  Bwtseqcontextiterator *bsci;
  Limdfsresult limdfsresult;
  GtIdxMatch match;

  gt_assert(child != NULL);
  bound = child->leftbound;
  bsci
    = gt_Bwtseqcontextiterator_new(limdfsresources->genericindex->packedindex,
                                   bound);
  initparentcopy(limdfsresources,adfst);
#ifdef SKDEBUG
  printf("retrieve context for bound = %lu\n",(unsigned long) bound);
#endif
  for (contextlength = 0; /* nothing */; contextlength++)
  {
    if (processinchar)
    {
      cc = child->inchar;
      processinchar = false;
    } else
    {
      cc = gt_Bwtseqcontextiterator_next(&bound,bsci);
    }
    if (cc != (GtUchar) SEPARATOR &&
        (!limdfsresources->nowildcards || cc != (GtUchar) WILDCARD))
    {
      Lcpintervalwithinfo *outstate;
#ifdef SKDEBUG
      printf("cc=%u\n",(unsigned int) cc);
#endif
      addpathchar(limdfsresources,
                  (unsigned long) (child->offset - 1 + contextlength),cc);
      outstate = expandsingleton(limdfsresources,
                                 &resetvalue,
                                 (contextlength > 0) ? true : false,
                                 cc,
                                 (unsigned long) (child->offset+contextlength),
                                 adfst);
      adfst->fullmatchLimdfsstate(&limdfsresult,
                                  limdfsresources->keepexpandedonstack
                                    ? outstate->aliasstate
                                    : limdfsresources->copyofcopyofparentstate,
                                  bound,
                                  bound+1,
                                  (unsigned long) 1,
                                  (unsigned long) (child->offset+contextlength),
                                  limdfsresources->dfsconstinfo);
      if (limdfsresult.status == Limdfsstop)
      {
        break;
      }
      if (limdfsresult.status == Limdfssuccess)
      {
        unsigned long startpos;

        startpos = gt_bwtseqfirstmatch(
                                    limdfsresources->genericindex->packedindex,
                                    child->leftbound);
        match.dbabsolute = true;
        match.dbstartpos = limdfsresources->genericindex->totallength -
                           (startpos + child->offset);
        match.dblen = child->offset + contextlength;
        match.dbsubstring = limdfsresources->currentpathspace;
        match.querylen = limdfsresult.pprefixlen;
        match.distance = limdfsresult.distance;
        if (limdfsresources->keepexpandedonstack)
        {
          tracethestackelems(&match,limdfsresources,limdfsresult.pprefixlen,
                             outstate);
        } else
        {
          match.querystartpos = 0;
          match.alignment = NULL;
        }
        /* call processmatch */
        limdfsresources->processmatch(limdfsresources->processmatchinfo,&match);
        limdfsresources->numberofmatches++;
        break;
      }
    } else
    {
      break;
    }
  }
  if (limdfsresources->keepexpandedonstack)
  {
    gt_assert(resetvalue > 0);
    limdfsresources->stack.nextfreeLcpintervalwithinfo = resetvalue;
  }
  gt_Bwtseqcontextiterator_delete(bsci);
  bsci = NULL;
}

static const Lcpintervalwithinfo *currentparent(const Limdfsresources
                                                *limdfsresources)
{
  const Lcpintervalwithinfo *ptr;

  if (limdfsresources->keepexpandedonstack)
  {
    gt_assert(limdfsresources->parentindex <
              limdfsresources->stack.nextfreeLcpintervalwithinfo);
  }
  ptr = limdfsresources->keepexpandedonstack
           ? (limdfsresources->stack.spaceLcpintervalwithinfo +
              limdfsresources->parentindex)
           : &limdfsresources->copyofparent;
  return ptr;
}

/* the following function does not contain any iteration */

static void pushandpossiblypop(Limdfsresources *limdfsresources,
                               const Indexbounds *child,
                               const AbstractDfstransformer *adfst)
{
  Limdfsresult limdfsresult;
  unsigned long width;
  Lcpintervalwithinfo *stackptr;

#ifdef SKDEBUG
  printf("(2) nextLimdfsstate(");
  adfst->showLimdfsstate(currentparent(limdfsresources)->aliasstate,
                         (unsigned long) (child->offset-1),
                         limdfsresources->dfsconstinfo);
  printf(",%u)=",(unsigned int) child->inchar);
#endif
  stackptr = allocateStackspace(limdfsresources,adfst);
  adfst->nextLimdfsstate(limdfsresources->dfsconstinfo,
                         stackptr->aliasstate,
                         (unsigned long) child->offset,
                         child->inchar,
                         currentparent(limdfsresources)->aliasstate);
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
  if (limdfsresult.status == Limdfscontinue)
  {
    stackptr->lcpitv = *child;
    if (limdfsresources->keepexpandedonstack)
    {
      stackptr->keeponstack = true;
      stackptr->previousstackelem = limdfsresources->parentindex;
    }
    return; /* no success, but still have the chance to find result */
  }
  if (limdfsresult.status == Limdfssuccess)
  {
    GtIdxMatch match;

    match.dbabsolute = true;
    match.querylen = limdfsresult.pprefixlen;
    match.distance = limdfsresult.distance;
    match.dblen = child->offset;
    match.dbsubstring = limdfsresources->currentpathspace;
    if (limdfsresources->keepexpandedonstack)
    {
      gt_assert(stackptr >= limdfsresources->stack.spaceLcpintervalwithinfo &&
                stackptr < limdfsresources->stack.spaceLcpintervalwithinfo +
                           limdfsresources->stack.nextfreeLcpintervalwithinfo);
      stackptr->lcpitv = *child;
      stackptr->keeponstack = true;
      stackptr->previousstackelem = limdfsresources->parentindex;
      tracethestackelems(&match,limdfsresources,limdfsresult.pprefixlen,
                         stackptr);
    } else
    {
      match.querystartpos = 0;
      match.alignment = NULL;
    }
    /* success with match of length pprefixlen */
    (limdfsresources->genericindex->withesa
         ? esa_overinterval
         : pck_overinterval) (limdfsresources,child,&match);
  }
  /* now status == Limdfssuccess || status == Limdfsstop */
  /* pop the element from the stack as there has been success or stop event */
  limdfsresources->stack.nextfreeLcpintervalwithinfo--;
}

static void processchildinterval(Limdfsresources *limdfsresources,
                                 const Indexbounds *child,
                                 const AbstractDfstransformer *adfst)
{
  if (child->leftbound + 1 < child->rightbound ||
      (limdfsresources->genericindex->withesa &&
       child->leftbound + 1 == child->rightbound))
  {
    pushandpossiblypop(limdfsresources, child, adfst);
  } else
  {
    if (limdfsresources->genericindex->withesa)
    {
      esa_overcontext(limdfsresources,child,adfst);
    } else
    {
      pck_overcontext(limdfsresources,child,adfst);
    }
  }
}

#ifdef SKDEBUG

static void showLCPinterval(bool withesa,const Indexbounds *itv)
{
  unsigned long width;

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
  unsigned long firstspecial;
  GtUchar extendchar;
  unsigned long idx;
  const Indexbounds *parent = &(currentparent(limdfsresources)->lcpitv);

  extendchar = gt_lcpintervalextendlcp(
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
    gt_lcpintervalsplitwithoutspecial(
                     &limdfsresources->bwci,
                     limdfsresources->genericindex->suffixarray->encseq,
                     limdfsresources->genericindex->suffixarray->readmode,
                     limdfsresources->genericindex->totallength,
                     limdfsresources->genericindex->suffixarray->suftab,
                     parent->offset,
                     parent->leftbound,
                     parent->rightbound);
  }
  firstspecial = parent->leftbound;
  for (idx = 0; idx < limdfsresources->bwci.nextfreeBoundswithchar; idx++)
  {
    Indexbounds child;

    /* reset as parentptr may have been moved */
    parent = &(currentparent(limdfsresources)->lcpitv);
    child.inchar = limdfsresources->bwci.spaceBoundswithchar[idx].inchar;
    child.offset = parent->offset+1;
    child.leftbound = limdfsresources->bwci.spaceBoundswithchar[idx].lbound;
    child.rightbound = limdfsresources->bwci.spaceBoundswithchar[idx].rbound;
    child.code = 0; /* not used, but we better define it */
    addpathchar(limdfsresources,(unsigned long) parent->offset,child.inchar);
#ifdef SKDEBUG
    printf("%u-child of ",(unsigned int) child.inchar);
    showLCPinterval(limdfsresources->genericindex->withesa,parent);
    printf(" is ");
    showLCPinterval(limdfsresources->genericindex->withesa,&child);
    printf("\n");
#endif
    processchildinterval(limdfsresources, &child, adfst);
    firstspecial = child.rightbound+1;
  }
  if (!limdfsresources->nowildcards)
  {
    Indexbounds child;
    unsigned long bound;

    child.inchar = (GtUchar) WILDCARD;
    child.offset = parent->offset+1;
    child.code = 0;  /* not used, but we better define it */
    for (bound = firstspecial; bound <= parent->rightbound; bound++)
    {
      child.leftbound = child.rightbound = bound;
      esa_overcontext(limdfsresources,&child,adfst);
    }
  }
}

static void smalldepthbwtrangesplitwithoutspecial(GtArrayBoundswithchar *bwci,
                                                  const Mbtab **mbtab,
                                                  GtUchar alphasize,
                                                  GtCodetype parentcode,
                                                  unsigned long childdepth)
{
  GtCodetype childcode;
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
        = (GtUchar) (mbptr - (mbtab[childdepth] + childcode));
      bwci->spaceBoundswithchar[bwci->nextfreeBoundswithchar].lbound
        = mbptr->lowerbound;
      bwci->spaceBoundswithchar[bwci->nextfreeBoundswithchar++].rbound
        = mbptr->upperbound;
    }
  }
}

static void pck_splitandprocess(Limdfsresources *limdfsresources,
                                const AbstractDfstransformer *adfst)
{
  unsigned long idx;
  unsigned long sumwidth = 0;
  const Indexbounds *parent = &(currentparent(limdfsresources)->lcpitv);
  GtCodetype startcode;

  if (parent->offset < (unsigned long) limdfsresources->genericindex->maxdepth)
  {
    smalldepthbwtrangesplitwithoutspecial(&limdfsresources->bwci,
                                          limdfsresources->genericindex->mbtab,
                                          limdfsresources->alphasize,
                                          parent->code,
                                          (unsigned long) (parent->offset + 1));
    startcode = parent->code * limdfsresources->alphasize;
  } else
  {
    gt_bwtrangesplitwithoutspecial(&limdfsresources->bwci,
                                limdfsresources->rangeOccs,
                                limdfsresources->genericindex->packedindex,
                                parent->leftbound,
                                parent->rightbound);
    startcode = 0;
  }
  for (idx = 0; idx < limdfsresources->bwci.nextfreeBoundswithchar; idx++)
  {
    Indexbounds child;

    /* reset as parentptr may have been moved */
    parent = &(currentparent(limdfsresources)->lcpitv);
    child.inchar = limdfsresources->bwci.spaceBoundswithchar[idx].inchar;
    child.offset = parent->offset+1;
    child.leftbound = limdfsresources->bwci.spaceBoundswithchar[idx].lbound;
    child.rightbound = limdfsresources->bwci.spaceBoundswithchar[idx].rbound;
    gt_assert(child.inchar < limdfsresources->alphasize);
    child.code = startcode + child.inchar;
    addpathchar(limdfsresources,(unsigned long) parent->offset,child.inchar);
    sumwidth += child.rightbound - child.leftbound;
#ifdef SKDEBUG
    printf("%u-child of ",(unsigned int) child.inchar);
    showLCPinterval(limdfsresources->genericindex->withesa,parent);
    printf(" is ");
    showLCPinterval(limdfsresources->genericindex->withesa,&child);
    printf("\n");
#endif
    processchildinterval(limdfsresources, &child, adfst);
  }
  if (!limdfsresources->nowildcards)
  {
    unsigned long bound;
    Indexbounds child;

    for (bound = parent->leftbound + sumwidth;
         bound < parent->rightbound; bound++)
    {
      GtUchar cc = gt_bwtseqgetsymbol(bound,
                                   limdfsresources->genericindex->packedindex);

      child.offset = parent->offset+1;
      child.code = 0;  /* not used, but we better define it */
      if (cc != (GtUchar) SEPARATOR)
      {
        child.leftbound = bound;
        child.inchar = cc;
        pck_overcontext(limdfsresources,&child,adfst);
      }
    }
  }
}

#ifdef SKDEBUG
#define SHOWSTACKTOP(STACKPTR)\
        printf("top=");\
        showLCPinterval(limdfsresources->genericindex->withesa,\
                        &(STACKPTR)->lcpitv);\
        adfst->showLimdfsstate((STACKPTR)->aliasstate,\
                               (unsigned long) (STACKPTR)->lcpitv.offset,\
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
            sizeof (limdfsresources->copyofparent.aliasstate));
  limdfsresources->numberofmatches = 0;
  initlcpinfostack(limdfsresources,adfst);
  while (limdfsresources->stack.nextfreeLcpintervalwithinfo > 0)
  {
    gt_assert(limdfsresources->stack.spaceLcpintervalwithinfo != NULL);
    /*
    printf("nextfreeLcpintervalwithinfo=%lu\n",
               limdfsresources->stack.nextfreeLcpintervalwithinfo);
    */
    stackptr = limdfsresources->stack.spaceLcpintervalwithinfo +
               limdfsresources->stack.nextfreeLcpintervalwithinfo - 1;
    SHOWSTACKTOP(stackptr);
    if (limdfsresources->keepexpandedonstack)
    {
      if (stackptr->keeponstack)
      {
        limdfsresources->parentindex
          = limdfsresources->stack.nextfreeLcpintervalwithinfo - 1;
      } else
      {
        limdfsresources->stack.nextfreeLcpintervalwithinfo--;
        continue;
      }
    } else
    {
      /* make a copy of the top most stack element to be used as source */
      if (adfst->copyLimdfsstate == NULL)
      {
        limdfsresources->copyofparent = *stackptr; /* make a copy */
      } else
      {
        limdfsresources->copyofparent.lcpitv = stackptr->lcpitv;
        adfst->copyLimdfsstate(limdfsresources->copyofparent.aliasstate,
                               stackptr->aliasstate,
                               limdfsresources->dfsconstinfo);
      }
      /* now parentptr always points to copyofparent */
    }
    if (currentparent(limdfsresources)->lcpitv.offset > 0)
    {
      addpathchar(limdfsresources,
                  (unsigned long)
                   (currentparent(limdfsresources)->lcpitv.offset-1),
                  currentparent(limdfsresources)->lcpitv.inchar);
    }
    gt_assert(limdfsresources->stack.nextfreeLcpintervalwithinfo > 0);
    /* now delete the top element from the stack as we have made a copy */
    if (limdfsresources->keepexpandedonstack)
    {
      stackptr->keeponstack = false;
    } else
    {
      limdfsresources->stack.nextfreeLcpintervalwithinfo--;
    }
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

typedef enum
{
  Popitv,
  Splititv,
  Processitv,
  Pushitv,
  Processcontext
} Runlimdfsstate;

bool gt_indexbasedapproxpatternmatching(Limdfsresources *limdfsresources,
                                     const GtUchar *pattern,
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

void gt_indexbasedmstats(Limdfsresources *limdfsresources,
                      const GtUchar *pattern,
                      unsigned long patternlength,
                      const AbstractDfstransformer *adfst)
{
  adfst->initdfsconstinfo(limdfsresources->dfsconstinfo,
                          (unsigned int) limdfsresources->alphasize,
                          pattern,
                          patternlength);
  runlimdfs(limdfsresources,adfst);
}

void gt_indexbasedspacedseeds(Limdfsresources *limdfsresources,
                           const GtUchar *pattern,
                           GtBitsequence seedbitvector,
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

void gt_indexbasedlocali(Limdfsresources *limdfsresources,
                      long matchscore,
                      long mismatchscore,
                      long gapstart,
                      long gapextend,
                      unsigned long threshold,
                      const GtUchar *query,
                      unsigned long querylength,
                      const AbstractDfstransformer *adfst)
{
  adfst->initdfsconstinfo(limdfsresources->dfsconstinfo,
                          (unsigned int) limdfsresources->alphasize,
                          matchscore,
                          mismatchscore,
                          gapstart,
                          gapextend,
                          threshold,
                          query,
                          querylength);
  runlimdfs(limdfsresources,adfst);
}

unsigned long genericmstats(const Limdfsresources *limdfsresources,
                            const GtUchar *qstart,
                            const GtUchar *qend)
{
  if (limdfsresources->genericindex->withesa)
  {
    return gt_suffixarraymstats (limdfsresources->genericindex->suffixarray,
                              0,
                              0,
                              limdfsresources->genericindex->totallength,
                              NULL,
                              qstart,
                              qend);
  }
  return gt_voidpackedindexmstatsforward(limdfsresources->genericindex->
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
                                     const GtUchar *pattern,
                                     unsigned long patternlength,
                                     GT_UNUSED const GtUchar *dbsubstring,
                                     ProcessIdxMatch processmatch,
                                     void *processmatchinfo)
{
  GtMMsearchiterator *mmsi;
  unsigned long dbstartpos,
         totallength = gt_encseq_total_length(suffixarray->encseq);
  bool nomatches;
  GtIdxMatch match;

  mmsi = gt_mmsearchiterator_new_complete_olain(suffixarray->encseq,
                                           suffixarray->suftab,
                                           0,  /* leftbound */
                                           totallength, /* rightbound */
                                           0, /* offset */
                                           suffixarray->readmode,
                                           pattern,
                                           patternlength);
  nomatches = gt_mmsearchiterator_isempty(mmsi);
  match.dbabsolute = true;
  match.dblen = (unsigned long) patternlength;
  match.dbsubstring = pattern;
  match.querystartpos = 0;
  match.querylen = patternlength;
  match.distance = 0;
  match.alignment = NULL;
  while (gt_mmsearchiterator_next(&dbstartpos,mmsi))
  {
    /* call processmatch */
    match.dbstartpos = dbstartpos;
    processmatch(processmatchinfo,&match);
  }
  gt_mmsearchiterator_delete(mmsi);
  return nomatches ? false : true;
}

bool gt_indexbasedexactpatternmatching(const Limdfsresources *limdfsresources,
                                    const GtUchar *pattern,
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
    return gt_pck_exactpatternmatching(
                                    limdfsresources->genericindex->packedindex,
                                    pattern,
                                    patternlength,
                                    limdfsresources->genericindex->totallength,
                                    limdfsresources->currentpathspace,
                                    limdfsresources->processmatch,
                                    limdfsresources->processmatchinfo);
  }
}

GtUchar gt_limdfs_getencodedchar(const Limdfsresources *limdfsresources,
                              unsigned long pos,
                              GtReadmode readmode)
{
  gt_assert(limdfsresources->encseq != NULL);

  return gt_encseq_get_encoded_char(limdfsresources->encseq,
                                           pos,
                                           readmode);
}

bool gt_intervalwidthleq(const Limdfsresources *limdfsresources,
                      unsigned long leftbound,unsigned long rightbound)
{
  unsigned long width;

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
  if (width > 0 && width <= (unsigned long) limdfsresources->maxintervalwidth)
  {
    return true;
  }
  return false;
}

const FMindex *genericindex_get_packedindex(const Genericindex *genericindex)
{
  gt_assert(genericindex->packedindex != NULL);
  return genericindex->packedindex;
}

unsigned long genericindex_get_totallength(const Genericindex *genericindex)
{
  gt_assert(genericindex && genericindex->totallength != 0);
  return genericindex->totallength;
}

unsigned long esa_exact_pattern_count(const Suffixarray *suffixarray,
                                      const GtUchar *pattern,
                                      unsigned long patternlength) {
  GtMMsearchiterator *mmsi;
  unsigned long count,
                totallength = gt_encseq_total_length(suffixarray->encseq);

  mmsi = gt_mmsearchiterator_new_complete_olain(suffixarray->encseq,
                                           suffixarray->suftab,
                                           0,  /* leftbound */
                                           totallength, /* rightbound */
                                           0, /* offset */
                                           suffixarray->readmode,
                                           pattern,
                                           patternlength);

  count = gt_mmsearchiterator_count(mmsi);
  gt_mmsearchiterator_delete(mmsi);
  return count;
}

unsigned long gt_indexbased_exact_pattern_count(
                                              const Genericindex *genericindex,
                                              const GtUchar *pattern,
                                              unsigned long patternlength) {
  if (genericindex->withesa) {
    return esa_exact_pattern_count(genericindex->suffixarray,
                                    pattern,
                                    patternlength);
  }
  else {
    return gt_pck_exact_pattern_count(genericindex->packedindex,
                                       pattern,
                                       patternlength);
  }
}
