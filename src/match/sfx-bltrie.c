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

#include "core/unused_api.h"
#include "core/arraydef.h"
#include "core/divmodmul.h"
#include "core/minmax.h"
#include "core/encseq.h"
#include "core/arraydef.h"

#include "lcpoverflow.h"
#include "sfx-bltrie.h"
#include "sfx-suffixgetset.h"

#define BLINDTRIECHAR_ISSPECIAL(X) ((X) >= (Blindtriesymbol) WILDCARD)

typedef unsigned long Blindtriesymbol;

typedef struct Blindtrienode
{
  struct Blindtrienode *rightsibling;
  union
  {
    struct
    {
      struct Blindtrienode *firstchild;
      unsigned long depth;
    } internalinfo;
    struct
    {
      unsigned long nodestartpos,
                    nodestoppos;
    } leafinfo;
  } either;
  Blindtriesymbol firstchar;
  bool isleaf;
} Blindtrienode;

typedef Blindtrienode * Blindtrienodeptr;

GT_DECLAREARRAYSTRUCT(Blindtrienodeptr);

struct Blindtrie
{
  /* The following belongs to the state and is initialized by blindtrie_new */
  unsigned long allocatedBlindtrienode,
                nextfreeBlindtrienode;
  Blindtrienode *spaceBlindtrienode;
  Blindtrienodeptr root;
  GtArrayBlindtrienodeptr stack;
  GtArrayGtUlong overflowsuffixes;

  /* The following needs to be supplied for each insertion */
  const GtEncseq *encseq;
  GtEncseqReader *esr1, *esr2;
  GtReadmode readmode;
  unsigned long totallength,
                maxdepthminusoffset;
  bool cmpcharbychar,
       has_twobitencoding_stoppos_support;
  GtSuffixsortspace *sssp;
};

static bool blindtrie_isleaf(const Blindtrienodeptr node)
{
  return node->isleaf;
}

static void blindtrie_setleaf(const Blindtrienodeptr node,bool isleaf)
{
  node->isleaf = isleaf;
}

static unsigned long blindtrie_getdepth(const Blindtrienodeptr node)
{
  gt_assert(!blindtrie_isleaf(node));
  return node->either.internalinfo.depth;
}

static void blindtrie_setdepth(const Blindtrienodeptr node,unsigned long depth)
{
  gt_assert(!blindtrie_isleaf(node));
  node->either.internalinfo.depth = depth;
}

static Blindtriesymbol blindtrie_firstchar_get(const Blindtrienodeptr node)
{
  return node->firstchar;
}

static void blindtrie_firstchar_set(Blindtrienodeptr node,bool isleaf,
                                    Blindtriesymbol firstchar)
{
  blindtrie_setleaf(node,isleaf);
  gt_assert(isleaf || !GT_ISUNIQUEINT(firstchar));
  node->firstchar = firstchar;
}

static bool blindtrie_isleftofboundary(const Blindtrie *blindtrie,
                                       unsigned long currentstartpos,
                                       unsigned long add)
{
  unsigned long endpos;

  if (blindtrie->maxdepthminusoffset == 0)
  {
    endpos = blindtrie->totallength;
  } else
  {
    endpos = currentstartpos + blindtrie->maxdepthminusoffset;
    if (endpos >= blindtrie->totallength)
    {
      endpos = blindtrie->totallength;
    }
  }
  return (currentstartpos + add < endpos) ? true : false;
}

static Blindtrienodeptr blindtrie_newnode(Blindtrie *blindtrie)
{
  gt_assert(blindtrie->nextfreeBlindtrienode <
            blindtrie->allocatedBlindtrienode);
  return blindtrie->spaceBlindtrienode + blindtrie->nextfreeBlindtrienode++;
}

static Blindtrienode *blindtrie_newleaf(Blindtrie *blindtrie,
                                        unsigned long currentstartpos,
                                        unsigned long
                                          currenttwobitencodingstoppos,
                                        Blindtriesymbol firstchar,
                                        struct Blindtrienode *rightsibling)
{
  Blindtrienode *newleaf;

  newleaf = blindtrie_newnode(blindtrie);
  blindtrie_firstchar_set(newleaf,true,firstchar);
  newleaf->either.leafinfo.nodestartpos = currentstartpos;
  newleaf->either.leafinfo.nodestoppos = currenttwobitencodingstoppos;
  newleaf->rightsibling = rightsibling;
  return newleaf;
}

static unsigned long blindtrie_currenttwobitencodingstoppos_get(
                                     const Blindtrie *blindtrie,
                                     unsigned long currentstartpos)
{
  if (!blindtrie->cmpcharbychar &&
      blindtrie->has_twobitencoding_stoppos_support)
  {
    gt_encseq_reader_reinit_with_readmode(blindtrie->esr1,blindtrie->encseq,
                                          blindtrie->readmode,
                                          currentstartpos);
    return gt_getnexttwobitencodingstoppos(GT_ISDIRREVERSE(blindtrie->readmode)
                                           ? false : true,blindtrie->esr1);
  }
  return GT_TWOBITENCODINGSTOPPOSUNDEF(blindtrie);
}

static Blindtrienodeptr blindtrie_makeroot(Blindtrie *blindtrie,
                                           unsigned long currentstartpos)
{
  Blindtrienode *root;
  Blindtriesymbol firstchar;
  unsigned long currenttwobitencodingstoppos;

  root = blindtrie_newnode(blindtrie);
  blindtrie_firstchar_set(root,false,0); /* firstchar of root will
                                            never be used */
  blindtrie_setdepth(root,0);
  root->rightsibling = NULL;
  if (blindtrie_isleftofboundary(blindtrie,currentstartpos,0))
  {
    /* Random access */
    firstchar = (Blindtriesymbol)
                gt_encseq_get_encoded_char(blindtrie->encseq,
                                           currentstartpos,
                                           blindtrie->readmode);
    if (BLINDTRIECHAR_ISSPECIAL(firstchar))
    {
      firstchar = GT_UNIQUEINT(currentstartpos);
    }
    currenttwobitencodingstoppos
      = blindtrie_currenttwobitencodingstoppos_get(blindtrie,currentstartpos);
  } else
  {
    firstchar = GT_UNIQUEINT(currentstartpos);
    currenttwobitencodingstoppos = GT_TWOBITENCODINGSTOPPOSUNDEF(blindtrie);
  }
  root->either.internalinfo.firstchild
    = blindtrie_newleaf(blindtrie,currentstartpos,
                        currenttwobitencodingstoppos,
                        firstchar,NULL);
  return root;
}

static inline Blindtrienodeptr blindtrie_extractleafnode(Blindtrienodeptr head)
{
  gt_assert(!blindtrie_isleaf(head));
  do
  {
    head = head->either.internalinfo.firstchild;
  } while (!blindtrie_isleaf(head));
  return head;
}

static inline int blindtrie_comparecharacters(Blindtriesymbol oldchar,
                                              Blindtriesymbol newchar)
{
  return (oldchar > newchar) ? 1 : ((oldchar < newchar) ? -1 : 0);
}

static Blindtrienodeptr blindtrie_findsucc(Blindtrienodeptr node,
                                           Blindtriesymbol newchar)
{
  int retval;

  for (;;)
  {
    retval = blindtrie_comparecharacters(blindtrie_firstchar_get(node),newchar);
    if (retval == 0)
    {              /* found branch corresponding to newchar */
      return node;
    }
    if (retval == 1)
    {               /* found sibling which is already greater than newchar */
      return NULL;
    }
    node = node->rightsibling;
    if (node == NULL) /* no more siblings: mismatch */
    {
      return NULL;
    }
  }
}

static Blindtrienodeptr blindtrie_findcompanion(
                                    Blindtrie *blindtrie,
                                    unsigned long currentstartpos,
                                    unsigned long currenttwobitencodingstoppos)
{
  Blindtriesymbol newchar;
  Blindtrienodeptr head, succ;
  unsigned long headdepth;

  blindtrie->stack.nextfreeBlindtrienodeptr = 0;
  head = blindtrie->root;
  while (!blindtrie_isleaf(head))
  {
    GT_STOREINARRAY (&blindtrie->stack, Blindtrienodeptr, 128, head);
    headdepth = blindtrie_getdepth(head);
    if (blindtrie_isleftofboundary(blindtrie,currentstartpos,headdepth))
    {
      /* Random access */
      if (blindtrie->has_twobitencoding_stoppos_support)
      {
        if ((GT_ISDIRREVERSE(blindtrie->readmode) &&
            GT_REVERSEPOS(blindtrie->totallength,currentstartpos+headdepth)
            >= currenttwobitencodingstoppos) ||
            (!GT_ISDIRREVERSE(blindtrie->readmode) &&
             currentstartpos + headdepth < currenttwobitencodingstoppos))
        {
          newchar = (Blindtriesymbol)
                    gt_encseq_extract_encoded_char(blindtrie->encseq,
                                                   currentstartpos + headdepth,
                                                   blindtrie->readmode);
        } else
        {
          newchar = GT_UNIQUEINT(currentstartpos + headdepth);
        }
      } else
      {
        newchar = (Blindtriesymbol)
                  gt_encseq_get_encoded_char(blindtrie->encseq,
                                             currentstartpos + headdepth,
                                             blindtrie->readmode);
        if (BLINDTRIECHAR_ISSPECIAL(newchar))
        {
          newchar = GT_UNIQUEINT(currentstartpos + headdepth);
        }
      }
    } else
    {
      newchar = GT_UNIQUEINT(currentstartpos + headdepth);
    }
    if (GT_ISUNIQUEINT(newchar))
    {
      return blindtrie_extractleafnode(head);
    }
    succ = blindtrie_findsucc(head->either.internalinfo.firstchild,newchar);
    if (succ == NULL)
    {
      return blindtrie_extractleafnode(head);
    }
    head = succ;
  }
  GT_STOREINARRAY (&blindtrie->stack, Blindtrienodeptr, 128, head);
  return head;
}

static void blindtrie_insertatsplitnode(Blindtrie *blindtrie,
                                        Blindtrienodeptr oldnode,
                                        Blindtriesymbol mm_oldsuffix,
                                        unsigned long lcp,
                                        Blindtriesymbol mm_newsuffix,
                                        unsigned long currentstartpos,
                                        unsigned long
                                           currenttwobitencodingstoppos)
{
  Blindtrienodeptr newleaf, newnode, previous, current;

  gt_assert(GT_ISUNIQUEINT(mm_oldsuffix) ||
            GT_ISUNIQUEINT(mm_newsuffix) ||
            mm_oldsuffix != mm_newsuffix ||
            blindtrie_isleaf(oldnode) ||
            blindtrie_getdepth(oldnode) == lcp);

  /* insert a new node before node oldnode if necessary */
  if (blindtrie_isleaf(oldnode))
  {
    gt_assert(lcp > 0);
  }
  if (blindtrie_isleaf(oldnode) || blindtrie_getdepth(oldnode) > lcp)
  {
    newnode = blindtrie_newnode(blindtrie);
    blindtrie_firstchar_set(newnode,blindtrie_isleaf(oldnode),mm_oldsuffix);
    if (!blindtrie_isleaf(oldnode))
    {
      blindtrie_setdepth(newnode,blindtrie_getdepth(oldnode));
      /* newnode inherits depth+children */
    }
    newnode->either = oldnode->either;
    newnode->rightsibling = NULL;
    blindtrie_setleaf(oldnode,false);
    gt_assert(!GT_ISUNIQUEINT(blindtrie_firstchar_get(oldnode)));
    gt_assert(lcp > 0);
    blindtrie_setdepth(oldnode,lcp);
    /* oldnode has newnode as only child*/
    oldnode->either.internalinfo.firstchild = newnode;
  }
  gt_assert(blindtrie_isleaf(oldnode) || blindtrie_getdepth(oldnode) == lcp);
  previous = NULL;
  current = oldnode->either.internalinfo.firstchild;
  while (current != NULL &&
         blindtrie_comparecharacters(blindtrie_firstchar_get(current),
                                     mm_newsuffix) < 0)
  {
    previous = current;
    current = current->rightsibling;
  }
  /* insert new leaf with current suffix */
  /* search S[lcp] among the offsprings */
  newleaf = blindtrie_newleaf(blindtrie,currentstartpos,
                              currenttwobitencodingstoppos,mm_newsuffix,
                              current);
  if (previous != NULL)
  {
    previous->rightsibling = newleaf;
  } else
  {
    oldnode->either.internalinfo.firstchild = newleaf;
  }
}

static unsigned long blindtrie_cmpcharbychar_getlcp(
                                 Blindtriesymbol *mm_oldsuffix,
                                 GT_UNUSED bool *mm_oldsuffixisseparator,
                                 Blindtriesymbol *mm_newsuffix,
                                 const Blindtrie *blindtrie,
                                 unsigned long leafpos,
                                 unsigned long currentstartpos)
{
  unsigned long lcp;
  Blindtriesymbol cc1, cc2;

  if (blindtrie_isleftofboundary(blindtrie,leafpos,0))
  {
    gt_encseq_reader_reinit_with_readmode(blindtrie->esr1,blindtrie->encseq,
                                          blindtrie->readmode,leafpos);
  }
  if (blindtrie_isleftofboundary(blindtrie,currentstartpos,0))
  {
    gt_encseq_reader_reinit_with_readmode(blindtrie->esr2,blindtrie->encseq,
                                          blindtrie->readmode,currentstartpos);
  }
  for (lcp = 0; /* Nothing */; lcp++)
  {
    if (blindtrie_isleftofboundary(blindtrie,leafpos,lcp))
    {
      cc1 = (Blindtriesymbol)
            gt_encseq_reader_next_encoded_char(blindtrie->esr1);
      if (BLINDTRIECHAR_ISSPECIAL(cc1))
      {
        cc1 = GT_UNIQUEINT(leafpos + lcp);
      }
    } else
    {
      cc1 = GT_UNIQUEINT(leafpos + lcp);
    }
    if (blindtrie_isleftofboundary(blindtrie,currentstartpos,lcp))
    {
      cc2 = (Blindtriesymbol)
            gt_encseq_reader_next_encoded_char(blindtrie->esr2);
      if (BLINDTRIECHAR_ISSPECIAL(cc2))
      {
        cc2 = GT_UNIQUEINT(currentstartpos + lcp);
      }
    } else
    {
      cc2 = GT_UNIQUEINT(currentstartpos + lcp);
    }
    if (blindtrie_comparecharacters(cc1,cc2) != 0)
    {
      *mm_oldsuffix = cc1;
      *mm_newsuffix = cc2;
      break;
    }
  }
  gt_assert(blindtrie->maxdepthminusoffset == 0 ||
            lcp <= blindtrie->maxdepthminusoffset);
  return lcp;
}

static unsigned long blindtrie_twobitencoding_getlcp(
                                 Blindtriesymbol *mm_oldsuffix,
                                 GT_UNUSED bool *mm_oldsuffixisseparator,
                                 Blindtriesymbol *mm_newsuffix,
                                 const Blindtrie *blindtrie,
                                 unsigned long leafpos,
                                 unsigned long leaftwobitencodingstoppos,
                                 unsigned long currentstartpos,
                                 unsigned long currenttwobitencodingstoppos)
{
  GtCommonunits commonunits;
  GtViatwobitkeyvalues vtk1, vtk2;
  const unsigned long depth = 0;

  gt_assert(leafpos != currentstartpos);
  gt_assignvittwobitkeyvalues(&vtk1,blindtrie->encseq,blindtrie->readmode,
                              NULL,leafpos,depth,
                              blindtrie->maxdepthminusoffset);
  gt_assignvittwobitkeyvalues(&vtk2,blindtrie->encseq,blindtrie->readmode,
                              NULL,currentstartpos,depth,
                              blindtrie->maxdepthminusoffset);
  vtk1.twobitencodingstoppos = leaftwobitencodingstoppos;
  vtk2.twobitencodingstoppos = currenttwobitencodingstoppos;
  (void) gt_encseq_process_viatwobitencoding(&commonunits,
                                             blindtrie->encseq,
                                             blindtrie->readmode,
                                             depth,
                                             blindtrie->maxdepthminusoffset,
                                             &vtk1,
                                             &vtk2);
  if (blindtrie_isleftofboundary(blindtrie,leafpos,commonunits.finaldepth) &&
      !commonunits.leftspecial)
  {
    *mm_oldsuffix = (Blindtriesymbol)
                    gt_encseq_extract_encoded_char(blindtrie->encseq,
                                                   leafpos +
                                                     commonunits.finaldepth,
                                                   blindtrie->readmode);
  } else
  {
    *mm_oldsuffix = GT_UNIQUEINT(leafpos + commonunits.finaldepth);
  }
  if (blindtrie_isleftofboundary(blindtrie,currentstartpos,
                                 commonunits.finaldepth) &&
      !commonunits.rightspecial)
  {
    *mm_newsuffix = (Blindtriesymbol)
                    gt_encseq_extract_encoded_char(blindtrie->encseq,
                                                   currentstartpos +
                                                   commonunits.finaldepth,
                                                   blindtrie->readmode);
  } else
  {
    *mm_newsuffix = GT_UNIQUEINT(currentstartpos + commonunits.finaldepth);
  }
  return commonunits.finaldepth;
}

static unsigned long blindtrie_getlcp(Blindtriesymbol *mm_oldsuffix,
                                      bool *mm_oldsuffixisseparator,
                                      Blindtriesymbol *mm_newsuffix,
                                      const Blindtrie *blindtrie,
                                      const Blindtrienodeptr lis,
                                      unsigned long currentstartpos,
                                      unsigned long
                                        currenttwobitencodingstoppos)
{
  if (blindtrie->cmpcharbychar)
  {
    return blindtrie_cmpcharbychar_getlcp (mm_oldsuffix,
                                           mm_oldsuffixisseparator,
                                           mm_newsuffix,
                                           blindtrie,
                                           lis->either.leafinfo.nodestartpos,
                                           currentstartpos);
  }
  return blindtrie_twobitencoding_getlcp(mm_oldsuffix,
                                         mm_oldsuffixisseparator,
                                         mm_newsuffix,
                                         blindtrie,
                                         lis->either.leafinfo.nodestartpos,
                                         lis->either.leafinfo.nodestoppos,
                                         currentstartpos,
                                         currenttwobitencodingstoppos);
}

static void blindtrie_suffixout(Blindtrie *blindtrie,
                                unsigned long subbucketleft,
                                unsigned long offset,
                                unsigned long nextfree,
                                unsigned long startpos)
{
  gt_suffixsortspace_set(blindtrie->sssp,subbucketleft,nextfree,
                         startpos - offset);
}

#define BLINDTRIE_SETCURRENTNODE(NODEPTR)\
        currentnodeisleaf = blindtrie_isleaf(NODEPTR) ? true : false;\
        currentnode = NODEPTR

static unsigned long blindtrie_enumeratetrieleaves (
                           Blindtrie *blindtrie,
                           unsigned long subbucketleft,
                           unsigned long offset,
                           unsigned long maxdepth,
                           unsigned long *lcpsubtab,
                           unsigned long *numoflargelcpvalues,
                           void *voiddcov,
                           Dc_processunsortedrange dc_processunsortedrange)
{
  bool readyforpop = false, currentnodeisleaf;
  Blindtrienodeptr currentnode, siblval, lcpnode = blindtrie->root;
  unsigned long nextfree = 0, equalsrangewidth = 0, lcpnodedepth,
                bucketleftidxplussubbucketleft;

  blindtrie->stack.nextfreeBlindtrienodeptr = 0;
  GT_STOREINARRAY (&blindtrie->stack, Blindtrienodeptr, 128, blindtrie->root);
  BLINDTRIE_SETCURRENTNODE(blindtrie->root->either.internalinfo.firstchild);
  gt_assert(maxdepth == 0 || dc_processunsortedrange != NULL);
  bucketleftidxplussubbucketleft
    = gt_suffixsortspace_bucketleftidx_get(blindtrie->sssp) + subbucketleft;
  for (;;)
  {
    lcpnodedepth = blindtrie_getdepth(lcpnode);
    if (currentnodeisleaf)
    {
      if (nextfree > 0)
      {
        if (lcpsubtab != NULL)
        {
          lcpsubtab[nextfree] = lcpnodedepth + offset;
          if (lcpnodedepth + offset >= (unsigned long) LCPOVERFLOW)
          {
            (*numoflargelcpvalues)++;
          }
        }
        if (maxdepth > 0)
        {
          if (lcpnodedepth + offset == maxdepth)
          {
            equalsrangewidth++;
          } else
          {
#ifndef NDEBUG
            if (lcpnodedepth + offset >= maxdepth)
            {
              fprintf(stderr,"lcpnode.depth=%lu,offset=%lu,maxdepth=%lu\n",
                              lcpnodedepth,
                              offset,
                              maxdepth);
              exit(EXIT_FAILURE);
            }
#endif
            gt_assert(lcpnodedepth + offset < maxdepth);
            if (equalsrangewidth > 0)
            {
              dc_processunsortedrange(
                               voiddcov,
                               bucketleftidxplussubbucketleft
                                 + nextfree - 1 - equalsrangewidth,
                               equalsrangewidth + 1,
                               maxdepth);
              equalsrangewidth = 0;
            }
          }
        }
      }
      blindtrie_suffixout(blindtrie,subbucketleft,offset,nextfree,
                          currentnode->either.leafinfo.nodestartpos);
      nextfree++;
      siblval = currentnode->rightsibling;
      if (siblval == NULL)
      {
        readyforpop = true;
        currentnodeisleaf = false; /* STATE 1 */
      } else
      {
        BLINDTRIE_SETCURRENTNODE(siblval);  /* current comes from brother */
        lcpnode = blindtrie->stack.spaceBlindtrienodeptr[
                             blindtrie->stack.nextfreeBlindtrienodeptr-1];
      }
    } else
    {
      if (readyforpop)
      {
        if (blindtrie->stack.nextfreeBlindtrienodeptr == 1UL)
        {
          break;
        }
        blindtrie->stack.nextfreeBlindtrienodeptr--;
        siblval = blindtrie->stack.spaceBlindtrienodeptr[
                       blindtrie->stack.nextfreeBlindtrienodeptr]->rightsibling;
        if (siblval != NULL)
        {
          BLINDTRIE_SETCURRENTNODE(siblval);   /* current comes from brother */
          lcpnode = blindtrie->stack.spaceBlindtrienodeptr[
                             blindtrie->stack.nextfreeBlindtrienodeptr - 1];
          readyforpop = false;
        }
      } else
      {
        GT_STOREINARRAY (&blindtrie->stack, Blindtrienodeptr, 128, currentnode);
        BLINDTRIE_SETCURRENTNODE(currentnode->either.internalinfo.firstchild);
      }
    }
  }
  if (nextfree > 0 && equalsrangewidth > 0)
  {
    dc_processunsortedrange(voiddcov,
                            bucketleftidxplussubbucketleft
                              + nextfree - 1 - equalsrangewidth,
                            equalsrangewidth + 1,
                            maxdepth);
    equalsrangewidth = 0;
  }
  return nextfree;
}

Blindtrie *gt_blindtrie_new(GtSuffixsortspace *suffixsortspace,
                            unsigned long numofsuffixes,
                            const GtEncseq *encseq,
                            bool cmpcharbychar,
                            GtEncseqReader *esr1,
                            GtEncseqReader *esr2,
                            GtReadmode readmode)
{
  Blindtrie *blindtrie;

  blindtrie = gt_malloc(sizeof (*blindtrie));
  blindtrie->allocatedBlindtrienode = GT_MULT2(numofsuffixes + 1) + 1;
  blindtrie->spaceBlindtrienode
    = gt_malloc(sizeof (*blindtrie->spaceBlindtrienode) *
                blindtrie->allocatedBlindtrienode);
  /*
  printf("# sizeof (blindtrie)=%lu\n",
            (unsigned long) (sizeof (Blindtrie) +
                             blindtrie->allocatedBlindtrienode *
                             sizeof (Blindtrienode)));
  */
  blindtrie->nextfreeBlindtrienode = 0;
  GT_INITARRAY (&blindtrie->overflowsuffixes, GtUlong);
  GT_INITARRAY (&blindtrie->stack, Blindtrienodeptr);

  blindtrie->encseq = encseq;
  blindtrie->has_twobitencoding_stoppos_support
    = gt_has_twobitencoding_stoppos_support(encseq);
  blindtrie->readmode = readmode;
  blindtrie->root = NULL;
  blindtrie->esr1 = esr1;
  blindtrie->esr2 = esr2;
  blindtrie->totallength = gt_encseq_total_length(encseq);
  blindtrie->cmpcharbychar = cmpcharbychar;
  blindtrie->sssp = suffixsortspace;
  return blindtrie;
}

void gt_blindtrie_delete(Blindtrie *blindtrie)
{
  if (blindtrie == NULL)
  {
    return;
  }
  gt_free(blindtrie->spaceBlindtrienode);
  GT_FREEARRAY(&blindtrie->overflowsuffixes, GtUlong);
  GT_FREEARRAY(&blindtrie->stack, Blindtrienodeptr);
  gt_free(blindtrie);
}

#undef SKDEBUG
#ifdef SKDEBUG

#define NODENUM(PTR)\
        ((PTR) == NULL\
           ? 99UL\
           : (unsigned long) ((PTR) - blindtrie->spaceBlindtrienode))

static void gt_blindtrie_showleaf(const Blindtrie *blindtrie,unsigned int level,
                                  Blindtrienodeptr current)
{
  printf("%*.*s",(int) (6 * level),(int) (6 * level)," ");
  gt_assert(current != NULL);
  printf("Leaf(add=%lu,firstchar=%u,startpos=%lu,rightsibling=%lu)\n",
         NODENUM(current),
         (unsigned int) current->firstchar,
         current->either.leafinfo.nodestartpos,
         NODENUM(current->rightsibling));
}

static void gt_blindtrie_showintern(const Blindtrie *blindtrie,
                                    unsigned int level,
                                    Blindtrienodeptr current)
{
  printf("%*.*s",(int) (6 * level),(int) (6 * level)," ");
  gt_assert(current != NULL);
  printf("Intern(add=%lu,firstchar=%u,depth=%lu"
         ",firstchild=%lu,rightsibling=%lu)\n",
          NODENUM(current),
          (unsigned int) current->firstchar,
          blindtrie_getdepth(current),
          NODENUM(current->either.internalinfo.firstchild),
          NODENUM(current->rightsibling));
}

static void gt_blindtrie_showrecursive(const Blindtrie *blindtrie,
                                       unsigned int level,
                                       Blindtrienodeptr node)
{
  Blindtrienodeptr current;

  for (current = node->either.internalinfo.firstchild;
       current != NULL;
       current = current->rightsibling)
  {
    if (blindtrie_isleaf(current))
    {
      gt_blindtrie_showleaf(blindtrie,level,current);
    } else
    {
      gt_blindtrie_showintern(blindtrie,level,current);
      gt_blindtrie_showrecursive(blindtrie,level+1,current);
    }
  }
}

static void gt_blindtrie_show(const Blindtrie *blindtrie)
{
  gt_blindtrie_showrecursive(blindtrie,0,blindtrie->root);
}

static void blindtrie_showstate(const Blindtrie *blindtrie,
                                unsigned long subbucketleft,
                                unsigned long numberofsuffixes,
                                unsigned long offset)
{
  unsigned long idx;

  printf("insert suffixes at offset %lu:\n",offset);
  for (idx=0; idx < numberofsuffixes; idx++)
  {
    printf("%lu ",
           gt_suffixsortspace_get(blindtrie->sssp,subbucketleft,idx) + offset);
  }
  printf("\nstep 0\n");
  gt_blindtrie_show(blindtrie);
}

#endif

static int blindtrie_compare_ascending(const void *a,const void *b)
{
  unsigned long *aptr = (unsigned long *) a;
  unsigned long *bptr = (unsigned long *) b;
  gt_assert(*aptr != *bptr);
  return *aptr < *bptr ? -1 : 1;
}

static void processoverflowsuffixes(Blindtrie *blindtrie,
                                    unsigned long offset,
                                    unsigned long *lcpsubtab,
                                    unsigned long subbucketleft,
                                    unsigned long nextsuffixtooutput)
{
  if (blindtrie->overflowsuffixes.nextfreeGtUlong > 0)
  {
    unsigned long idx;

    if (blindtrie->overflowsuffixes.nextfreeGtUlong > 1UL)
    {
      qsort(blindtrie->overflowsuffixes.spaceGtUlong,
            sizeof (*blindtrie->overflowsuffixes.spaceGtUlong),
            (size_t) blindtrie->overflowsuffixes.nextfreeGtUlong,
            blindtrie_compare_ascending);
    }
    for (idx = 0; idx < blindtrie->overflowsuffixes.nextfreeGtUlong; idx++)
    {
      if (lcpsubtab != NULL)
      {
        lcpsubtab[nextsuffixtooutput] = offset;
      }
      blindtrie_suffixout(blindtrie,subbucketleft,offset,nextsuffixtooutput,
                          blindtrie->overflowsuffixes.spaceGtUlong[idx]);
      nextsuffixtooutput++;
    }
  }
}

static unsigned long gt_blindtrie2sorting(Blindtrie *blindtrie,
                                          unsigned long subbucketleft,
                                          unsigned long *lcpsubtab,
                                          unsigned long offset,
                                          unsigned long maxdepth,
                                          void *voiddcov,
                                          Dc_processunsortedrange
                                            dc_processunsortedrange)
{
  unsigned long nextsuffixtosort, numoflargelcpvalues = 0;

  nextsuffixtosort
    = blindtrie_enumeratetrieleaves (blindtrie,subbucketleft,offset,maxdepth,
                                     lcpsubtab,&numoflargelcpvalues,
                                     voiddcov,dc_processunsortedrange);
  processoverflowsuffixes(blindtrie,
                          offset,
                          lcpsubtab,
                          subbucketleft,
                          nextsuffixtosort);
  if (lcpsubtab != NULL && blindtrie->overflowsuffixes.nextfreeGtUlong > 0 &&
      offset >= (unsigned long) LCPOVERFLOW)
  {
    numoflargelcpvalues += blindtrie->overflowsuffixes.nextfreeGtUlong;
  }
  return numoflargelcpvalues;
}

void gt_blindtrie_insertsuffix(Blindtrie *blindtrie,
                               unsigned long offset,
                               unsigned long maxdepth,
                               unsigned long currentstartpos)
{

  gt_assert(maxdepth == 0 || maxdepth > offset);
  if (blindtrie->nextfreeBlindtrienode == 0)
  { /* empty tree */
    if (maxdepth == 0)
    {
      blindtrie->maxdepthminusoffset = 0;
    } else
    {
      blindtrie->maxdepthminusoffset = maxdepth - offset;
    }
    blindtrie->overflowsuffixes.nextfreeGtUlong = 0;
    blindtrie->root = blindtrie_makeroot(blindtrie,currentstartpos);
  } else
  {
    unsigned long stackidx,
                  lcp,
                  currenttwobitencodingstoppos;
    Blindtrienodeptr leafinsubtrie, currentnode;
    Blindtriesymbol mm_oldsuffix, mm_newsuffix;

    if (blindtrie_isleftofboundary(blindtrie,currentstartpos,0))
    {
      currenttwobitencodingstoppos
        = blindtrie_currenttwobitencodingstoppos_get(blindtrie,currentstartpos);
      leafinsubtrie = blindtrie_findcompanion(blindtrie,currentstartpos,
                                              currenttwobitencodingstoppos);
      gt_assert(blindtrie_isleaf(leafinsubtrie));
      lcp = blindtrie_getlcp (&mm_oldsuffix,
                              NULL,
                              &mm_newsuffix,
                              blindtrie,
                              leafinsubtrie,
                              currentstartpos,
                              currenttwobitencodingstoppos);
      currentnode = blindtrie->root;
      for (stackidx=0;stackidx<blindtrie->stack.nextfreeBlindtrienodeptr;
           stackidx++)
      {
        currentnode = blindtrie->stack.spaceBlindtrienodeptr[stackidx];
        if (blindtrie_isleaf(currentnode) ||
            blindtrie_getdepth(currentnode) >= lcp)
        {
          break;
        }
      }
      blindtrie_insertatsplitnode(blindtrie,
                                  currentnode,
                                  mm_oldsuffix,
                                  lcp,
                                  mm_newsuffix,
                                  currentstartpos,
                                  currenttwobitencodingstoppos);
    } else
    {
      GT_STOREINARRAY(&blindtrie->overflowsuffixes,GtUlong,32,currentstartpos);
    }
  }
}

unsigned long gt_blindtrie_suffixsort(
                            Blindtrie *blindtrie,
                            unsigned long subbucketleft,
                            unsigned long *lcpsubtab,
                            unsigned long numberofsuffixes,
                            unsigned long offset,
                            unsigned long maxdepth,
                            void *voiddcov,
                            Dc_processunsortedrange dc_processunsortedrange)
{
  unsigned long idx, currentstartpos;

  blindtrie->nextfreeBlindtrienode = 0;
  for (idx=0; idx < numberofsuffixes; idx++)
  {
    currentstartpos = gt_suffixsortspace_get(blindtrie->sssp,subbucketleft,idx);
    gt_blindtrie_insertsuffix(blindtrie,
                              offset,
                              maxdepth,
                              currentstartpos + offset);
  }
  return gt_blindtrie2sorting(blindtrie,
                              subbucketleft,
                              lcpsubtab,
                              offset,
                              maxdepth,
                              voiddcov,
                              dc_processunsortedrange);
}
