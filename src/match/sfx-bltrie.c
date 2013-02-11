/*
  Copyright (c) 2007/2010 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007/2010 Center for Bioinformatics, University of Hamburg

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

#define GT_BLINDTRIECHAR_ISSPECIAL(X)     ((X) >= (GtBlindtriesymbol) WILDCARD)
#define GT_BLINDTRIE_REFNULL              0
#define GT_BLINDTRIE_BITSFORRIGHTSIBLING  31
#define GT_BLINDTRIE_ROOTIDX              0

typedef unsigned long GtBlindtriesymbol;
typedef unsigned int GtBlindtriesnodeptr;

/*
  For clarity, the two unions either1 and either2 should have been
  combined into one union with a struct for internal nodes and
  struct for leaves, like the following:
  union
  {
    struct
    {
      unsigned long depth;
      GtBlindtriesnodeptr firstchild;
    } internalinfo;
    struct
    {
      unsigned long nodestartpos;
      unsigned int nodestartstopposoffset;
    } leafinfo;
  } either;
  But as unions are also aligned to word boundaries, this means for
  -m64 compilations this union requires 16 bytes, which would result in
  32 bytes for a GtBlindtrienode. Instead, in our solution either1 requires
  8 bytes and either2 requires 4 bytes which leads to 24 bytes for the entire
  GtBlindtrienode.
*/

typedef struct
{
  GtBlindtriesymbol firstchar;
  union
  {
    unsigned long depth; /* for internal nodes */
    unsigned long nodestartpos; /* for leaves */
  } either1;
  union
  {
    GtBlindtriesnodeptr firstchild; /* for internal nodes */
    unsigned int nodestartstopposoffset; /* for leaves */
  } either2;
  unsigned int rightsibling:GT_BLINDTRIE_BITSFORRIGHTSIBLING;
  unsigned int isleafbit:1;
} GtBlindtrienode;

GT_DECLAREARRAYSTRUCT(GtBlindtriesnodeptr);

struct GtBlindtrie
{
  /* The following belongs to the state and is initialized by blindtrie_new */
  GtBlindtriesnodeptr allocatedBlindtrienode,
                   nextfreeBlindtrienode;
  GtBlindtrienode *spaceBlindtrienode;
  GtArrayGtBlindtriesnodeptr stack;
  GtArrayGtUlong overflowsuffixes;

  /* The following needs to be supplied for each insertion */
  const GtEncseq *encseq;
  GtEncseqReader *esr1, *esr2;
  GtReadmode readmode;
  unsigned long totallength,
                logicaltotallength,
                sortmaxdepthminusoffset;
  bool cmpcharbychar,
       hasmirror,
       has_twobitencoding_stoppos_support;
  unsigned int nodenumberincrement;
  GtViatwobitkeyvalues *vtk1, *vtk2;
  GtSuffixsortspace *sssp;
};

static bool blindtrie_isleaf(const GtBlindtrie *blindtrie,
                             const GtBlindtriesnodeptr node)
{
  gt_assert(node < blindtrie->nextfreeBlindtrienode);
  return blindtrie->spaceBlindtrienode[node].isleafbit ? true : false;
}

static void blindtrie_setleaf(GtBlindtrie *blindtrie,
                              const GtBlindtriesnodeptr node,bool isleaf)
{
  gt_assert(node < blindtrie->nextfreeBlindtrienode);
  blindtrie->spaceBlindtrienode[node].isleafbit = isleaf ? 1U : 0;
}

static unsigned long blindtrie_getdepth(const GtBlindtrie *blindtrie,
                                        const GtBlindtriesnodeptr node)
{
  gt_assert(!blindtrie_isleaf(blindtrie,node));
  return blindtrie->spaceBlindtrienode[node].either1.depth;
}

static void blindtrie_setdepth(const GtBlindtrie *blindtrie,
                               const GtBlindtriesnodeptr node,
                               unsigned long depth)
{
  gt_assert(!blindtrie_isleaf(blindtrie,node));
  blindtrie->spaceBlindtrienode[node].either1.depth = depth;
}

static GtBlindtriesymbol blindtrie_firstchar_get(const GtBlindtrie *blindtrie,
                                               const GtBlindtriesnodeptr node)
{
  gt_assert(node < blindtrie->nextfreeBlindtrienode);
  return blindtrie->spaceBlindtrienode[node].firstchar;
}

static void blindtrie_firstchar_set(GtBlindtrie *blindtrie,
                                    GtBlindtriesnodeptr node,bool isleaf,
                                    GtBlindtriesymbol firstchar)
{
  blindtrie_setleaf(blindtrie,node,isleaf);
  gt_assert(isleaf || !GT_ISUNIQUEINT(firstchar));
  blindtrie->spaceBlindtrienode[node].firstchar = firstchar;
}

static GtBlindtriesnodeptr blindtrie_rightsibling_get(
                        const GtBlindtrie *blindtrie,
                        const GtBlindtriesnodeptr node)
{
  gt_assert(node < blindtrie->nextfreeBlindtrienode);
  return blindtrie->spaceBlindtrienode[node].rightsibling;
}

static void blindtrie_rightsibling_set(const GtBlindtrie *blindtrie,
                                       GtBlindtriesnodeptr node,
                                       GtBlindtriesnodeptr rightsibling)
{
  gt_assert(node < blindtrie->nextfreeBlindtrienode);
  blindtrie->spaceBlindtrienode[node].rightsibling = rightsibling;
}

static GtBlindtriesnodeptr blindtrie_firstchild_get(
                     const GtBlindtrie *blindtrie,
                     const GtBlindtriesnodeptr node)
{
  gt_assert(!blindtrie_isleaf(blindtrie,node));
  return blindtrie->spaceBlindtrienode[node].either2.firstchild;
}

static void blindtrie_firstchild_set(const GtBlindtrie *blindtrie,
                                     GtBlindtriesnodeptr node,
                                     GtBlindtriesnodeptr firstchild)
{
  gt_assert(!blindtrie_isleaf(blindtrie,node));
  blindtrie->spaceBlindtrienode[node].either2.firstchild
    = firstchild;
}

#define GT_TWOBITENCODINGSTARTSTOPOFFSETUNDEF UINT_MAX
#define GT_BLTRIESTOPPOSUNDEF(BLTRIE)         ((BLTRIE)->totallength)

static void blindtrie_leafinfo_set(GtBlindtrie *blindtrie,
                                   GtBlindtriesnodeptr node,
                                   unsigned long currentstartpos,
                                   unsigned long currenttwobitencodingstoppos)
{
  GtBlindtrienode *leafptr = blindtrie->spaceBlindtrienode + node;

  gt_assert(node < blindtrie->nextfreeBlindtrienode);
  leafptr->either1.nodestartpos = currentstartpos;
  if (currenttwobitencodingstoppos != GT_BLTRIESTOPPOSUNDEF(blindtrie))
  {
    if (GT_ISDIRREVERSE(blindtrie->readmode))
    {
      gt_assert(GT_REVERSEPOS(blindtrie->totallength,currentstartpos) + 1 >=
                currenttwobitencodingstoppos);
      if (currenttwobitencodingstoppos == 0)
      {
        leafptr->either2.nodestartstopposoffset = 0;
      } else
      {
        unsigned long tmp = GT_REVERSEPOS(blindtrie->totallength,
                                          currentstartpos) + 2
                                          - currenttwobitencodingstoppos;
        gt_assert(tmp > 0 &&
                  tmp < (unsigned long) GT_TWOBITENCODINGSTARTSTOPOFFSETUNDEF);
        leafptr->either2.nodestartstopposoffset = (unsigned int) tmp;
      }
    } else
    {
      gt_assert(currentstartpos <= currenttwobitencodingstoppos);
      if (currenttwobitencodingstoppos == blindtrie->totallength)
      {
        leafptr->either2.nodestartstopposoffset = 0;
      } else
      {
        unsigned long tmp = currenttwobitencodingstoppos - currentstartpos + 1;

        gt_assert(tmp > 0 &&
                  tmp < (unsigned long) GT_TWOBITENCODINGSTARTSTOPOFFSETUNDEF);
        leafptr->either2.nodestartstopposoffset = (unsigned int) tmp;
      }
    }
  } else
  {
    leafptr->either2.nodestartstopposoffset
      = GT_TWOBITENCODINGSTARTSTOPOFFSETUNDEF;
  }
}

static unsigned long blindtrie_nodestartpos_get(const GtBlindtrie *blindtrie,
                                                GtBlindtriesnodeptr node)
{
  gt_assert(blindtrie_isleaf(blindtrie,node));
  return blindtrie->spaceBlindtrienode[node].either1.nodestartpos;
}

static unsigned long blindtrie_nodestoppos_get(const GtBlindtrie *blindtrie,
                                               GtBlindtriesnodeptr node)
{
  const GtBlindtrienode *leafptr = blindtrie->spaceBlindtrienode + node;

  gt_assert(blindtrie_isleaf(blindtrie,node));
  if (leafptr->either2.nodestartstopposoffset
      != GT_TWOBITENCODINGSTARTSTOPOFFSETUNDEF)
  {
    unsigned long nodestoppos;

    if (GT_ISDIRREVERSE(blindtrie->readmode))
    {
      if (leafptr->either2.nodestartstopposoffset == 0)
      {
        nodestoppos = 0;
      } else
      {
        gt_assert(GT_REVERSEPOS(blindtrie->totallength,
                                leafptr->either1.nodestartpos)+2
                  >= (unsigned long) leafptr->either2.nodestartstopposoffset);
        nodestoppos
          = GT_REVERSEPOS(blindtrie->totallength,leafptr->either1.nodestartpos)
            + 2 - (unsigned long) leafptr->either2.nodestartstopposoffset;
      }
    } else
    {
      if (leafptr->either2.nodestartstopposoffset == 0)
      {
        nodestoppos = blindtrie->totallength;
      } else
      {
        nodestoppos = leafptr->either1.nodestartpos +
                      (unsigned long) leafptr->either2.nodestartstopposoffset
                      - 1;
      }
    }
    return nodestoppos;
  }
  return GT_BLTRIESTOPPOSUNDEF(blindtrie);
}

static void blindtrie_copy_either(GtBlindtrie *blindtrie,
                                  GtBlindtriesnodeptr destnode,
                                  GtBlindtriesnodeptr srcnode)
{
  gt_assert(destnode < blindtrie->nextfreeBlindtrienode &&
            srcnode < blindtrie->nextfreeBlindtrienode);
  blindtrie->spaceBlindtrienode[destnode].either1
    = blindtrie->spaceBlindtrienode[srcnode].either1;
  blindtrie->spaceBlindtrienode[destnode].either2
    = blindtrie->spaceBlindtrienode[srcnode].either2;
}

static bool blindtrie_isleftofboundary(const GtBlindtrie *blindtrie,
                                       unsigned long currentstartpos,
                                       unsigned long add)
{
  unsigned long endpos;

  if (blindtrie->sortmaxdepthminusoffset == 0)
  {
    endpos = blindtrie->totallength;
  } else
  {
    endpos = currentstartpos + blindtrie->sortmaxdepthminusoffset;
    if (endpos >= blindtrie->totallength)
    {
      endpos = blindtrie->totallength;
    }
  }
  return (currentstartpos + add < endpos) ? true : false;
}

static GtBlindtriesnodeptr blindtrie_newnode(GtBlindtrie *blindtrie)
{
  if (blindtrie->nextfreeBlindtrienode >= blindtrie->allocatedBlindtrienode)
  {
    gt_assert(blindtrie->nodenumberincrement >= 1U);
    gt_assert(blindtrie->allocatedBlindtrienode +
              blindtrie->nodenumberincrement <=
              (1U << GT_BLINDTRIE_BITSFORRIGHTSIBLING) - 1);
    blindtrie->allocatedBlindtrienode += blindtrie->nodenumberincrement;
    blindtrie->spaceBlindtrienode
      = gt_realloc(blindtrie->spaceBlindtrienode,
                   sizeof (*blindtrie->spaceBlindtrienode) *
                          blindtrie->allocatedBlindtrienode);
  }
  return blindtrie->nextfreeBlindtrienode++;
}

static GtBlindtriesnodeptr blindtrie_newleaf(GtBlindtrie *blindtrie,
                                          unsigned long currentstartpos,
                                          unsigned long
                                            currenttwobitencodingstoppos,
                                          GtBlindtriesymbol firstchar,
                                          GtBlindtriesnodeptr rightsibling)
{
  GtBlindtriesnodeptr newleaf;

  newleaf = blindtrie_newnode(blindtrie);
  blindtrie_firstchar_set(blindtrie,newleaf,true,firstchar);
  blindtrie_leafinfo_set(blindtrie,
                         newleaf,
                         currentstartpos,
                         currenttwobitencodingstoppos);
  blindtrie_rightsibling_set(blindtrie,newleaf,rightsibling);
  return newleaf;
}

static unsigned long blindtrie_currenttwobitencodingstoppos_get(
                                     const GtBlindtrie *blindtrie,
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
  return GT_BLTRIESTOPPOSUNDEF(blindtrie);
}

static void blindtrie_makeroot(GtBlindtrie *blindtrie,
                               unsigned long currentstartpos,
                               unsigned long currenttwobitencodingstoppos)
{
  GtBlindtriesymbol firstchar;

  gt_assert(blindtrie->nextfreeBlindtrienode == GT_BLINDTRIE_ROOTIDX);
  (void) blindtrie_newnode(blindtrie);
  blindtrie_firstchar_set(blindtrie,GT_BLINDTRIE_ROOTIDX,
                          false,0); /* firstchar of root will never be used */
  blindtrie_setdepth(blindtrie,GT_BLINDTRIE_ROOTIDX,0);
  /*blindtrie_rightsibling_set(blindtrie,root,GT_BLINDTRIE_REFNULL); */
  if (blindtrie_isleftofboundary(blindtrie,currentstartpos,0))
  {
    /* Random access */
    firstchar = (GtBlindtriesymbol)
                gt_encseq_get_encoded_char(blindtrie->encseq,
                                           currentstartpos,
                                           blindtrie->readmode);
    if (GT_BLINDTRIECHAR_ISSPECIAL(firstchar))
    {
      firstchar = GT_UNIQUEINT(currentstartpos);
    }
    if (currenttwobitencodingstoppos == ULONG_MAX)
    {
      currenttwobitencodingstoppos
        = blindtrie_currenttwobitencodingstoppos_get(blindtrie,currentstartpos);
    }
  } else
  {
    firstchar = GT_UNIQUEINT(currentstartpos);
    currenttwobitencodingstoppos = GT_BLTRIESTOPPOSUNDEF(blindtrie);
  }
  blindtrie_firstchild_set(blindtrie,GT_BLINDTRIE_ROOTIDX,
                           blindtrie_newleaf(blindtrie,currentstartpos,
                                             currenttwobitencodingstoppos,
                                             firstchar,GT_BLINDTRIE_REFNULL));
}

static GtBlindtriesnodeptr blindtrie_extractleafnode(GtBlindtrie *blindtrie,
                                                     GtBlindtriesnodeptr head)
{
  gt_assert(!blindtrie_isleaf(blindtrie,head));
  do
  {
    head = blindtrie_firstchild_get(blindtrie,head);
  } while (!blindtrie_isleaf(blindtrie,head));
  return head;
}

static int blindtrie_comparecharacters(GtBlindtriesymbol oldchar,
                                              GtBlindtriesymbol newchar)
{
  return (oldchar > newchar) ? 1 : ((oldchar < newchar) ? -1 : 0);
}

static GtBlindtriesnodeptr blindtrie_findsucc(const GtBlindtrie *blindtrie,
                                           GtBlindtriesnodeptr node,
                                           GtBlindtriesymbol newchar)
{
  int retval;

  for (;;)
  {
    retval = blindtrie_comparecharacters(
                     blindtrie_firstchar_get(blindtrie,node),
                     newchar);
    if (retval == 0)
    {              /* found branch corresponding to newchar */
      return node;
    }
    if (retval == 1)
    {               /* found sibling which is already greater than newchar */
      return GT_BLINDTRIE_REFNULL;
    }
    node = blindtrie_rightsibling_get(blindtrie,node);
    if (node == GT_BLINDTRIE_REFNULL) /* no more siblings: mismatch */
    {
      return GT_BLINDTRIE_REFNULL;
    }
  }
}

static GtBlindtriesnodeptr blindtrie_findcompanion(
                                    GtBlindtrie *blindtrie,
                                    unsigned long currentstartpos,
                                    unsigned long currenttwobitencodingstoppos)
{
  GtBlindtriesymbol newchar;
  GtBlindtriesnodeptr head, succ;
  unsigned long headdepth;

  blindtrie->stack.nextfreeGtBlindtriesnodeptr = 0;
  head = GT_BLINDTRIE_ROOTIDX;
  while (!blindtrie_isleaf(blindtrie,head))
  {
    GT_STOREINARRAY (&blindtrie->stack, GtBlindtriesnodeptr, 128, head);
    headdepth = blindtrie_getdepth(blindtrie,head);
    if (blindtrie_isleftofboundary(blindtrie,currentstartpos,headdepth))
    {
      if (!blindtrie->cmpcharbychar &&
          blindtrie->has_twobitencoding_stoppos_support)
      {
        if ((GT_ISDIRREVERSE(blindtrie->readmode) &&
            GT_REVERSEPOS(blindtrie->totallength,currentstartpos+headdepth)
            >= currenttwobitencodingstoppos) ||
            (!GT_ISDIRREVERSE(blindtrie->readmode) &&
             currentstartpos + headdepth < currenttwobitencodingstoppos))
        {
          newchar = (GtBlindtriesymbol)
                    gt_encseq_get_encoded_char_nospecial(
                                                   blindtrie->encseq,
                                                   currentstartpos + headdepth,
                                                   blindtrie->readmode);
        } else
        {
          newchar = GT_UNIQUEINT(currentstartpos + headdepth);
        }
      } else
      {
        newchar = (GtBlindtriesymbol)
                  gt_encseq_get_encoded_char(blindtrie->encseq,
                                             currentstartpos + headdepth,
                                             blindtrie->readmode);
        if (GT_BLINDTRIECHAR_ISSPECIAL(newchar))
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
      return blindtrie_extractleafnode(blindtrie,head);
    }
    succ = blindtrie_findsucc(blindtrie,
                              blindtrie_firstchild_get(blindtrie,head),
                              newchar);
    if (succ == GT_BLINDTRIE_REFNULL)
    {
      return blindtrie_extractleafnode(blindtrie,head);
    }
    head = succ;
  }
  GT_STOREINARRAY (&blindtrie->stack, GtBlindtriesnodeptr, 128, head);
  return head;
}

static void blindtrie_insertatsplitnode(GtBlindtrie *blindtrie,
                                        GtBlindtriesnodeptr oldnode,
                                        GtBlindtriesymbol mm_oldsuffix,
                                        unsigned long lcp,
                                        GtBlindtriesymbol mm_newsuffix,
                                        unsigned long currentstartpos,
                                        unsigned long
                                           currenttwobitencodingstoppos)
{
  GtBlindtriesnodeptr newleaf, newnode, previousnode, currentnode;

  gt_assert(GT_ISUNIQUEINT(mm_oldsuffix) ||
            GT_ISUNIQUEINT(mm_newsuffix) ||
            mm_oldsuffix != mm_newsuffix ||
            blindtrie_isleaf(blindtrie,oldnode) ||
            blindtrie_getdepth(blindtrie,oldnode) == lcp);

  /* insert a new node before node oldnode if necessary */
  if (blindtrie_isleaf(blindtrie,oldnode))
  {
    gt_assert(lcp > 0);
  }
  if (blindtrie_isleaf(blindtrie,oldnode) ||
      blindtrie_getdepth(blindtrie,oldnode) > lcp)
  {
    newnode = blindtrie_newnode(blindtrie);
    blindtrie_firstchar_set(blindtrie,newnode,
                            blindtrie_isleaf(blindtrie,oldnode),mm_oldsuffix);
    if (!blindtrie_isleaf(blindtrie,oldnode))
    {
      blindtrie_setdepth(blindtrie,newnode,
                         blindtrie_getdepth(blindtrie,oldnode));
      /* newnode inherits depth+children */
    }
    blindtrie_copy_either(blindtrie,newnode,oldnode);
    blindtrie_rightsibling_set(blindtrie,newnode,GT_BLINDTRIE_REFNULL);
    blindtrie_setleaf(blindtrie,oldnode,false);
    gt_assert(!GT_ISUNIQUEINT(blindtrie_firstchar_get(blindtrie,oldnode)));
    gt_assert(lcp > 0);
    blindtrie_setdepth(blindtrie,oldnode,lcp);
    /* oldnode has newnode as only child*/
    blindtrie_firstchild_set(blindtrie,oldnode,newnode);
  }
  gt_assert(blindtrie_isleaf(blindtrie,oldnode) ||
            blindtrie_getdepth(blindtrie,oldnode) == lcp);
  previousnode = GT_BLINDTRIE_REFNULL;
  currentnode = blindtrie_firstchild_get(blindtrie,oldnode);
  while (currentnode != GT_BLINDTRIE_REFNULL &&
         blindtrie_comparecharacters(blindtrie_firstchar_get(blindtrie,
                                                             currentnode),
                                     mm_newsuffix) < 0)
  {
    previousnode = currentnode;
    currentnode = blindtrie_rightsibling_get(blindtrie,currentnode);
  }
  /* insert new leaf with current suffix */
  /* search S[lcp] among the offsprings */
  newleaf = blindtrie_newleaf(blindtrie,currentstartpos,
                              currenttwobitencodingstoppos,mm_newsuffix,
                              currentnode);
  if (previousnode != GT_BLINDTRIE_REFNULL)
  {
    blindtrie_rightsibling_set(blindtrie,previousnode,newleaf);
  } else
  {
    blindtrie_firstchild_set(blindtrie,oldnode,newleaf);
  }
}

static unsigned long blindtrie_cmpcharbychar_getlcp(
                                 GtBlindtriesymbol *mm_oldsuffix,
                                 bool *mm_oldsuffixisseparator,
                                 GtBlindtriesymbol *mm_newsuffix,
                                 const GtBlindtrie *blindtrie,
                                 unsigned long leafpos,
                                 unsigned long currentstartpos)
{
  unsigned long lcp;
  GtBlindtriesymbol cc1, cc2;

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
      cc1 = (GtBlindtriesymbol)
            gt_encseq_reader_next_encoded_char(blindtrie->esr1);
      if (GT_BLINDTRIECHAR_ISSPECIAL(cc1))
      {
        if (mm_oldsuffixisseparator != NULL)
        {
          *mm_oldsuffixisseparator = (cc1 == SEPARATOR) ? true : false;
        }
        cc1 = GT_UNIQUEINT(leafpos + lcp);
      }
    } else
    {
      if (mm_oldsuffixisseparator != NULL)
      {
        *mm_oldsuffixisseparator = true;
      }
      cc1 = GT_UNIQUEINT(leafpos + lcp);
    }
    if (blindtrie_isleftofboundary(blindtrie,currentstartpos,lcp))
    {
      cc2 = (GtBlindtriesymbol)
            gt_encseq_reader_next_encoded_char(blindtrie->esr2);
      if (GT_BLINDTRIECHAR_ISSPECIAL(cc2))
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
  gt_assert(blindtrie->sortmaxdepthminusoffset == 0 ||
            lcp <= blindtrie->sortmaxdepthminusoffset);
  return lcp;
}

static unsigned long blindtrie_twobitencoding_getlcp(
                                 GtBlindtriesymbol *mm_oldsuffix,
                                 bool *mm_oldsuffixisseparator,
                                 GtBlindtriesymbol *mm_newsuffix,
                                 const GtBlindtrie *blindtrie,
                                 unsigned long leafpos,
                                 unsigned long leaftwobitencodingstoppos,
                                 unsigned long currentstartpos,
                                 unsigned long currenttwobitencodingstoppos)
{
  GtCommonunits commonunits;
  const unsigned long depth = 0;

  gt_assert(leafpos != currentstartpos);
  gt_Viatwobitkeyvalues_reinit(blindtrie->vtk1,
                               blindtrie->encseq,
                               blindtrie->readmode,
                               NULL, /* esr */
                               leafpos,
                               depth,
                               blindtrie->sortmaxdepthminusoffset,
                               leaftwobitencodingstoppos);
  gt_Viatwobitkeyvalues_reinit(blindtrie->vtk2,
                               blindtrie->encseq,
                               blindtrie->readmode,
                               NULL, /* esr */
                               currentstartpos,
                               depth,
                               blindtrie->sortmaxdepthminusoffset,
                               currenttwobitencodingstoppos);
  (void) gt_encseq_twobitencoding_strcmp(&commonunits,
                                         blindtrie->encseq,
                                         blindtrie->encseq,
                                         blindtrie->readmode,
                                         depth,
                                         blindtrie->sortmaxdepthminusoffset,
                                         blindtrie->vtk1,
                                         blindtrie->vtk2);
  if (blindtrie_isleftofboundary(blindtrie,leafpos,commonunits.finaldepth) &&
      !commonunits.leftspecial)
  {
    *mm_oldsuffix = (GtBlindtriesymbol)
                    gt_encseq_get_encoded_char_nospecial(blindtrie->encseq,
                                                         leafpos +
                                                         commonunits.finaldepth,
                                                         blindtrie->readmode);
  } else
  {
    if (mm_oldsuffixisseparator != NULL)
    {
      if (blindtrie_isleftofboundary(blindtrie,leafpos,commonunits.finaldepth))
      {
        *mm_oldsuffixisseparator = true;
      } else
      {
        gt_assert(commonunits.leftspecial);
        /* be careful: if wildcards occur in the sequence,
           the following statement may be wrong. Please contact stefan
           if program is applied for such sequences. */
        *mm_oldsuffixisseparator = true;
      }
    }
    *mm_oldsuffix = GT_UNIQUEINT(leafpos + commonunits.finaldepth);
  }
  if (blindtrie_isleftofboundary(blindtrie,currentstartpos,
                                 commonunits.finaldepth) &&
      !commonunits.rightspecial)
  {
    *mm_newsuffix = (GtBlindtriesymbol)
                    gt_encseq_get_encoded_char_nospecial(blindtrie->encseq,
                                                         currentstartpos +
                                                         commonunits.finaldepth,
                                                         blindtrie->readmode);
  } else
  {
    *mm_newsuffix = GT_UNIQUEINT(currentstartpos + commonunits.finaldepth);
  }
  return commonunits.finaldepth;
}

static unsigned long blindtrie_getlcp(GtBlindtriesymbol *mm_oldsuffix,
                                      bool *mm_oldsuffixisseparator,
                                      GtBlindtriesymbol *mm_newsuffix,
                                      const GtBlindtrie *blindtrie,
                                      const GtBlindtriesnodeptr lis,
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
                                           blindtrie_nodestartpos_get(blindtrie,
                                                                      lis),
                                           currentstartpos);
  }
  return blindtrie_twobitencoding_getlcp(mm_oldsuffix,
                                         mm_oldsuffixisseparator,
                                         mm_newsuffix,
                                         blindtrie,
                                         blindtrie_nodestartpos_get(blindtrie,
                                                                    lis),
                                         blindtrie_nodestoppos_get(blindtrie,
                                                                   lis),
                                         currentstartpos,
                                         currenttwobitencodingstoppos);
}

static void blindtrie_suffixout(GtBlindtrie *blindtrie,
                                unsigned long subbucketleft,
                                unsigned long offset,
                                unsigned long nextfree,
                                unsigned long startpos)
{
  gt_suffixsortspace_set(blindtrie->sssp,subbucketleft,nextfree,
                         startpos - offset);
}

#define BLINDTRIE_SETCURRENTNODE(NODEPTR)\
        currentnodeisleaf = blindtrie_isleaf(blindtrie,NODEPTR) ? true : false;\
        currentnode = NODEPTR

static unsigned long blindtrie_enumeratetrieleaves (
                           GtBlindtrie *blindtrie,
                           unsigned long subbucketleft,
                           unsigned long offset,
                           unsigned long sortmaxdepth,
                           GtLcpvalues *tableoflcpvalues,
                           void *voiddcov,
                           GtProcessunsortedsuffixrange
                             processunsortedsuffixrange)
{
  bool readyforpop = false, currentnodeisleaf;
  GtBlindtriesnodeptr currentnode, siblval, lcpnode = GT_BLINDTRIE_ROOTIDX;
  unsigned long nextfree = 0, equalsrangewidth = 0, lcpnodedepth,
                bucketleftidxplussubbucketleft;

  blindtrie->stack.nextfreeGtBlindtriesnodeptr = 0;
  GT_STOREINARRAY (&blindtrie->stack, GtBlindtriesnodeptr, 128,
                  GT_BLINDTRIE_ROOTIDX);
  BLINDTRIE_SETCURRENTNODE(blindtrie_firstchild_get(blindtrie,
                           GT_BLINDTRIE_ROOTIDX));
  bucketleftidxplussubbucketleft
    = gt_suffixsortspace_bucketleftidx_get(blindtrie->sssp) + subbucketleft;
  for (;;)
  {
    lcpnodedepth = blindtrie_getdepth(blindtrie,lcpnode);
    if (currentnodeisleaf)
    {
      if (nextfree > 0)
      {
        if (tableoflcpvalues != NULL)
        {
          gt_lcptab_update(tableoflcpvalues,subbucketleft,nextfree,
                           lcpnodedepth + offset);
        }
        if (sortmaxdepth > 0)
        {
          if (lcpnodedepth + offset == sortmaxdepth)
          {
            equalsrangewidth++;
          } else
          {
            gt_assert(lcpnodedepth + offset < sortmaxdepth);
            if (equalsrangewidth > 0)
            {
              if (processunsortedsuffixrange != NULL)
              {
                processunsortedsuffixrange(
                               voiddcov,
                               bucketleftidxplussubbucketleft
                                 + nextfree - 1 - equalsrangewidth,
                               equalsrangewidth + 1,
                               sortmaxdepth);
              }
              equalsrangewidth = 0;
            }
          }
        }
      }
      blindtrie_suffixout(blindtrie,subbucketleft,offset,nextfree,
                          blindtrie_nodestartpos_get(blindtrie,currentnode));
      nextfree++;
      siblval = blindtrie_rightsibling_get(blindtrie,currentnode);
      if (siblval == GT_BLINDTRIE_REFNULL)
      {
        readyforpop = true;
        currentnodeisleaf = false; /* STATE 1 */
      } else
      {
        BLINDTRIE_SETCURRENTNODE(siblval);  /* current comes from brother */
        lcpnode = blindtrie->stack.spaceGtBlindtriesnodeptr[
                             blindtrie->stack.nextfreeGtBlindtriesnodeptr-1];
      }
    } else
    {
      if (readyforpop)
      {
        if (blindtrie->stack.nextfreeGtBlindtriesnodeptr == 1UL)
        {
          break;
        }
        blindtrie->stack.nextfreeGtBlindtriesnodeptr--;
        siblval = blindtrie_rightsibling_get(blindtrie,
                       blindtrie->stack.spaceGtBlindtriesnodeptr[
                       blindtrie->stack.nextfreeGtBlindtriesnodeptr]);
        if (siblval != GT_BLINDTRIE_REFNULL)
        {
          BLINDTRIE_SETCURRENTNODE(siblval);   /* current comes from brother */
          lcpnode = blindtrie->stack.spaceGtBlindtriesnodeptr[
                             blindtrie->stack.nextfreeGtBlindtriesnodeptr - 1];
          readyforpop = false;
        }
      } else
      {
        GT_STOREINARRAY (&blindtrie->stack, GtBlindtriesnodeptr, 128,
                         currentnode);
        BLINDTRIE_SETCURRENTNODE(blindtrie_firstchild_get(blindtrie,
                                                          currentnode));
      }
    }
  }
  if (nextfree > 0 && equalsrangewidth > 0)
  {
    if (processunsortedsuffixrange != NULL)
    {
      processunsortedsuffixrange(voiddcov,
                                 bucketleftidxplussubbucketleft
                                   + nextfree - 1 - equalsrangewidth,
                                 equalsrangewidth + 1,
                                 sortmaxdepth);
    }
    equalsrangewidth = 0;
  }
  return nextfree;
}

GtBlindtrie *gt_blindtrie_new(GtSuffixsortspace *suffixsortspace,
                              unsigned long maxnumofsuffixes,
                              unsigned int nodenumberincrement,
                              const GtEncseq *encseq,
                              bool cmpcharbychar,
                              GtEncseqReader *esr1,
                              GtEncseqReader *esr2,
                              GtReadmode readmode)
{
  GtBlindtrie *blindtrie;

  if (maxnumofsuffixes == 1UL)
  {
    return NULL;
  }
  blindtrie = gt_malloc(sizeof (*blindtrie));
  if (nodenumberincrement == 0)
  {
    gt_assert(maxnumofsuffixes >= 2UL);
    blindtrie->allocatedBlindtrienode
      = (GtBlindtriesnodeptr) GT_MULT2(maxnumofsuffixes + 1) + 1;
    blindtrie->nodenumberincrement = 0;
  } else
  {
    gt_assert(maxnumofsuffixes == 0);
    blindtrie->allocatedBlindtrienode = 0;
    blindtrie->nodenumberincrement = nodenumberincrement;
  }
  blindtrie->spaceBlindtrienode
    = gt_malloc(sizeof (*blindtrie->spaceBlindtrienode) *
                blindtrie->allocatedBlindtrienode);
  /*
  printf("# sizeof (blindtrie)=%lu\n",
            (unsigned long) (sizeof (GtBlindtrie) +
                             blindtrie->allocatedBlindtrienode *
                             sizeof (GtBlindtrienode)));
  */
  GT_INITARRAY (&blindtrie->overflowsuffixes, GtUlong);
  GT_INITARRAY (&blindtrie->stack, GtBlindtriesnodeptr);
  blindtrie->nextfreeBlindtrienode = 0;
  blindtrie->encseq = encseq;
  blindtrie->has_twobitencoding_stoppos_support
    = gt_encseq_has_twobitencoding_stoppos_support(encseq);
  blindtrie->readmode = readmode;
  blindtrie->esr1 = esr1;
  blindtrie->esr2 = esr2;
  blindtrie->totallength = gt_encseq_total_length(encseq);
  blindtrie->cmpcharbychar = cmpcharbychar;
  blindtrie->sssp = suffixsortspace;
  if (cmpcharbychar)
  {
    blindtrie->vtk1 = blindtrie->vtk2 = NULL;
  } else
  {
    blindtrie->vtk1 = gt_Viatwobitkeyvalues_new();
    blindtrie->vtk2 = gt_Viatwobitkeyvalues_new();
  }
  return blindtrie;
}

void gt_blindtrie_resize(GtBlindtrie *blindtrie,unsigned int maxnumofnodes)
{
  gt_assert(maxnumofnodes > 0);
  if (blindtrie->allocatedBlindtrienode > maxnumofnodes)
  {
    blindtrie->allocatedBlindtrienode = maxnumofnodes;
    blindtrie->spaceBlindtrienode
      = gt_realloc(blindtrie->spaceBlindtrienode,
                   sizeof (*blindtrie->spaceBlindtrienode) *
                   blindtrie->allocatedBlindtrienode);
  }
}

size_t gt_blindtrie_size(unsigned long maxnumofsuffixes)
{
  if (maxnumofsuffixes == 1UL)
  {
    return 0;
  } else
  {
    return sizeof (GtBlindtrie) +
           (GT_MULT2(maxnumofsuffixes + 1) + 1) * sizeof (GtBlindtrienode);
  }
}

size_t gt_blindtrie_current_size(const GtBlindtrie *blindtrie)
{
  return sizeof (GtBlindtrie) +
         sizeof (GtBlindtrienode) * blindtrie->allocatedBlindtrienode;
}

void gt_blindtrie_reset(GtBlindtrie *blindtrie)
{
  blindtrie->nextfreeBlindtrienode = 0;
  blindtrie->stack.nextfreeGtBlindtriesnodeptr = 0;
}

void gt_blindtrie_delete(GtBlindtrie *blindtrie)
{
  if (blindtrie == NULL)
  {
    return;
  }
  gt_Viatwobitkeyvalues_delete(blindtrie->vtk1);
  gt_Viatwobitkeyvalues_delete(blindtrie->vtk2);
  gt_free(blindtrie->spaceBlindtrienode);
  GT_FREEARRAY(&blindtrie->overflowsuffixes, GtUlong);
  GT_FREEARRAY(&blindtrie->stack, GtBlindtriesnodeptr);
  gt_free(blindtrie);
}

#undef SKDEBUG
#ifdef SKDEBUG

#define NODENUM(PTR) PTR

static void gt_blindtrie_showleaf(const GtBlindtrie *blindtrie,
                                  unsigned int level,
                                  GtBlindtriesnodeptr current)
{
  printf("%*.*s",(int) (6 * level),(int) (6 * level)," ");
  gt_assert(current != GT_BLINDTRIE_REFNULL);
  printf("Leaf(%u,firstchar=%u,startpos=%lu,rightsibling=%u)\n",
         NODENUM(current),
         (unsigned int) blindtrie_firstchar_get(blindtrie,current),
         blindtrie_nodestartpos_get(blindtrie,current),
         NODENUM(blindtrie_rightsibling_get(blindtrie,current)));
}

static void gt_blindtrie_showintern(const GtBlindtrie *blindtrie,
                                    unsigned int level,
                                    GtBlindtriesnodeptr current)
{
  printf("%*.*s",(int) (6 * level),(int) (6 * level)," ");
  gt_assert(current != GT_BLINDTRIE_REFNULL);
  printf("Intern(%u,firstchar=%lu,depth=%lu"
         ",firstchild=%u,rightsibling=%u)\n",
          NODENUM(current),
          blindtrie_firstchar_get(blindtrie,current),
          blindtrie_getdepth(blindtrie,current),
          NODENUM(blindtrie_firstchild_get(blindtrie,current)),
          NODENUM(blindtrie_rightsibling_get(blindtrie,current)));
}

static void gt_blindtrie_showrecursive(const GtBlindtrie *blindtrie,
                                       unsigned int level,
                                       GtBlindtriesnodeptr node)
{
  GtBlindtriesnodeptr current;

  for (current = blindtrie_firstchild_get(blindtrie,node);
       current != GT_BLINDTRIE_REFNULL;
       current = blindtrie_rightsibling_get(blindtrie,current))
  {
    if (blindtrie_isleaf(blindtrie,current))
    {
      gt_blindtrie_showleaf(blindtrie,level,current);
    } else
    {
      gt_blindtrie_showintern(blindtrie,level,current);
      gt_blindtrie_showrecursive(blindtrie,level+1,current);
    }
  }
}

static void gt_blindtrie_show(const GtBlindtrie *blindtrie)
{
  gt_blindtrie_showrecursive(blindtrie,0,GT_BLINDTRIE_ROOTIDX);
}

static void blindtrie_showstate(const GtBlindtrie *blindtrie,
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

static GtBlindtriesnodeptr blindtrie_findsplitnode(const GtBlindtrie *blindtrie,
                                                unsigned long lcp)
{
  GtBlindtriesnodeptr currentnode;
  unsigned long stackidx;

  currentnode = GT_BLINDTRIE_ROOTIDX;
  for (stackidx=0;stackidx<blindtrie->stack.nextfreeGtBlindtriesnodeptr;
       stackidx++)
  {
    currentnode = blindtrie->stack.spaceGtBlindtriesnodeptr[stackidx];
    if (blindtrie_isleaf(blindtrie,currentnode) ||
        blindtrie_getdepth(blindtrie,currentnode) >= lcp)
    {
      break;
    }
  }
  return currentnode;
}

static void gt_blindtrie_insertsuffix(GtBlindtrie *blindtrie,
                                      unsigned long offset,
                                      unsigned long sortmaxdepth,
                                      unsigned long currentstartpos)
{

  gt_assert(sortmaxdepth == 0 || sortmaxdepth > offset);
  if (blindtrie->nextfreeBlindtrienode == 0)
  { /* empty tree */
    if (sortmaxdepth == 0)
    {
      blindtrie->sortmaxdepthminusoffset = 0;
    } else
    {
      blindtrie->sortmaxdepthminusoffset = sortmaxdepth - offset;
    }
    blindtrie->overflowsuffixes.nextfreeGtUlong = 0;
    blindtrie_makeroot(blindtrie,currentstartpos,ULONG_MAX);
  } else
  {
    if (blindtrie_isleftofboundary(blindtrie,currentstartpos,0))
    {
      unsigned long lcp, currenttwobitencodingstoppos;
      GtBlindtriesnodeptr leafinsubtrie, splitnode;
      GtBlindtriesymbol mm_oldsuffix, mm_newsuffix;

      currenttwobitencodingstoppos
        = blindtrie_currenttwobitencodingstoppos_get(blindtrie,currentstartpos);
      leafinsubtrie = blindtrie_findcompanion(blindtrie,currentstartpos,
                                              currenttwobitencodingstoppos);
      gt_assert(blindtrie_isleaf(blindtrie,leafinsubtrie));
      lcp = blindtrie_getlcp(&mm_oldsuffix,
                             NULL,
                             &mm_newsuffix,
                             blindtrie,
                             leafinsubtrie,
                             currentstartpos,
                             currenttwobitencodingstoppos);
      splitnode = blindtrie_findsplitnode(blindtrie,lcp);
      blindtrie_insertatsplitnode(blindtrie,
                                  splitnode,
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

static int blindtrie_compare_ascending(const void *a,const void *b)
{
  unsigned long *aptr = (unsigned long *) a;
  unsigned long *bptr = (unsigned long *) b;
  gt_assert(*aptr != *bptr);
  return *aptr < *bptr ? -1 : 1;
}

static void processoverflowsuffixes(GtBlindtrie *blindtrie,
                                    unsigned long offset,
                                    GtLcpvalues *tableoflcpvalues,
                                    unsigned long subbucketleft,
                                    unsigned long nextsuffixtooutput)
{
  if (blindtrie->overflowsuffixes.nextfreeGtUlong > 0)
  {
    unsigned long idx;

    if (blindtrie->overflowsuffixes.nextfreeGtUlong > 1UL)
    {
      qsort(blindtrie->overflowsuffixes.spaceGtUlong,
            (size_t) blindtrie->overflowsuffixes.nextfreeGtUlong,
            sizeof (*blindtrie->overflowsuffixes.spaceGtUlong),
            blindtrie_compare_ascending);
    }
    for (idx = 0; idx < blindtrie->overflowsuffixes.nextfreeGtUlong; idx++)
    {
      if (tableoflcpvalues != NULL)
      {
        gt_lcptab_update(tableoflcpvalues,subbucketleft,nextsuffixtooutput,
                         offset);
      }
      blindtrie_suffixout(blindtrie,subbucketleft,offset,nextsuffixtooutput,
                          blindtrie->overflowsuffixes.spaceGtUlong[idx]);
      nextsuffixtooutput++;
    }
  }
}

static void gt_blindtrie2sorting(GtBlindtrie *blindtrie,
                                 unsigned long subbucketleft,
                                 GtLcpvalues *tableoflcpvalues,
                                 unsigned long offset,
                                 unsigned long sortmaxdepth,
                                 void *voiddcov,
                                 GtProcessunsortedsuffixrange
                                   processunsortedsuffixrange)
{
  unsigned long nextsuffixtosort;

  nextsuffixtosort
    = blindtrie_enumeratetrieleaves (blindtrie,subbucketleft,offset,
                                     sortmaxdepth,tableoflcpvalues,
                                     voiddcov,processunsortedsuffixrange);
  processoverflowsuffixes(blindtrie,
                          offset,
                          tableoflcpvalues,
                          subbucketleft,
                          nextsuffixtosort);
}

void gt_blindtrie_suffixsort(GtBlindtrie *blindtrie,
                             unsigned long subbucketleft,
                             GtLcpvalues *tableoflcpvalues,
                             unsigned long numberofsuffixes,
                             unsigned long offset,
                             unsigned long sortmaxdepth,
                             void *voiddcov,
                             GtProcessunsortedsuffixrange
                               processunsortedsuffixrange)
{
  unsigned long idx, currentstartpos;

  gt_blindtrie_reset(blindtrie);
  for (idx=0; idx < numberofsuffixes; idx++)
  {
    currentstartpos = gt_suffixsortspace_get(blindtrie->sssp,subbucketleft,idx);
    gt_blindtrie_insertsuffix(blindtrie,
                              offset,
                              sortmaxdepth,
                              currentstartpos + offset);
#ifdef SKDEBUG
    blindtrie_showstate(blindtrie,
                        subbucketleft,
                        idx+1,
                        offset);
#endif
  }
  gt_blindtrie2sorting(blindtrie,
                       subbucketleft,
                       tableoflcpvalues,
                       offset,
                       sortmaxdepth,
                       voiddcov,
                       processunsortedsuffixrange);
}

bool gt_blindtrie_retrieve(GtBlindtrie *blindtrie,
                           unsigned long currentstartpos,
                           unsigned long currenttwobitencodingstoppos)
{
  gt_assert(!blindtrie->cmpcharbychar);
  if (currenttwobitencodingstoppos != ULONG_MAX)
  {
    gt_assert(blindtrie->esr1 == NULL);
    gt_assert(blindtrie->esr2 == NULL);
  }
  if (blindtrie->nextfreeBlindtrienode == 0)
  {
    blindtrie->sortmaxdepthminusoffset = 0;
    blindtrie_makeroot(blindtrie,currentstartpos,currenttwobitencodingstoppos);
    return false;
  } else
  {
    GtBlindtriesnodeptr leafinsubtrie, splitnode;
    GtBlindtriesymbol mm_oldsuffix, mm_newsuffix;
    bool mm_oldsuffixisseparator;
    unsigned long lcp;

    if (currenttwobitencodingstoppos == ULONG_MAX)
    {
      currenttwobitencodingstoppos
        = blindtrie_currenttwobitencodingstoppos_get(blindtrie,currentstartpos);
    }
    leafinsubtrie = blindtrie_findcompanion(blindtrie,currentstartpos,
                                            currenttwobitencodingstoppos);
    gt_assert(blindtrie_isleaf(blindtrie,leafinsubtrie));
    lcp = blindtrie_getlcp(&mm_oldsuffix,
                           &mm_oldsuffixisseparator,
                           &mm_newsuffix,
                           blindtrie,
                           leafinsubtrie,
                           currentstartpos,
                           currenttwobitencodingstoppos);
    splitnode = blindtrie_findsplitnode(blindtrie,lcp);
    if (blindtrie_isleaf(blindtrie,splitnode) && GT_ISUNIQUEINT(mm_oldsuffix) &&
        mm_oldsuffixisseparator)
    {
      return true;
    }
    blindtrie_insertatsplitnode(blindtrie,
                                splitnode,
                                mm_oldsuffix,
                                lcp,
                                mm_newsuffix,
                                currentstartpos,
                                currenttwobitencodingstoppos);
    return false;
  }
}
