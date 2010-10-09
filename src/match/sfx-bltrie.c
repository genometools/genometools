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

#define BLINDTRIECHAR_ISSPECIAL(X)     ((X) >= (Blindtriesymbol) WILDCARD)
#define BLINDTRIE_REFNULL              0
#define BLINDTRIE_BITSFORRIGHTSIBLING  31

typedef unsigned long Blindtriesymbol;
typedef unsigned int Blindtrienodeptr;

/*
  For clarity, the two unions either1 and either2 should have been
  combined into one union with a struct for internal nodes and
  struct for leaves, like the following:
  union
  {
    struct
    {
      unsigned long depth;
      Blindtrienodeptr firstchild;
    } internalinfo;
    struct Leafinfo
    {
      unsigned long nodestartpos;
      unsigned int nodestartstopposoffset;
    } leafinfo;
  } either;
  But as unions are also aligned to work boundaries, this means for
  -m64 compilations this union requires 16 bytes, which would result in
  32 bytes for a Blindtrienode. Instead, in our solution either1 requires
  8 bytes and either2 requires 4 bytes which leads to 24 bytes for the entire
  Blindtrienode.
*/

typedef struct
{
  Blindtriesymbol firstchar;
  union
  {
    unsigned long depth; /* for internal nodes */
    unsigned long nodestartpos; /* for leaves */
  } either1;
  union
  {
    Blindtrienodeptr firstchild;         /* for internal nodes */
    unsigned int nodestartstopposoffset; /* for leaves */
  } either2;
  unsigned int rightsibling:BLINDTRIE_BITSFORRIGHTSIBLING;
  unsigned int isleafbit:1;
} Blindtrienode;

GT_DECLAREARRAYSTRUCT(Blindtrienodeptr);

struct Blindtrie
{
  /* The following belongs to the state and is initialized by blindtrie_new */
  Blindtrienodeptr allocatedBlindtrienode,
                   nextfreeBlindtrienode;
  Blindtrienode *spaceBlindtrienode;
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
  unsigned int nodenumberincrement;
  GtSuffixsortspace *sssp;
};

static bool blindtrie_isleaf(const Blindtrie *blindtrie,
                             const Blindtrienodeptr node)
{
  return blindtrie->spaceBlindtrienode[node].isleafbit ? true : false;
}

static void blindtrie_setleaf(Blindtrie *blindtrie,
                              const Blindtrienodeptr node,bool isleaf)
{
  blindtrie->spaceBlindtrienode[node].isleafbit = isleaf ? 1U : 0;
}

static unsigned long blindtrie_getdepth(const Blindtrie *blindtrie,
                                        const Blindtrienodeptr node)
{
  gt_assert(!blindtrie_isleaf(blindtrie,node));
  return blindtrie->spaceBlindtrienode[node].either1.depth;
}

static void blindtrie_setdepth(const Blindtrie *blindtrie,
                               const Blindtrienodeptr node,
                               unsigned long depth)
{
  gt_assert(!blindtrie_isleaf(blindtrie,node));
  blindtrie->spaceBlindtrienode[node].either1.depth = depth;
}

static Blindtriesymbol blindtrie_firstchar_get(const Blindtrie *blindtrie,
                                               const Blindtrienodeptr node)
{
  return blindtrie->spaceBlindtrienode[node].firstchar;
}

static void blindtrie_firstchar_set(Blindtrie *blindtrie,
                                    Blindtrienodeptr node,bool isleaf,
                                    Blindtriesymbol firstchar)
{
  blindtrie_setleaf(blindtrie,node,isleaf);
  gt_assert(isleaf || !GT_ISUNIQUEINT(firstchar));
  blindtrie->spaceBlindtrienode[node].firstchar = firstchar;
}

static Blindtrienodeptr blindtrie_rightsibling_get(const Blindtrie *blindtrie,
                                                   const Blindtrienodeptr node)
{
  return blindtrie->spaceBlindtrienode[node].rightsibling;
}

static void blindtrie_rightsibling_set(const Blindtrie *blindtrie,
                                       Blindtrienodeptr node,
                                       Blindtrienodeptr rightsibling)
{
  blindtrie->spaceBlindtrienode[node].rightsibling = rightsibling;
}

static Blindtrienodeptr blindtrie_firstchild_get(
                     const Blindtrie *blindtrie,
                     const Blindtrienodeptr node)
{
  gt_assert(!blindtrie_isleaf(blindtrie,node));
  return blindtrie->spaceBlindtrienode[node].either2.firstchild;
}

static void blindtrie_firstchild_set(const Blindtrie *blindtrie,
                                     Blindtrienodeptr node,
                                     Blindtrienodeptr firstchild)
{
  gt_assert(!blindtrie_isleaf(blindtrie,node));
  blindtrie->spaceBlindtrienode[node].either2.firstchild
    = firstchild;
}

#define GT_TWOBITENCODINGSTARTSTOPOFFSETUNDEF UINT_MAX

static void blindtrie_leafinfo_set(Blindtrie *blindtrie,
                                   Blindtrienodeptr node,
                                   unsigned long currentstartpos,
                                   unsigned long currenttwobitencodingstoppos)
{
  Blindtrienode *leafptr = blindtrie->spaceBlindtrienode + node;

  leafptr->either1.nodestartpos = currentstartpos;
  if (currenttwobitencodingstoppos != GT_TWOBITENCODINGSTOPPOSUNDEF(blindtrie))
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

static unsigned long blindtrie_nodestartpos_get(const Blindtrie *blindtrie,
                                                Blindtrienodeptr node)
{
  gt_assert(blindtrie_isleaf(blindtrie,node));
  return blindtrie->spaceBlindtrienode[node].either1.nodestartpos;
}

static unsigned long blindtrie_nodestoppos_get(const Blindtrie *blindtrie,
                                               Blindtrienodeptr node)
{
  const Blindtrienode *leafptr = blindtrie->spaceBlindtrienode + node;

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
  return GT_TWOBITENCODINGSTOPPOSUNDEF(blindtrie);
}

static void blindtrie_copy_either(Blindtrie *blindtrie,
                                  Blindtrienodeptr destnode,
                                  Blindtrienodeptr srcnode)
{
  blindtrie->spaceBlindtrienode[destnode].either1
    = blindtrie->spaceBlindtrienode[srcnode].either1;
  blindtrie->spaceBlindtrienode[destnode].either2
    = blindtrie->spaceBlindtrienode[srcnode].either2;
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
  if (blindtrie->nextfreeBlindtrienode >= blindtrie->allocatedBlindtrienode)
  {
    gt_assert(blindtrie->nodenumberincrement >= 1U);
    gt_assert(blindtrie->allocatedBlindtrienode +
              blindtrie->nodenumberincrement <=
              (1U << BLINDTRIE_BITSFORRIGHTSIBLING) - 1);
    blindtrie->allocatedBlindtrienode += blindtrie->nodenumberincrement;
    blindtrie->spaceBlindtrienode
      = gt_realloc(blindtrie->spaceBlindtrienode,
                   sizeof (*blindtrie->spaceBlindtrienode) *
                          blindtrie->allocatedBlindtrienode);
  }
  return blindtrie->nextfreeBlindtrienode++;
}

static Blindtrienodeptr blindtrie_newleaf(Blindtrie *blindtrie,
                                          unsigned long currentstartpos,
                                          unsigned long
                                            currenttwobitencodingstoppos,
                                          Blindtriesymbol firstchar,
                                          Blindtrienodeptr rightsibling)
{
  Blindtrienodeptr newleaf;

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

#define ROOTIDX 0

static void blindtrie_makeroot(Blindtrie *blindtrie,
                               unsigned long currentstartpos)
{
  Blindtriesymbol firstchar;
  unsigned long currenttwobitencodingstoppos;

  gt_assert(blindtrie->nextfreeBlindtrienode == ROOTIDX);
  (void) blindtrie_newnode(blindtrie);
  blindtrie_firstchar_set(blindtrie,ROOTIDX,false,0); /* firstchar of root will
                                                         never be used */
  blindtrie_setdepth(blindtrie,ROOTIDX,0);
  /*blindtrie_rightsibling_set(blindtrie,root,BLINDTRIE_REFNULL); */
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
  blindtrie_firstchild_set(blindtrie,ROOTIDX,
                           blindtrie_newleaf(blindtrie,currentstartpos,
                                             currenttwobitencodingstoppos,
                                             firstchar,BLINDTRIE_REFNULL));
}

static Blindtrienodeptr blindtrie_extractleafnode(Blindtrie *blindtrie,
                                                         Blindtrienodeptr head)
{
  gt_assert(!blindtrie_isleaf(blindtrie,head));
  do
  {
    head = blindtrie_firstchild_get(blindtrie,head);
  } while (!blindtrie_isleaf(blindtrie,head));
  return head;
}

static int blindtrie_comparecharacters(Blindtriesymbol oldchar,
                                              Blindtriesymbol newchar)
{
  return (oldchar > newchar) ? 1 : ((oldchar < newchar) ? -1 : 0);
}

static Blindtrienodeptr blindtrie_findsucc(const Blindtrie *blindtrie,
                                           Blindtrienodeptr node,
                                           Blindtriesymbol newchar)
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
      return BLINDTRIE_REFNULL;
    }
    node = blindtrie_rightsibling_get(blindtrie,node);
    if (node == BLINDTRIE_REFNULL) /* no more siblings: mismatch */
    {
      return BLINDTRIE_REFNULL;
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
  head = ROOTIDX;
  while (!blindtrie_isleaf(blindtrie,head))
  {
    GT_STOREINARRAY (&blindtrie->stack, Blindtrienodeptr, 128, head);
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
          newchar = (Blindtriesymbol)
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
      return blindtrie_extractleafnode(blindtrie,head);
    }
    succ = blindtrie_findsucc(blindtrie,
                              blindtrie_firstchild_get(blindtrie,head),
                              newchar);
    if (succ == BLINDTRIE_REFNULL)
    {
      return blindtrie_extractleafnode(blindtrie,head);
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
  Blindtrienodeptr newleaf, newnode, previousnode, currentnode;

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
    blindtrie_rightsibling_set(blindtrie,newnode,BLINDTRIE_REFNULL);
    blindtrie_setleaf(blindtrie,oldnode,false);
    gt_assert(!GT_ISUNIQUEINT(blindtrie_firstchar_get(blindtrie,oldnode)));
    gt_assert(lcp > 0);
    blindtrie_setdepth(blindtrie,oldnode,lcp);
    /* oldnode has newnode as only child*/
    blindtrie_firstchild_set(blindtrie,oldnode,newnode);
  }
  gt_assert(blindtrie_isleaf(blindtrie,oldnode) ||
            blindtrie_getdepth(blindtrie,oldnode) == lcp);
  previousnode = BLINDTRIE_REFNULL;
  currentnode = blindtrie_firstchild_get(blindtrie,oldnode);
  while (currentnode != BLINDTRIE_REFNULL &&
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
  if (previousnode != BLINDTRIE_REFNULL)
  {
    blindtrie_rightsibling_set(blindtrie,previousnode,newleaf);
  } else
  {
    blindtrie_firstchild_set(blindtrie,oldnode,newleaf);
  }
}

static unsigned long blindtrie_cmpcharbychar_getlcp(
                                 Blindtriesymbol *mm_oldsuffix,
                                 bool *mm_oldsuffixisseparator,
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
                                 bool *mm_oldsuffixisseparator,
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
    *mm_newsuffix = (Blindtriesymbol)
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
        currentnodeisleaf = blindtrie_isleaf(blindtrie,NODEPTR) ? true : false;\
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
  Blindtrienodeptr currentnode, siblval, lcpnode = ROOTIDX;
  unsigned long nextfree = 0, equalsrangewidth = 0, lcpnodedepth,
                bucketleftidxplussubbucketleft;

  blindtrie->stack.nextfreeBlindtrienodeptr = 0;
  GT_STOREINARRAY (&blindtrie->stack, Blindtrienodeptr, 128, ROOTIDX);
  BLINDTRIE_SETCURRENTNODE(blindtrie_firstchild_get(blindtrie,ROOTIDX));
  gt_assert(maxdepth == 0 || dc_processunsortedrange != NULL);
  bucketleftidxplussubbucketleft
    = gt_suffixsortspace_bucketleftidx_get(blindtrie->sssp) + subbucketleft;
  for (;;)
  {
    lcpnodedepth = blindtrie_getdepth(blindtrie,lcpnode);
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
                          blindtrie_nodestartpos_get(blindtrie,currentnode));
      nextfree++;
      siblval = blindtrie_rightsibling_get(blindtrie,currentnode);
      if (siblval == BLINDTRIE_REFNULL)
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
        siblval = blindtrie_rightsibling_get(blindtrie,
                       blindtrie->stack.spaceBlindtrienodeptr[
                       blindtrie->stack.nextfreeBlindtrienodeptr]);
        if (siblval != BLINDTRIE_REFNULL)
        {
          BLINDTRIE_SETCURRENTNODE(siblval);   /* current comes from brother */
          lcpnode = blindtrie->stack.spaceBlindtrienodeptr[
                             blindtrie->stack.nextfreeBlindtrienodeptr - 1];
          readyforpop = false;
        }
      } else
      {
        GT_STOREINARRAY (&blindtrie->stack, Blindtrienodeptr, 128, currentnode);
        BLINDTRIE_SETCURRENTNODE(blindtrie_firstchild_get(blindtrie,
                                                          currentnode));
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
                            unsigned long maxnumofsuffixes,
                            unsigned int nodenumberincrement,
                            const GtEncseq *encseq,
                            bool cmpcharbychar,
                            GtEncseqReader *esr1,
                            GtEncseqReader *esr2,
                            GtReadmode readmode)
{
  Blindtrie *blindtrie;

  blindtrie = gt_malloc(sizeof (*blindtrie));
  /*printf("sizeof(Blindtrienode)=%lu\n",
           (unsigned long) sizeof (Blindtrienode));*/
  if (nodenumberincrement == 0)
  {
    gt_assert(maxnumofsuffixes >= 2UL);
    blindtrie->allocatedBlindtrienode
      = (Blindtrienodeptr) GT_MULT2(maxnumofsuffixes + 1) + 1;
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
            (unsigned long) (sizeof (Blindtrie) +
                             blindtrie->allocatedBlindtrienode *
                             sizeof (Blindtrienode)));
  */
  GT_INITARRAY (&blindtrie->overflowsuffixes, GtUlong);
  GT_INITARRAY (&blindtrie->stack, Blindtrienodeptr);
  blindtrie->nextfreeBlindtrienode = 0;
  blindtrie->encseq = encseq;
  blindtrie->has_twobitencoding_stoppos_support
    = gt_has_twobitencoding_stoppos_support(encseq);
  blindtrie->readmode = readmode;
  blindtrie->esr1 = esr1;
  blindtrie->esr2 = esr2;
  blindtrie->totallength = gt_encseq_total_length(encseq);
  blindtrie->cmpcharbychar = cmpcharbychar;
  blindtrie->sssp = suffixsortspace;
  return blindtrie;
}

void gt_blindtrie_reset(Blindtrie *blindtrie)
{
  blindtrie->nextfreeBlindtrienode = 0;
  blindtrie->stack.nextfreeBlindtrienodeptr = 0;
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

#define NODENUM(PTR) PTR

static void gt_blindtrie_showleaf(const Blindtrie *blindtrie,unsigned int level,
                                  Blindtrienodeptr current)
{
  printf("%*.*s",(int) (6 * level),(int) (6 * level)," ");
  gt_assert(current != BLINDTRIE_REFNULL);
  printf("Leaf(%u,firstchar=%u,startpos=%lu,rightsibling=%u)\n",
         NODENUM(current),
         (unsigned int) blindtrie_firstchar_get(blindtrie,current),
         blindtrie_nodestartpos_get(blindtrie,current),
         NODENUM(blindtrie_rightsibling_get(blindtrie,current)));
}

static void gt_blindtrie_showintern(const Blindtrie *blindtrie,
                                    unsigned int level,
                                    Blindtrienodeptr current)
{
  printf("%*.*s",(int) (6 * level),(int) (6 * level)," ");
  gt_assert(current != BLINDTRIE_REFNULL);
  printf("Intern(%u,firstchar=%lu,depth=%lu"
         ",firstchild=%u,rightsibling=%u)\n",
          NODENUM(current),
          blindtrie_firstchar_get(blindtrie,current),
          blindtrie_getdepth(blindtrie,current),
          NODENUM(blindtrie_firstchild_get(blindtrie,current)),
          NODENUM(blindtrie_rightsibling_get(blindtrie,current)));
}

static void gt_blindtrie_showrecursive(const Blindtrie *blindtrie,
                                       unsigned int level,
                                       Blindtrienodeptr node)
{
  Blindtrienodeptr current;

  for (current = blindtrie_firstchild_get(blindtrie,node);
       current != BLINDTRIE_REFNULL;
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

static void gt_blindtrie_show(const Blindtrie *blindtrie)
{
  gt_blindtrie_showrecursive(blindtrie,0,ROOTIDX);
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

static Blindtrienodeptr blindtrie_findsplitnode(const Blindtrie *blindtrie,
                                                unsigned long lcp)
{
  Blindtrienodeptr currentnode;
  unsigned long stackidx;

  currentnode = ROOTIDX;
  for (stackidx=0;stackidx<blindtrie->stack.nextfreeBlindtrienodeptr;
       stackidx++)
  {
    currentnode = blindtrie->stack.spaceBlindtrienodeptr[stackidx];
    if (blindtrie_isleaf(blindtrie,currentnode) ||
        blindtrie_getdepth(blindtrie,currentnode) >= lcp)
    {
      break;
    }
  }
  return currentnode;
}

static void gt_blindtrie_insertsuffix(Blindtrie *blindtrie,
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
    blindtrie_makeroot(blindtrie,currentstartpos);
  } else
  {
    if (blindtrie_isleftofboundary(blindtrie,currentstartpos,0))
    {
      unsigned long lcp, currenttwobitencodingstoppos;
      Blindtrienodeptr leafinsubtrie, splitnode;
      Blindtriesymbol mm_oldsuffix, mm_newsuffix;

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

  gt_blindtrie_reset(blindtrie);
  for (idx=0; idx < numberofsuffixes; idx++)
  {
    currentstartpos = gt_suffixsortspace_get(blindtrie->sssp,subbucketleft,idx);
    gt_blindtrie_insertsuffix(blindtrie,
                              offset,
                              maxdepth,
                              currentstartpos + offset);
#ifdef SKDEBUG
    blindtrie_showstate(blindtrie,
                        subbucketleft,
                        idx+1,
                        offset);
#endif
  }
  return gt_blindtrie2sorting(blindtrie,
                              subbucketleft,
                              lcpsubtab,
                              offset,
                              maxdepth,
                              voiddcov,
                              dc_processunsortedrange);
}

bool gt_blindtrie_retrieve(Blindtrie *blindtrie,
                           unsigned long currentstartpos)
{
  if (blindtrie->nextfreeBlindtrienode == 0)
  {
    blindtrie->maxdepthminusoffset = 0;
    blindtrie_makeroot(blindtrie,currentstartpos);
    return false;
  } else
  {
    Blindtrienodeptr leafinsubtrie, splitnode;
    Blindtriesymbol mm_oldsuffix, mm_newsuffix;
    bool mm_oldsuffixisseparator;
    unsigned long currenttwobitencodingstoppos, lcp;

    currenttwobitencodingstoppos
      = blindtrie_currenttwobitencodingstoppos_get(blindtrie,currentstartpos);
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
