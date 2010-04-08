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

#include "core/arraydef.h"
#include "core/divmodmul.h"
#include "core/minmax.h"

#include "spacedef.h"
#include "core/encodedsequence.h"
#include "lcpoverflow.h"
#include "bltrie-ssort.h"

#undef SKDEBUG
#ifdef SKDEBUG

#define NODENUM(PTR)\
        ((PTR) == NULL\
           ? 99UL\
           : (unsigned long) ((PTR) - blindtrie->spaceBlindtrienode))
#endif

#define ISLEAF(NODE)      ((NODE) != blindtrie->root && (NODE)->depth == 0)
#define ISNOTLEAF(NODE)   ((NODE) == blindtrie->root || (NODE)->depth > 0)
#define SETLEAF(NODE,VAL) /* Nothing */

typedef struct Blindtrienode
{
  unsigned long depth;
  struct Blindtrienode *rightsibling;
  union
  {
    struct Blindtrienode *firstchild;
    unsigned long nodestartpos;
  } either;
  GtUchar firstchar;
} Blindtrienode;

typedef Blindtrienode * Nodeptr;

GT_DECLAREARRAYSTRUCT(Nodeptr);

struct Blindtrie
{
  const GtEncodedsequence *encseq;
  GtEncodedsequenceScanstate *esr1, *esr2;
  GtReadmode readmode;
  unsigned long totallength,
         offset,
         maxdepth,
         maxdepthminusoffset;
  Nodeptr root;
  bool cmpcharbychar;
  unsigned long allocatedBlindtrienode,
                nextfreeBlindtrienode;
  Blindtrienode *spaceBlindtrienode;
  GtArrayNodeptr stack;
};

static bool isleftofboundary(unsigned long currentstartpos,unsigned long add,
                             const Blindtrie *blindtrie)
{
  unsigned long endpos;

  gt_assert(currentstartpos >= blindtrie->offset);
  if (blindtrie->maxdepth == 0)
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

static Nodeptr newBlindtrienode(Blindtrie *blindtrie)
{
  gt_assert(blindtrie->nextfreeBlindtrienode <
            blindtrie->allocatedBlindtrienode);
  return blindtrie->spaceBlindtrienode + blindtrie->nextfreeBlindtrienode++;
}

static Blindtrienode *makenewleaf(Blindtrie *blindtrie,
                                  unsigned long currentstartpos,
                                  GtUchar firstchar)
{
  Blindtrienode *newleaf;

  newleaf = newBlindtrienode(blindtrie);
  newleaf->either.nodestartpos = currentstartpos;
  newleaf->depth = 0;
  SETLEAF(newleaf,true);
  newleaf->firstchar = firstchar;
  newleaf->rightsibling = NULL;
  return newleaf;
}

static Nodeptr makeroot(Blindtrie *blindtrie,unsigned long currentstartpos)
{
  Blindtrienode *root;
  GtUchar firstchar;

  root = newBlindtrienode(blindtrie);
  root->depth = 0;
  root->firstchar = 0; /* undefined */
  root->rightsibling = NULL;
  SETLEAF(root,false);
  if (isleftofboundary(currentstartpos,0,blindtrie))
  {
    /* Random access */
    firstchar = gt_encodedsequence_getencodedchar(blindtrie->encseq,
                                                  currentstartpos,
                                                  blindtrie->readmode);
    if (firstchar == (GtUchar) WILDCARD)
    {
      firstchar = (GtUchar) SEPARATOR;
    }
  } else
  {
    firstchar = (GtUchar) SEPARATOR;
  }
  root->either.firstchild = makenewleaf(blindtrie,currentstartpos,firstchar);
  return root;
}

static inline Nodeptr extractleafnode(const Blindtrie *blindtrie,Nodeptr head)
{
  gt_assert(ISNOTLEAF(head));
  do
  {
    head = head->either.firstchild;
  } while (ISNOTLEAF(head));
  return head;
}

static inline int comparecharacters(GtUchar oldchar,GtUchar newchar)
{
  return (oldchar > newchar)
           ? 1
           : ((oldchar < newchar || ISSPECIAL(oldchar))
                  ? -1
                  : 0);
}

static Nodeptr findsucc(Nodeptr node,GtUchar newchar)
{
  int retval;

  for (;;)
  {
    retval = comparecharacters(node->firstchar,newchar);
    if (retval == 0)
    {              /* found branch corresponding to newchar */
      return node;
    }
    if (retval == 1)
    {               /* found branch which is already greater than newchar */
      return NULL;
    }
    node = node->rightsibling;
    if (node == NULL) /* no other branches: mismatch */
    {
      return NULL;
    }
  }
}

static Nodeptr findcompanion(Blindtrie *blindtrie,unsigned long currentstartpos)
{
  GtUchar newchar;
  Nodeptr head, succ;

  blindtrie->stack.nextfreeNodeptr = 0;
  head = blindtrie->root;
  while (ISNOTLEAF(head))
  {
    GT_STOREINARRAY (&blindtrie->stack, Nodeptr, 128, head);
    if (isleftofboundary(currentstartpos,head->depth,blindtrie))
    {
      /* Random access */
      newchar = gt_encodedsequence_getencodedchar(blindtrie->encseq,
                                                  currentstartpos + head->depth,
                                                  blindtrie->readmode);
      if (newchar == (GtUchar) WILDCARD)
      {
        newchar = (GtUchar) SEPARATOR;
      }
    } else
    {
      newchar = (GtUchar) SEPARATOR;
    }
    if (ISSPECIAL(newchar))
    {
      return extractleafnode(blindtrie,head);
    }
    succ = findsucc(head->either.firstchild,newchar);
    if (succ == NULL)
    {
      return extractleafnode(blindtrie,head);
    }
    head = succ;
  }
  GT_STOREINARRAY (&blindtrie->stack, Nodeptr, 128, head);
  return head;
}

static void insertsuffixintoblindtrie(Blindtrie *blindtrie,
                                      Nodeptr oldnode,
                                      GtUchar mm_oldsuffix,
                                      unsigned long lcp,
                                      GtUchar mm_newsuffix,
                                      unsigned long currentstartpos)
{
  Nodeptr newleaf, newnode, previous, current;

  gt_assert(ISSPECIAL(mm_oldsuffix) || ISSPECIAL(mm_newsuffix) ||
            mm_oldsuffix != mm_newsuffix || ISLEAF(oldnode) ||
            oldnode->depth == lcp);

  /* insert a new node before node oldnode if necessary */
  if (oldnode->depth != lcp)
  {
    newnode = newBlindtrienode(blindtrie);
    newnode->firstchar = mm_oldsuffix;
    newnode->depth = oldnode->depth; /* newnode inherits depth+children */
    SETLEAF(newnode,oldnode->isleaf);
    newnode->either = oldnode->either;
    newnode->rightsibling = NULL;
    oldnode->depth = lcp;
    SETLEAF(oldnode,false);
    oldnode->either.firstchild = newnode; /* oldnode has newnode as only child*/
  }
  gt_assert(oldnode->depth == lcp);

  /* search S[lcp] among the offsprings */
  newleaf = newBlindtrienode(blindtrie);
  newleaf->depth = 0;
  SETLEAF(newleaf,true);
  newleaf->firstchar = mm_newsuffix;
  newleaf->either.nodestartpos = currentstartpos;
  previous = NULL;
  current = oldnode->either.firstchild;
  while (current != NULL &&
         comparecharacters(current->firstchar,mm_newsuffix) < 0)
  {
    previous = current;
    current = current->rightsibling;
  }
  /* insert new leaf with current suffix */
  if (previous != NULL)
  {
    previous->rightsibling = newleaf;
  } else
  {
    oldnode->either.firstchild = newleaf;
  }
  newleaf->rightsibling = current;
}

static unsigned long cmpcharbychargetlcp(GtUchar *mm_oldsuffix,
                                  GtUchar *mm_newsuffix,
                                  const Blindtrie *blindtrie,
                                  unsigned long leafpos,
                                  unsigned long currentstartpos)
{
  unsigned long lcp;
  GtUchar cc1, cc2;

  gt_encodedsequence_scanstate_init(blindtrie->esr1,blindtrie->encseq,
                               blindtrie->readmode,leafpos);
  gt_encodedsequence_scanstate_init(blindtrie->esr2,blindtrie->encseq,
                               blindtrie->readmode,currentstartpos);
  for (lcp = 0; /* Nothing */; lcp++)
  {
    if (isleftofboundary(leafpos,lcp,blindtrie))
    {
      cc1 = gt_encodedsequence_sequentialgetencodedchar(blindtrie->encseq,
                                                        blindtrie->esr1,
                                                        leafpos+lcp,
                                                        blindtrie->readmode);
      if (cc1 == (GtUchar) WILDCARD)
      {
        cc1 = (GtUchar) SEPARATOR;
      }
    } else
    {
      cc1 = (GtUchar) SEPARATOR;
    }
    if (isleftofboundary(currentstartpos,lcp,blindtrie))
    {
      cc2 = gt_encodedsequence_sequentialgetencodedchar(blindtrie->encseq,
                                                        blindtrie->esr2,
                                                        currentstartpos+lcp,
                                                        blindtrie->readmode);
      if (cc2 == (GtUchar) WILDCARD)
      {
        cc2 = (GtUchar) SEPARATOR;
      }
    } else
    {
      cc2 = (GtUchar) SEPARATOR;
    }
    if (comparecharacters(cc1,cc2) != 0)
    {
      *mm_oldsuffix = cc1;
      *mm_newsuffix = cc2;
      break;
    }
  }
  gt_assert(blindtrie->maxdepth == 0 || lcp <= blindtrie->maxdepthminusoffset);
  return lcp;
}

static unsigned long fastgetlcp(GtUchar *mm_oldsuffix,
                         GtUchar *mm_newsuffix,
                         const Blindtrie *blindtrie,
                         unsigned long leafpos,
                         unsigned long currentstartpos)
{
  GtCommonunits commonunits;

  if (blindtrie->maxdepth == 0)
  {
    (void) gt_encodedsequence_compare(blindtrie->encseq,
                                      &commonunits,
                                      GT_ISDIRREVERSE(blindtrie->readmode)
                                      ? false : true,
                                      GT_ISDIRCOMPLEMENT(blindtrie->readmode)
                                      ? true : false,
                                      blindtrie->esr1,
                                      blindtrie->esr2,
                                      leafpos,
                                      currentstartpos,
                                      0);
  } else
  {
    (void) gt_encodedsequence_compare_maxdepth(blindtrie->encseq,
                                               &commonunits,
                                               GT_ISDIRREVERSE(
                                                            blindtrie->readmode)
                                               ? false : true,
                                               GT_ISDIRCOMPLEMENT(
                                                            blindtrie->readmode)
                                               ? true : false,
                                               blindtrie->esr1,
                                               blindtrie->esr2,
                                               leafpos,
                                               currentstartpos,
                                               0,
                                               blindtrie->maxdepthminusoffset);
  }
  if (isleftofboundary(leafpos,commonunits.finaldepth,blindtrie) &&
      !commonunits.leftspecial)
  {
#ifdef OLDVERSION
    /*
    *mm_oldsuffix = gt_encodedsequence_getencodedchar(blindtrie->encseq,
                                   leafpos + commonunits.finaldepth,
                                   blindtrie->readmode);
    if (ISSPECIAL(*mm_oldsuffix))
    {
      gt_assert(commonunits.leftspecial);
    } else
    {
      GtUchar tmp = gt_encodedsequence_extractencodedchar(blindtrie->encseq,
                                                         leafpos +
                                                         commonunits.finaldepth,
                                                         blindtrie->readmode);
      gt_assert(tmp == *mm_oldsuffix);
    }
    if (*mm_oldsuffix == (GtUchar) WILDCARD)
    {
      *mm_oldsuffix = (GtUchar) SEPARATOR;
    }
    */
#endif
    *mm_oldsuffix = gt_encodedsequence_extractencodedchar(blindtrie->encseq,
                                                         leafpos +
                                                         commonunits.finaldepth,
                                                         blindtrie->readmode);
  } else
  {
    *mm_oldsuffix = (GtUchar) SEPARATOR;
  }
  if (isleftofboundary(currentstartpos,commonunits.finaldepth,blindtrie) &&
      !commonunits.rightspecial)
  {
#ifdef OLDVERSION
    /* Random access */
    *mm_newsuffix = gt_encodedsequence_getencodedchar(blindtrie->encseq,
                                       currentstartpos + commonunits.finaldepth,
                                       blindtrie->readmode);
    if (ISSPECIAL(*mm_newsuffix))
    {
      gt_assert(commonunits.rightspecial);
    } else
    {
      /* Random access */
      GtUchar tmp = gt_encodedsequence_extractencodedchar(blindtrie->encseq,
                                       currentstartpos + commonunits.finaldepth,
                                       blindtrie->readmode);
      gt_assert(tmp == *mm_newsuffix);
    }
    if (*mm_newsuffix == (GtUchar) WILDCARD)
    {
      *mm_newsuffix = (GtUchar) SEPARATOR;
    }
#endif
    *mm_newsuffix = gt_encodedsequence_extractencodedchar(blindtrie->encseq,
                                       currentstartpos +
                                       commonunits.finaldepth,
                                       blindtrie->readmode);
  } else
  {
    *mm_newsuffix = (GtUchar) SEPARATOR;
  }
  return commonunits.finaldepth;
}

#define SETCURRENT(VAL)\
        currentnodeisleaf = ISLEAF(VAL) ? true : false;\
        currentnode = VAL

static unsigned long enumeratetrieleaves (unsigned long *suffixtable,
                                          unsigned long *lcpsubtab,
                                          unsigned long *numoflargelcpvalues,
                                          Blindtrie *blindtrie,
                                          void *voiddcov,
                                          void (*dc_processunsortedrange)(
                                                   void *,
                                                   unsigned long *,
                                                   unsigned long *,
                                                   unsigned long))
{
  bool readyforpop = false, currentnodeisleaf;
  Nodeptr currentnode, siblval, lcpnode = blindtrie->root;
  unsigned long nextfree = 0, equalsrangewidth = 0;

  blindtrie->stack.nextfreeNodeptr = 0;
  GT_STOREINARRAY (&blindtrie->stack, Nodeptr, 128, blindtrie->root);
  SETCURRENT(blindtrie->root->either.firstchild);
  gt_assert(blindtrie->maxdepth == 0 || dc_processunsortedrange != NULL);
  for (;;)
  {
    if (currentnodeisleaf)
    {
      if (nextfree > 0)
      {
        if (lcpsubtab != NULL)
        {
          lcpsubtab[nextfree] = lcpnode->depth + blindtrie->offset;
          if (lcpnode->depth + blindtrie->offset >= (unsigned long) LCPOVERFLOW)
          {
            (*numoflargelcpvalues)++;
          }
        }
        if (blindtrie->maxdepth > 0)
        {
          if (lcpnode->depth + blindtrie->offset == blindtrie->maxdepth)
          {
            equalsrangewidth++;
          } else
          {
#ifndef NDEBUG
            if (lcpnode->depth + blindtrie->offset >= blindtrie->maxdepth)
            {
              fprintf(stderr,"lcpnode.depth=%lu,offset=%lu,maxdepth=%lu\n",
                              (unsigned long) lcpnode->depth,
                              (unsigned long) blindtrie->offset,
                              (unsigned long) blindtrie->maxdepth);
              exit(EXIT_FAILURE);
            }
#endif
            gt_assert(lcpnode->depth + blindtrie->offset < blindtrie->maxdepth);
            if (equalsrangewidth > 0)
            {
              dc_processunsortedrange(voiddcov,
                                      suffixtable + nextfree - 1
                                                  - equalsrangewidth,
                                      suffixtable + nextfree - 1,
                                      blindtrie->maxdepth);
              equalsrangewidth = 0;
            }
          }
        }
      }
      suffixtable[nextfree++]
        = currentnode->either.nodestartpos - blindtrie->offset;
      siblval = currentnode->rightsibling;
      if (siblval == NULL)
      {
        readyforpop = true;
        currentnodeisleaf = false; /* STATE 1 */
      } else
      {
        SETCURRENT (siblval);  /* current comes from brother */
        lcpnode
          = blindtrie->stack.spaceNodeptr[blindtrie->stack.nextfreeNodeptr-1];
      }
    } else
    {
      if (readyforpop)
      {
        if (blindtrie->stack.nextfreeNodeptr == 1UL)
        {
          break;
        }
        blindtrie->stack.nextfreeNodeptr--;
        siblval = blindtrie->stack.spaceNodeptr[
                           blindtrie->stack.nextfreeNodeptr]->rightsibling;
        if (siblval != NULL)
        {
          SETCURRENT (siblval);        /* current comes from brother */
          lcpnode = blindtrie->stack.spaceNodeptr[
                             blindtrie->stack.nextfreeNodeptr - 1];
          readyforpop = false;
        }
      } else
      {
        GT_STOREINARRAY (&blindtrie->stack, Nodeptr, 128, currentnode);
        SETCURRENT (currentnode->either.firstchild);
      }
    }
  }
  if (nextfree > 0 && equalsrangewidth > 0)
  {
    dc_processunsortedrange(voiddcov,
                            suffixtable + nextfree - 1 - equalsrangewidth,
                            suffixtable + nextfree - 1,
                            blindtrie->maxdepth);
    equalsrangewidth = 0;
  }
  return nextfree;
}

Blindtrie *gt_blindtrie_new(unsigned long numofsuffixes,
                         const GtEncodedsequence *encseq,
                         bool cmpcharbychar,
                         GtEncodedsequenceScanstate *esr1,
                         GtEncodedsequenceScanstate *esr2,
                         GtReadmode readmode)
{
  Blindtrie *blindtrie;

  ALLOCASSIGNSPACE(blindtrie,NULL,Blindtrie,1);
  blindtrie->allocatedBlindtrienode = GT_MULT2(numofsuffixes + 1) + 1;
  ALLOCASSIGNSPACE(blindtrie->spaceBlindtrienode,NULL,Blindtrienode,
                   blindtrie->allocatedBlindtrienode);
  /*
  printf("# sizeof (blindtrie)=%lu\n",
            (unsigned long) (sizeof (Blindtrie) +
                             blindtrie->allocatedBlindtrienode *
                             sizeof (Blindtrienode)));
  */
  blindtrie->nextfreeBlindtrienode = 0;
  blindtrie->encseq = encseq;
  blindtrie->readmode = readmode;
  blindtrie->root = NULL;
  blindtrie->esr1 = esr1;
  blindtrie->esr2 = esr2;
  blindtrie->totallength = gt_encodedsequence_total_length(encseq);
  blindtrie->cmpcharbychar = cmpcharbychar;
  GT_INITARRAY (&blindtrie->stack, Nodeptr);
  return blindtrie;
}

void gt_blindtrie_delete(Blindtrie **blindtrie)
{
  FREESPACE((*blindtrie)->spaceBlindtrienode);
  GT_FREEARRAY(&(*blindtrie)->stack, Nodeptr);
  FREESPACE(*blindtrie);
}

#ifdef SKDEBUG
static void checkcurrentblindtrie(Blindtrie *blindtrie)
{
  unsigned long suffixtable[6];
  unsigned long idx, numofsuffixes;
  unsigned long maxcommon;
  int retval;

  numofsuffixes = enumeratetrieleaves (&suffixtable[0], NULL, NULL, blindtrie);
  for (idx=1UL; idx < numofsuffixes; idx++)
  {
    maxcommon = 0;
    retval = gt_encodedsequence_comparetwostringsgeneric(blindtrie->encseq,
                                      GT_ISDIRREVERSE(blindtrie->readmode)
                                        ? false : true,
                                      GT_ISDIRCOMPLEMENT(blindtrie->readmode)
                                        ? true : false,
                                      &maxcommon,
                                      suffixtable[idx-1],
                                      suffixtable[idx],
                                      0);
    if (retval >= 0)
    {
      fprintf(stderr,"retval = %d, maxcommon = %u for idx = %lu\n",
              retval,maxcommon,idx);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

static void showleaf(const Blindtrie *blindtrie,unsigned int level,
                     Nodeptr current)
{
  printf("%*.*s",(int) (6 * level),(int) (6 * level)," ");
  gt_assert(current != NULL);
  printf("Leaf(add=%lu,firstchar=%u,startpos=" FormatSeqpos
         ",rightsibling=%lu)\n",
         NODENUM(current),
         (unsigned int) current->firstchar,
         PRINTSeqposcast(current->either.startpos),
         NODENUM(current->rightsibling));
}

static void showintern(const Blindtrie *blindtrie,unsigned int level,
                       Nodeptr current)
{
  printf("%*.*s",(int) (6 * level),(int) (6 * level)," ");
  gt_assert(current != NULL);
  printf("Intern(add=%lu,firstchar=%u,depth=" FormatSeqpos
         ",firstchild=%lu,rightsibling=%lu)\n",
          NODENUM(current),
          (unsigned int) current->firstchar,
          PRINTSeqposcast(current->depth),
          NODENUM(current->either.firstchild),
          NODENUM(current->rightsibling));
}

static void showblindtrie2(const Blindtrie *blindtrie,
                           unsigned int level,
                           Nodeptr node)
{
  Nodeptr current;

  for (current = node->either.firstchild;
       current != NULL;
       current = current->rightsibling)
  {
    if (ISLEAF(current))
    {
      showleaf(blindtrie,level,current);
    } else
    {
      showintern(blindtrie,level,current);
      showblindtrie2(blindtrie,level+1,current);
    }
  }
}

static void showblindtrie(const Blindtrie *blindtrie)
{
  showblindtrie2(blindtrie,0,blindtrie->root);
}
#endif

static int suffixcompare(const void *a, const void *b)
{
  gt_assert(*((unsigned long *) a) != *((unsigned long *) b));
  if (*((unsigned long *) a) < *((unsigned long *) b))
  {
    return -1;
  }
  return 1;
}

#ifndef NDEBUG

static void checksorting(bool ascending,
                         const unsigned long *suffixtable,
                         unsigned long numberofsuffixes)
{
  unsigned long idx;

  gt_assert(numberofsuffixes > 1UL);
  for (idx = 0; idx < numberofsuffixes - 1; idx++)
  {
    if ((ascending && suffixtable[idx] >= suffixtable[idx+1]) ||
        (!ascending && suffixtable[idx] <= suffixtable[idx+1]))
    {
      fprintf(stderr,"not %s: ",ascending ? "ascending" : "descending");
      fprintf(stderr,"suffixtable[%lu]=%lu vs %lu=suffixtable[%lu]\n",
                      idx,(unsigned long) suffixtable[idx],
                          (unsigned long) suffixtable[idx+1],idx+1);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

#endif

static void inplace_reverseSeqpos(unsigned long *tab,unsigned long len)
{
  unsigned long tmp, *frontptr, *backptr;

  for (frontptr = tab, backptr = tab + len - 1;
       frontptr < backptr; frontptr++, backptr--)
  {
    tmp = *frontptr;
    *frontptr = *backptr;
    *backptr = tmp;
  }
}

unsigned long gt_blindtrie_suffixsort(Blindtrie *blindtrie,
                            unsigned long *suffixtable,
                            unsigned long *lcpsubtab,
                            unsigned long numberofsuffixes,
                            unsigned long offset,
                            unsigned long maxdepth,
                            Ordertype ordertype,
                            void *voiddcov,
                            void (*dc_processunsortedrange)(void *,
                                                            unsigned long *,
                                                            unsigned long *,
                                                            unsigned long))
{
  unsigned long idx, stackidx;
  Nodeptr leafinsubtree, currentnode;
  unsigned long lcp, numoflargelcpvalues = 0;
  GtUchar mm_oldsuffix, mm_newsuffix;

  if (ordertype == Noorder)
  {
    qsort(suffixtable,
          (size_t) numberofsuffixes,
          sizeof (unsigned long),
          suffixcompare);
  } else
  {
    if (ordertype == Descending)
    {
#ifndef NDEBUG
      checksorting(false,suffixtable,numberofsuffixes);
#endif
      inplace_reverseSeqpos(suffixtable,numberofsuffixes);
    } else
    {
#ifndef NDEBUG
      checksorting(true,suffixtable,numberofsuffixes);
#endif
    }
  }
  gt_assert(maxdepth == 0 || maxdepth > offset);
  blindtrie->maxdepth = maxdepth;
  blindtrie->offset = offset;
  if (maxdepth > 0)
  {
    blindtrie->maxdepthminusoffset = maxdepth - offset;
  } else
  {
    blindtrie->maxdepthminusoffset = 0;
  }
  blindtrie->nextfreeBlindtrienode = 0;
  blindtrie->root = makeroot(blindtrie,suffixtable[0] + offset);
#ifdef SKDEBUG
  printf("insert suffixes at offset " FormatSeqpos ":\n",
          PRINTSeqposcast(offset));
  for (i=0; i < numberofsuffixes; i++)
  {
    printf(FormatSeqpos " ",PRINTSeqposcast(suffixtable[i] + offset));
  }
  printf("\nstep 0\n");
  showblindtrie(blindtrie);
#endif
  for (idx=1UL; idx < numberofsuffixes; idx++)
  {
    if (isleftofboundary(suffixtable[idx] + offset,0,blindtrie))
    {
      leafinsubtree = findcompanion(blindtrie,suffixtable[idx] + offset);
      gt_assert(ISLEAF(leafinsubtree));
      lcp = (blindtrie->cmpcharbychar ? cmpcharbychargetlcp : fastgetlcp)
                               (&mm_oldsuffix,
                                &mm_newsuffix,
                                blindtrie,
                                leafinsubtree->either.nodestartpos,
                                suffixtable[idx] + offset);
      currentnode = blindtrie->root;
      for (stackidx=0;stackidx<blindtrie->stack.nextfreeNodeptr;stackidx++)
      {
        currentnode = blindtrie->stack.spaceNodeptr[stackidx];
        if (ISLEAF(currentnode) || currentnode->depth >= lcp)
        {
          break;
        }
      }
      insertsuffixintoblindtrie(blindtrie,
                                currentnode,
                                mm_oldsuffix,
                                lcp,
                                mm_newsuffix,
                                suffixtable[idx] + offset);
#ifdef SKDEBUG
      printf("step %lu\n",i);
      showblindtrie(blindtrie);
      checkcurrentblindtrie(blindtrie);
#endif
    } else
    {
      break;
    }
  }
  (void) enumeratetrieleaves (suffixtable, lcpsubtab, &numoflargelcpvalues,
                              blindtrie,voiddcov,dc_processunsortedrange);
  if (lcpsubtab != NULL)
  {
    if (idx < numberofsuffixes && offset >= (unsigned long) LCPOVERFLOW)
    {
      numoflargelcpvalues += numberofsuffixes - idx;
    }
    while (idx < numberofsuffixes)
    {
      lcpsubtab[idx++] = offset;
    }
  }
  return numoflargelcpvalues;
}
