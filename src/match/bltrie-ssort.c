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
#include "divmodmul.h"
#include "seqpos-def.h"
#include "spacedef.h"
#include "encseq-def.h"
#include "bltrie-ssort.h"
#include "lcpoverflow.h"

#undef SKDEBUG
#ifdef SKDEBUG
#include "sfx-cmpsuf.pr"

#define NODENUM(PTR)\
        ((PTR) == NULL ? 99UL\
                       : (unsigned long) ((PTR) - trierep->spaceBlindtrienode))
#endif

#define ISLEAF(NODE)      ((NODE) != trierep->root && (NODE)->depth == 0)
#define ISNOTLEAF(NODE)   ((NODE) == trierep->root || (NODE)->depth > 0)
#define SETLEAF(NODE,VAL) /* Nothing */

typedef struct Blindtrienode
{
  Seqpos depth;
  struct Blindtrienode *rightsibling;
  union
  {
    struct Blindtrienode *firstchild;
    Seqpos nodestartpos;
  } either;
  Uchar firstchar;
} Blindtrienode;

typedef Blindtrienode * Nodeptr;

DECLAREARRAYSTRUCT(Nodeptr);

struct Blindtrierep
{
  const Encodedsequence *encseq;
  Encodedsequencescanstate *esr1, *esr2;
  Readmode readmode;
  Seqpos totallength,
         offset;
  Nodeptr root;
  bool cmpcharbychar;
  unsigned long allocatedBlindtrienode,
                nextfreeBlindtrienode;
  Blindtrienode *spaceBlindtrienode;
  ArrayNodeptr stack;
};

static bool isleftofboundary(Seqpos currentstartpos,
                             const Blindtrierep *trierep)
{
  return (currentstartpos < trierep->totallength) ? true : false;
}

static Nodeptr newBlindtrienode(Blindtrierep *trierep)
{
  gt_assert(trierep->nextfreeBlindtrienode < trierep->allocatedBlindtrienode);
  return trierep->spaceBlindtrienode + trierep->nextfreeBlindtrienode++;
}

static Blindtrienode *makenewleaf(Blindtrierep *trierep,
                                  Seqpos currentstartpos,
                                  Uchar firstchar)
{
  Blindtrienode *newleaf;

  newleaf = newBlindtrienode(trierep);
  newleaf->either.nodestartpos = currentstartpos;
  newleaf->depth = 0;
  SETLEAF(newleaf,true);
  newleaf->firstchar = firstchar;
  newleaf->rightsibling = NULL;
  return newleaf;
}

static Nodeptr makeroot(Blindtrierep *trierep,Seqpos currentstartpos)
{
  Blindtrienode *root;
  Uchar firstchar;

  root = newBlindtrienode(trierep);
  root->depth = 0;
  root->firstchar = 0; /* undefined */
  root->rightsibling = NULL;
  SETLEAF(root,false);
  if (isleftofboundary(currentstartpos,trierep))
  {
    firstchar = getencodedchar(trierep->encseq, /* Random access */
                               currentstartpos,
                               trierep->readmode);
    if (firstchar == (Uchar) WILDCARD)
    {
      firstchar = (Uchar) SEPARATOR;
    }
  } else
  {
    firstchar = (Uchar) SEPARATOR;
  }
  root->either.firstchild = makenewleaf(trierep,currentstartpos,firstchar);
  return root;
}

static Nodeptr extractleafnode(const Blindtrierep *trierep,Nodeptr head)
{
  gt_assert(ISNOTLEAF(head));
  do
  {
    head = head->either.firstchild;
  } while (ISNOTLEAF(head));
  return head;
}

static int comparecharacters(Uchar oldchar,Uchar newchar)
{
  if (oldchar > newchar)
  {
    return 1;
  }
  if (oldchar < newchar || ISSPECIAL(oldchar))
  {
    return -1;
  }
  return 0;
}

static Nodeptr findsucc(Nodeptr node,Uchar newchar)
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

static Nodeptr findcompanion(Blindtrierep *trierep,Seqpos currentstartpos)
{
  Uchar newchar;
  Nodeptr head, succ;

  trierep->stack.nextfreeNodeptr = 0;
  head = trierep->root;
  while (ISNOTLEAF(head))
  {
    STOREINARRAY (&trierep->stack, Nodeptr, 128, head);
    if (isleftofboundary(currentstartpos+head->depth,trierep))
    {
      newchar = getencodedchar(trierep->encseq, /* Random access */
                               currentstartpos + head->depth,
                               trierep->readmode);
      if (newchar == (Uchar) WILDCARD)
      {
        newchar = (Uchar) SEPARATOR;
      }
    } else
    {
      newchar = (Uchar) SEPARATOR;
    }
    if (ISSPECIAL(newchar))
    {
      return extractleafnode(trierep,head);
    }
    succ = findsucc(head->either.firstchild,newchar);
    if (succ == NULL)
    {
      return extractleafnode(trierep,head);
    }
    head = succ;
  }
  STOREINARRAY (&trierep->stack, Nodeptr, 128, head);
  return head;
}

static void insertsuffixintoblindtrie(Blindtrierep *trierep,
                                      Nodeptr oldnode,
                                      Uchar mm_oldsuffix,
                                      Seqpos lcp,
                                      Uchar mm_newsuffix,
                                      Seqpos currentstartpos)
{
  Nodeptr newleaf, newnode, previous, current;

  gt_assert(ISSPECIAL(mm_oldsuffix) || ISSPECIAL(mm_newsuffix) ||
            mm_oldsuffix != mm_newsuffix || ISLEAF(oldnode) ||
            oldnode->depth == lcp);

  /* insert a new node before node oldnode if necessary */
  if (oldnode->depth != lcp)
  {
    newnode = newBlindtrienode(trierep);
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
  newleaf = newBlindtrienode(trierep);
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

static Seqpos cmpcharbychargetlcp(Uchar *mm_oldsuffix,
                                  Uchar *mm_newsuffix,
                                  Blindtrierep *trierep,
                                  Seqpos leafpos,
                                  Seqpos currentstartpos)
{
  Seqpos idx1, idx2;
  Uchar cc1, cc2;

  initEncodedsequencescanstate(trierep->esr1,trierep->encseq,trierep->readmode,
                               leafpos);
  initEncodedsequencescanstate(trierep->esr2,trierep->encseq,trierep->readmode,
                               currentstartpos);
  for (idx1 = leafpos, idx2 = currentstartpos; /* Nothing */; idx1++, idx2++)
  {
    if (isleftofboundary(idx1,trierep))
    {
      cc1 = sequentialgetencodedchar(trierep->encseq,trierep->esr1,
                                     idx1,trierep->readmode);
      if (cc1 == (Uchar) WILDCARD)
      {
        cc1 = (Uchar) SEPARATOR;
      }
    } else
    {
      cc1 = (Uchar) SEPARATOR;
    }
    if (isleftofboundary(idx2,trierep))
    {
      cc2 = sequentialgetencodedchar(trierep->encseq,trierep->esr2,
                                     idx2,trierep->readmode);
      if (cc2 == (Uchar) WILDCARD)
      {
        cc2 = (Uchar) SEPARATOR;
      }
    } else
    {
      cc2 = (Uchar) SEPARATOR;
    }
    if (comparecharacters(cc1,cc2) != 0)
    {
      *mm_oldsuffix = cc1;
      *mm_newsuffix = cc2;
      break;
    }
  }
  return idx1 - leafpos;
}

static Seqpos fastgetlcp(Uchar *mm_oldsuffix,
                         Uchar *mm_newsuffix,
                         Blindtrierep *trierep,
                         Seqpos leafpos,
                         Seqpos currentstartpos)
{
  Seqpos lcp;

  (void) compareEncseqsequences(&lcp,
                                trierep->encseq,
                                ISDIRREVERSE(trierep->readmode) ? false : true,
                                ISDIRCOMPLEMENT(trierep->readmode) ? true
                                                                   : false,
                                trierep->esr1,
                                trierep->esr2,
                                leafpos,
                                currentstartpos,
                                0);
  if (isleftofboundary(leafpos+lcp,trierep))
  {
    *mm_oldsuffix = getencodedchar(trierep->encseq, /* Random access */
                                   leafpos + lcp,
                                   trierep->readmode);
    if (*mm_oldsuffix == (Uchar) WILDCARD)
    {
      *mm_oldsuffix = (Uchar) SEPARATOR;
    }
  } else
  {
    *mm_oldsuffix = (Uchar) SEPARATOR;
  }
  if (isleftofboundary(currentstartpos+lcp,trierep))
  {
    *mm_newsuffix = getencodedchar(trierep->encseq, /* Random access */
                                   currentstartpos + lcp,
                                   trierep->readmode);
    if (*mm_newsuffix == (Uchar) WILDCARD)
    {
      *mm_newsuffix = (Uchar) SEPARATOR;
    }
  } else
  {
    *mm_newsuffix = (Uchar) SEPARATOR;
  }
  return lcp;
}

#define SETCURRENT(VAL)\
        currentnodeisleaf = ISLEAF(VAL) ? true : false;\
        currentnode = VAL

static unsigned long enumeratetrieleaves (Seqpos *suffixtable,
                                          Seqpos *lcpsubtab,
                                          Seqpos *numoflargelcpvalues,
                                          Blindtrierep *trierep)
{
  bool readyforpop = false, currentnodeisleaf;
  Nodeptr currentnode, siblval, lcpnode = trierep->root;
  unsigned long nextfree = 0;

  trierep->stack.nextfreeNodeptr = 0;
  STOREINARRAY (&trierep->stack, Nodeptr, 128, trierep->root);
  SETCURRENT(trierep->root->either.firstchild);
  for (;;)
  {
    if (currentnodeisleaf)
    {
      if (lcpsubtab != NULL && nextfree > 0)
      {
        lcpsubtab[nextfree] = lcpnode->depth + trierep->offset;
        if (lcpnode->depth + trierep->offset >= (Seqpos) LCPOVERFLOW)
        {
          (*numoflargelcpvalues)++;
        }
      }
      suffixtable[nextfree++]
        = currentnode->either.nodestartpos - trierep->offset;
      siblval = currentnode->rightsibling;
      if (siblval == NULL)
      {
        readyforpop = true;
        currentnodeisleaf = false; /* STATE 1 */
      } else
      {
        SETCURRENT (siblval);  /* current comes from brother */
        lcpnode = trierep->stack.spaceNodeptr[trierep->stack.nextfreeNodeptr-1];
      }
    } else
    {
      if (readyforpop)
      {
        if (trierep->stack.nextfreeNodeptr == 1UL)
        {
          break;
        }
        trierep->stack.nextfreeNodeptr--;
        siblval = trierep->stack.spaceNodeptr[
                           trierep->stack.nextfreeNodeptr]->rightsibling;
        if (siblval != NULL)
        {
          SETCURRENT (siblval);        /* current comes from brother */
          lcpnode = trierep->stack.spaceNodeptr[
                             trierep->stack.nextfreeNodeptr - 1];
          readyforpop = false;
        }
      } else
      {
        STOREINARRAY (&trierep->stack, Nodeptr, 128, currentnode);
        SETCURRENT (currentnode->either.firstchild);
      }
    }
  }
  return nextfree;
}

Blindtrierep *newBlindtrierep(unsigned long numofsuffixes,
                              const Encodedsequence *encseq,
                              bool cmpcharbychar,
                              Readmode readmode)
{
  Blindtrierep *trierep;

  ALLOCASSIGNSPACE(trierep,NULL,Blindtrierep,1);
  trierep->allocatedBlindtrienode = MULT2(numofsuffixes + 1) + 1;
  ALLOCASSIGNSPACE(trierep->spaceBlindtrienode,NULL,Blindtrienode,
                   trierep->allocatedBlindtrienode);
  printf("# sizeof (blindtrie)=%lu\n",
            (unsigned long) (sizeof (Blindtrierep) +
                             trierep->allocatedBlindtrienode *
                             sizeof (Blindtrienode)));
  trierep->nextfreeBlindtrienode = 0;
  trierep->encseq = encseq;
  trierep->readmode = readmode;
  trierep->root = NULL;
  trierep->esr1 = newEncodedsequencescanstate();
  trierep->esr2 = newEncodedsequencescanstate();
  trierep->totallength = getencseqtotallength(encseq);
  trierep->cmpcharbychar = cmpcharbychar;
  INITARRAY (&trierep->stack, Nodeptr);
  return trierep;
}

void freeBlindtrierep(Blindtrierep **trierep)
{
  FREESPACE((*trierep)->spaceBlindtrienode);
  FREEARRAY(&(*trierep)->stack, Nodeptr);
  freeEncodedsequencescanstate(&((*trierep)->esr1));
  freeEncodedsequencescanstate(&((*trierep)->esr2));
  FREESPACE(*trierep);
}

#ifdef SKDEBUG
static void checkcurrentblindtrie(Blindtrierep *trierep)
{
  Seqpos suffixtable[6];
  unsigned long idx, numofsuffixes;
  Seqpos maxcommon;
  int retval;

  numofsuffixes = enumeratetrieleaves (&suffixtable[0], NULL, NULL, trierep);
  for (idx=1UL; idx < numofsuffixes; idx++)
  {
    maxcommon = 0;
    retval = comparetwostringsgeneric(trierep->encseq,
                                      ISDIRREVERSE(trierep->readmode)
                                        ? false : true,
                                      ISDIRCOMPLEMENT(trierep->readmode)
                                        ? true : false,
                                      &maxcommon,
                                      suffixtable[idx-1],
                                      suffixtable[idx],
                                      0);
    if (retval >= 0)
    {
      fprintf(stderr,"retval = %d, maxcommon = %u for idx = %lu\n",
              retval,maxcommon,idx);
      gt_assert(retval < 0);
    }
  }
}

static void showleaf(const Blindtrierep *trierep,unsigned int level,
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

static void showintern(const Blindtrierep *trierep,unsigned int level,
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

static void showblindtrie2(const Blindtrierep *trierep,
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
      showleaf(trierep,level,current);
    } else
    {
      showintern(trierep,level,current);
      showblindtrie2(trierep,level+1,current);
    }
  }
}

static void showblindtrie(const Blindtrierep *trierep)
{
  showblindtrie2(trierep,0,trierep->root);
}
#endif

static int suffixcompare(const void *a, const void *b)
{
  gt_assert(*((Seqpos *) a) != *((Seqpos *) b));
  if (*((Seqpos *) a) < *((Seqpos *) b))
  {
    return -1;
  }
  return 1;
}

#ifndef NDEBUG

static void checksorting(bool ascending,
                         const Seqpos *suffixtable,
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
      exit(EXIT_FAILURE);
    }
  }
}

#endif

static void inplace_reverseSeqpos(Seqpos *tab,unsigned long len)
{
  Seqpos tmp, *frontptr, *backptr;

  for (frontptr = tab, backptr = tab + len - 1;
       frontptr < backptr; frontptr++, backptr--)
  {
    tmp = *frontptr;
    *frontptr = *backptr;
    *backptr = tmp;
  }
}

Seqpos blindtriesuffixsort(Blindtrierep *trierep,
                           Seqpos *suffixtable,
                           Seqpos *lcpsubtab,
                           unsigned long numberofsuffixes,
                           Seqpos offset,
                           Ordertype ordertype)
{
  unsigned long idx, stackidx;
  Nodeptr leafinsubtree, currentnode;
  Seqpos lcp, numoflargelcpvalues = 0;
  Uchar mm_oldsuffix, mm_newsuffix;

  if (ordertype == Noorder)
  {
    qsort(suffixtable,(size_t) numberofsuffixes,sizeof (Seqpos), suffixcompare);
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
  trierep->offset = offset;
  trierep->nextfreeBlindtrienode = 0;
  trierep->root = makeroot(trierep,suffixtable[0] + offset);
#ifdef SKDEBUG
  printf("insert suffixes at offset " FormatSeqpos ":\n",
          PRINTSeqposcast(offset));
  for (i=0; i < numberofsuffixes; i++)
  {
    printf(FormatSeqpos " ",PRINTSeqposcast(suffixtable[i] + offset));
  }
  printf("\nstep 0\n");
  showblindtrie(trierep);
#endif
  for (idx=1UL; idx < numberofsuffixes; idx++)
  {
    if (isleftofboundary(suffixtable[idx] + offset,trierep))
    {
      leafinsubtree = findcompanion(trierep,suffixtable[idx] + offset);
      gt_assert(ISLEAF(leafinsubtree));
      lcp = (trierep->cmpcharbychar ? cmpcharbychargetlcp : fastgetlcp)
                               (&mm_oldsuffix,
                                &mm_newsuffix,
                                trierep,
                                leafinsubtree->either.nodestartpos,
                                suffixtable[idx] + offset);
      currentnode = trierep->root;
      for (stackidx=0;stackidx<trierep->stack.nextfreeNodeptr;stackidx++)
      {
        currentnode = trierep->stack.spaceNodeptr[stackidx];
        if (ISLEAF(currentnode) || currentnode->depth >= lcp)
        {
          break;
        }
      }
      insertsuffixintoblindtrie(trierep,
                                currentnode,
                                mm_oldsuffix,
                                lcp,
                                mm_newsuffix,
                                suffixtable[idx] + offset);
#ifdef SKDEBUG
      printf("step %lu\n",i);
      showblindtrie(trierep);
      checkcurrentblindtrie(trierep);
#endif
    } else
    {
      break;
    }
  }
  (void) enumeratetrieleaves (suffixtable, lcpsubtab, &numoflargelcpvalues,
                              trierep);
  if (lcpsubtab != NULL)
  {
    if (idx < numberofsuffixes && offset >= (Seqpos) LCPOVERFLOW)
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
