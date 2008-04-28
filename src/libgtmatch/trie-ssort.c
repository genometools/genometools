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

#include "libgtcore/arraydef.h"
#include "divmodmul.h"
#include "seqpos-def.h"
#include "spacedef.h"
#include "encseq-def.h"
#include "trie-ssort.h"

#define ISLEAF(NODE) ((NODE)->firstchild == NULL)

typedef struct Trienode
{
  Seqpos startpos,
         depth;
  struct Trienode *firstchild,
                  *rightsibling;
} Trienode;

typedef Trienode * Trienodeptr;

DECLAREARRAYSTRUCT(Trienodeptr);

struct Trierep
{
  const Encodedsequence *encseq;
  Seqpos totallength;
  Readmode readmode;
  Trienode *nodetable,
           *root;
  unsigned long allocatedTrienode,
                nextfreeTrienode;
  ArrayTrienodeptr stack;
};

typedef struct
{
  Trienode *previous,
           *current;
} Nodepair;

static Uchar getfirstedgechar(const Trierep *trierep,
                              const Trienode *node,
                              Seqpos prevdepth)
{
  if (ISLEAF(node) && node->startpos + prevdepth >= trierep->totallength)
  {
    return (Uchar) SEPARATOR;
  }
  return getencodedchar(trierep->encseq, /* Random access */
                        node->startpos + prevdepth,
                        trierep->readmode);
}

static int comparecharacters(Uchar cc1,Seqpos idx1,
                             Uchar cc2,Seqpos idx2)
{
  if (ISSPECIAL(cc1))
  {
    if (ISSPECIAL(cc2))
    {
      if (idx1 <= idx2)
      {
        return -1;  /* cc1 < cc2 */
      } else
      {
        return 1;  /* cc1 > cc2 */
      }
    } else
    {
      return 1; /* cc1 > cc2 */
    }
  } else
  {
    if (ISSPECIAL(cc2))
    {
      return -1;  /* cc1 < cc2 */
    } else
    {
      if (cc1 < cc2)
      {
        return -1;  /* cc1 < cc2 */
      } else
      {
        if (cc1 > cc2)
        {
          return 1;  /* cc1 > cc2 */
        } else
        {
          return 0; /* cc1 == cc2 */
        }
      }
    }
  }
}

static Trienode *newTrienode(Trierep *trierep)
{
  assert(trierep->nextfreeTrienode < trierep->allocatedTrienode);
  return trierep->nodetable + trierep->nextfreeTrienode++;
}

static Trienode *makenewleaf(Trierep *trierep,Seqpos startpos)
{
  Trienode *newleaf;

  newleaf = newTrienode(trierep);
  newleaf->startpos = startpos;
  newleaf->firstchild = NULL;
  newleaf->rightsibling = NULL;
  return newleaf;
}

static Trienode *makeroot(Trierep *trierep,Seqpos startpos)
{
  Trienode *root, *newleaf;

  root = newTrienode(trierep);
  root->startpos = startpos;
  root->depth = 0;
  root->rightsibling = NULL;
  newleaf = makenewleaf(trierep,startpos);
  root->firstchild = newleaf;
  return root;
}

static void makesuccs(Trienode *newbranch,Trienode *first,
                      Trienode *second)
{
  second->rightsibling = NULL;
  first->rightsibling = second;
  newbranch->firstchild = first;
}

static Trienode *makenewbranch(Trierep *trierep,
                               Seqpos startpos,
                               Seqpos currentdepth,
                               Trienode *oldnode)
{
  Trienode *newbranch, *newleaf;
  Uchar cc1, cc2;

  newbranch = newTrienode(trierep);
  newbranch->startpos = startpos;
  newbranch->rightsibling = oldnode->rightsibling;
  newleaf = makenewleaf(trierep,startpos);
  cc1 = getfirstedgechar(trierep,oldnode,currentdepth);
  if (startpos + currentdepth >= trierep->totallength)
  {
    cc2 = (Uchar) SEPARATOR;
  } else
  {
    cc2 = getencodedchar(trierep->encseq, /* Random access */
                         startpos + currentdepth,
                         trierep->readmode);
  }
  if (comparecharacters(cc1,oldnode->startpos,cc2,startpos) <= 0)
  {
    makesuccs(newbranch,oldnode,newleaf);
  } else
  {
    makesuccs(newbranch,newleaf,oldnode);
  }
  newbranch->depth = currentdepth;
  return newbranch;
}

static Seqpos getlcp(const Encodedsequence *encseq,
                     Readmode readmode,
                     Seqpos start1,Seqpos end1,
                     Seqpos start2,Seqpos end2)
{
  Seqpos idx1, idx2;
  Uchar cc1, cc2;

  for (idx1=start1, idx2=start2; idx1 <= end1 && idx2 <= end2; idx1++, idx2++)
  {
    cc1 = getencodedchar(/*XXX*/ encseq,idx1,readmode);
    if (ISSPECIAL(cc1))
    {
      break;
    }
    cc2 = getencodedchar(/*XXX*/ encseq,idx2,readmode);
    if (ISSPECIAL(cc2) || cc1 != cc2)
    {
      break;
    }
  }
  return idx1 - start1;
}

static bool hasccsuccessor(const Trierep *trierep,
                           Nodepair *np,
                           const Trienode *node,
                           Uchar cc2,
                           Seqpos startpos2)
{
  Uchar cc1;
  int cmpresult;

  for (np->previous = NULL, np->current = node->firstchild;
       np->current != NULL;
       np->current = np->current->rightsibling)
  {
    cc1 = getfirstedgechar(trierep,np->current,node->depth);
    cmpresult = comparecharacters(cc1,np->current->startpos,cc2,startpos2);
    if (cmpresult == 1)
    {
      return false;
    }
    if (cmpresult == 0)
    {
      return true;
    }
    np->previous = np->current;
  }
  return false;
}

static void insertsuffixintotrie(Trierep *trierep,Seqpos startpos)
{
  if (trierep->root != NULL)
  {
    Seqpos lcpvalue;
    Trienode *currentnode, *newleaf, *newbranch, *succ;
    Nodepair np;
    Uchar cc;

    currentnode = trierep->root;
    while (true)
    {
      if (startpos + currentnode->depth >= trierep->totallength)
      {
        cc = (Uchar) SEPARATOR;
      } else
      {
        cc = getencodedchar(trierep->encseq, /* Random access */
                            startpos + currentnode->depth,
                            trierep->readmode);
      }
      assert(currentnode != NULL);
      assert(!ISLEAF(currentnode));
      if (!hasccsuccessor(trierep,&np,currentnode,cc,startpos))
      {
        newleaf = makenewleaf(trierep,startpos);
        newleaf->rightsibling = np.current;
        if (np.previous == NULL)
        {
          currentnode->firstchild = newleaf;
        } else
        {
          np.previous->rightsibling = newleaf;
        }
        return;
      }
      succ = np.current;
      if (ISLEAF(succ))
      {
        lcpvalue = getlcp(trierep->encseq,
                          trierep->readmode,
                          startpos + currentnode->depth + 1,
                          trierep->totallength - 1,
                          succ->startpos + currentnode->depth + 1,
                          trierep->totallength - 1);
        newbranch = makenewbranch(trierep,
                                  startpos,
                                  currentnode->depth + lcpvalue + 1,
                                  succ);
        if (np.previous == NULL)
        {
          currentnode->firstchild = newbranch;
        } else
        {
          np.previous->rightsibling = newbranch;
        }
        return;
      }
      lcpvalue = getlcp(trierep->encseq,
                        trierep->readmode,
                        startpos + currentnode->depth + 1,
                        trierep->totallength - 1,
                        succ->startpos + currentnode->depth + 1,
                        succ->startpos + succ->depth - 1);
      if (currentnode->depth + lcpvalue + 1 < succ->depth)
      {
        newbranch = makenewbranch(trierep,
                                  startpos,
                                  currentnode->depth + lcpvalue + 1,
                                  succ);
        if (np.previous == NULL)
        {
          currentnode->firstchild = newbranch;
        } else
        {
          np.previous->rightsibling = newbranch;
        }
        return;
      }
      currentnode = succ;
    }
  } else
  {
    trierep->root = makeroot(trierep,startpos);
  }
}

Trierep *newTrierep(unsigned long numofsuffixes,
                    const Encodedsequence *encseq,
                    Readmode readmode)
{
  Trierep *trierep;

  ALLOCASSIGNSPACE(trierep,NULL,Trierep,1);
  trierep->allocatedTrienode = MULT2(numofsuffixes + 1) + 1;
  ALLOCASSIGNSPACE(trierep->nodetable,NULL,Trienode,
                   trierep->allocatedTrienode);
  trierep->nextfreeTrienode = 0;
  trierep->root = NULL;
  trierep->encseq = encseq;
  trierep->readmode = readmode;
  trierep->totallength = getencseqtotallength(encseq);
  INITARRAY (&trierep->stack, Trienodeptr);
  return trierep;
}

/*
static unsigned long rekenumeratetrieleaves(Seqpos *suffixtable,
                                            unsigned long nextfree,
                                            const Trierep *trierep,
                                            const Trienode *node,
                                            Seqpos depth)
{
  Trienode *current;

  for (current = node->firstchild;
       current != NULL;
       current = current->rightsibling)
  {
    if (ISLEAF(current))
    {
      assert(current->startpos >= depth);
      suffixtable[nextfree++] = current->startpos - depth;
    } else
    {
      nextfree = rekenumeratetrieleaves(suffixtable,nextfree,
                                        trierep,current,depth);
    }
  }
  return nextfree;
}
*/

#define SETCURRENT(VAL)\
        currentnode = VAL;\
        currentnodeisleaf = ISLEAF(VAL) ? true : false

static void enumeratetrieleaves (Seqpos *suffixtable,
                                 Seqpos *lcpsubtab,
                                 Trierep *trierep,
                                 Seqpos depth)
{
  bool readyforpop = false, currentnodeisleaf;
  Trienode *siblval, *child, *lcpnode = NULL, *currentnode;
  unsigned long nextfree = 0;

  trierep->stack.nextfreeTrienodeptr = 0;
  STOREINARRAY (&trierep->stack, Trienodeptr, 128, trierep->root);
  assert (!ISLEAF (trierep->root));
  SETCURRENT (trierep->root->firstchild);
  for (;;)
  {
    if (currentnodeisleaf)
    {
      assert (currentnode->startpos >= depth);
      if (lcpsubtab != NULL && nextfree > 0)
      {
        assert(lcpnode != NULL);
        lcpsubtab[nextfree] = lcpnode->depth + depth;
      }
      suffixtable[nextfree++] = currentnode->startpos - depth;
      siblval = currentnode->rightsibling;
      if (siblval == NULL)
      {
        readyforpop = true;
        currentnodeisleaf = false;
      } else
      {
        SETCURRENT (siblval);          /* current comes from brother */
        lcpnode = trierep->stack.spaceTrienodeptr[
                           trierep->stack.nextfreeTrienodeptr - 1];
      }
    } else
    {
      if (readyforpop)
      {
        if (trierep->stack.nextfreeTrienodeptr == 1UL)
        {
          break;
        }
        trierep->stack.nextfreeTrienodeptr--;
        siblval = trierep->stack.spaceTrienodeptr[
                           trierep->stack.nextfreeTrienodeptr]->
                                              rightsibling;
        if (siblval != NULL)
        {
          SETCURRENT (siblval);        /* current comes from brother */
          lcpnode = trierep->stack.spaceTrienodeptr[
                             trierep->stack.nextfreeTrienodeptr - 1];
          readyforpop = false;
        }
      } else
      {
        STOREINARRAY (&trierep->stack, Trienodeptr, 128, currentnode);
        child = currentnode->firstchild;
        SETCURRENT (child);            /* current comes from child */
      }
    }
  }
}

void triesuffixsort(Trierep *trierep,
                    Seqpos *suffixtable,
                    Seqpos *lcpsubtab,
                    unsigned long numberofsuffixes,
                    Seqpos depth)
{
  unsigned long idx;

  trierep->nextfreeTrienode = 0;
  trierep->root = NULL;
  for (idx=0; idx<numberofsuffixes; idx++)
  {
    insertsuffixintotrie(trierep,suffixtable[idx]+depth);
  }
  enumeratetrieleaves(suffixtable,lcpsubtab,trierep,depth);
}

void freeTrierep(Trierep **trierep)
{
  FREESPACE((*trierep)->nodetable);
  FREEARRAY (&(*trierep)->stack, Trienodeptr);
  FREESPACE(*trierep);
}
