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
#include "bltrie-ssort.h"

/*
#define NODE(N)         (trierep->spaceBlindtrienode[N])
#define DEPTH(N)        (NODE(N).depth)
#define EITHER(N)       (NODE(N).either)
#define FIRSTCHILD(N)   (EITHER(N).firstchild)
#define STARTPOS(N)     (EITHER(N).startpos)
#define RIGHTSIBLING(N) (NODE(N).rightsibling)
#define ISLEAF(N)       (NODE(N).isleaf)
#define FIRSTCHAR(N)    (NODE(N).firstchar)
#define ROOTNODE        trierep->root
#define BlindtrieNULL   (trierep->allocatedBlindtrienode)
*/

#define DEPTH(N)        ((N)->depth)
#define EITHER(N)       ((N)->either)
#define FIRSTCHILD(N)   (EITHER(N).firstchild)
#define STARTPOS(N)     (EITHER(N).startpos)
#define RIGHTSIBLING(N) ((N)->rightsibling)
#define ISLEAF(N)       ((N)->isleaf)
#define FIRSTCHAR(N)    ((N)->firstchar)
#define ROOTNODE        trierep->root
#define BlindtrieNULL   NULL

typedef struct Blindtrienode
{
  Seqpos depth;
  struct Blindtrienode *rightsibling;
  union
  {
    struct Blindtrienode *firstchild;
    Seqpos startpos;
  } either;
  bool isleaf;
  Uchar firstchar;
} Blindtrienode;

typedef Blindtrienode * Nodeptr;

DECLAREARRAYSTRUCT(Nodeptr);

struct Blindtrierep
{
  const Encodedsequence *encseq;
  Readmode readmode;
  Seqpos totallength;
  Blindtrienode *spaceBlindtrienode;
  Nodeptr root;
  unsigned long allocatedBlindtrienode,
                nextfreeBlindtrienode;
  ArrayNodeptr stack;
};

static Nodeptr newBlindtrienode(Blindtrierep *trierep)
{
  assert(trierep->nextfreeBlindtrienode < trierep->allocatedBlindtrienode);
  return trierep->spaceBlindtrienode + trierep->nextfreeBlindtrienode++;
}

static Blindtrienode *makenewleaf(Blindtrierep *trierep,Seqpos startpos,
                                  Seqpos depth)
{
  Blindtrienode *newleaf;

  newleaf = newBlindtrienode(trierep);
  STARTPOS(newleaf) = startpos;
  if (startpos + depth >= trierep->totallength)
  {
    FIRSTCHAR(newleaf) = (Uchar) SEPARATOR;
  } else
  {
    FIRSTCHAR(newleaf) = getencodedchar(trierep->encseq, /* Random access */
                                        startpos + depth,
                                        trierep->readmode);
  }
  DEPTH(newleaf) = 0;
  ISLEAF(newleaf) = true;
  RIGHTSIBLING(newleaf) = NULL;
  return newleaf;
}

static Nodeptr makeroot(Blindtrierep *trierep,Seqpos startpos,
                               Seqpos depth)
{
  Blindtrienode *root, *newleaf;

  root = newBlindtrienode(trierep);
  DEPTH(root) = 0;
  FIRSTCHAR(root) = 0;
  RIGHTSIBLING(root) = NULL;
  ISLEAF(root) = false;
  newleaf = makenewleaf(trierep,startpos,depth);
  FIRSTCHILD(root) = newleaf;
  return root;
}

static Nodeptr extractleafnode(Nodeptr head)
{
  assert(!ISLEAF(head));
  do
  {
    head = FIRSTCHILD(head);
  } while (!ISLEAF(head));
  return head;
}

static Nodeptr findcompanion(Blindtrierep *trierep,Seqpos startpos)
{
  Uchar cc;
  Nodeptr head, currentnodenum;

  trierep->stack.nextfreeNodeptr = 0;
  head = ROOTNODE;
  while (!ISLEAF(head))
  {
    STOREINARRAY (&trierep->stack, Nodeptr, 128, head);
    if (startpos + DEPTH(head) >= trierep->totallength)
    {
      cc = (Uchar) SEPARATOR;
    } else
    {
      cc = getencodedchar(trierep->encseq, /* Random access */
                          startpos + DEPTH(head),
                          trierep->readmode);
    }
    if (ISSPECIAL(cc))
    {
      return extractleafnode(head);
    }
    currentnodenum = FIRSTCHILD(head);
    for (;;)
    {
      if (cc == FIRSTCHAR(currentnodenum))
      {              /* found branch corresponding to cc */
        head = FIRSTCHILD(currentnodenum);
        break;
      }
      if (cc < FIRSTCHAR(currentnodenum)) /* no branch corresponding to c */
      {
        return extractleafnode(head);
      }
      currentnodenum = RIGHTSIBLING(currentnodenum);
      if (currentnodenum == BlindtrieNULL) /* no other branches: mismatch */
      {
        return extractleafnode(head);
      }
    }
  }
  STOREINARRAY (&trierep->stack, Nodeptr, 128, head);
  return head;
}

static void insertcurrentsuffix(Blindtrierep *trierep,Seqpos startpos,
                                Seqpos lcp, Uchar mismatchchar,
                                Nodeptr currentnode)
{
  Uchar treecc;
  Nodeptr p, *pp;

  treecc = FIRSTCHAR(currentnode);
  assert(ISSPECIAL(treecc) || ISSPECIAL(mismatchchar) ||
         treecc != mismatchchar || ISLEAF(currentnode) ||
         DEPTH(currentnode) == lcp);

  /* insert a new node before node currentnode if necessary */
  if (DEPTH(currentnode) != lcp)
  {
    p = newBlindtrienode(trierep);
    FIRSTCHAR(p) = mismatchchar;
    DEPTH(p) = DEPTH(currentnode); /* p inherits depth+children of currentnode*/
    ISLEAF(p) = ISLEAF(currentnode);
    EITHER(p) = EITHER(currentnode);
    RIGHTSIBLING(p) = BlindtrieNULL;
    DEPTH(currentnode) = lcp;
    ISLEAF(currentnode) = false;
    FIRSTCHILD(currentnode) = p;        /* now *h has p as the only child */
  }
  assert(DEPTH(currentnode) == lcp);

  /* search the position of s[n] among *h offsprings */
  pp = &FIRSTCHILD(currentnode);
  while ((*pp) != BlindtrieNULL)
  {
    if (FIRSTCHAR(*pp) >= treecc)
    {
      break;
    }
    pp = &RIGHTSIBLING(*pp);
  }
  /* insert new node containing suf */
  p = newBlindtrienode(trierep);
  DEPTH(p) = 0;
  ISLEAF(p) = true;
  FIRSTCHAR(p) = treecc;
  RIGHTSIBLING(p) = *pp;
  *pp = p;
  STARTPOS(p) = startpos;
}

static Seqpos getlcp(Uchar *mismatchchar,
                     const Encodedsequence *encseq,
                     Readmode readmode,
                     Seqpos totallength,
                     Seqpos currentstartpos,
                     Seqpos leafpos)
{
  Seqpos idx1, idx2;
  Uchar cc1, cc2;

  for (idx1 = currentstartpos, idx2=leafpos;
       /* Nothing */;
       idx1++, idx2++)
  {
    if (idx1 == totallength)
    {
      *mismatchchar = (Uchar) SEPARATOR;
      break;
    }
    cc1 = getencodedchar(/*XXX*/ encseq,idx1,readmode);
    if (ISSPECIAL(cc1))
    {
      *mismatchchar = cc1;
      break;
    }
    if (idx2 == totallength)
    {
      *mismatchchar = cc1;
      break;
    }
    cc2 = getencodedchar(/*XXX*/ encseq,idx2,readmode);
    if (ISSPECIAL(cc2) || cc1 != cc2)
    {
      *mismatchchar = cc1;
      break;
    }
  }
  return idx1 - currentstartpos;
}

#define SETCURRENT(VAL)\
        currentnode = VAL;\
        currentnodeisleaf = ISLEAF(VAL) ? true : false

static void enumeratetrieleaves (Seqpos *suffixtable,
                                 Seqpos *lcpsubtab,
                                 Blindtrierep *trierep,
                                 Seqpos depth)
{
  bool readyforpop = false, currentnodeisleaf;
  Nodeptr currentnode, siblval, child, lcpnode = ROOTNODE;
  unsigned long nextfree = 0;

  trierep->stack.nextfreeNodeptr = 0;
  STOREINARRAY (&trierep->stack, Nodeptr, 128, ROOTNODE);
  SETCURRENT(FIRSTCHILD(ROOTNODE));
  for (;;)
  {
    if (currentnodeisleaf)
    {
      assert (STARTPOS(currentnode) >= depth);
      if (lcpsubtab != NULL && nextfree > 0)
      {
        lcpsubtab[nextfree] = DEPTH(lcpnode) + depth;
      }
      suffixtable[nextfree++] = STARTPOS(currentnode) - depth;
      siblval = RIGHTSIBLING(currentnode);
      if (siblval == BlindtrieNULL)
      {
        readyforpop = true;
        currentnodeisleaf = false;
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
        siblval = RIGHTSIBLING(trierep->stack.spaceNodeptr[
                               trierep->stack.nextfreeNodeptr]);
        if (siblval != BlindtrieNULL)
        {
          SETCURRENT (siblval);        /* current comes from brother */
          lcpnode = trierep->stack.spaceNodeptr[
                             trierep->stack.nextfreeNodeptr - 1];
          readyforpop = false;
        }
      } else
      {
        STOREINARRAY (&trierep->stack, Nodeptr, 128, currentnode);
        child = FIRSTCHILD(currentnode);
        SETCURRENT (child);            /* current comes from child */
      }
    }
  }
}

Blindtrierep *newBlindtrierep(unsigned long numofsuffixes,
                              const Encodedsequence *encseq,
                              Readmode readmode)
{
  Blindtrierep *trierep;

  ALLOCASSIGNSPACE(trierep,NULL,Blindtrierep,1);
  trierep->allocatedBlindtrienode = MULT2(numofsuffixes + 1) + 1;
  ALLOCASSIGNSPACE(trierep->spaceBlindtrienode,NULL,Blindtrienode,
                   trierep->allocatedBlindtrienode);
  trierep->nextfreeBlindtrienode = 0;
  trierep->encseq = encseq;
  trierep->readmode = readmode;
  trierep->root = NULL;
  trierep->totallength = getencseqtotallength(encseq);
  INITARRAY (&trierep->stack, Nodeptr);
  return trierep;
}

void freeBlindtrierep(Blindtrierep **trierep)
{
  FREESPACE((*trierep)->spaceBlindtrienode);
  FREEARRAY (&(*trierep)->stack, Nodeptr);
  FREESPACE(*trierep);
}

#define NODENUM(PTR) (unsigned long) (RIGHTSIBLING(PTR) - \
                                      trierep->spaceBlindtrienode)

static void showleaf(const Blindtrierep *trierep,unsigned int level,
                     Nodeptr current)
{
  printf("%*.*s",(int) (6 * level),(int) (6 * level)," ");
  printf("Leaf(add=%lu,firstchar=%u,startpos=" FormatSeqpos
         ",rightsibling=%lu)\n",
         NODENUM(current),
         (unsigned int) FIRSTCHAR(current),
         PRINTSeqposcast(STARTPOS(current)),
         NODENUM(RIGHTSIBLING(current)));
}

static void showintern(const Blindtrierep *trierep,unsigned int level,
                       Nodeptr current)
{
  printf("%*.*s",(int) (6 * level),(int) (6 * level)," ");
  printf("Intern(add=%lu,firstchar=%u,depth=" FormatSeqpos
         ",firstchild=%lu,rightsibling=%lu)\n",
          NODENUM(current),
          (unsigned int) FIRSTCHAR(current),
          PRINTSeqposcast(DEPTH(current)),
          NODENUM(FIRSTCHILD(current)),
          NODENUM(RIGHTSIBLING(current)));
}

static void showblindtrie2(const Blindtrierep *trierep,
                           unsigned int level,
                           Nodeptr node)
{
  Nodeptr current;

  for (current = FIRSTCHILD(node);
       current != BlindtrieNULL;
       current = RIGHTSIBLING(current))
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
  showblindtrie2(trierep,0,ROOTNODE);
}

void blindtriesuffixsort(Blindtrierep *trierep,
                         Seqpos *suffixtable,
                         Seqpos *lcpsubtab,
                         unsigned long numberofsuffixes,
                         Seqpos depth)
{
  unsigned long i, t;
  Nodeptr leafinsubtree, currentnode;
  Seqpos lcp;
  Uchar mismatchchar;

  trierep->nextfreeBlindtrienode = 0;
  trierep->root = makeroot(trierep,suffixtable[0],depth);
  printf("insert suffixes at depth " FormatSeqpos ":\n",PRINTSeqposcast(depth));
  for (i=0; i < numberofsuffixes; i++)
  {
    printf(FormatSeqpos " ",PRINTSeqposcast(suffixtable[i]));
  }
  printf("\nstep 0\n");
  showblindtrie(trierep);
  for (i=1UL; i < numberofsuffixes; i++)
  {
    leafinsubtree = findcompanion(trierep, suffixtable[i]);
    assert(ISLEAF(leafinsubtree));
    lcp = depth + getlcp(&mismatchchar,trierep->encseq,trierep->readmode,
                         trierep->totallength,
                         suffixtable[i] + depth,
                         STARTPOS(leafinsubtree) + depth);
    currentnode = ROOTNODE;
    for (t=0;t<trierep->stack.nextfreeNodeptr;t++)
    {
      currentnode = trierep->stack.spaceNodeptr[t];
      if (ISLEAF(currentnode) || DEPTH(currentnode) >= lcp)
      {
        break;
      }
    }
    insertcurrentsuffix(trierep,suffixtable[i],lcp,mismatchchar,currentnode);
    printf("step %lu\n",i);
    showblindtrie(trierep);
  }
  enumeratetrieleaves (suffixtable, lcpsubtab, trierep, depth);
}
