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

#define NODE(N)         (trierep->spaceBlindtrienode[N])
#define DEPTH(N)        (NODE(N).depth)
#define EITHER(N)       (NODE(N).either)
#define FIRSTCHILD(N)   (EITHER(N).firstchild)
#define STARTPOS(N)     (EITHER(N).startpos)
#define RIGHTSIBLING(N) (NODE(N).rightsibling)
#define ISLEAF(N)       (NODE(N).isleaf)
#define FIRSTCHAR(N)    (NODE(N).firstchar)
#define ROOTNODE        0
#define BlindtrieNULL   (trierep->allocatedBlindtrienode)

typedef unsigned int Nodenum;

typedef struct Blindtrienode
{
  Seqpos depth;
  Nodenum rightsibling;
  union
  {
    Nodenum firstchild;
    Seqpos startpos;
  } either;
  bool isleaf;
  Uchar firstchar;
} Blindtrienode;

DECLAREARRAYSTRUCT(Nodenum);

struct Blindtrierep
{
  const Encodedsequence *encseq;
  Readmode readmode;
  Seqpos totallength;
  Blindtrienode *spaceBlindtrienode;
  Nodenum allocatedBlindtrienode,
          nextfreeBlindtrienode;
  ArrayNodenum stack;
};

static Nodenum newBlindtrienode(Blindtrierep *trierep)
{
  assert(trierep->nextfreeBlindtrienode < trierep->allocatedBlindtrienode);
  return trierep->nextfreeBlindtrienode++;
}

static void makefirstleaf(Blindtrierep *trierep,Seqpos startpos,Seqpos depth)
{
  Nodenum newleaf;

  newleaf = newBlindtrienode(trierep);
  ISLEAF(newleaf) = true;
  DEPTH(newleaf) = 0;
  if (startpos + depth >= trierep->totallength)
  {
    FIRSTCHAR(newleaf) = (Uchar) SEPARATOR;
  } else
  {
    FIRSTCHAR(newleaf) = getencodedchar(trierep->encseq, /* Random access */
                                        startpos + depth,
                                        trierep->readmode);
  }
  STARTPOS(newleaf) = startpos;
  RIGHTSIBLING(newleaf) = BlindtrieNULL;
}

static Nodenum extractleafnode(const Blindtrierep *trierep,Nodenum head)
{
  assert(!ISLEAF(head));
  do
  {
    head = FIRSTCHILD(head);
  } while (!ISLEAF(head));
  return head;
}

static Nodenum findcompanion(Blindtrierep *trierep,Seqpos startpos)
{
  Uchar cc;
  Nodenum head, currentnodenum;

  trierep->stack.nextfreeNodenum = 0;
  head = ROOTNODE;
  while (!ISLEAF(head))
  {
    STOREINARRAY (&trierep->stack, Nodenum, 128, head);
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
      return extractleafnode(trierep,head);
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
        return extractleafnode(trierep,head);
      }
      currentnodenum = RIGHTSIBLING(currentnodenum);
      if (currentnodenum == BlindtrieNULL) /* no other branches: mismatch */
      {
        return extractleafnode(trierep,head);
      }
    }
  }
  STOREINARRAY (&trierep->stack, Nodenum, 128, head);
  return head;
}

static void insertcurrentsuffix(Blindtrierep *trierep,Seqpos startpos,
                                Seqpos lcp, Uchar mmchar)
{
  unsigned long t;
  Uchar cc;
  Nodenum currentnode = ROOTNODE, p, *pp;

  for (t=0;t<trierep->stack.nextfreeNodenum;t++)
  {
    currentnode = trierep->stack.spaceNodenum[t];
    if (ISLEAF(currentnode) || DEPTH(currentnode) >= lcp)
    {
      break;
    }
  }
  if (startpos + lcp >= trierep->totallength)
  {
    cc = (Uchar) SEPARATOR;
  } else
  {
    cc = getencodedchar(trierep->encseq, /* Random access */
                        startpos + lcp,
                        trierep->readmode);
  }
  assert(ISSPECIAL(cc) || ISSPECIAL(mmchar) ||
         cc != mmchar || ISLEAF(currentnode) || DEPTH(currentnode) == lcp);

  /* insert a new node before node currentnode if necessary */
  if (DEPTH(currentnode) != lcp)
  {
    p = newBlindtrienode(trierep);
    FIRSTCHAR(p) = mmchar;
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
    if (FIRSTCHAR(*pp) >= cc)
    {
      break;
    }
    pp = &RIGHTSIBLING(*pp);
  }
  /* insert new node containing suf */
  p = newBlindtrienode(trierep);
  DEPTH(p) = 0;
  ISLEAF(p) = true;
  FIRSTCHAR(p) = cc;
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
  Nodenum currentnode, siblval, child, lcpnode = 0;
  unsigned long nextfree = 0;

  trierep->stack.nextfreeNodenum = 0;
  STOREINARRAY (&trierep->stack, Nodenum, 128, ROOTNODE);
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
        lcpnode = trierep->stack.spaceNodenum[trierep->stack.nextfreeNodenum-1];
      }
    } else
    {
      if (readyforpop)
      {
        if (trierep->stack.nextfreeNodenum == 1UL)
        {
          break;
        }
        trierep->stack.nextfreeNodenum--;
        siblval = RIGHTSIBLING(trierep->stack.spaceNodenum[
                               trierep->stack.nextfreeNodenum]);
        if (siblval != BlindtrieNULL)
        {
          SETCURRENT (siblval);        /* current comes from brother */
          lcpnode = trierep->stack.spaceNodenum[
                             trierep->stack.nextfreeNodenum - 1];
          readyforpop = false;
        }
      } else
      {
        STOREINARRAY (&trierep->stack, Nodenum, 128, currentnode);
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
  trierep->allocatedBlindtrienode = (Nodenum) MULT2(numofsuffixes + 1) + 1;
  ALLOCASSIGNSPACE(trierep->spaceBlindtrienode,NULL,Blindtrienode,
                   trierep->allocatedBlindtrienode);
  trierep->nextfreeBlindtrienode = 0;
  trierep->encseq = encseq;
  trierep->readmode = readmode;
  trierep->totallength = getencseqtotallength(encseq);
  INITARRAY (&trierep->stack, Nodenum);
  return trierep;
}

void freeBlindtrierep(Blindtrierep **trierep)
{
  FREESPACE((*trierep)->spaceBlindtrienode);
  FREEARRAY (&(*trierep)->stack, Nodenum);
  FREESPACE(*trierep);
}

static void showleaf(const Blindtrierep *trierep,unsigned int level,
                     Nodenum current)
{
  printf("%*.*s",(int) (6 * level),(int) (6 * level)," ");
  printf("Leaf(add=%u,firstchar=%u,startpos=" FormatSeqpos
         ",rightsibling=%u)\n",
         current,
         (unsigned int) FIRSTCHAR(current),
         PRINTSeqposcast(STARTPOS(current)),
         RIGHTSIBLING(current));
}

static void showintern(const Blindtrierep *trierep,unsigned int level,
                       Nodenum current)
{
  printf("%*.*s",(int) (6 * level),(int) (6 * level)," ");
  printf("Intern(add=%u,firstchar=%u,depth=" FormatSeqpos
         ",rightsibling=%u)\n",
          current,
          (unsigned int) FIRSTCHAR(current),
          PRINTSeqposcast(DEPTH(current)),
          RIGHTSIBLING(current));
}

static void showblindtrie2(const Blindtrierep *trierep,
                           unsigned int level,
                           Nodenum node)
{
  Nodenum current;

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
  if (trierep->nextfreeBlindtrienode == (Nodenum) 1)
  {
    assert(ISLEAF(ROOTNODE));
    showleaf(trierep,1U,ROOTNODE);
  } else
  {
    assert(!ISLEAF(ROOTNODE));
    showintern(trierep,0,ROOTNODE);
    showblindtrie2(trierep,1U,ROOTNODE);
  }
}

void blindtriesuffixsort(Blindtrierep *trierep,
                         Seqpos *suffixtable,
                         Seqpos *lcpsubtab,
                         unsigned long numberofsuffixes,
                         Seqpos depth)
{
  unsigned long i;
  Nodenum leafinsubtree;
  Seqpos lcp;
  Uchar mismatchchar;

  trierep->nextfreeBlindtrienode = 0;
  makefirstleaf(trierep,suffixtable[0],depth);
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
    insertcurrentsuffix(trierep,suffixtable[i],lcp,mismatchchar);
    printf("step %lu\n",i);
    showblindtrie(trierep);
  }
  enumeratetrieleaves (suffixtable, lcpsubtab, trierep, depth);
}
