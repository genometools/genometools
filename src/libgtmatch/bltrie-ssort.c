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

typedef struct
{
  Blindtrienode *previous,
                *current;
} Nodepair;

static Nodeptr newBlindtrienode(Blindtrierep *trierep)
{
  assert(trierep->nextfreeBlindtrienode < trierep->allocatedBlindtrienode);
  return trierep->spaceBlindtrienode + trierep->nextfreeBlindtrienode++;
}

static Blindtrienode *makenewleaf(Blindtrierep *trierep,
                                  Seqpos startpos,Seqpos firstcharpos)
{
  Blindtrienode *newleaf;

  newleaf = newBlindtrienode(trierep);
  STARTPOS(newleaf) = startpos;
  if (startpos + firstcharpos >= trierep->totallength)
  {
    FIRSTCHAR(newleaf) = (Uchar) SEPARATOR;
  } else
  {
    FIRSTCHAR(newleaf) = getencodedchar(trierep->encseq, /* Random access */
                                        startpos + firstcharpos,
                                        trierep->readmode);
  }
  DEPTH(newleaf) = 0;
  ISLEAF(newleaf) = true;
  RIGHTSIBLING(newleaf) = NULL;
  return newleaf;
}

static Nodeptr makeroot(Blindtrierep *trierep,Seqpos startpos)
{
  Blindtrienode *root, *newleaf;

  root = newBlindtrienode(trierep);
  DEPTH(root) = 0;
  FIRSTCHAR(root) = 0; /* undefined */
  RIGHTSIBLING(root) = BlindtrieNULL;
  ISLEAF(root) = false;
  newleaf = makenewleaf(trierep,startpos,0);
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

static void makesuccs(Blindtrienode *newbranch,Blindtrienode *first,
                      Blindtrienode *second)
{
  RIGHTSIBLING(second) = BlindtrieNULL;
  RIGHTSIBLING(first) = second;
  FIRSTCHILD(newbranch) = first;
}

static void findsplit(Nodepair *np,
                      const Blindtrienode *node,
                      Uchar cc_newsuffix)
{
  Uchar edgechar;

  assert(!ISLEAF(node));
  for (np->previous = NULL, np->current = FIRSTCHILD(node);
       np->current != NULL;
       np->current = RIGHTSIBLING(np->current))
  {
    edgechar = FIRSTCHAR(np->current);
    assert(edgechar != cc_newsuffix || ISSPECIAL(cc_newsuffix));
    if (cc_newsuffix < edgechar)
    {
      break;
    }
    np->previous = np->current;
  }
}

static Blindtrienode *makenewbranch(Blindtrierep *trierep,
                                    Seqpos currentstartpos,
                                    Seqpos newdepth,
                                    Uchar mm_oldsuffix,
                                    Uchar mm_newsuffix,
                                    Blindtrienode *oldnode)
{
  Blindtrienode *newbranch, *newleaf;

  assert(mm_oldsuffix != mm_newsuffix || ISSPECIAL(mm_oldsuffix));
  newbranch = newBlindtrienode(trierep);
  RIGHTSIBLING(newbranch) = RIGHTSIBLING(oldnode);
  newleaf = makenewleaf(trierep,currentstartpos,newdepth);
  if (mm_oldsuffix < mm_newsuffix || ISSPECIAL(mm_newsuffix))
  {
    makesuccs(newbranch,oldnode,newleaf);
  } else
  {
    makesuccs(newbranch,newleaf,oldnode);
  }
  DEPTH(newbranch) = newdepth;
  return newbranch;
}

static void insertcurrentsuffix(Blindtrierep *trierep,Seqpos startpos,
                                Seqpos lcp, Uchar mm_newsuffix,
                                Nodeptr currentnode)
{
  Uchar treecc;
  Nodeptr p, *pp;

  treecc = FIRSTCHAR(currentnode);
  assert(ISSPECIAL(treecc) || ISSPECIAL(mm_newsuffix) ||
         treecc != mm_newsuffix || ISLEAF(currentnode) ||
         DEPTH(currentnode) == lcp);

  /* insert a new node before node currentnode if necessary */
  if (DEPTH(currentnode) != lcp)
  {
    p = newBlindtrienode(trierep);
    FIRSTCHAR(p) = mm_newsuffix;
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

static Seqpos getlcp(Uchar *mm_oldsuffix,
                     Uchar *mm_newsuffix,
                     const Encodedsequence *encseq,
                     Readmode readmode,
                     Seqpos totallength,
                     Seqpos leafpos,
                     Seqpos currentstartpos)
{
  Seqpos idx1, idx2;
  Uchar cc1, cc2;

  for (idx1 = currentstartpos, idx2=leafpos;
       /* Nothing */;
       idx1++, idx2++)
  {
    if (idx1 < totallength)
    {
      cc1 = getencodedchar(/*XXX*/ encseq,idx1,readmode);
    } else
    {
      cc1 = (Uchar) SEPARATOR;
    }
    if (idx2 < totallength)
    {
      cc2 = getencodedchar(/*XXX*/ encseq,idx2,readmode);
    } else
    {
      cc2 = (Uchar) SEPARATOR;
    }
    if (cc1 != cc2 || ISSPECIAL(cc1))
    {
      *mm_oldsuffix = cc1;
      *mm_newsuffix = cc2;
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
                         Seqpos offset)
{
  unsigned long i, t;
  Nodeptr leafinsubtree, currentnode;
  Seqpos lcp;
  Uchar mm_oldsuffix, mm_newsuffix;

  trierep->nextfreeBlindtrienode = 0;
  trierep->root = makeroot(trierep,suffixtable[0] + offset);
  printf("insert suffixes at offset " FormatSeqpos ":\n",
          PRINTSeqposcast(offset));
  for (i=0; i < numberofsuffixes; i++)
  {
    printf(FormatSeqpos " ",PRINTSeqposcast(suffixtable[i] + offset));
  }
  printf("\nstep 0\n");
  showblindtrie(trierep);
  for (i=1UL; i < numberofsuffixes; i++)
  {
    leafinsubtree = findcompanion(trierep, suffixtable[i] + offset);
    assert(ISLEAF(leafinsubtree));
    lcp = getlcp(&mm_oldsuffix,&mm_newsuffix,
                 trierep->encseq,trierep->readmode,
                 trierep->totallength,
                 STARTPOS(leafinsubtree),
                 suffixtable[i] + offset);
    currentnode = ROOTNODE;
    for (t=0;t<trierep->stack.nextfreeNodeptr;t++)
    {
      currentnode = trierep->stack.spaceNodeptr[t];
      if (ISLEAF(currentnode) || DEPTH(currentnode) >= lcp)
      {
        break;
      }
    }
    insertcurrentsuffix(trierep,suffixtable[i]+offset,
                        lcp,mm_newsuffix,currentnode);
    printf("step %lu\n",i);
    showblindtrie(trierep);
  }
  enumeratetrieleaves (suffixtable, lcpsubtab, trierep, offset);
}
