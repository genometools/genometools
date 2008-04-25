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
#define SUFFIXPOS(N)    (EITHER(N).startpos)
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

static Nodenum makenewleaf(Blindtrierep *trierep,Seqpos startpos)
{
  Nodenum newleaf;

  newleaf = newBlindtrienode(trierep);
  ISLEAF(newleaf) = true;
  DEPTH(newleaf) = 0;
  SUFFIXPOS(newleaf) = startpos;
  RIGHTSIBLING(newleaf) = BlindtrieNULL;
  return newleaf;
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
      if (cc < FIRSTCHAR(currentnodenum))/* no branch corresponding to c */
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
  Nodenum h = ROOTNODE, p, *pp;

  for (t=0;t<trierep->stack.nextfreeNodenum;t++)
  {
    h = trierep->stack.spaceNodenum[t];
    if (ISLEAF(h) || DEPTH(h) >= lcp)
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
  assert(cc != mmchar || ISLEAF(h) || DEPTH(h) == lcp);

  /* insert a new node before node *h if necessary */
  if (DEPTH(h) != lcp)
  {
    p = newBlindtrienode(trierep);
    FIRSTCHAR(p) = mmchar;
    DEPTH(p) = DEPTH(h);  /* p inherits skip and children of *h */
    ISLEAF(p) = ISLEAF(h);
    EITHER(p) = EITHER(h);
    RIGHTSIBLING(p) = BlindtrieNULL;
    DEPTH(h) = lcp;
    ISLEAF(h) = false;
    FIRSTCHILD(h) = p;        /* now *h has p as the only child */
  }
  assert(DEPTH(h) == lcp);

  /* -------- search the position of s[n] among *h offsprings */
  pp = &FIRSTCHILD(h);
  while ((*pp) != BlindtrieNULL)
  {
    if (FIRSTCHAR(*pp) >= cc)
    {
      break;
    }
    pp = &RIGHTSIBLING(*pp);
  }
  /* ------- insert new node containing suf */
  p = newBlindtrienode(trierep);
  DEPTH(p) = 0;
  ISLEAF(p) = true;
  FIRSTCHAR(p) = cc;
  RIGHTSIBLING(p) = *pp;
  *pp = p;
  SUFFIXPOS(p) = startpos;
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

void blindtriesuffixsort(Blindtrierep *trierep,
                         Seqpos *suffixtable,
                         unsigned long numberofsuffixes,
                         Seqpos depth)
{
  unsigned long i;
  Nodenum h;
  Seqpos lcp, currentstartpos;
  Uchar cc;

  trierep->stack.nextfreeNodenum = 0;
  (void) makenewleaf(trierep,suffixtable[0] + depth);
  for (i=1UL; i < numberofsuffixes; i++)
  {
    currentstartpos = suffixtable[i] + depth;
    h = findcompanion(trierep, currentstartpos);
    assert(ISLEAF(h));
    lcp = getlcp(trierep->encseq,trierep->readmode,
                 currentstartpos,trierep->totallength-1,
                 SUFFIXPOS(h),trierep->totallength-1);
    cc = getencodedchar(trierep->encseq, /* Random access */
                        currentstartpos + lcp,
                        trierep->readmode);
    insertcurrentsuffix(trierep,currentstartpos,lcp,cc);
  }
}
