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

#include "seqpos-def.h"
#include "encseq-def.h"

typedef struct Blindtrienode
{
  Seqpos startpos, depth;
  struct Blindtrienode *firstchild,
                       *rightsibling,
                       *parent;
} Blindtrienode;

typedef struct
{
  Encodedsequence *encseq;
  Readmode readmode;
  Blindtrienode *nodetable,
                *root,
                **unusedBlindtrienodes;
  unsigned int nextunused,
               allocatedBlindtrienode,
               nextfreeBlindtrienode;
} Blindtrierep;

typedef struct
{
  Blindtrienode *previous,
                *current;
} Nodepair;

#define ISLEAF(NODE) ((NODE)->firstchild == NULL)
#define SETFIRSTCHILD(NODE,VALUE)\
        (NODE)->firstchild = VALUE;\
        if ((VALUE) != NULL)\
        {\
          (VALUE)->parent = NODE;\
        }

#define SETFIRSTCHILDNULL(NODE)\
        (NODE)->firstchild = NULL

static Uchar getfirstedgechar(const Blindtrierep *trierep,
                              const Blindtrienode *node,
                              Seqpos prevdepth)
{
  if (ISLEAF(node) &&
      node->startpos + prevdepth >= getencseqtotallength(trierep->encseq))
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

static Blindtrienode *newBlindtrienode(Blindtrierep *trierep)
{
  if (trierep->nextfreeBlindtrienode >= trierep->allocatedBlindtrienode)
  {
    assert(trierep->nextunused > 0);
    trierep->nextunused--;
    return trierep->unusedBlindtrienodes[trierep->nextunused];
  }
  return trierep->nodetable + trierep->nextfreeBlindtrienode++;
}

static Blindtrienode *makenewleaf(Blindtrierep *trierep,Seqpos startpos)
{
  Blindtrienode *newleaf;

  newleaf = newBlindtrienode(trierep);
  newleaf->startpos = startpos;
  SETFIRSTCHILDNULL(newleaf);
  newleaf->rightsibling = NULL;
  return newleaf;
}

static Blindtrienode *makeroot(Blindtrierep *trierep,Seqpos startpos)
{
  Blindtrienode *root, *newleaf;

  root = newBlindtrienode(trierep);
  root->parent = NULL;
  root->startpos = startpos;
  root->depth = 0;
  root->rightsibling = NULL;
  newleaf = makenewleaf(trierep,startpos);
  SETFIRSTCHILD(root,newleaf);
  return root;
}

static void makesuccs(Blindtrienode *newbranch,Blindtrienode *first,
                      Blindtrienode *second)
{
  second->rightsibling = NULL;
  first->rightsibling = second;
  SETFIRSTCHILD(newbranch,first);
}

static Blindtrienode *makenewbranch(Blindtrierep *trierep,
                                    Seqpos startpos,
                                    Seqpos currentdepth,
                                    Blindtrienode *oldnode)
{
  Blindtrienode *newbranch, *newleaf;
  Uchar cc1, cc2;

  newbranch = newBlindtrienode(trierep);
  newbranch->startpos = startpos;
  newbranch->rightsibling = oldnode->rightsibling;
  cc1 = getfirstedgechar(trierep,oldnode,currentdepth);
  if (startpos + currentdepth >= getencseqtotallength(trierep->encseq))
  {
    cc2 = (Uchar) SEPARATOR;
  } else
  {
    cc2 = getencodedchar(trierep->encseq, /* Random access */
                         startpos + currentdepth,
                         trierep->readmode);
  }
  newleaf = makenewleaf(trierep,startpos);
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

static Seqpos getlcp(const Encodedsequence *encseq,Readmode readmode,
                     Seqpos start1,Seqpos end1,
                     Seqpos start2,Seqpos end2)
{
  Seqpos i1, i2;
  Uchar cc1;

  for (i1=start1, i2=start2; i1 <= end1 && i2 <= end2; i1++, i2++)
  {
    cc1 = getencodedchar(/*XXX*/ encseq,i1,readmode);
    if (cc1 != getencodedchar(/*XXX*/ encseq,i2,readmode) || ISSPECIAL(cc1))
    {
      break;
    }
  }
  return i1 - start1;
}

static bool hassuccessor(const Blindtrierep *trierep,
                         Nodepair *np,
                         Seqpos prevdepth,
                         const Blindtrienode *node,
                         Uchar cc2,
                         Seqpos startpos2)
{
  Uchar cc1;
  int cmpresult;

  for (np->previous = NULL, np->current = node->firstchild;
       np->current != NULL;
       np->current = np->current->rightsibling)
  {
    cc1 = getfirstedgechar(trierep,np->current,prevdepth);
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

void insertsuffixintoblindtrie(Blindtrierep *trierep,
                               Blindtrienode *node,
                               Seqpos startpos)
{
  if (trierep->root == NULL)
  {
    trierep->root = makeroot(trierep,startpos);
  } else
  {
    Seqpos currentdepth, lcpvalue, totallength;
    Blindtrienode *currentnode, *newleaf, *newbranch, *succ;
    Nodepair np;
    Uchar cc;

    assert(!ISLEAF(node));
    currentnode = node;
    currentdepth = node->depth;
    totallength = getencseqtotallength(trierep->encseq);
    while (true)
    {
      if (startpos + currentdepth >= totallength)
      {
        cc = (Uchar) SEPARATOR;
      } else
      {
        cc = getencodedchar(trierep->encseq, /* Random access */
                            startpos + currentdepth,
                            trierep->readmode);
      }
      assert(currentnode != NULL);
      assert(!ISLEAF(currentnode));
      if (!hassuccessor(trierep,&np,currentdepth,currentnode,cc,startpos))
      {
        newleaf = makenewleaf(trierep,startpos);
        newleaf->rightsibling = np.current;
        if (np.previous == NULL)
        {
          SETFIRSTCHILD(currentnode,newleaf);
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
                          startpos + currentdepth + 1,
                          totallength - 1,
                          succ->startpos + currentdepth + 1,
                          totallength - 1);
        newbranch = makenewbranch(trierep,
                                  startpos,
                                  currentdepth + lcpvalue + 1,
                                  succ);
        if (np.previous == NULL)
        {
          SETFIRSTCHILD(currentnode,newbranch);
        } else
        {
          np.previous->rightsibling = newbranch;
        }
        return;
      }
      lcpvalue = getlcp(trierep->encseq,
                        trierep->readmode,
                        startpos + currentdepth + 1,
                        totallength - 1,
                        succ->startpos + currentdepth + 1,
                        succ->startpos + succ->depth - 1);
      if (currentdepth + lcpvalue + 1 < succ->depth)
      {
        newbranch = makenewbranch(trierep,
                                  startpos,
                                  currentdepth + lcpvalue + 1,
                                  succ);
        if (np.previous == NULL)
        {
          SETFIRSTCHILD(currentnode,newbranch);
        } else
        {
          np.previous->rightsibling = newbranch;
        }
        return;
      }
      currentnode = succ;
      currentdepth = currentnode->depth;
    }
  }
}
