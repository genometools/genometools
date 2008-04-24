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

#define NODE(N)       (trierep->spaceBlindtrienode[N])
#define ISLEAF(N)     (NODE(N).isleaf)
#define ROOTNODE      0
#define BlindtrieNULL (trierep->allocatedBlindtrienode)

typedef struct Blindtrienode
{
  Seqpos depth,
         firstchild;
  unsigned int rightsibling;
  bool isleaf;
  Uchar firstchar;
} Blindtrienode;

DECLAREARRAYSTRUCT(Seqpos);

struct Blindtrierep
{
  const Encodedsequence *encseq;
  Readmode readmode;
  Seqpos totallength;
  Blindtrienode *spaceBlindtrienode;
  unsigned long allocatedBlindtrienode,
                nextfreeBlindtrienode;
  ArraySeqpos stack;
};

typedef struct
{
  Blindtrienode *previous,
                *current;
} Nodepair;

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
  assert(trierep->nextfreeBlindtrienode < trierep->allocatedBlindtrienode);
  return trierep->spaceBlindtrienode + trierep->nextfreeBlindtrienode++;
}

static Blindtrienode *makenewleaf(Blindtrierep *trierep,Seqpos startpos)
{
  Blindtrienode *newleaf;

  newleaf = newBlindtrienode(trierep);
  newleaf->isleaf = true;
  newleaf->depth = 0;
  newleaf->firstchild = startpos;
  newleaf->rightsibling = BlindtrieNULL;
  return newleaf;
}

static Seqpos extractleafnode(const Blindtrierep *trierep,Seqpos head)
{
  assert(!ISLEAF(head));
  do 
  {
    head = NODE(head).firstchild;
  } while(!ISLEAF(head));
  return head;
}

static Seqpos findcompanion(Blindtrierep *trierep,Seqpos startpos)
{
  Uchar cc;
  Seqpos head, currentnodenum;

  trierep->stack.nextfreeSeqpos = 0;
  head = ROOTNODE; 
  while(!ISLEAF(head)) 
  {
    STOREINARRAY (&trierep->stack, Seqpos, 128, head);
    cc = getencodedchar(trierep->encseq, /* Random access */
                        startpos + NODE(head).depth,
                        trierep->readmode);
    if (ISSPECIAL(cc))
    {
      return extractleafnode(trierep,head);
    }
    currentnodenum = NODE(head).firstchild;
    for(;;)
    {
      if(cc == NODE(currentnodenum).firstchar) 
      {              /* found branch corresponding to cc */
        head = NODE(currentnodenum).firstchild;
        break;
      }
      if(cc < NODE(currentnodenum).firstchar)/* no branch corresponding to c */
      {
        return extractleafnode(trierep,head);
      }
      currentnodenum = NODE(currentnodenum).rightsibling;
      if(currentnodenum == BlindtrieNULL) /* no other branches: mismatch */
      {
        return extractleafnode(trierep,head);
      }
    }
  }
  STOREINARRAY (&trierep->stack, Seqpos, 128, head);
  return head;
}
