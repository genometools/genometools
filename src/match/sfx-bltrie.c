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

#include "core/unused_api.h"
#include "core/arraydef.h"
#include "core/divmodmul.h"
#include "core/minmax.h"
#include "core/encseq.h"

#include "lcpoverflow.h"
#include "sfx-bltrie.h"
#include "sfx-suffixgetset.h"

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
  const GtEncseq *encseq;
  GtEncseqReader *esr1, *esr2;
  GtReadmode readmode;
  unsigned long totallength,
                offset,
                maxdepth,
                maxdepthminusoffset,
                allocatedBlindtrienode,
                nextfreeBlindtrienode,
                subbucketleft;
  Nodeptr root;
  bool cmpcharbychar;
  Blindtrienode *spaceBlindtrienode;
  GtArrayNodeptr stack;
  GtSuffixsortspace *sssp;
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
    firstchar = gt_encseq_get_encoded_char(blindtrie->encseq,
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
      newchar = gt_encseq_get_encoded_char(blindtrie->encseq,
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

  gt_encseq_reader_reinit_with_readmode(blindtrie->esr1,blindtrie->encseq,
                                        blindtrie->readmode,leafpos);
  gt_encseq_reader_reinit_with_readmode(blindtrie->esr2,blindtrie->encseq,
                                        blindtrie->readmode,currentstartpos);
  for (lcp = 0; /* Nothing */; lcp++)
  {
    if (isleftofboundary(leafpos,lcp,blindtrie))
    {
      cc1 = gt_encseq_reader_next_encoded_char(blindtrie->esr1);
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
      cc2 = gt_encseq_reader_next_encoded_char(blindtrie->esr2);
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

  (void) gt_encseq_compare_viatwobitencoding(&commonunits,
                                             blindtrie->encseq,
                                             blindtrie->readmode,
                                             blindtrie->esr1,
                                             blindtrie->esr2,
                                             leafpos,
                                             currentstartpos,
                                             0,
                                             (blindtrie->maxdepth == 0)
                                                ? 0
                                                : blindtrie->maxdepthminusoffset
                                             );
  if (isleftofboundary(leafpos,commonunits.finaldepth,blindtrie) &&
      !commonunits.leftspecial)
  {
    *mm_oldsuffix = gt_encseq_extract_encoded_char(blindtrie->encseq,
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
    *mm_newsuffix = gt_encseq_extract_encoded_char(blindtrie->encseq,
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

static unsigned long enumeratetrieleaves (Blindtrie *blindtrie,
                                          unsigned long subbucketleft,
                                          unsigned long *lcpsubtab,
                                          unsigned long *numoflargelcpvalues,
                                          void *voiddcov,
                                          Dc_processunsortedrange
                                            dc_processunsortedrange)
{
  bool readyforpop = false, currentnodeisleaf;
  Nodeptr currentnode, siblval, lcpnode = blindtrie->root;
  unsigned long nextfree = 0, equalsrangewidth = 0,
                bucketleftidxplussubbucketleft;

  blindtrie->stack.nextfreeNodeptr = 0;
  GT_STOREINARRAY (&blindtrie->stack, Nodeptr, 128, blindtrie->root);
  SETCURRENT(blindtrie->root->either.firstchild);
  gt_assert(blindtrie->maxdepth == 0 || dc_processunsortedrange != NULL);
  bucketleftidxplussubbucketleft
    = gt_suffixsortspace_bucketleftidx_get(blindtrie->sssp) + subbucketleft;
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
              dc_processunsortedrange(
                               voiddcov,
                               bucketleftidxplussubbucketleft
                                 + nextfree - 1 - equalsrangewidth,
                               equalsrangewidth + 1,
                               blindtrie->maxdepth);
              equalsrangewidth = 0;
            }
          }
        }
      }
      gt_suffixsortspace_set(blindtrie->sssp,subbucketleft,nextfree,
                             currentnode->either.nodestartpos
                               - blindtrie->offset);
      nextfree++;
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
                            bucketleftidxplussubbucketleft
                              + nextfree - 1 - equalsrangewidth,
                            equalsrangewidth + 1,
                            blindtrie->maxdepth);
    equalsrangewidth = 0;
  }
  return nextfree;
}

Blindtrie *gt_blindtrie_new(GtSuffixsortspace *suffixsortspace,
                            unsigned long numofsuffixes,
                            const GtEncseq *encseq,
                            bool cmpcharbychar,
                            GtEncseqReader *esr1,
                            GtEncseqReader *esr2,
                            GtReadmode readmode)
{
  Blindtrie *blindtrie;

  blindtrie = gt_malloc(sizeof (*blindtrie));
  blindtrie->allocatedBlindtrienode = GT_MULT2(numofsuffixes + 1) + 1;
  blindtrie->spaceBlindtrienode
    = gt_malloc(sizeof (*blindtrie->spaceBlindtrienode) *
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
  blindtrie->sssp = suffixsortspace;
  blindtrie->root = NULL;
  blindtrie->esr1 = esr1;
  blindtrie->esr2 = esr2;
  blindtrie->totallength = gt_encseq_total_length(encseq);
  blindtrie->cmpcharbychar = cmpcharbychar;
  GT_INITARRAY (&blindtrie->stack, Nodeptr);
  return blindtrie;
}

void gt_blindtrie_delete(Blindtrie *blindtrie)
{
  if (blindtrie == NULL)
  {
    return;
  }
  gt_free(blindtrie->spaceBlindtrienode);
  GT_FREEARRAY(&blindtrie->stack, Nodeptr);
  gt_free(blindtrie);
}

#ifdef SKDEBUG
static void gt_blindtrie_showleaf(const Blindtrie *blindtrie,unsigned int level,
                                  Nodeptr current)
{
  printf("%*.*s",(int) (6 * level),(int) (6 * level)," ");
  gt_assert(current != NULL);
  printf("Leaf(add=%lu,firstchar=%u,startpos=%lu,rightsibling=%lu)\n",
         NODENUM(current),
         (unsigned int) current->firstchar,
         current->either.nodestartpos,
         NODENUM(current->rightsibling));
}

static void gt_blindtrie_showintern(const Blindtrie *blindtrie,
                                    unsigned int level,
                                    Nodeptr current)
{
  printf("%*.*s",(int) (6 * level),(int) (6 * level)," ");
  gt_assert(current != NULL);
  printf("Intern(add=%lu,firstchar=%u,depth=%lu"
         ",firstchild=%lu,rightsibling=%lu)\n",
          NODENUM(current),
          (unsigned int) current->firstchar,
          current->depth,
          NODENUM(current->either.firstchild),
          NODENUM(current->rightsibling));
}

static void gt_blindtrie_showrecursive(const Blindtrie *blindtrie,
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
  gt_blindtrie_showrecursive(blindtrie,0,blindtrie->root);
}
#endif

#ifndef NDEBUG

static void checksorting(const Blindtrie *blindtrie,
                         unsigned long subbucketleft,
                         unsigned long numberofsuffixes,
                         bool ascending)
{
  unsigned long idx, pos1, pos2;

  gt_assert(numberofsuffixes > 1UL);
  for (idx = 0; idx < numberofsuffixes - 1; idx++)
  {
    pos1 = gt_suffixsortspace_get(blindtrie->sssp,subbucketleft,idx);
    pos2 = gt_suffixsortspace_get(blindtrie->sssp,subbucketleft,idx+1);
    if ((ascending && pos1 >= pos2) ||
        (!ascending && pos1 <= pos2))
    {
      fprintf(stderr,"not %s: ",ascending ? "ascending" : "descending");
      fprintf(stderr,"subbucket[%lu]=%lu vs %lu=subbucket[%lu]\n",
                      idx,pos1,pos2,idx+1);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

#endif

static void inplace_reverseSuffixptr(const Blindtrie *blindtrie,
                                     unsigned long subbucketleft,
                                     unsigned long len)
{
  unsigned long tmp, i, j;

  gt_assert(len > 0);
  for (i = 0, j = len - 1; i < j; i++, j--)
  {
    tmp = gt_suffixsortspace_get(blindtrie->sssp,subbucketleft,i);
    gt_suffixsortspace_set(blindtrie->sssp,subbucketleft,i,
                           gt_suffixsortspace_get(blindtrie->sssp,
                                                  subbucketleft,j));
    gt_suffixsortspace_set(blindtrie->sssp,subbucketleft,j,tmp);
  }
}

#ifdef  QSORTNAME
#undef  QSORTNAME
#endif

#define QSORTNAME(NAME) bltrie_##NAME

#ifdef QSORT_ARRAY_DECLARE
#undef QSORT_ARRAY_DECLARE
#endif

#define QSORT_ARRAY_DECLARE\
        Blindtrie *blindtrie = (Blindtrie *) data

#ifdef QSORT_ARRAY_GET
#undef QSORT_ARRAY_GET
#endif

#define QSORT_ARRAY_GET(ARR,RELIDX)\
        gt_suffixsortspace_get(blindtrie->sssp,blindtrie->subbucketleft,RELIDX)

#ifdef QSORT_ARRAY_SET
#undef QSORT_ARRAY_SET
#endif

#define QSORT_ARRAY_SET(ARR,RELIDX,VALUE)\
        gt_suffixsortspace_set(blindtrie->sssp,blindtrie->subbucketleft,RELIDX,\
                               VALUE)

static int QSORTNAME(qsortcmparr) (
                  GT_UNUSED const void *subbucket,
                  unsigned long a,
                  unsigned long b,
                  const void *data)
{
  const Blindtrie *blindtrie = (const Blindtrie *) data;
  unsigned long start1, start2;

  start1 = QSORT_ARRAY_GET(NULL,a);
  start2 = QSORT_ARRAY_GET(NULL,b);
  gt_assert(start1 != start2);
  if (start1 < start2)
  {
    return -1;
  }
  return 1;
}

typedef void * QSORTNAME(Sorttype);

#include "qsort-array.gen"

unsigned long gt_blindtrie_suffixsort(
                            Blindtrie *blindtrie,
                            unsigned long subbucketleft,
                            unsigned long *lcpsubtab,
                            unsigned long numberofsuffixes,
                            unsigned long offset,
                            unsigned long maxdepth,
                            Ordertype ordertype,
                            void *voiddcov,
                            Dc_processunsortedrange dc_processunsortedrange)
{
  unsigned long idx, stackidx;
  Nodeptr leafinsubtree, currentnode;
  unsigned long pos, lcp, numoflargelcpvalues = 0;
  GtUchar mm_oldsuffix, mm_newsuffix;

  if (ordertype == Noorder)
  {
    blindtrie->subbucketleft = subbucketleft;
    QSORTNAME(gt_inlinedarr_qsort_r) (NULL,numberofsuffixes,
                                      (void *) blindtrie);
  } else
  {
    if (ordertype == Descending)
    {
#ifndef NDEBUG
      checksorting(blindtrie,subbucketleft,numberofsuffixes,false);
#endif
      inplace_reverseSuffixptr(blindtrie,subbucketleft,numberofsuffixes);
    } else
    {
#ifndef NDEBUG
      checksorting(blindtrie,subbucketleft,numberofsuffixes,true);
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
  pos = gt_suffixsortspace_get(blindtrie->sssp,subbucketleft,0) + offset;
  blindtrie->root = makeroot(blindtrie,pos);
#ifdef SKDEBUG
  printf("insert suffixes at offset %lu:\n",offset);
  for (idx=0; idx < numberofsuffixes; idx++)
  {
    printf("%lu ",
           gt_suffixsortspace_get(blindtrie->sssp,subbucketleft,idx) + offset);
  }
  printf("\nstep 0\n");
  gt_blindtrie_show(blindtrie);
#endif
  for (idx=1UL; idx < numberofsuffixes; idx++)
  {
    pos = gt_suffixsortspace_get(blindtrie->sssp,subbucketleft,idx) + offset;
    if (isleftofboundary(pos,0,blindtrie))
    {
      leafinsubtree = findcompanion(blindtrie,pos);
      gt_assert(ISLEAF(leafinsubtree));
      lcp = (blindtrie->cmpcharbychar ? cmpcharbychargetlcp : fastgetlcp)
                               (&mm_oldsuffix,
                                &mm_newsuffix,
                                blindtrie,
                                leafinsubtree->either.nodestartpos,
                                pos);
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
                                pos);
#ifdef SKDEBUG
      printf("step %lu\n",idx);
      gt_blindtrie_show(blindtrie);
#endif
    } else
    {
      break;
    }
  }
  (void) enumeratetrieleaves (blindtrie, subbucketleft,lcpsubtab,
                              &numoflargelcpvalues,
                              voiddcov,dc_processunsortedrange);
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
