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
#define SETLEAF(NODE,VAL) /* Nothing: used for documentation purpose */

typedef struct Blindtrienode
{
  unsigned long depth;
  struct Blindtrienode *rightsibling;
  union
  {
    struct Blindtrienode *firstchild;
    struct
    {
      unsigned long nodestartpos,
                    nodestoppos;
    } leafinfo;
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

static bool blindtrie_isleftofboundary(const Blindtrie *blindtrie,
                                       unsigned long currentstartpos,
                                       unsigned long add)
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

static Nodeptr blindtrie_newnode(Blindtrie *blindtrie)
{
  gt_assert(blindtrie->nextfreeBlindtrienode <
            blindtrie->allocatedBlindtrienode);
  return blindtrie->spaceBlindtrienode + blindtrie->nextfreeBlindtrienode++;
}

static Blindtrienode *blindtrie_newleaf(Blindtrie *blindtrie,
                                        unsigned long currentstartpos,
                                        unsigned long currentstoppos,
                                        GtUchar firstchar,
                                        struct Blindtrienode *rightsibling)
{
  Blindtrienode *newleaf;

  newleaf = blindtrie_newnode(blindtrie);
  newleaf->depth = 0;
  newleaf->either.leafinfo.nodestartpos = currentstartpos;
  newleaf->either.leafinfo.nodestoppos = currentstoppos;
  SETLEAF(newleaf,true);
  newleaf->firstchar = firstchar;
  newleaf->rightsibling = rightsibling;
  return newleaf;
}

static Nodeptr blindtrie_makeroot(Blindtrie *blindtrie,
                                  unsigned long currentstartpos)
{
  Blindtrienode *root;
  GtUchar firstchar;

  root = blindtrie_newnode(blindtrie);
  root->depth = 0;
  root->firstchar = 0; /* undefined */
  root->rightsibling = NULL;
  SETLEAF(root,false);
  if (blindtrie_isleftofboundary(blindtrie,currentstartpos,0))
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
  root->either.firstchild = blindtrie_newleaf(blindtrie,currentstartpos,
                                              /* DO first fwd/rev getnext */ 0,
                                              firstchar,NULL);
  return root;
}

static inline Nodeptr blindtrie_extractleafnode(const Blindtrie *blindtrie,
                                                Nodeptr head)
{
  gt_assert(ISNOTLEAF(head));
  do
  {
    head = head->either.firstchild;
  } while (ISNOTLEAF(head));
  return head;
}

static inline int blindtrie_comparecharacters(GtUchar oldchar,GtUchar newchar)
{
  return (oldchar > newchar)
           ? 1
           : ((oldchar < newchar || ISSPECIAL(oldchar))
                  ? -1
                  : 0);
}

static Nodeptr blindtrie_findsucc(Nodeptr node,GtUchar newchar)
{
  int retval;

  for (;;)
  {
    retval = blindtrie_comparecharacters(node->firstchar,newchar);
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

static Nodeptr blindtrie_findcompanion(Blindtrie *blindtrie,
                                       unsigned long currentstartpos)
{
  GtUchar newchar;
  Nodeptr head, succ;

  blindtrie->stack.nextfreeNodeptr = 0;
  head = blindtrie->root;
  while (ISNOTLEAF(head))
  {
    GT_STOREINARRAY (&blindtrie->stack, Nodeptr, 128, head);
    if (blindtrie_isleftofboundary(blindtrie,currentstartpos,head->depth))
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
      return blindtrie_extractleafnode(blindtrie,head);
    }
    succ = blindtrie_findsucc(head->either.firstchild,newchar);
    if (succ == NULL)
    {
      return blindtrie_extractleafnode(blindtrie,head);
    }
    head = succ;
  }
  GT_STOREINARRAY (&blindtrie->stack, Nodeptr, 128, head);
  return head;
}

static void blindtrie_insertsuffix(Blindtrie *blindtrie,
                                   Nodeptr oldnode,
                                   GtUchar mm_oldsuffix,
                                   unsigned long lcp,
                                   GtUchar mm_newsuffix,
                                   unsigned long currentstartpos,
                                   unsigned long currentstoppos)
{
  Nodeptr newleaf, newnode, previous, current;

  gt_assert(ISSPECIAL(mm_oldsuffix) || ISSPECIAL(mm_newsuffix) ||
            mm_oldsuffix != mm_newsuffix || ISLEAF(oldnode) ||
            oldnode->depth == lcp);

  /* insert a new node before node oldnode if necessary */
  if (oldnode->depth != lcp)
  {
    newnode = blindtrie_newnode(blindtrie);
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
  previous = NULL;
  current = oldnode->either.firstchild;
  while (current != NULL &&
         blindtrie_comparecharacters(current->firstchar,mm_newsuffix) < 0)
  {
    previous = current;
    current = current->rightsibling;
  }
  /* insert new leaf with current suffix */
  /* search S[lcp] among the offsprings */
  newleaf = blindtrie_newleaf(blindtrie,currentstartpos,
                              currentstoppos,mm_newsuffix,current);
  if (previous != NULL)
  {
    previous->rightsibling = newleaf;
  } else
  {
    oldnode->either.firstchild = newleaf;
  }
}

static unsigned long blindtrie_cmpcharbychargetlcp(
                                 GtUchar *mm_oldsuffix,
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
    if (blindtrie_isleftofboundary(blindtrie,leafpos,lcp))
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
    if (blindtrie_isleftofboundary(blindtrie,currentstartpos,lcp))
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
    if (blindtrie_comparecharacters(cc1,cc2) != 0)
    {
      *mm_oldsuffix = cc1;
      *mm_newsuffix = cc2;
      break;
    }
  }
  gt_assert(blindtrie->maxdepth == 0 || lcp <= blindtrie->maxdepthminusoffset);
  return lcp;
}

static unsigned long blindtrie_twobitencodinggetlcp(
                                          GtUchar *mm_oldsuffix,
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
                                               : blindtrie->maxdepthminusoffset,
                                             NULL,
                                             NULL);
  if (blindtrie_isleftofboundary(blindtrie,leafpos,commonunits.finaldepth) &&
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
  if (blindtrie_isleftofboundary(blindtrie,currentstartpos,
                                 commonunits.finaldepth) &&
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

static unsigned long blindtrie_enumeratetrieleaves (
                           Blindtrie *blindtrie,
                           unsigned long subbucketleft,
                           unsigned long *lcpsubtab,
                           unsigned long *numoflargelcpvalues,
                           void *voiddcov,
                           Dc_processunsortedrange dc_processunsortedrange)
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
                             currentnode->either.leafinfo.nodestartpos
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
         current->either.leafinfo.nodestartpos,
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
  unsigned long idx, stackidx, currentstartpos, lcp, numoflargelcpvalues = 0;
  Nodeptr leafinsubtree, currentnode;
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
  currentstartpos = gt_suffixsortspace_get(blindtrie->sssp,subbucketleft,0)
                    + offset;
  blindtrie->root = blindtrie_makeroot(blindtrie,currentstartpos);
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
    currentstartpos = gt_suffixsortspace_get(blindtrie->sssp,subbucketleft,idx)
                      + offset;
    if (blindtrie_isleftofboundary(blindtrie,currentstartpos,0))
    {
      leafinsubtree = blindtrie_findcompanion(blindtrie,currentstartpos);
      gt_assert(ISLEAF(leafinsubtree));
      if (blindtrie->cmpcharbychar)
      {
        lcp = blindtrie_cmpcharbychargetlcp
                               (&mm_oldsuffix,
                                &mm_newsuffix,
                                blindtrie,
                                leafinsubtree->either.leafinfo.nodestartpos,
                                currentstartpos);
      } else
      {
        lcp = blindtrie_twobitencodinggetlcp
                               (&mm_oldsuffix,
                                &mm_newsuffix,
                                blindtrie,
                                leafinsubtree->either.leafinfo.nodestartpos,
                                currentstartpos);
      }
      currentnode = blindtrie->root;
      for (stackidx=0;stackidx<blindtrie->stack.nextfreeNodeptr;stackidx++)
      {
        currentnode = blindtrie->stack.spaceNodeptr[stackidx];
        if (ISLEAF(currentnode) || currentnode->depth >= lcp)
        {
          break;
        }
      }
      blindtrie_insertsuffix(blindtrie,
                             currentnode,
                             mm_oldsuffix,
                             lcp,
                             mm_newsuffix,
                             currentstartpos,
                             /* XXX compute stoppos for pos */ 0);
#ifdef SKDEBUG
      printf("step %lu\n",idx);
      gt_blindtrie_show(blindtrie);
#endif
    } else
    {
      break;
    }
  }
  (void) blindtrie_enumeratetrieleaves (blindtrie, subbucketleft,lcpsubtab,
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
