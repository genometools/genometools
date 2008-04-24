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

#include <assert.h>
#include "libgtcore/arraydef.h"
#include "libgtcore/chardef.h"
#include "format64.h"
#include "intbits-tab.h"
#include "divmodmul.h"
#include "encseq-def.h"

#include "merger-trie.h"

#define ISLEAF(NODE) ((NODE)->firstchild == NULL)
#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
#define SHOWNODERELATIONS(NODE)\
        shownoderelations(__LINE__,#NODE,NODE)
#else
#define SHOWNODERELATIONS(NODE) /* Nothing */
#endif
#else
#define SHOWNODERELATIONS(NODE) /* Nothing */
#endif

#define SETFIRSTCHILD(NODE,VALUE)\
        (NODE)->firstchild = VALUE;\
        if ((VALUE) != NULL)\
        {\
          (VALUE)->parent = NODE;\
        }

#define SETFIRSTCHILDNULL(NODE)\
        (NODE)->firstchild = NULL

typedef struct
{
  Mergertrienode *previous,
           *current;
} Nodepair;

static Uchar getfirstedgechar(const Mergertrierep *trierep,
                              const Mergertrienode *node,
                              Seqpos prevdepth)
{
  Encseqreadinfo *eri = trierep->encseqreadinfo + node->suffixinfo.idx;

  if (ISLEAF(node) &&
      node->suffixinfo.startpos + prevdepth >=
      getencseqtotallength(eri->encseqptr))
  {
    return (Uchar) SEPARATOR;
  }
  return getencodedchar(eri->encseqptr, /* Random access */
                        node->suffixinfo.startpos + prevdepth,
                        eri->readmode);
}

static int comparecharacters(Uchar cc1,unsigned int idx1,
                             Uchar cc2,unsigned int idx2)
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

#ifdef WITHTRIEIDENT
static void showmergertrie2(const Mergertrierep *trierep,
                            const Uchar *characters,
                            unsigned int level,
                            const Mergertrienode *node)
{
  Uchar cc = 0;
  Seqpos pos, endpos;
  Mergertrienode *current;

  for (current = node->firstchild;
       current != NULL;
       current = current->rightsibling)
  {
    printf("%*.*s",(int) (6 * level),(int) (6 * level)," ");
    if (ISLEAF(current))
    {
      endpos = getencseqtotallength(trierep->encseqtable[current->suffixinfo.
                                                         idx]);
    } else
    {
      endpos = current->suffixinfo.startpos + current->depth;
    }
    for (pos = current->suffixinfo.startpos + node->depth;
         pos < endpos; pos++)
    {
      cc = getencodedchar( /* just for testing */
              trierep->enseqreadinfo[current->suffixinfo.idx].
              encseqptr,
              pos,
              trierep->enseqreadinfo[current->suffixinfo.idx].
              readmode);
      if (ISSPECIAL(cc))
      {
        printf("#\n");
        break;
      }
      printf("%c",characters[(int) cc]);
    }
    if (ISLEAF(current))
    {
      if (!ISSPECIAL(cc))
      {
        printf("~\n");
      }
    } else
    {
      printf(" d=" FormatSeqpos ",i=" Formatuint64_t "\n",
            PRINTSeqposcast(current->depth),
            PRINTuint64_tcast(current->suffixinfo.ident));
      showmergertrie2(trierep,characters,level+1,current);
    }
  }
}

void showmergertrie(const Mergertrierep *trierep,
                    const Uchar *characters)
{
  if (trierep->root != NULL)
  {
    showmergertrie2(trierep,characters,0,trierep->root);
  }
}

/*
   Check the following:
   (1) for each branching node there exist at least 2 DONE
   (2) for each branching node the list of successors is strictly ordered
       according to the first character of the edge label DONE
   (3) there are no empty edge labels DONE
   (4) there are \(n+1\) leaves and for each leaf there is exactly one
       incoming edge DONE
*/

static void checkmergertrie2(Mergertrierep *trierep,
                             Mergertrienode *node,
                             Mergertrienode *father,
                             Bitstring *leafused,
                             unsigned int *numberofbitsset)
{
  Mergertrienode *current, *previous;

  if (ISLEAF(node))
  {
    Seqpos start = node->suffixinfo.startpos;
    if (ISIBITSET(leafused,start))
    {
      fprintf(stderr,"leaf " FormatSeqpos " already found\n",
              PRINTSeqposcast(start));
      exit(EXIT_FAILURE); /* programming error */
    }
    SETIBIT(leafused,start);
    (*numberofbitsset)++;
  } else
  {
    assert(node->depth == 0 || node->firstchild->rightsibling != NULL);
    if (father != NULL)
    {
      assert(!ISLEAF(father));
      if (father->depth >= node->depth)
      {
        fprintf(stderr,"father.depth = " FormatSeqpos " >= " FormatSeqpos
                       " = node.depth\n",
                       PRINTSeqposcast(father->depth),
                       PRINTSeqposcast(node->depth));
        exit(EXIT_FAILURE); /* programming error */
      }
    }
    previous = NULL;
    for (current = node->firstchild; current != NULL;
        current = current->rightsibling)
    {
      if (previous != NULL)
      {
        if (comparecharacters(
              getfirstedgechar(trierep,previous,node->depth),
              previous->suffixinfo.idx,
              getfirstedgechar(trierep,current,node->depth),
              current->suffixinfo.idx) >= 0)
        {
          fprintf(stderr,"nodes not correctly ordered\n");
          exit(EXIT_FAILURE); /* programming error */
        }
      }
      checkmergertrie2(trierep,current,node,leafused,numberofbitsset);
      previous = current;
    }
  }
}

void checkmergertrie(Mergertrierep *trierep,unsigned int numberofleaves,
                     unsigned int maxleafnum,Error *err)
{
  error_check(err);
  if (trierep->root != NULL)
  {
    Bitstring *leafused;
    unsigned int numberofbitsset = 0;

    INITBITTAB(leafused,maxleafnum+1);
    checkmergertrie2(trierep,trierep->root,NULL,leafused,&numberofbitsset);
    if (numberofbitsset != numberofleaves)
    {
      fprintf(stderr,"numberofbitsset = %u != %u = numberofleaves\n",
                      numberofbitsset,
                      numberofleaves);
      exit(EXIT_FAILURE); /* programming error */
    }
    FREESPACE(leafused);
  }
}

#ifdef WITHTRIESHOW
static void shownode(const Mergertrienode *node)
{
  if (node == NULL)
  {
    printf("NULL");
  } else
  {
    printf("%s " Formatuint64_t,ISLEAF(node) ? "leaf" : "branch",
                   PRINTuint64_tcast(node->suffixinfo.ident));
  }
}

static void showsimplenoderelations(const Mergertrienode *node)
{
  shownode(node);
  printf(".firstchild=");
  shownode(node->firstchild);
  printf("; ");
  shownode(node);
  printf(".rightsibling=");
  shownode(node->rightsibling);
  printf("\n");
}

static void shownoderelations(int line,char *nodestring,
                              const Mergertrienode *node)
{
  printf("l. %d: %s: ",line,nodestring);
  showsimplenoderelations(node);
}

void showallnoderelations(const Mergertrienode *node)
{
  Mergertrienode *tmp;

  showsimplenoderelations(node);
  for (tmp = node->firstchild; tmp != NULL; tmp = tmp->rightsibling)
  {
    if (tmp->firstchild == NULL)
    {
      showsimplenoderelations(tmp);
    } else
    {
      showallnoderelations(tmp);
    }
  }
}
#endif
#endif

static Mergertrienode *newMergertrienode(Mergertrierep *trierep)
{
#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
  printf("# available trie nodes: %u; ",
          trierep->allocatedMergertrienode - trierep->nextfreeMergertrienode);
  printf("unused trie nodes: %u\n",trierep->nextunused);
#endif
#endif
  if (trierep->nextfreeMergertrienode >= trierep->allocatedMergertrienode)
  {
    assert(trierep->nextunused > 0);
    trierep->nextunused--;
    return trierep->unusedMergertrienodes[trierep->nextunused];
  }
  return trierep->nodetable + trierep->nextfreeMergertrienode++;
}

static Mergertrienode *makenewleaf(Mergertrierep *trierep,
                                   Suffixinfo *suffixinfo)
{
  Mergertrienode *newleaf;

#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
  printf("makenewleaf(" Formatuint64_t ")\n",
         PRINTuint64_tcast(suffixinfo->ident));
#endif
#endif
  newleaf = newMergertrienode(trierep);
  newleaf->suffixinfo = *suffixinfo;
  SETFIRSTCHILDNULL(newleaf);
  newleaf->rightsibling = NULL;
  SHOWNODERELATIONS(newleaf);
  return newleaf;
}

static Mergertrienode *makeroot(Mergertrierep *trierep,Suffixinfo *suffixinfo)
{
  Mergertrienode *root, *newleaf;

#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
  printf("makeroot(" Formatuint64_t ")\n",PRINTuint64_tcast(suffixinfo->ident));
#endif
#endif
  root = newMergertrienode(trierep);
  root->parent = NULL;
  root->suffixinfo = *suffixinfo;
  root->depth = 0;
  root->rightsibling = NULL;
  newleaf = makenewleaf(trierep,suffixinfo);
  SETFIRSTCHILD(root,newleaf);
  SHOWNODERELATIONS(root);
  return root;
}

static void makesuccs(Mergertrienode *newbranch,Mergertrienode *first,
                      Mergertrienode *second)
{
  second->rightsibling = NULL;
  first->rightsibling = second;
  SETFIRSTCHILD(newbranch,first);
  SHOWNODERELATIONS(second);
  SHOWNODERELATIONS(first);
  SHOWNODERELATIONS(newbranch);
}

static Mergertrienode *makenewbranch(Mergertrierep *trierep,
                                     Suffixinfo *suffixinfo,
                                     Seqpos currentdepth,
                                     Mergertrienode *oldnode)
{
  Mergertrienode *newbranch, *newleaf;
  Uchar cc1, cc2;
  Encseqreadinfo *eri = trierep->encseqreadinfo + suffixinfo->idx;

#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
  printf("makenewbranch(ident=" Formatuint64_t ")\n",
          PRINTuint64_tcast(suffixinfo->ident));
#endif
#endif
  newbranch = newMergertrienode(trierep);
  newbranch->suffixinfo = *suffixinfo;
  newbranch->rightsibling = oldnode->rightsibling;
  cc1 = getfirstedgechar(trierep,oldnode,currentdepth);
  if (suffixinfo->startpos + currentdepth >=
      getencseqtotallength(eri->encseqptr))
  {
    cc2 = (Uchar) SEPARATOR;
  } else
  {
    cc2 = getencodedchar(eri->encseqptr, /* Random access */
                         suffixinfo->startpos + currentdepth,
                         eri->readmode);
  }
  newleaf = makenewleaf(trierep,suffixinfo);
  if (comparecharacters(cc1,oldnode->suffixinfo.idx,
                        cc2,suffixinfo->idx) <= 0)
  {
    makesuccs(newbranch,oldnode,newleaf);
  } else
  {
    makesuccs(newbranch,newleaf,oldnode);
  }
  newbranch->depth = currentdepth;
  return newbranch;
}

static Seqpos getlcp(const Encodedsequence *encseq1,Readmode readmode1,
                     Seqpos start1,Seqpos end1,
                     const Encodedsequence *encseq2,Readmode readmode2,
                     Seqpos start2,Seqpos end2)
{
  Seqpos i1, i2;
  Uchar cc1;

  for (i1=start1, i2=start2; i1 <= end1 && i2 <= end2; i1++, i2++)
  {
    cc1 = getencodedchar(/*XXX*/ encseq1,i1,readmode1);
    if (cc1 != getencodedchar(/*XXX*/ encseq2,i2,readmode2) || ISSPECIAL(cc1))
    {
      break;
    }
  }
  return i1 - start1;
}

static bool hassuccessor(const Mergertrierep *trierep,
                         Nodepair *np,
                         Seqpos prevdepth,
                         const Mergertrienode *node,
                         Uchar cc2,
                         unsigned int idx2)
{
  Uchar cc1;
  int cmpresult;

  for (np->previous = NULL, np->current = node->firstchild;
      np->current != NULL;
      np->current = np->current->rightsibling)
  {
    cc1 = getfirstedgechar(trierep,np->current,prevdepth);
    cmpresult = comparecharacters(cc1,np->current->suffixinfo.idx,cc2,idx2);
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

void insertsuffixintomergertrie(Mergertrierep *trierep,
                                Mergertrienode *node,
                                Suffixinfo *suffixinfo)
{
  if (trierep->root == NULL)
  {
    trierep->root = makeroot(trierep,suffixinfo);
  } else
  {
    Seqpos currentdepth, lcpvalue, totallength;
    Mergertrienode *currentnode, *newleaf, *newbranch, *succ;
    Nodepair np;
    Uchar cc;
    Encseqreadinfo *eri = trierep->encseqreadinfo + suffixinfo->idx;

    assert(!ISLEAF(node));
    currentnode = node;
    currentdepth = node->depth;
    totallength = getencseqtotallength(eri->encseqptr);
    while (true)
    {
      if (suffixinfo->startpos + currentdepth >= totallength)
      {
        cc = (Uchar) SEPARATOR;
      } else
      {
        cc = getencodedchar(eri->encseqptr, /* Random access */
                            suffixinfo->startpos + currentdepth,
                            eri->readmode);
      }
      assert(currentnode != NULL);
      assert(!ISLEAF(currentnode));
      if (!hassuccessor(trierep,&np,currentdepth,currentnode,cc,
                        suffixinfo->idx))
      {
        newleaf = makenewleaf(trierep,suffixinfo);
        newleaf->rightsibling = np.current;
        SHOWNODERELATIONS(newleaf);
        if (np.previous == NULL)
        {
          SETFIRSTCHILD(currentnode,newleaf);
          SHOWNODERELATIONS(currentnode);
        } else
        {
          np.previous->rightsibling = newleaf;
          SHOWNODERELATIONS(np.previous);
        }
        return;
      }
      succ = np.current;
      if (ISLEAF(succ))
      {
        lcpvalue = getlcp(eri->encseqptr,
                          eri->readmode,
                          suffixinfo->startpos + currentdepth + 1,
                          getencseqtotallength(eri->encseqptr) - 1,
                          trierep->encseqreadinfo[succ->suffixinfo.idx].
                                encseqptr,
                          trierep->encseqreadinfo[succ->suffixinfo.idx].
                                readmode,
                          succ->suffixinfo.startpos + currentdepth + 1,
                          getencseqtotallength(
                              trierep->encseqreadinfo[succ->suffixinfo.idx].
                                        encseqptr) - 1);
        newbranch = makenewbranch(trierep,
                                  suffixinfo,
                                  currentdepth + lcpvalue + 1,
                                  succ);
        if (np.previous == NULL)
        {
          SETFIRSTCHILD(currentnode,newbranch);
          SHOWNODERELATIONS(currentnode);
        } else
        {
          np.previous->rightsibling = newbranch;
          SHOWNODERELATIONS(np.previous);
        }
        return;
      }
      lcpvalue = getlcp(eri->encseqptr,
                        eri->readmode,
                        suffixinfo->startpos + currentdepth + 1,
                        getencseqtotallength(eri->encseqptr) - 1,
                        trierep->encseqreadinfo[succ->suffixinfo.idx].encseqptr,
                        trierep->encseqreadinfo[succ->suffixinfo.idx].readmode,
                        succ->suffixinfo.startpos + currentdepth + 1,
                        succ->suffixinfo.startpos + succ->depth - 1);
      if (currentdepth + lcpvalue + 1 < succ->depth)
      {
        newbranch = makenewbranch(trierep,
                                  suffixinfo,
                                  currentdepth + lcpvalue + 1,
                                  succ);
        if (np.previous == NULL)
        {
          SETFIRSTCHILD(currentnode,newbranch);
          SHOWNODERELATIONS(currentnode);
        } else
        {
          np.previous->rightsibling = newbranch;
          SHOWNODERELATIONS(np.previous);
        }
        return;
      }
      currentnode = succ;
      currentdepth = currentnode->depth;
    }
  }
}

Mergertrienode *findsmallestnodeintrie(const Mergertrierep *trierep)
{
  Mergertrienode *node;

  assert(trierep->root != NULL);
  for (node = trierep->root; node->firstchild != NULL; node = node->firstchild)
    /* Nothing */ ;
  return node;
}

void deletesmallestpath(Mergertrienode *smallest,Mergertrierep *trierep)
{
  Mergertrienode *father, *son;

  for (son = smallest; son->parent != NULL; son = son->parent)
  {
    father = son->parent;
    if (son->firstchild == NULL)
    {
      SETFIRSTCHILD(father,son->rightsibling);
      SHOWNODERELATIONS(father);
      son->rightsibling = NULL;
    } else
    {
      if (son->firstchild->rightsibling != NULL)
      {
        break;
      }
      son->firstchild->rightsibling = father->firstchild->rightsibling;
      SETFIRSTCHILD(father,son->firstchild);
      SHOWNODERELATIONS(father);
      son->rightsibling = NULL;
      SETFIRSTCHILDNULL(son);
    }
#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
    printf("delete %s " Formatuint64_t "\n",
             ISLEAF(son) ? "leaf" : "branch",
             PRINTuint64_tcast(son->suffixinfo.ident));
#endif
#endif
    trierep->unusedMergertrienodes[trierep->nextunused++] = son;
  }
  if (trierep->root->firstchild == NULL)
  {
    trierep->unusedMergertrienodes[trierep->nextunused++] = trierep->root;
    trierep->root = NULL;
  }
}

void initmergertrienodetable(Mergertrierep *trierep,Seqpos numofsuffixes,
                             unsigned int numofindexes)
{
  trierep->numofindexes = numofindexes;
  trierep->allocatedMergertrienode
    = (unsigned int) MULT2(numofsuffixes + 1) + 1;
  ALLOCASSIGNSPACE(trierep->nodetable,NULL,Mergertrienode,
                   trierep->allocatedMergertrienode);
  trierep->nextfreeMergertrienode = 0;
  trierep->root = NULL;
  trierep->nextunused = 0;
  ALLOCASSIGNSPACE(trierep->unusedMergertrienodes,NULL,Mergertrienodeptr,
                   trierep->allocatedMergertrienode);
}

void freemergertrierep(Mergertrierep *trierep)
{
  FREESPACE(trierep->nodetable);
  FREESPACE(trierep->unusedMergertrienodes);
  FREESPACE(trierep->encseqreadinfo);
}
