/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "format64.h"
#include "chardef.h"
#include "intbits-tab.h"
#include "divmodmul.h"
#include "arraydef.h"
#include "encseq-def.h"

#include "trieins-def.h"

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
        if((VALUE) != NULL)\
        {\
          (VALUE)->parent = NODE;\
        }

#define SETFIRSTCHILDNULL(NODE)\
        (NODE)->firstchild = NULL

typedef struct
{
  Trienode *previous,
           *current;
} Nodepair;

static Uchar getfirstedgechar(const Trierep *trierep,
                              const Trienode *node,
                              Seqpos prevdepth)
{
  Encseqreadinfo *eri = trierep->encseqreadinfo + node->suffixinfo.idx;

  if(ISLEAF(node) &&
     node->suffixinfo.startpos + prevdepth >=
     getencseqtotallength(eri->encseqptr))
  {
    return (Uchar) SEPARATOR;
  }
  return getencodedchar(eri->encseqptr,
                        node->suffixinfo.startpos + prevdepth,
                        eri->readmode);
}

static int comparecharacters(Uchar cc1,Seqpos idx1,Uchar cc2,Seqpos idx2)
{
  if(ISSPECIAL(cc1))
  {
    if(ISSPECIAL(cc2))
    {
      if(idx1 <= idx2)
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
    if(ISSPECIAL(cc2))
    {
      return -1;  /* cc1 < cc2 */
    } else
    {
      if(cc1 < cc2)
      {
        return -1;  /* cc1 < cc2 */
      } else
      {
        if(cc1 > cc2)
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
static void showtrie2(const Trierep *trierep,
                      const Uchar *characters,
                      uint32_t level,
                      const Trienode *node)
{
  Uchar cc = 0;
  Seqpos pos, endpos;
  Trienode *current;

  for(current = node->firstchild;
      current != NULL;
      current = current->rightsibling)
  {
    printf("%*.*s",(int) (6 * level),(int) (6 * level)," ");
    if(ISLEAF(current))
    {
      endpos = getencseqtotallength(
                  trierep->encseqtable[current->suffixinfo.idx]);
    } else
    {
      endpos = current->suffixinfo.startpos + current->depth;
    }
    for(pos = current->suffixinfo.startpos + node->depth;
        pos < endpos; pos++)
    {
      cc = getencodedchar(trierep->enseqreadinfo[current->suffixinfo.idx].
                          encseqptr,
                          pos,
                          trierep->enseqreadinfo[current->suffixinfo.idx].
                          readmode);
      if(ISSPECIAL(cc))
      {
        printf("#\n");
        break;
      }
      printf("%c",characters[(int) cc]);
    }
    if(ISLEAF(current))
    {
      if(!ISSPECIAL(cc))
      {
        printf("~\n");
      }
    } else
    {
      printf(" d=" FormatSeqpos ",i=" Formatuint64_t "\n",
            PRINTSeqposcast(current->depth),
            PRINTuint64_tcast(current->suffixinfo.ident));
      showtrie2(trierep,characters,level+1,current);
    }
  }
}

void showtrie(const Trierep *trierep,
              const Uchar *characters)
{
  if(trierep->root != NULL)
  {
    showtrie2(trierep,characters,0,trierep->root);
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

static void checktrie2(Trierep *trierep,
                       Trienode *node,
                       Trienode *father,
                       Bitstring *leafused,
                       uint32_t *numberofbitsset)
{
  Trienode *current, *previous;

  if(ISLEAF(node))
  {
    Seqpos start = node->suffixinfo.startpos;
    if(ISIBITSET(leafused,start))
    {
      fprintf(stderr,"leaf " FormatSeqpos " already found\n",
              PRINTSeqposcast(start));
      exit(EXIT_FAILURE);
    }
    SETIBIT(leafused,start);
    (*numberofbitsset)++;
  } else
  {
    if(node->depth > 0 && node->firstchild->rightsibling == NULL)
    {
      fprintf(stderr,"Node has less than two successors\n");
      exit(EXIT_FAILURE);
    }
    if(father != NULL)
    {
      if(ISLEAF(father))
      {
        fprintf(stderr,"father of branching node is a leaf\n");
        exit(EXIT_FAILURE);
      }
      if(father->depth >= node->depth)
      {
        fprintf(stderr,"father.depth = " FormatSeqpos " >= " FormatSeqpos
                       " = node.depth\n",
                       PRINTSeqposcast(father->depth),
                       PRINTSeqposcast(node->depth));
        exit(EXIT_FAILURE);
      }
    }
    previous = NULL;
    for(current = node->firstchild; current != NULL;
        current = current->rightsibling)
    {
      if(previous != NULL)
      {
        if(comparecharacters(
              getfirstedgechar(trierep,previous,node->depth),
              previous->suffixinfo.idx,
              getfirstedgechar(trierep,current,node->depth),
              current->suffixinfo.idx) >= 0)
        {
          fprintf(stderr,"nodes not correctly ordered\n");
          exit(EXIT_FAILURE);
        }
      }
      checktrie2(trierep,current,node,leafused,numberofbitsset);
      previous = current;
    }
  }
}

void checktrie(Trierep *trierep,uint32_t numberofleaves,
               uint32_t maxleafnum,Env *env)
{
  env_error_check(env);
  if(trierep->root != NULL)
  {
    Bitstring *leafused;
    uint32_t numberofbitsset = 0;

    INITBITTAB(leafused,maxleafnum+1);
    checktrie2(trierep,trierep->root,NULL,leafused,&numberofbitsset);
    if(numberofbitsset != numberofleaves)
    {
      fprintf(stderr,"numberofbitsset = %u != %u = numberofleaves\n",
                      (unsigned int) numberofbitsset,
                      (unsigned int) numberofleaves);
      exit(EXIT_FAILURE);
    }
    FREESPACE(leafused);
  }
}

#ifdef WITHTRIESHOW
static void shownode(const Trienode *node)
{
  if(node == NULL)
  {
    printf("NULL");
  } else
  {
    printf("%s " Formatuint64_t,ISLEAF(node) ? "leaf" : "branch",
                   PRINTuint64_tcast(node->suffixinfo.ident));
  }
}

static void showsimplenoderelations(const Trienode *node)
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

static void shownoderelations(int line,char *nodestring,const Trienode *node)
{
  printf("l. %d: %s: ",line,nodestring);
  showsimplenoderelations(node);
}

void showallnoderelations(const Trienode *node)
{
  Trienode *tmp;

  showsimplenoderelations(node);
  for(tmp = node->firstchild; tmp != NULL; tmp = tmp->rightsibling)
  {
    if(tmp->firstchild == NULL)
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

static Trienode *newTrienode(Trierep *trierep)
{
#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
  printf("# available trie nodes: %u; ",
          trierep->allocatedTrienode - trierep->nextfreeTrienode);
  printf("unused trie nodes: %u\n",trierep->nextunused);
#endif
#endif
  if(trierep->nextfreeTrienode >= trierep->allocatedTrienode)
  {
    if(trierep->nextunused == 0)
    {
      fprintf(stderr,"not enough nodes have been allocated\n");
      exit(EXIT_FAILURE);
    }
    trierep->nextunused--;
    return trierep->unusedTrienodes[trierep->nextunused];
  }
  return trierep->nodetable + trierep->nextfreeTrienode++;
}

static Trienode *makenewleaf(Trierep *trierep,Suffixinfo *suffixinfo)
{
  Trienode *newleaf;

#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
  printf("makenewleaf(" Formatuint64_t ")\n",
         PRINTuint64_tcast(suffixinfo->ident));
#endif
#endif
  newleaf = newTrienode(trierep);
  newleaf->suffixinfo = *suffixinfo;
  SETFIRSTCHILDNULL(newleaf);
  newleaf->rightsibling = NULL;
  SHOWNODERELATIONS(newleaf);
  return newleaf;
}

static Trienode *makeroot(Trierep *trierep,Suffixinfo *suffixinfo)
{
  Trienode *root, *newleaf;

#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
  printf("makeroot(" Formatuint64_t ")\n",PRINTuint64_tcast(suffixinfo->ident));
#endif
#endif
  root = newTrienode(trierep);
  root->parent = NULL;
  root->suffixinfo = *suffixinfo;
  root->depth = 0;
  root->rightsibling = NULL;
  newleaf = makenewleaf(trierep,suffixinfo);
  SETFIRSTCHILD(root,newleaf);
  SHOWNODERELATIONS(root);
  return root;
}

static void makesuccs(Trienode *newbranch,Trienode *first,Trienode *second)
{
  second->rightsibling = NULL;
  first->rightsibling = second;
  SETFIRSTCHILD(newbranch,first);
  SHOWNODERELATIONS(second);
  SHOWNODERELATIONS(first);
  SHOWNODERELATIONS(newbranch);
}

static Trienode *makenewbranch(Trierep *trierep,
                               Suffixinfo *suffixinfo,
                               Seqpos currentdepth,
                               Trienode *oldnode)
{
  Trienode *newbranch, *newleaf;
  Uchar cc1, cc2;
  Encseqreadinfo *eri = trierep->encseqreadinfo + suffixinfo->idx;

#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
  printf("makenewbranch(ident=" Formatuint64_t ")\n",
          PRINTuint64_tcast(suffixinfo->ident));
#endif
#endif
  newbranch = newTrienode(trierep);
  newbranch->suffixinfo = *suffixinfo;
  newbranch->rightsibling = oldnode->rightsibling;
  cc1 = getfirstedgechar(trierep,oldnode,currentdepth);
  if(suffixinfo->startpos + currentdepth >=
     getencseqtotallength(eri->encseqptr))
  {
    cc2 = (Uchar) SEPARATOR;
  } else
  {
    cc2 = getencodedchar(eri->encseqptr,
                         suffixinfo->startpos + currentdepth,
                         eri->readmode);
  }
  newleaf = makenewleaf(trierep,suffixinfo);
  if(comparecharacters(cc1,oldnode->suffixinfo.idx,
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

  for(i1=start1, i2=start2; i1 <= end1 && i2 <= end2; i1++, i2++)
  {
    cc1 = getencodedchar(encseq1,i1,readmode1);
    if(cc1 != getencodedchar(encseq2,i2,readmode2) || ISSPECIAL(cc1))
    {
      break;
    }
  }
  return i1 - start1;
}

static bool hassuccessor(const Trierep *trierep,
                         Nodepair *np,
                         Seqpos prevdepth,
                         const Trienode *node,
                         Uchar cc2,
                         Seqpos idx2)
{
  Uchar cc1;
  int cmpresult;

  for(np->previous = NULL, np->current = node->firstchild;
      np->current != NULL;
      np->current = np->current->rightsibling)
  {
    cc1 = getfirstedgechar(trierep,np->current,prevdepth);
    cmpresult = comparecharacters(cc1,np->current->suffixinfo.idx,cc2,idx2);
    if(cmpresult == 1)
    {
      return false;
    }
    if(cmpresult == 0)
    {
      return true;
    }
    np->previous = np->current;
  }
  return false;
}

void insertsuffixintotrie(Trierep *trierep,
                          Trienode *node,
                          Suffixinfo *suffixinfo)
{
  if(trierep->root == NULL)
  {
    trierep->root = makeroot(trierep,suffixinfo);
  } else
  {
    Seqpos currentdepth, lcpvalue;
    Trienode *currentnode, *newleaf, *newbranch, *succ;
    Nodepair np;
    Uchar cc;
    Encseqreadinfo *eri = trierep->encseqreadinfo + suffixinfo->idx;

    assert(!ISLEAF(node));
    currentnode = node;
    currentdepth = node->depth;
    while(true)
    {
      if(suffixinfo->startpos + currentdepth >=
         getencseqtotallength(eri->encseqptr))
      {
	cc = (Uchar) SEPARATOR;
      } else
      {
	cc = getencodedchar(eri->encseqptr,
                            suffixinfo->startpos + currentdepth,
                            eri->readmode);
      }
      assert(currentnode != NULL);
      assert(!ISLEAF(currentnode));
      if(!hassuccessor(trierep,&np,currentdepth,currentnode,cc,suffixinfo->idx))
      {
	newleaf = makenewleaf(trierep,suffixinfo);
	newleaf->rightsibling = np.current;
	SHOWNODERELATIONS(newleaf);
	if(np.previous == NULL)
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
      if(ISLEAF(succ))
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
	if(np.previous == NULL)
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
      if(currentdepth + lcpvalue + 1 < succ->depth)
      {
	newbranch = makenewbranch(trierep,
				  suffixinfo,
				  currentdepth + lcpvalue + 1,
				  succ);
	if(np.previous == NULL)
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

Trienode *findsmallestnodeintrie(const Trierep *trierep)
{
  Trienode *node;

  assert(trierep->root != NULL);
  for(node = trierep->root; node->firstchild != NULL; node = node->firstchild)
    /* Nothing */ ;
  return node;
}

void deletesmallestpath(Trienode *smallest,Trierep *trierep)
{
  Trienode *father, *son;

  for(son = smallest; son->parent != NULL; son = son->parent)
  {
    father = son->parent;
    if(son->firstchild == NULL)
    {
      SETFIRSTCHILD(father,son->rightsibling);
      SHOWNODERELATIONS(father);
      son->rightsibling = NULL;
    } else
    {
      if(son->firstchild->rightsibling != NULL)
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
    trierep->unusedTrienodes[trierep->nextunused++] = son;
  }
  if(trierep->root->firstchild == NULL)
  {
    trierep->unusedTrienodes[trierep->nextunused++] = trierep->root;
    trierep->root = NULL;
  }
}

void inittrienodetable(Trierep *trierep,Seqpos numofsuffixes,
                       uint32_t numofindexes,Env *env)
{
  env_error_check(env);
  trierep->numofindexes = numofindexes;
  trierep->allocatedTrienode = MULT2(numofsuffixes + 1) + 1;
  ALLOCASSIGNSPACE(trierep->nodetable,NULL,Trienode,trierep->allocatedTrienode);
  trierep->nextfreeTrienode = 0;
  trierep->root = NULL;
  trierep->nextunused = 0;
  ALLOCASSIGNSPACE(trierep->unusedTrienodes,NULL,Trienodeptr,
                   trierep->allocatedTrienode);
}

void freetrierep(Trierep *trierep,Env *env)
{
  FREESPACE(trierep->nodetable);
  FREESPACE(trierep->unusedTrienodes);
  FREESPACE(trierep->encseqreadinfo);
}
