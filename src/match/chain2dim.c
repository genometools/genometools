/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#include <stdbool.h>
#include "core/minmax.h"
#include "core/unused_api.h"
#include "extended/redblack.h"
#include "verbose-def.h"
#include "chaindef.h"

/*
  The basic information required for each fragment is stored
  in a structure of the following type. The user has to specify
  those components which are tagged `user defined'. The chaining
  algorithms computes the remaining components tagged 'compute'
  in chain.
*/

typedef struct
{
  GtChainpostype startpos[2], /* start of fragments in the 2 dimensions, 
                                 userdef */
                 endpos[2];  /* end of fragments in the 2 dimensions, userdef */
  unsigned long firstinchain,   /* first element in chain, compute */
                previousinchain;  /* previous index in chain, compute */
  GtChainscoretype
         weight, /* weight of fragment, user defined */
         initialgap, /* gap to start of sequences, user defined */
         terminalgap, /* gap to last positions of fragment, user defined */
         score; /* score of highest scoreing chain ending here, compute */
} Fragmentinfo;

struct GtFragmentinfotable
{
   Fragmentinfo *fragments;
   unsigned long nextfree, allocated;
};

GtFragmentinfotable *fragmentinfotable_new(unsigned long numberoffragments)
{
  GtFragmentinfotable *fragmentinfotable
    = gt_malloc(sizeof (*fragmentinfotable));
  fragmentinfotable->fragments
    = gt_malloc(sizeof (*fragmentinfotable->fragments) * numberoffragments);
  fragmentinfotable->nextfree = 0;
  fragmentinfotable->allocated = numberoffragments;
  return fragmentinfotable;
}

void fragmentinfotable_delete(GtFragmentinfotable *fragmentinfotable)
{
  gt_free(fragmentinfotable->fragments);
  gt_free(fragmentinfotable);
}

void fragmentinfotable_add(GtFragmentinfotable *fragmentinfotable,
                           GtChainpostype start1,
                           GtChainpostype end1,
                           GtChainpostype start2,
                           GtChainpostype end2,
                           GtChainscoretype initialgap,
                           GtChainscoretype terminalgap,
                           GtChainscoretype weight)
{
  Fragmentinfo *frag;

  gt_assert(fragmentinfotable->nextfree < fragmentinfotable->allocated);
  frag = fragmentinfotable->fragments + fragmentinfotable->nextfree;
  frag->startpos[0] = start1;
  frag->startpos[1] = start2;
  frag->endpos[0] = end1;
  frag->endpos[1] = end2;
  frag->initialgap = initialgap;
  frag->terminalgap = terminalgap;
  frag->weight = weight;
}

#define MAKEENDPOINT(FID)       (FID)
#define FRAGIDENT(FRAG)         ((FRAG)->fpident)

#define UNDEFPREVIOUS           numofmatches

#define GETSTOREDSTARTPOINT(DIM,IDX)\
        fragmentinfo[IDX].startpos[DIM]
#define GETSTOREDENDPOINT(DIM,IDX)\
        fragmentinfo[IDX].endpos[DIM]
#define INITIALGAP(IDX)\
        fragmentinfo[IDX].initialgap
#define TERMINALGAP(IDX)\
        fragmentinfo[IDX].terminalgap

#define CHECKCHAINSPACE\
        if (lengthofchain >= chain->chainedfragments.allocatedGtChainref)\
        {\
         chain->chainedfragments.spaceGtChainref = \
           gt_realloc(chain->chainedfragments.spaceGtChainref,\
                      sizeof (*chain->chainedfragments.spaceGtChainref) *\
                              lengthofchain);\
          chain->chainedfragments.allocatedGtChainref = lengthofchain;\
          chain->chainedfragments.nextfreeGtChainref = 0;\
        }

typedef GtChainscoretype (*GtChaingapcostfunction)(const Fragmentinfo *,
                                                   unsigned long,unsigned long);

typedef struct
{
  unsigned long fpident;
  GtChainpostype fpposition;
} Fragpoint;

typedef struct
{
  GtRBTnode *dictroot;
  unsigned long *endpointperm;
} Fragmentstore;

typedef struct
{
  GtChainscoretype maxscore;
  unsigned long maxfragnum;
  bool defined;
} Maxfragvalue;

/*
  The component isavailable is used to
  (1) indicate that some score is already stored (when generating the classes)
  (2) indicate that the class representative has not yet been processed
      further (after generation)
*/

typedef struct
{
  bool isavailable;
  GtChainscoretype score;
} Bestofclass;

static bool overlappingfragments(const Fragmentinfo *fragmentinfo,
                                 unsigned long i,
                                 unsigned long j)
{
  return (GETSTOREDENDPOINT(0,i) >= GETSTOREDSTARTPOINT(0,j) ||
          GETSTOREDENDPOINT(1,i) >= GETSTOREDSTARTPOINT(1,j)) ? true : false;
}

static bool colinearfragments(const Fragmentinfo *fragmentinfo,
                              unsigned long i,
                              unsigned long j)
{
  return (GETSTOREDSTARTPOINT(0, i) < GETSTOREDSTARTPOINT(0, j) &&
          GETSTOREDENDPOINT(0, i)   < GETSTOREDENDPOINT(0, j)   &&
          GETSTOREDSTARTPOINT(1, i) < GETSTOREDSTARTPOINT(1, j) &&
          GETSTOREDENDPOINT(1, i)   < GETSTOREDENDPOINT(1, j)) ? true : false;
}

static GtChainscoretype gapcostL1(const Fragmentinfo *fragmentinfo,
                                  unsigned long i,
                                  unsigned long j)
{
  return (GtChainscoretype)
         ((GETSTOREDSTARTPOINT(0,j) - GETSTOREDENDPOINT(0,i)) +
          (GETSTOREDSTARTPOINT(1,j) - GETSTOREDENDPOINT(1,i)));
}

static GtChainscoretype overlapcost(const Fragmentinfo *fragmentinfo,
                                    unsigned long i,
                                    unsigned long j)
{
  GtChainpostype overlaplength = 0;

  /* add overlap in first dimension */
  if (GETSTOREDSTARTPOINT(0, j) <= GETSTOREDENDPOINT(0, i))
  {
    overlaplength += GETSTOREDENDPOINT(0, i) - GETSTOREDSTARTPOINT(0, j) + 1;
  }

  /* add overlap in second dimension */
  if (GETSTOREDSTARTPOINT(1, j) <= GETSTOREDENDPOINT(1, i))
  {
    overlaplength += GETSTOREDENDPOINT(1, i) - GETSTOREDSTARTPOINT(1, j) + 1;
  }
  return (GtChainscoretype) overlaplength;
}

/*
  The Chebychev distance for two qhits h=(i,j) and h'=(i',j')
  is defined as:

                      max{|i'-i-q|,|j'-j-q|},

  whereas i+q-1 is the end point in i-dimension of fragment 1 and
  j+q-1 is the end point in j-dimension of fragment 1.
  In using the match specific end points, other than fixed values for q,
  following function generalizes to MEMs as well.
*/

static GtChainscoretype gapcostCc(const Fragmentinfo *fragmentinfo,
                                  unsigned long i,unsigned long j)
{
  GtChainpostype value1, value2;

  gt_assert(GETSTOREDSTARTPOINT(0,j) > GETSTOREDENDPOINT(0,i) &&
            GETSTOREDSTARTPOINT(1,j) > GETSTOREDENDPOINT(1,i));
  value1 = GETSTOREDSTARTPOINT(0,j) - GETSTOREDENDPOINT(0,i) - 1,
  value2 = GETSTOREDSTARTPOINT(1,j) - GETSTOREDENDPOINT(1,i) - 1;
  return (GtChainscoretype) MAX(value1,value2);
}

static void chainingboundarycases(const GtChainmode *chainmode,
                                  GtChain *chain,
                                  const Fragmentinfo *fragmentinfo,
                                  unsigned long numofmatches)
{
  if (numofmatches == 0)
  {
    chain->scoreofchain = 0;
    chain->chainedfragments.nextfreeGtChainref = 0;
  } else
  {
    if (numofmatches == 1UL)
    {
      unsigned long lengthofchain;

      chain->scoreofchain = fragmentinfo[0].weight;
      if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
      {
        chain->scoreofchain -= (INITIALGAP(0) + TERMINALGAP(0));
      }
      lengthofchain = 1UL;
      CHECKCHAINSPACE;
      chain->chainedfragments.spaceGtChainref[0] = 0;
      chain->chainedfragments.nextfreeGtChainref = 1UL;
    }
  }
}

static void retracepreviousinchain(GtChain *chain,
                                   const Fragmentinfo *fragmentinfo,
                                   unsigned long numofmatches,
                                   unsigned long retracestart)
{
  unsigned long fragnum, idx, lengthofchain;

  for (lengthofchain = 0, fragnum = retracestart;
       fragnum != UNDEFPREVIOUS; lengthofchain++)
  {
    fragnum = fragmentinfo[fragnum].previousinchain;
  }
  CHECKCHAINSPACE
  fragnum = retracestart;
  idx = lengthofchain;
  while (fragnum != UNDEFPREVIOUS)
  {
    gt_assert(idx > 0);
    idx--;
    chain->chainedfragments.spaceGtChainref[idx] = fragnum;
    fragnum = fragmentinfo[fragnum].previousinchain;
  }
  gt_assert(idx == 0);
  chain->chainedfragments.nextfreeGtChainref = lengthofchain;
}

static bool checkmaxgapwidth(const Fragmentinfo *fragmentinfo,
                             GtChainpostype maxgapwidth,
                             unsigned long leftfrag,
                             unsigned long rightfrag)
{
  GtChainpostype gapwidth, startpoint, endpoint;

  startpoint = GETSTOREDSTARTPOINT(0,rightfrag);
  endpoint = GETSTOREDENDPOINT(0,leftfrag);
  if (startpoint <= endpoint)
  {
    gapwidth = 0;
  } else
  {
    gapwidth = startpoint - endpoint - 1;
  }
  if (gapwidth > maxgapwidth)
  {
    return false;
  }
  startpoint = GETSTOREDSTARTPOINT(1,rightfrag);
  endpoint = GETSTOREDENDPOINT(1,leftfrag);
  if (startpoint <= endpoint)
  {
    gapwidth = 0;
  } else
  {
    gapwidth = startpoint - endpoint - 1;
  }
  if (gapwidth > maxgapwidth)
  {
    return false;
  }
  return true;
}

static void bruteforcechainingscores(const GtChainmode *chainmode,
                                     Fragmentinfo *fragmentinfo,
                                     unsigned long numofmatches,
                                     GtChaingapcostfunction
                                       chaingapcostfunction)
{
  unsigned long previous, leftfrag, rightfrag;
  GtChainscoretype weightright, score;
  Maxfragvalue localmaxfrag;
  bool combinable;

  if (numofmatches > 1UL)
  {
    fragmentinfo[0].firstinchain = 0;
    fragmentinfo[0].previousinchain = UNDEFPREVIOUS;
    fragmentinfo[0].score = fragmentinfo[0].weight;
    if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
    {
      fragmentinfo[0].score -= (INITIALGAP(0) + TERMINALGAP(0));
    }
    for (rightfrag=1Ul; rightfrag<numofmatches; rightfrag++)
    {
      weightright = fragmentinfo[rightfrag].weight;
      localmaxfrag.defined = false;
      localmaxfrag.maxscore = 0;
      localmaxfrag.maxfragnum = 0;
      for (leftfrag=0; leftfrag<rightfrag; leftfrag++)
      {
        if (chainmode->maxgapwidth != 0 &&
           !checkmaxgapwidth(fragmentinfo,
                             chainmode->maxgapwidth,
                             leftfrag,
                             rightfrag))
        {
          combinable = false;
        } else
        {
          if (chainmode->chainkind == GLOBALCHAININGWITHOVERLAPS)
          {
            combinable = colinearfragments(fragmentinfo,leftfrag,rightfrag);
          } else
          {
            if (overlappingfragments(fragmentinfo,leftfrag,rightfrag))
            {
              combinable = false;
            } else
            {
              combinable = true;
            }
          }
        }
        if (combinable)
        {
          score = fragmentinfo[leftfrag].score;
          if (chainmode->chainkind == GLOBALCHAINING)
          {
            /* process chainkinds without gap costs */
            score += weightright;
            previous = leftfrag;
          } else
          {
            score -= chaingapcostfunction(fragmentinfo,leftfrag,rightfrag);
            if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
            {
              score += (weightright + TERMINALGAP(leftfrag)
                                    - TERMINALGAP(rightfrag));
              previous = leftfrag;
            } else
            {
              if (score > 0)
              {
                score += weightright;
                previous = leftfrag;
              } else
              {
                score = weightright;
                previous = UNDEFPREVIOUS;
              }
            }
          }
          if (!localmaxfrag.defined || localmaxfrag.maxscore < score)
          {
            localmaxfrag.maxscore = score;
            localmaxfrag.maxfragnum = previous;
            localmaxfrag.defined = true;
          }
        }
      }
      if (localmaxfrag.defined)
      {
        fragmentinfo[rightfrag].previousinchain = localmaxfrag.maxfragnum;
        if (localmaxfrag.maxfragnum == UNDEFPREVIOUS)
        {
          fragmentinfo[rightfrag].firstinchain = rightfrag;
        } else
        {
          fragmentinfo[rightfrag].firstinchain
            = fragmentinfo[localmaxfrag.maxfragnum].firstinchain;
        }
        fragmentinfo[rightfrag].score = localmaxfrag.maxscore;
      } else
      {
        fragmentinfo[rightfrag].previousinchain = UNDEFPREVIOUS;
        fragmentinfo[rightfrag].firstinchain = rightfrag;
        fragmentinfo[rightfrag].score = weightright;
        if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
        {
          fragmentinfo[rightfrag].score
            -= (INITIALGAP(rightfrag) + TERMINALGAP(rightfrag));
        }
      }
    }
  }
}

/*
  The following function compares fragment points. These are
  end points of fragments in dimension 2.
  If the fragment points are identical, then the order is undefined.
*/

static int cmpendFragpoint2(const GtKeytype keya,
                            const GtKeytype keyb,
                            GT_UNUSED void *info)
{
  if (((Fragpoint *) keya)->fpposition < ((Fragpoint *) keyb)->fpposition)
  {
    return -1;
  }
  if (((Fragpoint *) keya)->fpposition > ((Fragpoint *) keyb)->fpposition)
  {
    return 1;
  }
  if (FRAGIDENT((Fragpoint *) keya) < FRAGIDENT((Fragpoint *) keyb))
  {
    return -1;
  }
  if (FRAGIDENT((Fragpoint *) keya) > FRAGIDENT((Fragpoint *) keyb))
  {
    return 1;
  }
  return 0;
}

static GtChainscoretype evalpriority(bool addterminal,
                                     const Fragmentinfo *fragmentinfo,
                                     unsigned long fragnum)
{
  if (addterminal)
  {
    return fragmentinfo[fragnum].score - TERMINALGAP(fragnum);
  }
  return fragmentinfo[fragnum].score;
}

static void insertintodict(bool addterminal,
                           const Fragmentinfo *fragmentinfo,
                           Fragmentstore *fragmentstore,
                           Fragpoint *qfrag2)
{
  Fragpoint *retval2;
  bool nodecreated;

  retval2 = (Fragpoint *) gt_rbt_search ((const GtKeytype) qfrag2,
                                           &nodecreated,
                                           &fragmentstore->dictroot,
                                           cmpendFragpoint2,
                                           NULL);
  gt_assert(retval2 != NULL);
  if (!nodecreated)
  {
    if (evalpriority(addterminal,fragmentinfo,FRAGIDENT(retval2)) <
        evalpriority(addterminal,fragmentinfo,FRAGIDENT(qfrag2)))
    {
      gt_assert(retval2->fpposition == qfrag2->fpposition);
      retval2->fpident = qfrag2->fpident;
    }
    gt_free(qfrag2);
  }
}

static void activatefragpoint(bool addterminal,
                              const Fragmentinfo *fragmentinfo,
                              Fragmentstore *fragmentstore,
                              Fragpoint *qfrag2)
{
  Fragpoint *tmp2;
  GtChainscoretype qpriority;

  qpriority = evalpriority(addterminal,fragmentinfo,FRAGIDENT(qfrag2));
  tmp2 = (Fragpoint *) gt_rbt_previousequalkey ((const GtKeytype) qfrag2,
                                                fragmentstore->dictroot,
                                                cmpendFragpoint2,
                                                NULL);
  if (tmp2 == NULL ||
      qpriority > evalpriority(addterminal,fragmentinfo,FRAGIDENT(tmp2)))
  {
    insertintodict(addterminal,fragmentinfo,fragmentstore,qfrag2);
    while (true)
    {
      tmp2 = (Fragpoint *) gt_rbt_nextkey ((const GtKeytype) qfrag2,
                                                fragmentstore->dictroot,
                                                cmpendFragpoint2,
                                                NULL);
      if (tmp2 == NULL || qpriority <= evalpriority(addterminal,
                                                    fragmentinfo,
                                                    FRAGIDENT(tmp2)))
      {
        break;
      }
      if (gt_rbt_delete ((const GtKeytype) tmp2,
                        &fragmentstore->dictroot,
                        cmpendFragpoint2,
                        NULL) != 0)
      {
        fprintf(stderr,"cannot delete successor node\n");
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
      gt_free(tmp2);
    }
  } else
  {
    gt_free(qfrag2);
  }
}

static void evalfragmentscore(const GtChainmode *chainmode,
                              Fragmentinfo *fragmentinfo,
                              unsigned long numofmatches,
                              Fragmentstore *fragmentstore,
                              bool gapsL1,
                              unsigned long fragpointident,
                              unsigned int presortdim)
{
  unsigned long previous;
  GtChainpostype startpos2;
  Fragpoint *qfrag2;
  GtChainscoretype score;

  startpos2 = GETSTOREDSTARTPOINT(1-presortdim,fragpointident);
  if (startpos2 == 0)
  {
    qfrag2 = NULL;
  } else
  {
    Fragpoint keyfrag2;

    keyfrag2.fpposition = startpos2 - 1;  /* it is a start position */
    keyfrag2.fpident = MAKEENDPOINT(fragpointident);
                       /* but considered as endpoint */
    qfrag2 = (Fragpoint *) gt_rbt_previousequalkey(
                                  (const GtKeytype) &keyfrag2,
                                  fragmentstore->dictroot,
                                  cmpendFragpoint2,
                                  NULL);
    if (qfrag2 != NULL)
    {
      if (chainmode->maxgapwidth != 0 &&
         !checkmaxgapwidth(fragmentinfo,
                           chainmode->maxgapwidth,
                           FRAGIDENT(qfrag2),
                           fragpointident))
      {
        qfrag2 = NULL;
      }
    }
  }
  if (qfrag2 == NULL)
  {
    score = fragmentinfo[fragpointident].weight;
    if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
    {
      score -= INITIALGAP(fragpointident);
    }
    previous = UNDEFPREVIOUS;
  } else
  {
    score = fragmentinfo[FRAGIDENT(qfrag2)].score;
    if (chainmode->chainkind == GLOBALCHAINING)
    {
      score += fragmentinfo[fragpointident].weight;
      previous = FRAGIDENT(qfrag2);
    } else
    {
      GtChainscoretype tmpgc;

      if (gapsL1)
      {
        tmpgc = gapcostL1(fragmentinfo,FRAGIDENT(qfrag2),fragpointident);
      } else
      {
        tmpgc = gapcostCc(fragmentinfo,FRAGIDENT(qfrag2),fragpointident);
      }
      if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST || score > tmpgc)
      {
        score += (fragmentinfo[fragpointident].weight - tmpgc);
        previous = FRAGIDENT(qfrag2);
      } else
      {
        score = fragmentinfo[fragpointident].weight;
        previous = UNDEFPREVIOUS;
      }
    }
  }
  fragmentinfo[fragpointident].score = score;
  fragmentinfo[fragpointident].previousinchain = previous;
  if (previous == UNDEFPREVIOUS)
  {
    fragmentinfo[fragpointident].firstinchain = fragpointident;
  } else
  {
    fragmentinfo[fragpointident].firstinchain
      = fragmentinfo[previous].firstinchain;
  }
}

static bool isrightmaximallocalchain(const Fragmentinfo *fragmentinfo,
                                     unsigned long numofmatches,
                                     unsigned long currentfrag)
{
  if (currentfrag == numofmatches - 1)
  {
    return true;
  }
  if (fragmentinfo[currentfrag+1].previousinchain != currentfrag)
  {
    return true;
  }
  if (fragmentinfo[currentfrag+1].score < fragmentinfo[currentfrag].score)
  {
    return true;
  }
  return false;
}

static void determineequivreps(Bestofclass *chainequivalenceclasses,
                               const Fragmentinfo *fragmentinfo,
                               unsigned long numofmatches)
{
  unsigned long matchnum;
  Bestofclass *classptr, *classrep;

  for (classptr = chainequivalenceclasses;
       classptr < chainequivalenceclasses + numofmatches;
       classptr++)
  {
    classptr->isavailable = false;
  }
  for (matchnum=0; matchnum<numofmatches; matchnum++)
  {
    if (isrightmaximallocalchain(fragmentinfo,numofmatches,matchnum))
    {
      classrep = chainequivalenceclasses +
                 fragmentinfo[matchnum].firstinchain;
      if (!classrep->isavailable ||
          classrep->score < fragmentinfo[matchnum].score)
      {
        classrep->score = fragmentinfo[matchnum].score;
        classrep->isavailable = true;
      }
    }
  }
}

static bool retrievemaximalscore(GtChainscoretype *maxscore,
                                 const GtChainmode *chainmode,
                                 const Fragmentinfo *fragmentinfo,
                                 unsigned long numofmatches)
{
  unsigned long matchnum;
  GtChainscoretype tgap;
  bool maxscoredefined = false;

  *maxscore = 0;
  for (matchnum=0; matchnum<numofmatches; matchnum++)
  {
    if (isrightmaximallocalchain(fragmentinfo,numofmatches,matchnum))
    {
      if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
      {
        tgap = TERMINALGAP(matchnum);
      } else
      {
        tgap = 0;
      }
      if (!maxscoredefined || *maxscore < fragmentinfo[matchnum].score - tgap)
      {
        *maxscore = fragmentinfo[matchnum].score - tgap;
        maxscoredefined = true;
      }
    }
  }
  return maxscoredefined;
}

static int comparescores(const GtKeytype key1,
                         const GtKeytype key2,
                         GT_UNUSED void *info)
{
  if (*((GtChainscoretype *) key1) < *(((GtChainscoretype *) key2)))
  {
    return -1;
  }
  if (*((GtChainscoretype *) key1) > *(((GtChainscoretype *) key2)))
  {
    return 1;
  }
  return 0;
}

typedef struct
{
  unsigned long currentdictsize,     /* current size of the dictionary */
                maxdictsize;         /* maximal size of the dictionary */
  GtRBTnode *worstelement,           /* reference to worst key */
            *root,                   /* root of tree */
            *lastcallinsertedelem,   /* element inserted in last call */
            *lastcalldeletedelem;    /* element deleted in last call */
} Dictmaxsize;

static Dictmaxsize *dictmaxsize_new(unsigned long maxsize)
{
  Dictmaxsize *dict;

  dict = gt_malloc(sizeof (*dict));
  dict->currentdictsize = 0;
  dict->maxdictsize = maxsize;
  dict->root = NULL;
  dict->lastcallinsertedelem = NULL;
  dict->lastcalldeletedelem = NULL;
  dict->worstelement = NULL;
  return dict;
}

static void dictmaxsize_delete(Dictmaxsize *dict)
{
  gt_rbt_destroy (false,NULL,NULL,dict->root);
  gt_free(dict);
}

typedef int (*Dictcomparefunction)(const GtKeytype,const GtKeytype,void *);
typedef void (*Freekeyfunction)(const GtKeytype,void *);
typedef bool (*Comparewithkey)(const GtKeytype,void *);

static void insertDictmaxsize(Dictmaxsize *dict,
                              Dictcomparefunction comparefunction,
                              void *cmpinfo,
                              void *elemin)
{
  bool nodecreated;

  if (dict->currentdictsize < dict->maxdictsize)
  {
    if (dict->currentdictsize == 0 ||
        comparefunction(elemin,dict->worstelement,cmpinfo) < 0)
    {
      dict->worstelement = elemin;
    }
    (void) gt_rbt_search(elemin,
                         &nodecreated,
                         &dict->root,
                         comparefunction,
                         cmpinfo);
    if (nodecreated)
    {
      dict->currentdictsize++;
      dict->lastcallinsertedelem = elemin;
    }
  } else
  {
/*
  new element is not as worse as worst element, so insert it and
  and delete the worst element
*/
    if (comparefunction(dict->worstelement,elemin,cmpinfo) < 0)
    {
      (void) gt_rbt_search(elemin,
                           &nodecreated,
                           &dict->root,
                           comparefunction,
                           cmpinfo);
      if (nodecreated)
      {
        dict->lastcallinsertedelem = elemin;
        if (gt_rbt_delete(dict->worstelement,
                          &dict->root,
                          comparefunction,
                          cmpinfo) != 0)
        {
          fprintf(stderr,"insertDictmaxsize: deletion failed\n");
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
        dict->lastcalldeletedelem = dict->worstelement;
        dict->worstelement = gt_rbt_minimumkey(dict->root);
      }
    }
  }
}

static void retrievechainbestscores(bool *minscoredefined,
                                    GtChainscoretype *minscore,
                                    const Fragmentinfo *fragmentinfo,
                                    unsigned long numofmatches,
                                    unsigned long howmanybest)
{
  unsigned long matchnum, fragnum = 0;
  GtChainscoretype *scores;
  Dictmaxsize *dictbestfragments;
  void *minkey;

  scores = gt_malloc(sizeof (*scores) * numofmatches);
  dictbestfragments = dictmaxsize_new(howmanybest);
  for (matchnum=0; matchnum<numofmatches; matchnum++)
  {
    if (isrightmaximallocalchain(fragmentinfo,numofmatches,matchnum))
    {
      scores[fragnum] = fragmentinfo[matchnum].score;
      insertDictmaxsize(dictbestfragments,
                            comparescores,
                            NULL,
                            (void *) (scores + fragnum));
      fragnum++;
    }
  }
  if (fragnum == 0)
  {
    *minscoredefined = false;
  } else
  {
    minkey = gt_rbt_minimumkey(dictbestfragments->root);
    gt_assert(minkey != NULL);
    *minscore = *((GtChainscoretype *) minkey);
    *minscoredefined = true;
  }
  dictmaxsize_delete(dictbestfragments);
  gt_free(scores);
}

static int retrievechainthreshold(const GtChainmode *chainmode,
                                  Fragmentinfo *fragmentinfo,
                                  unsigned long numofmatches,
                                  GtChain *chain,
                                  GtChainscoretype minscore,
                                  Bestofclass *chainequivalenceclasses,
                                  GtChainprocessor chainprocessor,
                                  void *cpinfo,
                                  GtError *err)
{
  unsigned long matchnum;
  GtChainscoretype tgap;
  Bestofclass *classrep;

  for (matchnum=0; matchnum < numofmatches; matchnum++)
  {
    if (isrightmaximallocalchain(fragmentinfo,numofmatches,matchnum))
    {
      if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
      {
        tgap = TERMINALGAP(matchnum);
      } else
      {
        tgap = 0;
      }
      if (fragmentinfo[matchnum].score - tgap >= minscore)
      {
        if (chainequivalenceclasses != NULL)
        {
          classrep = chainequivalenceclasses +
                     fragmentinfo[matchnum].firstinchain;
          gt_assert(classrep != NULL);
          if (classrep->isavailable &&
              classrep->score == fragmentinfo[matchnum].score - tgap)
          {
            chain->scoreofchain = classrep->score;
            classrep->isavailable = false;
            retracepreviousinchain(chain,
                                   fragmentinfo,
                                   numofmatches,
                                   matchnum);
            if (chainprocessor(cpinfo,chain,err) != 0)
            {
              return -1;
            }
          }
        } else
        {
          chain->scoreofchain = fragmentinfo[matchnum].score - tgap;
          retracepreviousinchain(chain,
                                 fragmentinfo,
                                 numofmatches,
                                 matchnum);
          if (chainprocessor(cpinfo,chain,err) != 0)
          {
            return -1;
          }
        }
      }
    }
  }
  return 0;
}

static int comparestartandend(const Fragmentinfo *sortedstartpoints,
                              const Fragmentinfo *sortedendpoints,
                              unsigned int presortdim)
{
  if (sortedstartpoints->startpos[presortdim] <
      sortedendpoints->endpos[presortdim])
  {
    return -1;
  }
  if (sortedstartpoints->startpos[presortdim] >
      sortedendpoints->endpos[presortdim])
  {
    return 1;
  }
  return -1;
}

static Fragpoint *makeactivationpoint(const Fragmentinfo *fragmentinfo,
                                      unsigned long fpident,
                                      unsigned int postsortdim)
{
  Fragpoint *fpptr;

  fpptr = gt_malloc(sizeof (*fpptr));
  gt_assert(fpptr != NULL);
  fpptr->fpident = MAKEENDPOINT(fpident);
  fpptr->fpposition = GETSTOREDENDPOINT(postsortdim,fpident);
  return fpptr;
}

static void mergestartandendpoints(const GtChainmode *chainmode,
                                   Fragmentinfo *fragmentinfo,
                                   unsigned long numofmatches,
                                   Fragmentstore *fragmentstore,
                                   bool gapsL1,
                                   unsigned int presortdim)
{
  unsigned long xidx, startcount, endcount;
  bool addterminal;
  unsigned int postsortdim = 1U - presortdim;

  addterminal = (chainmode->chainkind == GLOBALCHAINING) ? false : true;
  fragmentstore->dictroot = NULL;
  for (xidx = 0, startcount = 0, endcount = 0;
       startcount < numofmatches && endcount < numofmatches;
       xidx++)
  {
    if (comparestartandend(fragmentinfo + startcount,
                          fragmentinfo + fragmentstore->endpointperm[endcount],
                          presortdim) < 0)
    {
      evalfragmentscore(chainmode,
                        fragmentinfo,
                        numofmatches,
                        fragmentstore,
                        gapsL1,
                        startcount,
                        presortdim);
      startcount++;
    } else
    {
      activatefragpoint(addterminal,fragmentinfo,fragmentstore,
                        makeactivationpoint(fragmentinfo,
                                            fragmentstore->
                                            endpointperm[endcount],
                                            postsortdim));
      endcount++;
    }
  }
  while (startcount < numofmatches)
  {
    evalfragmentscore(chainmode,
                      fragmentinfo,
                      numofmatches,
                      fragmentstore,
                      gapsL1,
                      startcount,
                      presortdim);
    startcount++;
    xidx++;
  }
  while (endcount < numofmatches)
  {
    activatefragpoint(addterminal,fragmentinfo,fragmentstore,
                      makeactivationpoint(fragmentinfo,
                                          fragmentstore->
                                          endpointperm[endcount],
                                          postsortdim));
    endcount++;
    xidx++;
  }
}

static int findmaximalscores(const GtChainmode *chainmode,
                             GtChain *chain,
                             Fragmentinfo *fragmentinfo,
                             unsigned long numofmatches,
                             Fragmentstore *fragmentstore,
                             GtChainprocessor chainprocessor,
                             bool withequivclasses,
                             void *cpinfo,
                             Verboseinfo *verboseinfo,
                             GtError *err)
{
  unsigned long fragnum;
  GtChainscoretype minscore = 0;
  Fragpoint *maxpoint;
  bool minscoredefined = false, haserr = false;
  Bestofclass *chainequivalenceclasses;
  int retval;

  if (withequivclasses)
  {
    if (chainmode->chainkind == LOCALCHAININGMAX ||
        chainmode->chainkind == LOCALCHAININGTHRESHOLD ||
        chainmode->chainkind == LOCALCHAININGBEST ||
        chainmode->chainkind == LOCALCHAININGPERCENTAWAY)
    {
      chainequivalenceclasses = gt_malloc(sizeof (*chainequivalenceclasses) *
                                          numofmatches);
      determineequivreps(chainequivalenceclasses, fragmentinfo, numofmatches);
    } else
    {
      chainequivalenceclasses = NULL;
    }
  } else
  {
    chainequivalenceclasses = NULL;
  }
  switch (chainmode->chainkind)
  {
    case GLOBALCHAINING:
      maxpoint = gt_rbt_maximumkey(fragmentstore->dictroot);
      gt_assert(maxpoint != NULL);
      fragnum = FRAGIDENT(maxpoint);
      minscore = fragmentinfo[fragnum].score;
      minscoredefined = true;
      break;
    case GLOBALCHAININGWITHGAPCOST:
    case GLOBALCHAININGWITHOVERLAPS:
    case LOCALCHAININGMAX:
      minscoredefined
        = retrievemaximalscore(&minscore,chainmode,fragmentinfo,numofmatches);
      break;
    case LOCALCHAININGTHRESHOLD:
      minscore = chainmode->minimumscore;
      minscoredefined = true;
      break;
    case LOCALCHAININGBEST:
      retrievechainbestscores(&minscoredefined,
                              &minscore,
                              fragmentinfo,
                              numofmatches,
                              chainmode->howmanybest);
      break;
    case LOCALCHAININGPERCENTAWAY:
      minscoredefined
        = retrievemaximalscore(&minscore,chainmode,fragmentinfo,numofmatches);
      if (minscoredefined)
      {
        minscore = (GtChainscoretype)
                   ((double) minscore *
                    (1.0 - (double) chainmode->percentawayfrombest/100.0));
      }
      break;
    default:
      fprintf(stderr,"chainkind = %d not valid\n",
               (int) chainmode->chainkind);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  if (minscoredefined)
  {
    showverbose(verboseinfo,
               "compute optimal %s chains with score >= %ld",
               (chainmode->chainkind == GLOBALCHAINING ||
                chainmode->chainkind == GLOBALCHAININGWITHGAPCOST ||
                chainmode->chainkind == GLOBALCHAININGWITHOVERLAPS)
                ? "global"
                : "local",
               minscore);
    if (retrievechainthreshold(chainmode,
                               fragmentinfo,
                               numofmatches,
                               chain,
                               minscore,
                               chainequivalenceclasses,
                               chainprocessor,
                               cpinfo,
                               err) != 0)
    {
      haserr = true;
    }
    retval = 0;
  } else
  {
    retval = 1;
  }
  if (chainequivalenceclasses != NULL)
  {
    gt_free(chainequivalenceclasses);
  }
  return haserr ? -1 : retval;
}

static void makesortedendpointpermutation(unsigned long *perm,
                                          Fragmentinfo *fragmentinfo,
                                          unsigned long numofmatches,
                                          unsigned int presortdim)
{
  unsigned long temp, *iptr, *jptr, i, moves = 0;

  for (i = 0; i < numofmatches; i++)
  {
    perm[i] = i;
  }
  for (iptr = perm + 1UL; iptr < perm + numofmatches; iptr++)
  {
    for (jptr = iptr; jptr > perm; jptr--)
    {
      if (GETSTOREDENDPOINT(presortdim,*(jptr-1)) <=
          GETSTOREDENDPOINT(presortdim,*jptr))
      {
        break;
      }
      temp = *(jptr-1);
      *(jptr-1) = *jptr;
      *jptr = temp;
      moves++;
    }
  }
}

static void fastchainingscores(const GtChainmode *chainmode,
                               Fragmentinfo *fragmentinfo,
                               unsigned long numofmatches,
                               Fragmentstore *fragmentstore,
                               unsigned int presortdim,
                               bool gapsL1)
{
  fragmentstore->endpointperm
    = gt_malloc(sizeof (*fragmentstore->endpointperm) * numofmatches);
  makesortedendpointpermutation(fragmentstore->endpointperm,
                                fragmentinfo,
                                numofmatches,
                                presortdim);
  mergestartandendpoints(chainmode,
                         fragmentinfo,
                         numofmatches,
                         fragmentstore,
                         gapsL1,
                         presortdim);
  gt_free(fragmentstore->endpointperm);
}

static GtChaingapcostfunction assignchaingapcostfunction(GtChainkind chainkind,
                                                         bool gapsL1)
{
  if (chainkind == GLOBALCHAININGWITHOVERLAPS)
  {
    return overlapcost;
  }
  if (chainkind != GLOBALCHAINING)
  {
    if (gapsL1)
    {
      return gapcostL1;
    }
    return gapcostCc;
  }
  return NULL;
}

/*EE
  The following function implements the different kinds of chaining
  algorithms for global and local chaining. The mode is specified
  by \texttt{chainmode}. There are \texttt{numofmatches} many matches,
  each specified in an array \texttt{fragmentinfo} points to.
  \texttt{chain} is the array in which a chain is stored. However,
  do not further process chain. Just use
  INITARRAY(&chain.chainedfragments,Uint) to initialize it before calling
  \textttt{fastchaining} and
  FREEARRAY(&chain.chainedfragments,Uint) to free the space after calling
  \textttt{fastchaining}. The function \texttt{chainprocessor} processes
  each chain found, and \texttt{cpinfo} is used as a the first argument
  in each call to \texttt{chainprocessor}. The function returns
  0 upon success and otherwise, if an error occurs.
  Finally \texttt{showverbose} is applied to each status message generated
  during the execution of \texttt{fastchaining}. If \texttt{showverbose}
  is \texttt{NULL}, then nothing is generated and shown.
*/

int fastchaining(const GtChainmode *chainmode,
                 GtChain *chain,
                 Fragmentinfo *fragmentinfo,
                 unsigned long numofmatches,
                 bool gapsL1,
                 unsigned int presortdim,
                 bool withequivclasses,
                 GtChainprocessor chainprocessor,
                 void *cpinfo,
                 Verboseinfo *verboseinfo,
                 GtError *err)
{
  int retcode;
  GtChaingapcostfunction chaingapcostfunction;
  bool haserr = false;

  gt_assert(presortdim <= 1U);
  chaingapcostfunction
    = assignchaingapcostfunction(chainmode->chainkind,gapsL1);
  if (numofmatches > 1UL)
  {
    Fragmentstore fragmentstore;

    showverbose(verboseinfo,"compute chain scores");
    if (chainmode->chainkind == GLOBALCHAININGWITHOVERLAPS)
    {
      bruteforcechainingscores(chainmode,
                               fragmentinfo,
                               numofmatches,
                               chaingapcostfunction);
    } else
    {
      fastchainingscores(chainmode,
                         fragmentinfo,
                         numofmatches,
                         &fragmentstore,
                         presortdim,
                         gapsL1);
    }
    showverbose(verboseinfo,"retrieve optimal chains");
    retcode = findmaximalscores(chainmode,
                                chain,
                                fragmentinfo,
                                numofmatches,
                                &fragmentstore,
                                chainprocessor,
                                withequivclasses,
                                cpinfo,
                                verboseinfo,
                                err);
    if (chainmode->chainkind != GLOBALCHAININGWITHOVERLAPS)
    {
      gt_rbt_destroy (true,NULL,NULL,fragmentstore.dictroot);
    }
    if (retcode < 0)
    {
      haserr = true;
    }
  } else
  {
    chainingboundarycases(chainmode,
                          chain,
                          fragmentinfo,
                          numofmatches);
    if (chainprocessor(cpinfo,chain,err))
    {
      haserr = true;
    }
    retcode = 0;
  }
  return haserr ? -1 : retcode;
}
