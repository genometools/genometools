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
#include <string.h>
#include "core/minmax.h"
#include "core/mathsupport.h"
#include "core/str_api.h"
#include "core/str_array_api.h"
#include "core/unused_api.h"
#include "core/arraydef.h"
#include "stamp.h"
#include "extended/redblack.h"
#include "verbose-def.h"
#include "chain2dim.h"
#include "prsqualint.h"

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
   GtChainscoretype largestdim1,
                    largestdim2;
   unsigned long nextfree,
                 allocated;
};

/*
  The following type defines the possible kinds of chaining.
  The mode can be one of the two following values.
*/

typedef enum
{
  GLOBALCHAINING,              /* global chaining without gap costs */
  GLOBALCHAININGWITHGAPCOST,   /* global chaining with L1 gap costs */
  GLOBALCHAININGWITHOVERLAPS,  /* chaining with overlaps */
  LOCALCHAININGMAX,            /* local chaining; one maximum is reported */
  LOCALCHAININGTHRESHOLD,      /* local chaining; all chains >= minscore */
  LOCALCHAININGBEST,           /* local chaining; k best local chains */
  LOCALCHAININGPERCENTAWAY     /* local chaining; percent away from best */
} GtChainkind;

/*
  The following type defines the chain mode consisting of a chainkind.
  If chainkind = LOCALCHAININGTHRESHOLD, then an additional
  component minimumscore is used.
  If chainkind = LOCALCHAININGBEST, then  an additional
  component howmanybest is used.
  If chainkind = LOCALCHAININGPERCENTAWAY, then  an additional
  component percentawayfrombest is defined
*/

struct GtChainmode
{
  GtChainkind chainkind;
  GtChainpostype maxgapwidth;  /* 0 if undefined or
                                  otherwise maximal width of gap */
  GtChainscoretype minimumscore; /* only defined if
                                  chainkind = LOCALCHAININGTHRESHOLD */
  unsigned long howmanybest,   /* only defined if
                                  chainkind = LOCALCHAININGBEST */
                percentawayfrombest;  /* only defined if
                                         chainkind = LOCALCHAININGPERCENTAWAY */
};

typedef unsigned long GtChainref;

GT_DECLAREARRAYSTRUCT(GtChainref);

/*
  A chain consists of an array of integers. These refer to the array of
  fragment informations.
*/

struct GtChain
{
  GtArrayGtChainref chainedfragments;
  GtChainscoretype scoreofchain;
};

GtFragmentinfotable *gt_chain_fragmentinfotable_new(
                               unsigned long numberoffragments)
{
  GtFragmentinfotable *fragmentinfotable
    = gt_malloc(sizeof (*fragmentinfotable));
  fragmentinfotable->fragments
    = gt_malloc(sizeof (*fragmentinfotable->fragments) * numberoffragments);
  fragmentinfotable->nextfree = 0;
  fragmentinfotable->allocated = numberoffragments;
  fragmentinfotable->largestdim1 = fragmentinfotable->largestdim2 = 0;
  return fragmentinfotable;
}

void gt_chain_fragmentinfotable_delete(GtFragmentinfotable *fragmentinfotable)
{
  if (fragmentinfotable != NULL)
  {
    gt_free(fragmentinfotable->fragments);
    gt_free(fragmentinfotable);
  }
}

void gt_chain_fragmentinfotable_empty(GtFragmentinfotable *fragmentinfotable)
{
  gt_assert(fragmentinfotable != NULL);
  fragmentinfotable->largestdim1 = fragmentinfotable->largestdim2 = 0;
  fragmentinfotable->nextfree = 0;
}

void gt_chain_fragmentinfotable_add(GtFragmentinfotable *fragmentinfotable,
                                    GtChainpostype start1,
                                    GtChainpostype end1,
                                    GtChainpostype start2,
                                    GtChainpostype end2,
                                    GtChainscoretype weight)
{
  Fragmentinfo *frag;

  gt_assert(fragmentinfotable->nextfree < fragmentinfotable->allocated);
  frag = fragmentinfotable->fragments + fragmentinfotable->nextfree++;
  frag->startpos[0] = start1;
  frag->startpos[1] = start2;
  frag->endpos[0] = end1;
  frag->endpos[1] = end2;
  frag->weight = weight;
  if (fragmentinfotable->largestdim1 < (GtChainscoretype) end1)
  {
    fragmentinfotable->largestdim1 = (GtChainscoretype) end1;
  }
  if (fragmentinfotable->largestdim2 < (GtChainscoretype) end2)
  {
    fragmentinfotable->largestdim2 = (GtChainscoretype) end2;
  }
}

void gt_chain_fillthegapvalues(GtFragmentinfotable *fragmentinfotable)
{
  Fragmentinfo *fiptr;

  for (fiptr = fragmentinfotable->fragments;
       fiptr < fragmentinfotable->fragments + fragmentinfotable->nextfree;
       fiptr++)
  {
    fiptr->initialgap
      = (GtChainscoretype) (fiptr->startpos[0] + fiptr->startpos[1]);
    fiptr->terminalgap
      = (GtChainscoretype) (fragmentinfotable->largestdim1 - fiptr->endpos[0] +
                            fragmentinfotable->largestdim2 - fiptr->endpos[1]);
  }
}

#define MAKEENDPOINT(FID)       (FID)
#define FRAGIDENT(FRAG)         ((FRAG)->fpident)

#define UNDEFPREVIOUS           fragmentinfotable->nextfree

#define GETSTOREDSTARTPOINT(DIM,IDX)\
        fragmentinfotable->fragments[IDX].startpos[DIM]
#define GETSTOREDENDPOINT(DIM,IDX)\
        fragmentinfotable->fragments[IDX].endpos[DIM]
#define INITIALGAP(IDX)\
        fragmentinfotable->fragments[IDX].initialgap
#define TERMINALGAP(IDX)\
        fragmentinfotable->fragments[IDX].terminalgap

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

typedef GtChainscoretype (*GtChaingapcostfunction)(const GtFragmentinfotable *,
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

static bool overlappingfragments(const GtFragmentinfotable *fragmentinfotable,
                                 unsigned long i,
                                 unsigned long j)
{
  return (GETSTOREDENDPOINT(0,i) >= GETSTOREDSTARTPOINT(0,j) ||
          GETSTOREDENDPOINT(1,i) >= GETSTOREDSTARTPOINT(1,j)) ? true : false;
}

static bool colinearfragments(const GtFragmentinfotable *fragmentinfotable,
                              unsigned long i,
                              unsigned long j)
{
  return (GETSTOREDSTARTPOINT(0, i) < GETSTOREDSTARTPOINT(0, j) &&
          GETSTOREDENDPOINT(0, i)   < GETSTOREDENDPOINT(0, j)   &&
          GETSTOREDSTARTPOINT(1, i) < GETSTOREDSTARTPOINT(1, j) &&
          GETSTOREDENDPOINT(1, i)   < GETSTOREDENDPOINT(1, j)) ? true : false;
}

static GtChainscoretype gapcostL1(const GtFragmentinfotable *fragmentinfotable,
                                  unsigned long i,
                                  unsigned long j)
{
  return (GtChainscoretype)
         ((GETSTOREDSTARTPOINT(0,j) - GETSTOREDENDPOINT(0,i)) +
          (GETSTOREDSTARTPOINT(1,j) - GETSTOREDENDPOINT(1,i)));
}

static GtChainscoretype overlapcost(
                            const GtFragmentinfotable *fragmentinfotable,
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

static GtChainscoretype gapcostCc(const GtFragmentinfotable *fragmentinfotable,
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
                                  const GtFragmentinfotable *fragmentinfotable)
{
  if (fragmentinfotable->nextfree == 0)
  {
    chain->scoreofchain = 0;
    chain->chainedfragments.nextfreeGtChainref = 0;
  } else
  {
    if (fragmentinfotable->nextfree == 1UL)
    {
      unsigned long lengthofchain;

      chain->scoreofchain = fragmentinfotable->fragments[0].weight;
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
                                   const GtFragmentinfotable *fragmentinfotable,
                                   unsigned long retracestart)
{
  unsigned long fragnum, idx, lengthofchain;

  for (lengthofchain = 0, fragnum = retracestart;
       fragnum != UNDEFPREVIOUS; lengthofchain++)
  {
    fragnum = fragmentinfotable->fragments[fragnum].previousinchain;
  }
  CHECKCHAINSPACE
  fragnum = retracestart;
  idx = lengthofchain;
  while (fragnum != UNDEFPREVIOUS)
  {
    gt_assert(idx > 0);
    idx--;
    chain->chainedfragments.spaceGtChainref[idx] = fragnum;
    fragnum = fragmentinfotable->fragments[fragnum].previousinchain;
  }
  gt_assert(idx == 0);
  chain->chainedfragments.nextfreeGtChainref = lengthofchain;
}

static bool checkmaxgapwidth(const GtFragmentinfotable *fragmentinfotable,
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
                                     GtFragmentinfotable *fragmentinfotable,
                                     GtChaingapcostfunction
                                       chaingapcostfunction)
{
  unsigned long previous, leftfrag, rightfrag;
  GtChainscoretype weightright, score;
  Maxfragvalue localmaxfrag;
  bool combinable;

  if (fragmentinfotable->nextfree > 1UL)
  {
    fragmentinfotable->fragments[0].firstinchain = 0;
    fragmentinfotable->fragments[0].previousinchain = UNDEFPREVIOUS;
    fragmentinfotable->fragments[0].score
      = fragmentinfotable->fragments[0].weight;
    if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
    {
      fragmentinfotable->fragments[0].score -= (INITIALGAP(0) + TERMINALGAP(0));
    }
    for (rightfrag=1Ul; rightfrag<fragmentinfotable->nextfree; rightfrag++)
    {
      weightright = fragmentinfotable->fragments[rightfrag].weight;
      localmaxfrag.defined = false;
      localmaxfrag.maxscore = 0;
      localmaxfrag.maxfragnum = 0;
      for (leftfrag=0; leftfrag<rightfrag; leftfrag++)
      {
        if (chainmode->maxgapwidth != 0 &&
           !checkmaxgapwidth(fragmentinfotable,
                             chainmode->maxgapwidth,
                             leftfrag,
                             rightfrag))
        {
          combinable = false;
        } else
        {
          if (chainmode->chainkind == GLOBALCHAININGWITHOVERLAPS)
          {
            combinable = colinearfragments(fragmentinfotable,leftfrag,
                                           rightfrag);
          } else
          {
            if (overlappingfragments(fragmentinfotable,leftfrag,rightfrag))
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
          score = fragmentinfotable->fragments[leftfrag].score;
          if (chainmode->chainkind == GLOBALCHAINING)
          {
            /* process chainkinds without gap costs */
            score += weightright;
            previous = leftfrag;
          } else
          {
            score -= chaingapcostfunction(fragmentinfotable,leftfrag,rightfrag);
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
        fragmentinfotable->fragments[rightfrag].previousinchain
          = localmaxfrag.maxfragnum;
        if (localmaxfrag.maxfragnum == UNDEFPREVIOUS)
        {
          fragmentinfotable->fragments[rightfrag].firstinchain = rightfrag;
        } else
        {
          fragmentinfotable->fragments[rightfrag].firstinchain
            = fragmentinfotable->fragments[localmaxfrag.maxfragnum].
              firstinchain;
        }
        fragmentinfotable->fragments[rightfrag].score = localmaxfrag.maxscore;
      } else
      {
        fragmentinfotable->fragments[rightfrag].previousinchain = UNDEFPREVIOUS;
        fragmentinfotable->fragments[rightfrag].firstinchain = rightfrag;
        fragmentinfotable->fragments[rightfrag].score = weightright;
        if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
        {
          fragmentinfotable->fragments[rightfrag].score
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

static GtChainscoretype evalpriority(
                               bool addterminal,
                               const GtFragmentinfotable *fragmentinfotable,
                               unsigned long fragnum)
{
  if (addterminal)
  {
    return fragmentinfotable->fragments[fragnum].score - TERMINALGAP(fragnum);
  }
  return fragmentinfotable->fragments[fragnum].score;
}

static void insertintodict(bool addterminal,
                           const GtFragmentinfotable *fragmentinfotable,
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
    if (evalpriority(addterminal,fragmentinfotable,FRAGIDENT(retval2)) <
        evalpriority(addterminal,fragmentinfotable,FRAGIDENT(qfrag2)))
    {
      gt_assert(retval2->fpposition == qfrag2->fpposition);
      retval2->fpident = qfrag2->fpident;
    }
    gt_free(qfrag2);
  }
}

static void activatefragpoint(bool addterminal,
                              const GtFragmentinfotable *fragmentinfotable,
                              Fragmentstore *fragmentstore,
                              Fragpoint *qfrag2)
{
  Fragpoint *tmp2;
  GtChainscoretype qpriority;

  qpriority = evalpriority(addterminal,fragmentinfotable,FRAGIDENT(qfrag2));
  tmp2 = (Fragpoint *) gt_rbt_previousequalkey ((const GtKeytype) qfrag2,
                                                fragmentstore->dictroot,
                                                cmpendFragpoint2,
                                                NULL);
  if (tmp2 == NULL ||
      qpriority > evalpriority(addterminal,fragmentinfotable,FRAGIDENT(tmp2)))
  {
    insertintodict(addterminal,fragmentinfotable,fragmentstore,qfrag2);
    while (true)
    {
      tmp2 = (Fragpoint *) gt_rbt_nextkey ((const GtKeytype) qfrag2,
                                           fragmentstore->dictroot,
                                           cmpendFragpoint2,
                                           NULL);
      if (tmp2 == NULL || qpriority <= evalpriority(addterminal,
                                                    fragmentinfotable,
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
                              GtFragmentinfotable *fragmentinfotable,
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
         !checkmaxgapwidth(fragmentinfotable,
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
    score = fragmentinfotable->fragments[fragpointident].weight;
    if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
    {
      score -= INITIALGAP(fragpointident);
    }
    previous = UNDEFPREVIOUS;
  } else
  {
    score = fragmentinfotable->fragments[FRAGIDENT(qfrag2)].score;
    if (chainmode->chainkind == GLOBALCHAINING)
    {
      score += fragmentinfotable->fragments[fragpointident].weight;
      previous = FRAGIDENT(qfrag2);
    } else
    {
      GtChainscoretype tmpgc;

      if (gapsL1)
      {
        tmpgc = gapcostL1(fragmentinfotable,FRAGIDENT(qfrag2),fragpointident);
      } else
      {
        tmpgc = gapcostCc(fragmentinfotable,FRAGIDENT(qfrag2),fragpointident);
      }
      if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST || score > tmpgc)
      {
        score += (fragmentinfotable->fragments[fragpointident].weight - tmpgc);
        previous = FRAGIDENT(qfrag2);
      } else
      {
        score = fragmentinfotable->fragments[fragpointident].weight;
        previous = UNDEFPREVIOUS;
      }
    }
  }
  fragmentinfotable->fragments[fragpointident].score = score;
  fragmentinfotable->fragments[fragpointident].previousinchain = previous;
  if (previous == UNDEFPREVIOUS)
  {
    fragmentinfotable->fragments[fragpointident].firstinchain = fragpointident;
  } else
  {
    fragmentinfotable->fragments[fragpointident].firstinchain
      = fragmentinfotable->fragments[previous].firstinchain;
  }
}

static bool isrightmaximallocalchain(
                        const GtFragmentinfotable *fragmentinfotable,
                        unsigned long currentfrag)
{
  if (currentfrag == fragmentinfotable->nextfree - 1)
  {
    return true;
  }
  if (fragmentinfotable->fragments[currentfrag+1].previousinchain !=
      currentfrag)
  {
    return true;
  }
  if (fragmentinfotable->fragments[currentfrag+1].score <
      fragmentinfotable->fragments[currentfrag].score)
  {
    return true;
  }
  return false;
}

static void determineequivreps(Bestofclass *chainequivalenceclasses,
                               const GtFragmentinfotable *fragmentinfotable)
{
  unsigned long matchnum;
  Bestofclass *classptr, *classrep;

  for (classptr = chainequivalenceclasses;
       classptr < chainequivalenceclasses + fragmentinfotable->nextfree;
       classptr++)
  {
    classptr->isavailable = false;
  }
  for (matchnum=0; matchnum<fragmentinfotable->nextfree; matchnum++)
  {
    if (isrightmaximallocalchain(fragmentinfotable,matchnum))
    {
      classrep = chainequivalenceclasses +
                 fragmentinfotable->fragments[matchnum].firstinchain;
      if (!classrep->isavailable ||
          classrep->score < fragmentinfotable->fragments[matchnum].score)
      {
        classrep->score = fragmentinfotable->fragments[matchnum].score;
        classrep->isavailable = true;
      }
    }
  }
}

static bool retrievemaximalscore(GtChainscoretype *maxscore,
                                 const GtChainmode *chainmode,
                                 const GtFragmentinfotable *fragmentinfotable)
{
  unsigned long matchnum;
  GtChainscoretype tgap;
  bool maxscoredefined = false;

  *maxscore = 0;
  for (matchnum=0; matchnum<fragmentinfotable->nextfree; matchnum++)
  {
    if (isrightmaximallocalchain(fragmentinfotable,matchnum))
    {
      if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
      {
        tgap = TERMINALGAP(matchnum);
      } else
      {
        tgap = 0;
      }
      if (!maxscoredefined ||
          *maxscore < fragmentinfotable->fragments[matchnum].score - tgap)
      {
        *maxscore = fragmentinfotable->fragments[matchnum].score - tgap;
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

static void retrievechainbestscores(
                          bool *minscoredefined,
                          GtChainscoretype *minscore,
                          const GtFragmentinfotable *fragmentinfotable,
                          unsigned long howmanybest)
{
  unsigned long matchnum, fragnum = 0;
  GtChainscoretype *scores;
  Dictmaxsize *dictbestfragments;
  void *minkey;

  scores = gt_malloc(sizeof (*scores) * fragmentinfotable->nextfree);
  dictbestfragments = dictmaxsize_new(howmanybest);
  for (matchnum=0; matchnum<fragmentinfotable->nextfree; matchnum++)
  {
    if (isrightmaximallocalchain(fragmentinfotable,matchnum))
    {
      scores[fragnum] = fragmentinfotable->fragments[matchnum].score;
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
                                  GtFragmentinfotable *fragmentinfotable,
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

  for (matchnum=0; matchnum < fragmentinfotable->nextfree; matchnum++)
  {
    if (isrightmaximallocalchain(fragmentinfotable,matchnum))
    {
      if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
      {
        tgap = TERMINALGAP(matchnum);
      } else
      {
        tgap = 0;
      }
      if (fragmentinfotable->fragments[matchnum].score - tgap >= minscore)
      {
        if (chainequivalenceclasses != NULL)
        {
          classrep = chainequivalenceclasses +
                     fragmentinfotable->fragments[matchnum].firstinchain;
          gt_assert(classrep != NULL);
          if (classrep->isavailable &&
              classrep->score == fragmentinfotable->fragments[matchnum].score -
                                 tgap)
          {
            chain->scoreofchain = classrep->score;
            classrep->isavailable = false;
            retracepreviousinchain(chain, fragmentinfotable, matchnum);
            if (chainprocessor(cpinfo,fragmentinfotable,chain,err) != 0)
            {
              return -1;
            }
          }
        } else
        {
          chain->scoreofchain = fragmentinfotable->fragments[matchnum].score -
                                tgap;
          retracepreviousinchain(chain, fragmentinfotable, matchnum);
          if (chainprocessor(cpinfo,fragmentinfotable,chain,err) != 0)
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

static Fragpoint *makeactivationpoint(
                            const GtFragmentinfotable *fragmentinfotable,
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
                                   GtFragmentinfotable *fragmentinfotable,
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
       startcount < fragmentinfotable->nextfree &&
       endcount < fragmentinfotable->nextfree;
       xidx++)
  {
    if (comparestartandend(fragmentinfotable->fragments + startcount,
                           fragmentinfotable->fragments +
                           fragmentstore->endpointperm[endcount],
                           presortdim) < 0)
    {
      evalfragmentscore(chainmode,
                        fragmentinfotable,
                        fragmentstore,
                        gapsL1,
                        startcount,
                        presortdim);
      startcount++;
    } else
    {
      activatefragpoint(addterminal,fragmentinfotable,fragmentstore,
                        makeactivationpoint(fragmentinfotable,
                                            fragmentstore->
                                            endpointperm[endcount],
                                            postsortdim));
      endcount++;
    }
  }
  while (startcount < fragmentinfotable->nextfree)
  {
    evalfragmentscore(chainmode,
                      fragmentinfotable,
                      fragmentstore,
                      gapsL1,
                      startcount,
                      presortdim);
    startcount++;
    xidx++;
  }
  while (endcount < fragmentinfotable->nextfree)
  {
    activatefragpoint(addterminal,fragmentinfotable,fragmentstore,
                      makeactivationpoint(fragmentinfotable,
                                          fragmentstore->
                                          endpointperm[endcount],
                                          postsortdim));
    endcount++;
    xidx++;
  }
}

static int findmaximalscores(const GtChainmode *chainmode,
                             GtChain *chain,
                             GtFragmentinfotable *fragmentinfotable,
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
                                          fragmentinfotable->nextfree);
      determineequivreps(chainequivalenceclasses, fragmentinfotable);
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
      minscore = fragmentinfotable->fragments[fragnum].score;
      minscoredefined = true;
      break;
    case GLOBALCHAININGWITHGAPCOST:
    case GLOBALCHAININGWITHOVERLAPS:
    case LOCALCHAININGMAX:
      minscoredefined
        = retrievemaximalscore(&minscore,chainmode,fragmentinfotable);
      break;
    case LOCALCHAININGTHRESHOLD:
      minscore = chainmode->minimumscore;
      minscoredefined = true;
      break;
    case LOCALCHAININGBEST:
      retrievechainbestscores(&minscoredefined,
                              &minscore,
                              fragmentinfotable,
                              chainmode->howmanybest);
      break;
    case LOCALCHAININGPERCENTAWAY:
      minscoredefined
        = retrievemaximalscore(&minscore,chainmode,fragmentinfotable);
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
                               fragmentinfotable,
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

static void makesortedendpointpermutation(
                                     unsigned long *perm,
                                     GtFragmentinfotable *fragmentinfotable,
                                     unsigned int presortdim)
{
  unsigned long temp, *iptr, *jptr, i, moves = 0;

  for (i = 0; i < fragmentinfotable->nextfree; i++)
  {
    perm[i] = i;
  }
  for (iptr = perm + 1UL; iptr < perm + fragmentinfotable->nextfree; iptr++)
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
                               GtFragmentinfotable *fragmentinfotable,
                               Fragmentstore *fragmentstore,
                               unsigned int presortdim,
                               bool gapsL1)
{
  fragmentstore->endpointperm
    = gt_malloc(sizeof (*fragmentstore->endpointperm) *
                fragmentinfotable->nextfree);
  makesortedendpointpermutation(fragmentstore->endpointperm,
                                fragmentinfotable,
                                presortdim);
  mergestartandendpoints(chainmode,
                         fragmentinfotable,
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

int gt_chain_fastchaining(const GtChainmode *chainmode,
                          GtChain *chain,
                          GtFragmentinfotable *fragmentinfotable,
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
  if (fragmentinfotable->nextfree > 1UL)
  {
    Fragmentstore fragmentstore;

    showverbose(verboseinfo,"compute chain scores");
    if (chainmode->chainkind == GLOBALCHAININGWITHOVERLAPS)
    {
      bruteforcechainingscores(chainmode,fragmentinfotable,
                               chaingapcostfunction);
    } else
    {
      fastchainingscores(chainmode,
                         fragmentinfotable,
                         &fragmentstore,
                         presortdim,
                         gapsL1);
    }
    showverbose(verboseinfo,"retrieve optimal chains");
    retcode = findmaximalscores(chainmode,
                                chain,
                                fragmentinfotable,
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
    chainingboundarycases(chainmode, chain, fragmentinfotable);
    if (chainprocessor(cpinfo,fragmentinfotable,chain,err))
    {
      haserr = true;
    }
    retcode = 0;
  }
  return haserr ? -1 : retcode;
}

static int cmpFragmentinfo0(const void *keya,const void *keyb)
{
  if (((const Fragmentinfo *) keya)->startpos[0] <
      ((const Fragmentinfo *) keyb)->startpos[0])
  {
    return -1;
  }
  if (((const Fragmentinfo *) keya)->startpos[0] >
      ((const Fragmentinfo *) keyb)->startpos[0])
  {
    return 1;
  }
  return 0;
}

static int cmpFragmentinfo1(const void *keya,const void *keyb)
{
  if (((const Fragmentinfo *) keya)->startpos[1] <
      ((const Fragmentinfo *) keyb)->startpos[1])
  {
    return -1;
  }
  if (((const Fragmentinfo *) keya)->startpos[1] >
      ((const Fragmentinfo *) keyb)->startpos[1])
  {
    return 1;
  }
  return 0;
}

typedef int (*Qsortcomparefunction)(const void *,const void *);

void gt_chain_possiblysortopenformatfragments(
                             Verboseinfo *verboseinfo,
                             GtFragmentinfotable *fragmentinfotable,
                             unsigned int presortdim)
{
  if (fragmentinfotable->nextfree > 1UL)
  {
    Qsortcomparefunction qsortcomparefunction;
    Fragmentinfo *fptr;
    bool fragmentsaresorted = true;

    gt_assert(presortdim <= 1U);
    qsortcomparefunction
      = (presortdim == 0) ? cmpFragmentinfo0 : cmpFragmentinfo1;
    for (fptr = fragmentinfotable->fragments;
         fptr < fragmentinfotable->fragments + fragmentinfotable->nextfree - 1;
         fptr++)
    {
      if (qsortcomparefunction((const void *) fptr,
                               (const void *) (fptr+1)) == 1)
      {
        fragmentsaresorted = false;
        break;
      }
    }
    if (!fragmentsaresorted)
    {
      showverbose(verboseinfo,"input fragments are not yet sorted => "
                              "sort them");
    }
    qsort(fragmentinfotable->fragments,(size_t) fragmentinfotable->nextfree,
          sizeof (Fragmentinfo),qsortcomparefunction);
  }
}

static int parselocalchainingparameter(GtChainmode *gtchainmode,
                                       const char *option,
                                       const char *lparam,
                                       GtError *err)
{
  Qualifiedinteger *qualint;

  qualint = parsequalifiedinteger(option,lparam,err);
  if (qualint == NULL)
  {
    return -1;
  }
  switch (qualint->qualtag)
  {
    case Qualbestof:
      gtchainmode->chainkind = LOCALCHAININGBEST;
      gtchainmode->howmanybest = qualint->integervalue;
      break;
    case Qualpercentaway:
      gtchainmode->chainkind = LOCALCHAININGPERCENTAWAY;
      gtchainmode->percentawayfrombest = qualint->integervalue;
      break;
    case Qualabsolute:
      gtchainmode->chainkind = LOCALCHAININGTHRESHOLD;
      gtchainmode->minimumscore = (GtChainscoretype) qualint->integervalue;
      break;
  }
  gt_free(qualint);
  return 0;
}

/*
  The following string is used to trigger the usage of gap costs
  for global chaining.
*/

#define GAPCOSTSWITCH        "gc"

/*
  The following string is used to trigger the use of a chaining algorithm
  allowing for overlaps between the hits.
*/

#define OVERLAPSWITCH        "ov"

static int parseglobalchainingparameter(GtChainmode *chainmode,
                                        const char *option,
                                        const char *gparam,
                                        GtError *err)
{
  if (strcmp(gparam,GAPCOSTSWITCH) == 0)
  {
    chainmode->chainkind = GLOBALCHAININGWITHGAPCOST;
    return 0;
  }
  if (strcmp(gparam,OVERLAPSWITCH) == 0)
  {
    chainmode->chainkind = GLOBALCHAININGWITHOVERLAPS;
    return 0;
  }
  gt_error_set(err,"argument of option -%s must be %s or %s: ",
                    option,
                    GAPCOSTSWITCH,
                    OVERLAPSWITCH);
  return -1;
}

GtChainmode *gt_chain_chainmode_new(bool weightfactorset,
                                    unsigned long maxgap,
                                    bool globalset,
                                    const GtStrArray *globalargs,
                                    bool localset,
                                    const GtStrArray *localargs,
                                    GtError *err)
{
  GtChainmode *gtchainmode = gt_malloc(sizeof( GtChainmode));
  bool haserr = false;

  gt_assert(!(globalset && localset));
  gtchainmode->chainkind = GLOBALCHAINING;
  gtchainmode->maxgapwidth = (GtChainpostype) maxgap;
  if (localset)
  {
    unsigned long localargsnum = gt_str_array_size(localargs);
    if (localargsnum == 0)
    {
      gtchainmode->chainkind = LOCALCHAININGMAX;
    } else
    {
      if (localargsnum == 1UL)
      {
        if (parselocalchainingparameter(gtchainmode,
                                        "local",
                                        gt_str_array_get(localargs,0),
                                        err) != 0)
        {
          haserr = true;
        }
      } else
      {
        gt_error_set(err,"option -local can only have one optional argument");
        haserr = true;
      }
    }
  }
  if (globalset)
  {
    unsigned long globalargsnum = gt_str_array_size(globalargs);
    if (globalargsnum == 0)
    {
      gtchainmode->chainkind = GLOBALCHAINING;
    } else
    {
      if (globalargsnum == 1UL)
      {
        if (parseglobalchainingparameter(gtchainmode,
                                         "global",
                                          gt_str_array_get(globalargs,0),
                                          err) != 0)
        {
          haserr = true;
        }
      } else
      {
        gt_error_set(err,"option -global can only have one optional argument");
        haserr = true;
      }
    }
  }
  if (!haserr && weightfactorset)
  {
    if (!localset &&
        gtchainmode->chainkind != GLOBALCHAININGWITHGAPCOST &&
        gtchainmode->chainkind != GLOBALCHAININGWITHOVERLAPS)
    {
      gt_error_set(err,
                   "option wf requires either option -local or option -global "
                   "with argument %s or %s",
                   GAPCOSTSWITCH,
                   OVERLAPSWITCH);
      haserr = true;
    }
  }
  if (haserr)
  {
    gt_free(gtchainmode);
    return NULL;
  }
  return gtchainmode;
}

void gt_chain_chainmode_delete(GtChainmode *gtchainmode)
{
  if (gtchainmode != NULL)
  {
    gt_free(gtchainmode);
  }
}

GtChain *gt_chain_chain_new(void)
{
  GtChain *chain;

  chain = gt_malloc(sizeof (GtChain));
  GT_INITARRAY(&chain->chainedfragments,GtChainref);
  return chain;
}

void gt_chain_chain_delete(GtChain *chain)
{
  if (chain != NULL)
  {
    GT_FREEARRAY(&chain->chainedfragments,GtChainref);
    gt_free(chain);
  }
}

static int gt_outputformatchaingeneric(
                                bool silent,
                                GT_UNUSED void *data,
                                const GtFragmentinfotable *fragmentinfotable,
                                const GtChain *chain,
                                GT_UNUSED GtError *err)
{
  const Fragmentinfo *fiptr;
  unsigned long chaincounter = 0, i;

  printf("# chain %lu: length %lu score %ld\n",
         chaincounter,
         chain->chainedfragments.nextfreeGtChainref,
         chain->scoreofchain);
  if (!silent)
  {
    for (i=0; i < chain->chainedfragments.nextfreeGtChainref; i++)
    {
      fiptr = fragmentinfotable->fragments +
              chain->chainedfragments.spaceGtChainref[i];
      printf(FormatSeqpos " " FormatSeqpos " " FormatSeqpos " " FormatSeqpos
             " %ld\n",PRINTSeqposcast(fiptr->startpos[0]),
                      PRINTSeqposcast(fiptr->endpos[0]),
                      PRINTSeqposcast(fiptr->startpos[1]),
                      PRINTSeqposcast(fiptr->endpos[1]),
                      fiptr->weight);
    }
  }
  chaincounter++;
  return 0;
}

int gt_outputformatchainsilent(void *data,
                               const GtFragmentinfotable *fragmentinfotable,
                               const GtChain *chain,
                               GT_UNUSED GtError *err)
{
  return gt_outputformatchaingeneric(true,
                                     data,
                                     fragmentinfotable,
                                     chain,
                                     err);
}

int gt_outputformatchain(void *data,
                         const GtFragmentinfotable *fragmentinfotable,
                         const GtChain *chain,
                         GT_UNUSED GtError *err)
{
  return gt_outputformatchaingeneric(false,
                                     data,
                                     fragmentinfotable,
                                     chain,
                                     err);
}
