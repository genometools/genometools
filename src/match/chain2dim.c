#include <stdbool.h>
#include "core/minmax.h"
#include "core/unused_api.h"
#include "extended/redblack.h"
#include "chaindef.h"

/*
  The basic information required for each fragment is stored
  in a structure of the following type. The user has to specify
  those components which a tagged `user defined'. The chaining
  algorithms computes the remaining components score and previous
  in chain.
*/

typedef struct
{
  Seqpos startpos[2],  /* start of fragments in the 2 dimensions, userdef */
         endpos[2];    /* end of fragments in the 2 dimensions, userdef */
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
    = gt_malloc(sizeof(*fragmentinfotable));
  fragmentinfotable->fragments
    = gt_malloc(sizeof(*fragmentinfotable->fragments) * numberoffragments);
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
                           Seqpos start1,Seqpos end1,
                           Seqpos start2,Seqpos end2,
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
                      sizeof (*chain->chainedfragments.spaceGtChainref *\
                              lengthofchain));\
          chain->chainedfragments.allocatedGtChainref = lengthofchain;\
          chain->chainedfragments.nextfreeGtChainref = 0;\
        }

typedef GtChainscoretype (*GtChaingapcostfunction)(Fragmentinfo *,
                                                   unsigned long,unsigned long);

typedef struct
{
  unsigned long fpident;
  Seqpos fpposition;
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

static bool overlappingfragments(Fragmentinfo *fragmentinfo,unsigned long i,
                                 unsigned long j)
{
  return (GETSTOREDENDPOINT(0,i) >= GETSTOREDSTARTPOINT(0,j) ||
          GETSTOREDENDPOINT(1,i) >= GETSTOREDSTARTPOINT(1,j)) ? true : false;
}

static bool colinearfragments(Fragmentinfo *fragmentinfo,unsigned long i,
                              unsigned long j)
{
  return (GETSTOREDSTARTPOINT(0, i) < GETSTOREDSTARTPOINT(0, j) &&
          GETSTOREDENDPOINT(0, i)   < GETSTOREDENDPOINT(0, j)   &&
          GETSTOREDSTARTPOINT(1, i) < GETSTOREDSTARTPOINT(1, j) &&
          GETSTOREDENDPOINT(1, i)   < GETSTOREDENDPOINT(1, j)) ? true : false;
}

static GtChainscoretype gapcostL1(Fragmentinfo *fragmentinfo,unsigned long i,
                                  unsigned long j)
{
  return (GtChainscoretype)
         ((GETSTOREDSTARTPOINT(0,j) - GETSTOREDENDPOINT(0,i)) +
          (GETSTOREDSTARTPOINT(1,j) - GETSTOREDENDPOINT(1,i)));
}

static GtChainscoretype overlapcost(Fragmentinfo *fragmentinfo,
                                    unsigned long i,
                                    unsigned long j)
{
  Seqpos overlaplength = 0;

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

static GtChainscoretype gapcostCc(Fragmentinfo *fragmentinfo,
                                  unsigned long i,unsigned long j)
{
  Seqpos value1, value2;

  gt_assert(GETSTOREDSTARTPOINT(0,j) > GETSTOREDENDPOINT(0,i) &&
            GETSTOREDSTARTPOINT(1,j) > GETSTOREDENDPOINT(1,i));
  value1 = GETSTOREDSTARTPOINT(0,j) - GETSTOREDENDPOINT(0,i) - 1,
  value2 = GETSTOREDSTARTPOINT(1,j) - GETSTOREDENDPOINT(1,i) - 1;
  return (GtChainscoretype) MAX(value1,value2);
}

static void chainingboundarycases(GtChainmode *chainmode,
                                  GtChain *chain,
                                  Fragmentinfo *fragmentinfo,
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
                                   Fragmentinfo *fragmentinfo,
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

static bool checkmaxgapwidth(Fragmentinfo *fragmentinfo,
                             Seqpos maxgapwidth,
                             unsigned long leftfrag,
                             unsigned long rightfrag)
{
  Seqpos gapwidth, startpoint, endpoint;

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

static void bruteforcechainingscores(GtChainmode *chainmode,
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
                           Fragmentinfo *fragmentinfo,
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
                              Fragmentinfo *fragmentinfo,
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

static void evalfragmentscore(GtChainmode *chainmode,
                              Fragmentinfo *fragmentinfo,
                              unsigned long numofmatches,
                              Fragmentstore *fragmentstore,
                              bool gapsL1,
                              unsigned long fragpointident,
                              int presortdim)
{
  unsigned long previous;
  Seqpos startpos2;
  Fragpoint *qfrag2;
  GtChainscoretype score;

  gt_assert(presortdim == 0 || presortdim == 1);
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

static bool isrightmaximallocalchain(Fragmentinfo *fragmentinfo,
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
