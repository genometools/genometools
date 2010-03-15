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
#include "extended/redblack.h"
#include "core/logger.h"
#include "chain2dim.h"
#include "prsqualint.h"

/*
  The basic information required for each match is stored
  in a structure of the following type. The user has to specify
  those components which are tagged `user defined'. The chaining
  algorithms computes the remaining components tagged 'compute'
  in chain.
*/

typedef struct
{
  GtChainpostype startpos[2], /* start of matches in the 2 dimensions,
                                 userdef */
                 endpos[2];  /* end of matches in the 2 dimensions, userdef */
  unsigned long firstinchain,   /* first element in chain, compute */
                previousinchain;  /* previous index in chain, compute */
  GtChainscoretype
         weight, /* weight of match, user defined */
         initialgap, /* gap to start of sequences, user defined */
         terminalgap, /* gap to last positions of match, user defined */
         score; /* score of highest scoreing chain ending here, compute */
} Matchchaininfo;

struct GtChainmatchtable
{
   Matchchaininfo *matches;
   GtChainscoretype largestdim0,
                    largestdim1;
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
  A chain consists of an array of references to chained matches.
  These refer to an array of matchchaininfos.
*/

struct GtChain
{
  GtArrayGtChainref chainedmatches;
  GtChainscoretype scoreofchain;
};

GtChainmatchtable *gt_chain_matchtable_new(unsigned long numberofmatches)
{
  GtChainmatchtable *matchtable
    = gt_malloc(sizeof (*matchtable));
  matchtable->matches
    = gt_malloc(sizeof (*matchtable->matches) * numberofmatches);
  matchtable->nextfree = 0;
  matchtable->allocated = numberofmatches;
  matchtable->largestdim0 = matchtable->largestdim1 = 0;
  return matchtable;
}

void gt_chain_matchtable_delete(GtChainmatchtable *matchtable)
{
  if (matchtable != NULL)
  {
    gt_free(matchtable->matches);
    gt_free(matchtable);
  }
}

void gt_chain_matchtable_empty(GtChainmatchtable *matchtable)
{
  gt_assert(matchtable != NULL);
  matchtable->largestdim0 = matchtable->largestdim1 = 0;
  matchtable->nextfree = 0;
}

void gt_chain_matchtable_add(GtChainmatchtable *matchtable,
                             const GtChainmatchvalues *inmatch)
{
  Matchchaininfo *matchchaininfo;

  gt_assert(matchtable->nextfree < matchtable->allocated);
  gt_assert(inmatch->startpos[0] <= inmatch->endpos[0]);
  gt_assert(inmatch->startpos[1] <= inmatch->endpos[1]);
  matchchaininfo = matchtable->matches + matchtable->nextfree++;
  matchchaininfo->startpos[0] = inmatch->startpos[0];
  matchchaininfo->startpos[1] = inmatch->startpos[1];
  matchchaininfo->endpos[0] = inmatch->endpos[0];
  matchchaininfo->endpos[1] = inmatch->endpos[1];
  matchchaininfo->weight = inmatch->weight;
  if (matchtable->largestdim0 < (GtChainscoretype) inmatch->endpos[0])
  {
    matchtable->largestdim0 = (GtChainscoretype) inmatch->endpos[0];
  }
  if (matchtable->largestdim1 < (GtChainscoretype) inmatch->endpos[1])
  {
    matchtable->largestdim1 = (GtChainscoretype) inmatch->endpos[1];
  }
}

void gt_chain_fillthegapvalues(GtChainmatchtable *matchtable)
{
  Matchchaininfo *fiptr;

  for (fiptr = matchtable->matches;
       fiptr < matchtable->matches + matchtable->nextfree;
       fiptr++)
  {
    fiptr->initialgap
      = (GtChainscoretype) (fiptr->startpos[0] + fiptr->startpos[1]);
    fiptr->terminalgap
      = (GtChainscoretype) (matchtable->largestdim0 - fiptr->endpos[0] +
                            matchtable->largestdim1 - fiptr->endpos[1]);
  }
}

void gt_chain_applyweight(double weightfactor,GtChainmatchtable *matchtable)
{
  if (!gt_double_equals_double(weightfactor, 1.0))
  {
    Matchchaininfo *fiptr;

    for (fiptr = matchtable->matches;
         fiptr < matchtable->matches + matchtable->nextfree;
         fiptr++)
    {
      fiptr->weight *= weightfactor;
    }
  }
}

#define MAKEENDPOINT(FID)       (FID)
#define FRAGIDENT(FRAG)         ((FRAG)->fpident)

#define UNDEFPREVIOUS           matchtable->nextfree

#define GETSTOREDSTARTPOINT(DIM,IDX)\
        matchtable->matches[IDX].startpos[DIM]
#define GETSTOREDENDPOINT(DIM,IDX)\
        matchtable->matches[IDX].endpos[DIM]
#define INITIALGAP(IDX)\
        matchtable->matches[IDX].initialgap
#define TERMINALGAP(IDX)\
        matchtable->matches[IDX].terminalgap

#define CHECKCHAINSPACE\
        if (lengthofchain >= chain->chainedmatches.allocatedGtChainref)\
        {\
         chain->chainedmatches.spaceGtChainref = \
           gt_realloc(chain->chainedmatches.spaceGtChainref,\
                      sizeof (*chain->chainedmatches.spaceGtChainref) *\
                              lengthofchain);\
          chain->chainedmatches.allocatedGtChainref = lengthofchain;\
          chain->chainedmatches.nextfreeGtChainref = 0;\
        }

typedef GtChainscoretype (*GtChaingapcostfunction)(const GtChainmatchtable *,
                                                   unsigned long,unsigned long);

typedef struct
{
  unsigned long fpident;
  GtChainpostype fpposition;
} Matchpoint;

typedef struct
{
  GtRBTnode *dictroot;
  unsigned long *endpointperm;
} Matchstore;

typedef struct
{
  GtChainscoretype maxscore;
  unsigned long maxmatchnum;
  bool defined;
} Maxmatchvalue;

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

static bool overlappingmatches(const GtChainmatchtable *matchtable,
                               unsigned long i,
                               unsigned long j)
{
  return (GETSTOREDENDPOINT(0,i) >= GETSTOREDSTARTPOINT(0,j) ||
          GETSTOREDENDPOINT(1,i) >= GETSTOREDSTARTPOINT(1,j)) ? true : false;
}

static bool colinearmatches(const GtChainmatchtable *matchtable,
                            unsigned long i,
                            unsigned long j)
{
  return (GETSTOREDSTARTPOINT(0, i) < GETSTOREDSTARTPOINT(0, j) &&
          GETSTOREDENDPOINT(0, i)   < GETSTOREDENDPOINT(0, j)   &&
          GETSTOREDSTARTPOINT(1, i) < GETSTOREDSTARTPOINT(1, j) &&
          GETSTOREDENDPOINT(1, i)   < GETSTOREDENDPOINT(1, j)) ? true : false;
}

static GtChainscoretype gapcostL1(const GtChainmatchtable *matchtable,
                                  unsigned long i,
                                  unsigned long j)
{
  return (GtChainscoretype)
         ((GETSTOREDSTARTPOINT(0,j) - GETSTOREDENDPOINT(0,i)) +
          (GETSTOREDSTARTPOINT(1,j) - GETSTOREDENDPOINT(1,i)));
}

static GtChainscoretype overlapcost(const GtChainmatchtable *matchtable,
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

  whereas i+q-1 is the end point in i-dimension of match 1 and
  j+q-1 is the end point in j-dimension of match 1.
  In using the match specific end points, other than fixed values for q,
  following function generalizes to MEMs as well.
*/

static GtChainscoretype gapcostCc(const GtChainmatchtable *matchtable,
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
                                  const GtChainmatchtable *matchtable)
{
  if (matchtable->nextfree == 0)
  {
    chain->scoreofchain = 0;
    chain->chainedmatches.nextfreeGtChainref = 0;
  } else
  {
    if (matchtable->nextfree == 1UL)
    {
      unsigned long lengthofchain;

      chain->scoreofchain = matchtable->matches[0].weight;
      if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
      {
        chain->scoreofchain -= (INITIALGAP(0) + TERMINALGAP(0));
      }
      lengthofchain = 1UL;
      CHECKCHAINSPACE;
      chain->chainedmatches.spaceGtChainref[0] = 0;
      chain->chainedmatches.nextfreeGtChainref = 1UL;
    }
  }
}

static void retracepreviousinchain(GtChain *chain,
                                   const GtChainmatchtable *matchtable,
                                   unsigned long retracestart)
{
  unsigned long matchnum, idx, lengthofchain;

  for (lengthofchain = 0, matchnum = retracestart;
       matchnum != UNDEFPREVIOUS; lengthofchain++)
  {
    matchnum = matchtable->matches[matchnum].previousinchain;
  }
  CHECKCHAINSPACE
  matchnum = retracestart;
  idx = lengthofchain;
  while (matchnum != UNDEFPREVIOUS)
  {
    gt_assert(idx > 0);
    idx--;
    chain->chainedmatches.spaceGtChainref[idx] = matchnum;
    matchnum = matchtable->matches[matchnum].previousinchain;
  }
  gt_assert(idx == 0);
  chain->chainedmatches.nextfreeGtChainref = lengthofchain;
}

static bool checkmaxgapwidth(const GtChainmatchtable *matchtable,
                             GtChainpostype maxgapwidth,
                             unsigned long leftmatch,
                             unsigned long rightmatch)
{
  GtChainpostype gapwidth, startpoint, endpoint;

  startpoint = GETSTOREDSTARTPOINT(0,rightmatch);
  endpoint = GETSTOREDENDPOINT(0,leftmatch);
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
  startpoint = GETSTOREDSTARTPOINT(1,rightmatch);
  endpoint = GETSTOREDENDPOINT(1,leftmatch);
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
                                     GtChainmatchtable *matchtable,
                                     GtChaingapcostfunction
                                       chaingapcostfunction)
{
  unsigned long previous, leftmatch, rightmatch;
  GtChainscoretype weightright, score;
  Maxmatchvalue localmaxmatch;
  bool combinable;

  if (matchtable->nextfree > 1UL)
  {
    matchtable->matches[0].firstinchain = 0;
    matchtable->matches[0].previousinchain = UNDEFPREVIOUS;
    matchtable->matches[0].score
      = matchtable->matches[0].weight;
    if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
    {
      matchtable->matches[0].score -= (INITIALGAP(0) + TERMINALGAP(0));
    }
    for (rightmatch=1Ul; rightmatch<matchtable->nextfree; rightmatch++)
    {
      weightright = matchtable->matches[rightmatch].weight;
      localmaxmatch.defined = false;
      localmaxmatch.maxscore = 0;
      localmaxmatch.maxmatchnum = 0;
      for (leftmatch=0; leftmatch<rightmatch; leftmatch++)
      {
        if (chainmode->maxgapwidth != 0 &&
           !checkmaxgapwidth(matchtable,
                             chainmode->maxgapwidth,
                             leftmatch,
                             rightmatch))
        {
          combinable = false;
        } else
        {
          if (chainmode->chainkind == GLOBALCHAININGWITHOVERLAPS)
          {
            combinable = colinearmatches(matchtable,leftmatch,
                                           rightmatch);
          } else
          {
            if (overlappingmatches(matchtable,leftmatch,rightmatch))
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
          score = matchtable->matches[leftmatch].score;
          if (chainmode->chainkind == GLOBALCHAINING)
          {
            /* process chainkinds without gap costs */
            score += weightright;
            previous = leftmatch;
          } else
          {
            score -= chaingapcostfunction(matchtable,leftmatch,rightmatch);
            if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
            {
              score += (weightright + TERMINALGAP(leftmatch)
                                    - TERMINALGAP(rightmatch));
              previous = leftmatch;
            } else
            {
              if (score > 0)
              {
                score += weightright;
                previous = leftmatch;
              } else
              {
                score = weightright;
                previous = UNDEFPREVIOUS;
              }
            }
          }
          if (!localmaxmatch.defined || localmaxmatch.maxscore < score)
          {
            localmaxmatch.maxscore = score;
            localmaxmatch.maxmatchnum = previous;
            localmaxmatch.defined = true;
          }
        }
      }
      if (localmaxmatch.defined)
      {
        matchtable->matches[rightmatch].previousinchain
          = localmaxmatch.maxmatchnum;
        if (localmaxmatch.maxmatchnum == UNDEFPREVIOUS)
        {
          matchtable->matches[rightmatch].firstinchain = rightmatch;
        } else
        {
          matchtable->matches[rightmatch].firstinchain
            = matchtable->matches[localmaxmatch.maxmatchnum].
              firstinchain;
        }
        matchtable->matches[rightmatch].score = localmaxmatch.maxscore;
      } else
      {
        matchtable->matches[rightmatch].previousinchain = UNDEFPREVIOUS;
        matchtable->matches[rightmatch].firstinchain = rightmatch;
        matchtable->matches[rightmatch].score = weightright;
        if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
        {
          matchtable->matches[rightmatch].score
            -= (INITIALGAP(rightmatch) + TERMINALGAP(rightmatch));
        }
      }
    }
  }
}

/*
  The following function compares match points. These are
  end points of matches in dimension 2.
  If the match points are identical, then the order is undefined.
*/

static int cmpendMatchpoint2(const GtKeytype keya,
                             const GtKeytype keyb,
                             GT_UNUSED void *info)
{
  if (((Matchpoint *) keya)->fpposition < ((Matchpoint *) keyb)->fpposition)
  {
    return -1;
  }
  if (((Matchpoint *) keya)->fpposition > ((Matchpoint *) keyb)->fpposition)
  {
    return 1;
  }
  if (FRAGIDENT((Matchpoint *) keya) < FRAGIDENT((Matchpoint *) keyb))
  {
    return -1;
  }
  if (FRAGIDENT((Matchpoint *) keya) > FRAGIDENT((Matchpoint *) keyb))
  {
    return 1;
  }
  return 0;
}

static GtChainscoretype evalpriority(bool addterminal,
                                     const GtChainmatchtable *matchtable,
                                     unsigned long matchnum)
{
  if (addterminal)
  {
    return matchtable->matches[matchnum].score - TERMINALGAP(matchnum);
  }
  return matchtable->matches[matchnum].score;
}

static void insertintodict(bool addterminal,
                           const GtChainmatchtable *matchtable,
                           Matchstore *matchstore,
                           Matchpoint *qmatch2)
{
  Matchpoint *retval2;
  bool nodecreated;

  retval2 = (Matchpoint *) gt_rbt_search ((const GtKeytype) qmatch2,
                                          &nodecreated,
                                          &matchstore->dictroot,
                                          cmpendMatchpoint2,
                                          NULL);
  gt_assert(retval2 != NULL);
  if (!nodecreated)
  {
    if (evalpriority(addterminal,matchtable,FRAGIDENT(retval2)) <
        evalpriority(addterminal,matchtable,FRAGIDENT(qmatch2)))
    {
      gt_assert(retval2->fpposition == qmatch2->fpposition);
      retval2->fpident = qmatch2->fpident;
    }
    gt_free(qmatch2);
  }
}

static void activatematchpoint(bool addterminal,
                               const GtChainmatchtable *matchtable,
                               Matchstore *matchstore,
                               Matchpoint *qmatch2)
{
  Matchpoint *tmp2;
  GtChainscoretype qpriority;

  qpriority = evalpriority(addterminal,matchtable,FRAGIDENT(qmatch2));
  tmp2 = (Matchpoint *) gt_rbt_previousequalkey ((const GtKeytype) qmatch2,
                                                matchstore->dictroot,
                                                cmpendMatchpoint2,
                                                NULL);
  if (tmp2 == NULL ||
      qpriority > evalpriority(addterminal,matchtable,FRAGIDENT(tmp2)))
  {
    insertintodict(addterminal,matchtable,matchstore,qmatch2);
    while (true)
    {
      tmp2 = (Matchpoint *) gt_rbt_nextkey ((const GtKeytype) qmatch2,
                                           matchstore->dictroot,
                                           cmpendMatchpoint2,
                                           NULL);
      if (tmp2 == NULL || qpriority <= evalpriority(addterminal,
                                                    matchtable,
                                                    FRAGIDENT(tmp2)))
      {
        break;
      }
      if (gt_rbt_delete ((const GtKeytype) tmp2,
                        &matchstore->dictroot,
                        cmpendMatchpoint2,
                        NULL) != 0)
      {
        fprintf(stderr,"cannot delete successor node\n");
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
      gt_free(tmp2);
    }
  } else
  {
    gt_free(qmatch2);
  }
}

static void evalmatchscore(const GtChainmode *chainmode,
                           GtChainmatchtable *matchtable,
                           Matchstore *matchstore,
                           bool gapsL1,
                           unsigned long matchpointident,
                           unsigned int presortdim)
{
  unsigned long previous;
  GtChainpostype startpos2;
  Matchpoint *qmatch2;
  GtChainscoretype score;

  startpos2 = GETSTOREDSTARTPOINT(1-presortdim,matchpointident);
  if (startpos2 == 0)
  {
    qmatch2 = NULL;
  } else
  {
    Matchpoint keymatch2;

    keymatch2.fpposition = startpos2 - 1;  /* it is a start position */
    keymatch2.fpident = MAKEENDPOINT(matchpointident);
                       /* but considered as endpoint */
    qmatch2 = (Matchpoint *) gt_rbt_previousequalkey(
                                  (const GtKeytype) &keymatch2,
                                  matchstore->dictroot,
                                  cmpendMatchpoint2,
                                  NULL);
    if (qmatch2 != NULL)
    {
      if (chainmode->maxgapwidth != 0 &&
         !checkmaxgapwidth(matchtable,
                           chainmode->maxgapwidth,
                           FRAGIDENT(qmatch2),
                           matchpointident))
      {
        qmatch2 = NULL;
      }
    }
  }
  if (qmatch2 == NULL)
  {
    score = matchtable->matches[matchpointident].weight;
    if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
    {
      score -= INITIALGAP(matchpointident);
    }
    previous = UNDEFPREVIOUS;
  } else
  {
    score = matchtable->matches[FRAGIDENT(qmatch2)].score;
    if (chainmode->chainkind == GLOBALCHAINING)
    {
      score += matchtable->matches[matchpointident].weight;
      previous = FRAGIDENT(qmatch2);
    } else
    {
      GtChainscoretype tmpgc;

      if (gapsL1)
      {
        tmpgc = gapcostL1(matchtable,FRAGIDENT(qmatch2),matchpointident);
      } else
      {
        tmpgc = gapcostCc(matchtable,FRAGIDENT(qmatch2),matchpointident);
      }
      if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST || score > tmpgc)
      {
        score += (matchtable->matches[matchpointident].weight - tmpgc);
        previous = FRAGIDENT(qmatch2);
      } else
      {
        score = matchtable->matches[matchpointident].weight;
        previous = UNDEFPREVIOUS;
      }
    }
  }
  matchtable->matches[matchpointident].score = score;
  matchtable->matches[matchpointident].previousinchain = previous;
  if (previous == UNDEFPREVIOUS)
  {
    matchtable->matches[matchpointident].firstinchain = matchpointident;
  } else
  {
    matchtable->matches[matchpointident].firstinchain
      = matchtable->matches[previous].firstinchain;
  }
}

static bool isrightmaximallocalchain(const GtChainmatchtable *matchtable,
                                     unsigned long currentmatch)
{
  if (currentmatch == matchtable->nextfree - 1)
  {
    return true;
  }
  if (matchtable->matches[currentmatch+1].previousinchain !=
      currentmatch)
  {
    return true;
  }
  if (matchtable->matches[currentmatch+1].score <
      matchtable->matches[currentmatch].score)
  {
    return true;
  }
  return false;
}

static void determineequivreps(Bestofclass *chainequivalenceclasses,
                               const GtChainmatchtable *matchtable)
{
  unsigned long matchnum;
  Bestofclass *classptr, *classrep;

  for (classptr = chainequivalenceclasses;
       classptr < chainequivalenceclasses + matchtable->nextfree;
       classptr++)
  {
    classptr->isavailable = false;
  }
  for (matchnum=0; matchnum<matchtable->nextfree; matchnum++)
  {
    if (isrightmaximallocalchain(matchtable,matchnum))
    {
      classrep = chainequivalenceclasses +
                 matchtable->matches[matchnum].firstinchain;
      if (!classrep->isavailable ||
          classrep->score < matchtable->matches[matchnum].score)
      {
        classrep->score = matchtable->matches[matchnum].score;
        classrep->isavailable = true;
      }
    }
  }
}

static bool retrievemaximalscore(GtChainscoretype *maxscore,
                                 const GtChainmode *chainmode,
                                 const GtChainmatchtable *matchtable)
{
  unsigned long matchnum;
  GtChainscoretype tgap;
  bool maxscoredefined = false;

  *maxscore = 0;
  for (matchnum=0; matchnum<matchtable->nextfree; matchnum++)
  {
    if (isrightmaximallocalchain(matchtable,matchnum))
    {
      if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
      {
        tgap = TERMINALGAP(matchnum);
      } else
      {
        tgap = 0;
      }
      if (!maxscoredefined ||
          *maxscore < matchtable->matches[matchnum].score - tgap)
      {
        *maxscore = matchtable->matches[matchnum].score - tgap;
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
                                    const GtChainmatchtable *matchtable,
                                    unsigned long howmanybest)
{
  unsigned long idx, matchnum = 0;
  GtChainscoretype *scores;
  Dictmaxsize *dictbestmatches;
  void *minkey;

  scores = gt_malloc(sizeof (*scores) * matchtable->nextfree);
  dictbestmatches = dictmaxsize_new(howmanybest);
  for (idx=0; idx<matchtable->nextfree; idx++)
  {
    if (isrightmaximallocalchain(matchtable,idx))
    {
      scores[matchnum] = matchtable->matches[idx].score;
      insertDictmaxsize(dictbestmatches,
                        comparescores,
                        NULL,
                        (void *) (scores + matchnum));
      matchnum++;
    }
  }
  if (matchnum == 0)
  {
    *minscoredefined = false;
  } else
  {
    minkey = gt_rbt_minimumkey(dictbestmatches->root);
    gt_assert(minkey != NULL);
    *minscore = *((GtChainscoretype *) minkey);
    *minscoredefined = true;
  }
  dictmaxsize_delete(dictbestmatches);
  gt_free(scores);
}

static void retrievechainthreshold(const GtChainmode *chainmode,
                                   GtChainmatchtable *matchtable,
                                   GtChain *chain,
                                   GtChainscoretype minscore,
                                   Bestofclass *chainequivalenceclasses,
                                   GtChainprocessor chainprocessor,
                                   void *cpinfo)
{
  unsigned long matchnum;
  GtChainscoretype tgap;
  Bestofclass *classrep;

  for (matchnum=0; matchnum < matchtable->nextfree; matchnum++)
  {
    if (isrightmaximallocalchain(matchtable,matchnum))
    {
      if (chainmode->chainkind == GLOBALCHAININGWITHGAPCOST)
      {
        tgap = TERMINALGAP(matchnum);
      } else
      {
        tgap = 0;
      }
      if (matchtable->matches[matchnum].score - tgap >= minscore)
      {
        if (chainequivalenceclasses != NULL)
        {
          classrep = chainequivalenceclasses +
                     matchtable->matches[matchnum].firstinchain;
          gt_assert(classrep != NULL);
          if (classrep->isavailable &&
              classrep->score == matchtable->matches[matchnum].score - tgap)
          {
            chain->scoreofchain = classrep->score;
            classrep->isavailable = false;
            retracepreviousinchain(chain, matchtable, matchnum);
            chainprocessor(cpinfo,matchtable,chain);
          }
        } else
        {
          chain->scoreofchain = matchtable->matches[matchnum].score - tgap;
          retracepreviousinchain(chain, matchtable, matchnum);
          chainprocessor(cpinfo,matchtable,chain);
        }
      }
    }
  }
}

static int comparestartandend(const Matchchaininfo *sortedstartpoints,
                              const Matchchaininfo *sortedendpoints,
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

static Matchpoint *makeactivationpoint(const GtChainmatchtable *matchtable,
                                       unsigned long fpident,
                                       unsigned int postsortdim)
{
  Matchpoint *fpptr;

  fpptr = gt_malloc(sizeof (*fpptr));
  gt_assert(fpptr != NULL);
  fpptr->fpident = MAKEENDPOINT(fpident);
  fpptr->fpposition = GETSTOREDENDPOINT(postsortdim,fpident);
  return fpptr;
}

static void mergestartandendpoints(const GtChainmode *chainmode,
                                   GtChainmatchtable *matchtable,
                                   Matchstore *matchstore,
                                   bool gapsL1,
                                   unsigned int presortdim)
{
  unsigned long xidx, startcount, endcount;
  bool addterminal;
  unsigned int postsortdim = 1U - presortdim;

  addterminal = (chainmode->chainkind == GLOBALCHAINING) ? false : true;
  matchstore->dictroot = NULL;
  for (xidx = 0, startcount = 0, endcount = 0;
       startcount < matchtable->nextfree &&
       endcount < matchtable->nextfree;
       xidx++)
  {
    if (comparestartandend(matchtable->matches + startcount,
                           matchtable->matches +
                           matchstore->endpointperm[endcount],
                           presortdim) < 0)
    {
      evalmatchscore(chainmode,
                     matchtable,
                     matchstore,
                     gapsL1,
                     startcount,
                     presortdim);
      startcount++;
    } else
    {
      activatematchpoint(addterminal,matchtable,matchstore,
                        makeactivationpoint(matchtable,
                                            matchstore->
                                            endpointperm[endcount],
                                            postsortdim));
      endcount++;
    }
  }
  while (startcount < matchtable->nextfree)
  {
    evalmatchscore(chainmode,
                   matchtable,
                   matchstore,
                   gapsL1,
                   startcount,
                   presortdim);
    startcount++;
    xidx++;
  }
  while (endcount < matchtable->nextfree)
  {
    activatematchpoint(addterminal,matchtable,matchstore,
                      makeactivationpoint(matchtable,
                                          matchstore->
                                          endpointperm[endcount],
                                          postsortdim));
    endcount++;
    xidx++;
  }
}

static unsigned int findmaximalscores(const GtChainmode *chainmode,
                                      GtChain *chain,
                                      GtChainmatchtable *matchtable,
                                      Matchstore *matchstore,
                                      GtChainprocessor chainprocessor,
                                      bool withequivclasses,
                                      void *cpinfo,
                                      GtLogger *logger)
{
  unsigned long matchnum;
  GtChainscoretype minscore = 0;
  Matchpoint *maxpoint;
  bool minscoredefined = false;
  Bestofclass *chainequivalenceclasses;
  unsigned int retval;

  if (withequivclasses)
  {
    if (chainmode->chainkind == LOCALCHAININGMAX ||
        chainmode->chainkind == LOCALCHAININGTHRESHOLD ||
        chainmode->chainkind == LOCALCHAININGBEST ||
        chainmode->chainkind == LOCALCHAININGPERCENTAWAY)
    {
      chainequivalenceclasses = gt_malloc(sizeof (*chainequivalenceclasses) *
                                          matchtable->nextfree);
      determineequivreps(chainequivalenceclasses, matchtable);
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
      maxpoint = gt_rbt_maximumkey(matchstore->dictroot);
      gt_assert(maxpoint != NULL);
      matchnum = FRAGIDENT(maxpoint);
      minscore = matchtable->matches[matchnum].score;
      minscoredefined = true;
      break;
    case GLOBALCHAININGWITHGAPCOST:
    case GLOBALCHAININGWITHOVERLAPS:
    case LOCALCHAININGMAX:
      minscoredefined
        = retrievemaximalscore(&minscore,chainmode,matchtable);
      break;
    case LOCALCHAININGTHRESHOLD:
      minscore = chainmode->minimumscore;
      minscoredefined = true;
      break;
    case LOCALCHAININGBEST:
      retrievechainbestscores(&minscoredefined,
                              &minscore,
                              matchtable,
                              chainmode->howmanybest);
      break;
    case LOCALCHAININGPERCENTAWAY:
      minscoredefined
        = retrievemaximalscore(&minscore,chainmode,matchtable);
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
    gt_logger_log(logger,
               "compute optimal %s chains with score >= %ld",
               (chainmode->chainkind == GLOBALCHAINING ||
                chainmode->chainkind == GLOBALCHAININGWITHGAPCOST ||
                chainmode->chainkind == GLOBALCHAININGWITHOVERLAPS)
                ? "global"
                : "local",
               minscore);
    retrievechainthreshold(chainmode,
                           matchtable,
                           chain,
                           minscore,
                           chainequivalenceclasses,
                           chainprocessor,
                           cpinfo);
    retval = 0;
  } else
  {
    retval = 1U;
  }
  if (chainequivalenceclasses != NULL)
  {
    gt_free(chainequivalenceclasses);
  }
  return retval;
}

static void makesortedendpointpermutation(unsigned long *perm,
                                          GtChainmatchtable *matchtable,
                                          unsigned int presortdim)
{
  unsigned long temp, *iptr, *jptr, i, moves = 0;

  for (i = 0; i < matchtable->nextfree; i++)
  {
    perm[i] = i;
  }
  for (iptr = perm + 1UL; iptr < perm + matchtable->nextfree; iptr++)
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
                               GtChainmatchtable *matchtable,
                               Matchstore *matchstore,
                               unsigned int presortdim,
                               bool gapsL1)
{
  matchstore->endpointperm
    = gt_malloc(sizeof (*matchstore->endpointperm) *
                matchtable->nextfree);
  makesortedendpointpermutation(matchstore->endpointperm,
                                matchtable,
                                presortdim);
  mergestartandendpoints(chainmode,
                         matchtable,
                         matchstore,
                         gapsL1,
                         presortdim);
  gt_free(matchstore->endpointperm);
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
  each specified in an array \texttt{matchtable} points to.
  \texttt{chain} is the array in which a chain is stored. However,
  do not further process chain. Just use
  INITARRAY(&chain.chainedmatches,Uint) to initialize it before calling
  \textttt{fastchaining} and
  FREEARRAY(&chain.chainedmatches,Uint) to free the space after calling
  \textttt{fastchaining}. The function \texttt{chainprocessor} processes
  each chain found, and \texttt{cpinfo} is used as a the first argument
  in each call to \texttt{chainprocessor}. The function returns
  0 upon success and otherwise, if an error occurs.
  Finally \texttt{gt_logger_log} is applied to each status message generated
  during the execution of \texttt{fastchaining}. If \texttt{gt_logger_log}
  is \texttt{NULL}, then nothing is generated and shown.
*/

void gt_chain_fastchaining(const GtChainmode *chainmode,
                           GtChain *chain,
                           GtChainmatchtable *matchtable,
                           bool gapsL1,
                           unsigned int presortdim,
                           bool withequivclasses,
                           GtChainprocessor chainprocessor,
                           void *cpinfo,
                           GtLogger *logger)
{
  unsigned int retval;
  GtChaingapcostfunction chaingapcostfunction;

  gt_assert(presortdim <= 1U);
  chaingapcostfunction
    = assignchaingapcostfunction(chainmode->chainkind,gapsL1);
  if (matchtable->nextfree > 1UL)
  {
    Matchstore matchstore;

    gt_logger_log(logger,"compute chain scores");
    if (chainmode->chainkind == GLOBALCHAININGWITHOVERLAPS)
    {
      bruteforcechainingscores(chainmode,matchtable,
                               chaingapcostfunction);
    } else
    {
      fastchainingscores(chainmode,
                         matchtable,
                         &matchstore,
                         presortdim,
                         gapsL1);
    }
    gt_logger_log(logger,"retrieve optimal chains");
    retval = findmaximalscores(chainmode,
                               chain,
                               matchtable,
                               &matchstore,
                               chainprocessor,
                               withequivclasses,
                               cpinfo,
                               logger);
    if (chainmode->chainkind != GLOBALCHAININGWITHOVERLAPS)
    {
      gt_rbt_destroy (true,NULL,NULL,matchstore.dictroot);
    }
  } else
  {
    chainingboundarycases(chainmode, chain, matchtable);
    chainprocessor(cpinfo,matchtable,chain);
    retval = 0;
  }
  /* retval is not reported. */
}

static int cmpMatchchaininfo0(const void *keya,const void *keyb)
{
  if (((const Matchchaininfo *) keya)->startpos[0] <
      ((const Matchchaininfo *) keyb)->startpos[0])
  {
    return -1;
  }
  if (((const Matchchaininfo *) keya)->startpos[0] >
      ((const Matchchaininfo *) keyb)->startpos[0])
  {
    return 1;
  }
  return 0;
}

static int cmpMatchchaininfo1(const void *keya,const void *keyb)
{
  if (((const Matchchaininfo *) keya)->startpos[1] <
      ((const Matchchaininfo *) keyb)->startpos[1])
  {
    return -1;
  }
  if (((const Matchchaininfo *) keya)->startpos[1] >
      ((const Matchchaininfo *) keyb)->startpos[1])
  {
    return 1;
  }
  return 0;
}

typedef int (*Qsortcomparefunction)(const void *,const void *);

void gt_chain_possiblysortmatches(GtLogger *logger,
                                    GtChainmatchtable *matchtable,
                                    unsigned int presortdim)
{
  if (matchtable->nextfree > 1UL)
  {
    Qsortcomparefunction qsortcomparefunction;
    Matchchaininfo *fptr;
    bool matchesaresorted = true;

    gt_assert(presortdim <= 1U);
    qsortcomparefunction
      = (presortdim == 0) ? cmpMatchchaininfo0 : cmpMatchchaininfo1;
    for (fptr = matchtable->matches;
         fptr < matchtable->matches + matchtable->nextfree - 1;
         fptr++)
    {
      if (qsortcomparefunction((const void *) fptr,
                               (const void *) (fptr+1)) == 1)
      {
        matchesaresorted = false;
        break;
      }
    }
    if (!matchesaresorted)
    {
      gt_logger_log(logger,"input matches are not yet sorted => "
                              "sort them");
    }
    qsort(matchtable->matches,(size_t) matchtable->nextfree,
          sizeof (*matchtable->matches),qsortcomparefunction);
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
  if (err != NULL)
  {
    gt_error_set(err,"argument of option -%s must be %s or %s: ",
                      option,
                      GAPCOSTSWITCH,
                      OVERLAPSWITCH);
  } else
  {
    fprintf(stderr,"argument of option -%s must be %s or %s: ",
                   option,
                   GAPCOSTSWITCH,
                   OVERLAPSWITCH);
  }
  return -1;
}

GtChainmode *gt_chain_chainmode_new(unsigned long maxgap,
                                    bool globalset,
                                    const char *globalargs,
                                    bool localset,
                                    const char *localargs,
                                    GtError *err)
{
  GtChainmode *gtchainmode;
  bool haserr = false;

  gt_assert(!(globalset && localset));
  gtchainmode = gt_malloc(sizeof (*gtchainmode));
  gtchainmode->chainkind = GLOBALCHAINING;
  gtchainmode->maxgapwidth = (GtChainpostype) maxgap;
  if (localset)
  {
    if (localargs == NULL)
    {
      gtchainmode->chainkind = LOCALCHAININGMAX;
    } else
    {
      if (parselocalchainingparameter(gtchainmode,
                                      "local",
                                      localargs,
                                      err) != 0)
      {
        haserr = true;
      }
    }
  }
  if (globalset)
  {
    if (globalargs == NULL)
    {
      gtchainmode->chainkind = GLOBALCHAINING;
    } else
    {
      if (parseglobalchainingparameter(gtchainmode,
                                       "global",
                                        globalargs,
                                        err) != 0)
      {
        haserr = true;
      }
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

  chain = gt_malloc(sizeof (*chain));
  GT_INITARRAY(&chain->chainedmatches,GtChainref);
  return chain;
}

void gt_chain_chain_delete(GtChain *chain)
{
  if (chain != NULL)
  {
    GT_FREEARRAY(&chain->chainedmatches,GtChainref);
    gt_free(chain);
  }
}

GtChainscoretype gt_chain_chainscore(const GtChain *chain)
{
  return chain->scoreofchain;
}

unsigned long gt_chain_chainlength(const GtChain *chain)
{
  return chain->chainedmatches.nextfreeGtChainref;
}

void gt_chain_extractchainelem(GtChainmatchvalues *value,
                               const GtChainmatchtable *matchtable,
                               const GtChain *chain,unsigned long idx)
{
  const Matchchaininfo *fiptr;

  gt_assert(idx <  gt_chain_chainlength(chain));
  fiptr = matchtable->matches +
          chain->chainedmatches.spaceGtChainref[idx];
  value->startpos[0] = fiptr->startpos[0];
  value->startpos[1] = fiptr->startpos[1];
  value->startpos[2] = fiptr->startpos[2];
  value->startpos[3] = fiptr->startpos[3];
  value->weight = fiptr->weight;
}

void gt_chain_printchainelem(FILE *outfp,const GtChainmatchvalues *value)
{
  fprintf(outfp,
          FormatSeqpos " " FormatSeqpos " " FormatSeqpos " " FormatSeqpos
          " %ld\n",PRINTSeqposcast(value->startpos[0]),
                   PRINTSeqposcast(value->endpos[0]),
                   PRINTSeqposcast(value->startpos[1]),
                   PRINTSeqposcast(value->endpos[1]),
                   value->weight);
}
