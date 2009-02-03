#include <stdarg.h>
#include "core/ma_api.h"
#include "core/symboldef.h"
#include "core/chardef.h"
#include "core/unused_api.h"
#include "core/assert_api.h"
#include "absdfstrans-imp.h"

#define MINUSINFTY (-1L)
#define REPLACEMENTSCORE(A,B) (((A) != (B) || ISSPECIAL(A))\
                                   ? scorevalues->mismatchscore\
                                   : scorevalues->matchscore)

typedef struct
{
  long matchscore,   /* must be positive */
       mismatchscore,/* must be negative */
       gapstart,     /* must be negative */
       gapextend;    /* must be negative */
} Scorevalues;

typedef struct
{
  Scorevalues scorevalues;
  const Uchar *query;
  unsigned long querylength,
                threshold;
} Limdfsconstinfo;

typedef struct
{
  long repcell,
       inscell,
       delcell,
       bestcell;
} Matrixvalue;

typedef struct
{
  Matrixvalue *colvalues;
  unsigned long pprefixlen,
                maxvalue;
} Column;

typedef struct
{
  Limdfsstatus status;
  unsigned long qseqendpos;
  long alignmentscore;
} Idxlocaliresult;

static inline long max2 (long a,long b)
{
  return (a < b) ? b : a;
}

static inline long max3 (long a,long b,long c)
{
  long temp;

  temp = (a < b) ? b : a;
  return (temp < c) ? c : temp;
}

#ifdef SKDEBUG
static void showscorecolumn(const Column *column,
                            unsigned long querylength,
                            unsigned long currentdepth)
{
  printf("at depth %lu\n",currentdepth);
  if (column->colvalues == NULL)
  {
    printf("empty column\n");
  } else
  {
    unsigned long idx;

    for (idx = 0; idx <= querylength; idx++)
    {
      if (column->colvalues[idx].bestcell > 0)
      {
        printf("(%lu,%ld) ",idx,column->colvalues[idx].bestcell);
      }
    }
    printf("max=%lu\n",column->maxvalue);
  }
}

void locali_showLimdfsstate(const DECLAREPTRDFSSTATE(aliasstate),
                            unsigned long currentdepth,
                            const void *dfsconstinfo)
{
  const Limdfsconstinfo *lci = (const Limdfsconstinfo *) dfsconstinfo;

  showscorecolumn((const Column *) aliasstate,lci->querylength,currentdepth);
}
#endif

static void secondcolumn (Column *column,
                          const Scorevalues *scorevalues,
                          const Uchar *query,
                          unsigned long querylength,
                          Uchar dbchar)
{
  unsigned long i;

  if (column->colvalues == NULL)
  {
    column->colvalues = gt_malloc (sizeof (Matrixvalue) * (querylength + 1));
  }
  column->colvalues[0].repcell = MINUSINFTY;
  column->colvalues[0].inscell = scorevalues->gapstart +
                                 scorevalues->gapextend;
  column->colvalues[0].delcell = MINUSINFTY;
  column->colvalues[0].bestcell = MINUSINFTY;
  column->maxvalue = 0;
  column->pprefixlen = 0;
  for (i = 1UL; i <= querylength; i++)
  {
    column->colvalues[i].delcell = MINUSINFTY;
    column->colvalues[i].inscell = MINUSINFTY;
    column->colvalues[i].repcell = REPLACEMENTSCORE(dbchar,query[i-1]);
    column->colvalues[i].bestcell = max2(column->colvalues[i].delcell,
                                         column->colvalues[i].repcell);
    if (column->colvalues[i].bestcell > 0 &&
        column->colvalues[i].bestcell > (long) column->maxvalue)
    {
      column->maxvalue = (unsigned long) column->colvalues[i].bestcell;
      column->pprefixlen = i;
    }
  }
}

static void nextcolumn (Column *outcol,
                        const Scorevalues *scorevalues,
                        const Uchar dbchar,
                        const Uchar *query,
                        unsigned long querylength,
                        const Column *incol)
{
  unsigned long i;

  gt_assert(outcol != incol);
  gt_assert(outcol->colvalues != incol->colvalues);
  if (outcol->colvalues == NULL)
  {
    outcol->colvalues = gt_malloc (sizeof (Matrixvalue) * (querylength + 1));
  }
  outcol->colvalues[0].repcell = outcol->colvalues[0].delcell = MINUSINFTY;
  if (incol->colvalues[0].inscell > 0)
  {
    if (incol->colvalues[0].bestcell > 0)
    {
      outcol->colvalues[0].inscell
        = max2 (incol->colvalues[0].inscell + scorevalues->gapextend,
                incol->colvalues[0].bestcell + scorevalues->gapstart +
                                               scorevalues->gapextend);
    } else
    {
      outcol->colvalues[0].inscell = incol->colvalues[0].inscell +
                                     scorevalues->gapextend;
    }
  } else
  {
    if (incol->colvalues[0].bestcell > 0)
    {
      outcol->colvalues[0].inscell = incol->colvalues[0].bestcell +
                                     scorevalues->gapstart +
                                     scorevalues->gapextend;
    } else
    {
      outcol->colvalues[0].inscell = MINUSINFTY;
    }
  }
  outcol->colvalues[0].bestcell = max3 (outcol->colvalues[0].repcell,
                                        outcol->colvalues[0].inscell,
                                        outcol->colvalues[0].delcell);
  outcol->maxvalue = (unsigned long) max2 (0,outcol->colvalues[0].bestcell);
  outcol->pprefixlen = 0;
  for (i = 1UL; i <= querylength; i++)
  {
    if (incol->colvalues[i-1].bestcell > 0)
    {
      outcol->colvalues[i].repcell = incol->colvalues[i-1].bestcell +
                                     REPLACEMENTSCORE(dbchar,query[i-1]);
    } else
    {
      outcol->colvalues[i].repcell = MINUSINFTY;
    }
    if (incol->colvalues[i].inscell > 0)
    {
      if (incol->colvalues[i].bestcell > 0)
      {
        outcol->colvalues[i].inscell
          = max2 (incol->colvalues[i].inscell + scorevalues->gapextend,
                  incol->colvalues[i].bestcell + scorevalues->gapstart +
                                                 scorevalues->gapextend);
      } else
      {
        outcol->colvalues[i].inscell = incol->colvalues[i].inscell +
                                       scorevalues->gapextend;
      }
    } else
    {
      if (incol->colvalues[i].bestcell > 0)
      {
        outcol->colvalues[i].inscell = incol->colvalues[i].bestcell +
                                       scorevalues->gapstart +
                                       scorevalues->gapextend;
      } else
      {
        outcol->colvalues[i].inscell = MINUSINFTY;
      }
    }
    if (outcol->colvalues[i-1].delcell > 0)
    {
      if (outcol->colvalues[i-1].bestcell > 0)
      {
        outcol->colvalues[i].delcell
          = max2 (outcol->colvalues[i-1].delcell + scorevalues->gapextend,
                  outcol->colvalues[i-1].bestcell + scorevalues->gapstart +
                                                    scorevalues->gapextend);
      } else
      {
        outcol->colvalues[i].delcell = outcol->colvalues[i-1].delcell +
                                       scorevalues->gapextend;
      }
    } else
    {
      if (outcol->colvalues[i-1].bestcell > 0)
      {
        outcol->colvalues[i].delcell = outcol->colvalues[i-1].bestcell +
                                       scorevalues->gapstart +
                                       scorevalues->gapextend;
      } else
      {
        outcol->colvalues[i].delcell = MINUSINFTY;
      }
    }
    outcol->colvalues[i].bestcell = max3 (outcol->colvalues[i].repcell,
                                          outcol->colvalues[i].inscell,
                                          outcol->colvalues[i].delcell);
    if (outcol->colvalues[i].bestcell > 0 &&
        outcol->colvalues[i].bestcell > (long) outcol->maxvalue)
    {
      outcol->maxvalue = (unsigned long) outcol->colvalues[i].bestcell;
      outcol->pprefixlen = i;
    }
  }
}

static void inplacenextcolumn (const Scorevalues *scorevalues,
                               const Uchar dbchar,
                               const Uchar *query,
                               unsigned long querylength,
                               Column *column)
{
  unsigned long i;
  Matrixvalue nw, west;

  column->colvalues[0].repcell = column->colvalues[0].delcell = MINUSINFTY;
  if (column->colvalues[0].inscell > 0)
  {
    if (column->colvalues[0].bestcell > 0)
    {
      column->colvalues[0].inscell
        = max2 (column->colvalues[0].inscell + scorevalues->gapextend,
                column->colvalues[0].bestcell + scorevalues->gapstart +
                                                scorevalues->gapextend);
    } else
    {
      column->colvalues[0].inscell
        = column->colvalues[0].inscell + scorevalues->gapextend;
    }
  } else
  {
    if (column->colvalues[0].bestcell > 0)
    {
      column->colvalues[0].inscell = column->colvalues[0].bestcell +
                                     scorevalues->gapstart +
                                     scorevalues->gapextend;
    } else
    {
      column->colvalues[0].inscell = MINUSINFTY;
    }
  }
  column->colvalues[0].bestcell = max3 (column->colvalues[0].repcell,
                                        column->colvalues[0].inscell,
                                        column->colvalues[0].delcell);
  column->maxvalue = (unsigned long) max2(0,column->colvalues[0].bestcell);
  column->pprefixlen = 0;
  nw = column->colvalues[0];
  for (i = 1UL; i <= querylength; i++)
  {
    west = column->colvalues[i];
    if (nw.bestcell > 0)
    {
      column->colvalues[i].repcell = nw.bestcell +
                                     REPLACEMENTSCORE(dbchar,query[i-1]);
    } else
    {
      column->colvalues[i].repcell = MINUSINFTY;
    }
    if (west.inscell > 0)
    {
      if (west.bestcell > 0)
      {
        column->colvalues[i].inscell
          = max2 (west.inscell + scorevalues->gapextend,
                  west.bestcell + scorevalues->gapstart
                                + scorevalues->gapextend);
      } else
      {
        column->colvalues[i].inscell = west.inscell + scorevalues->gapextend;
      }
    } else
    {
      if (west.bestcell > 0)
      {
        column->colvalues[i].inscell = west.bestcell + scorevalues->gapstart +
                                                       scorevalues->gapextend;
      } else
      {
        column->colvalues[i].inscell = MINUSINFTY;
      }
    }
    if (column->colvalues[i-1].delcell > 0)
    {
      if (column->colvalues[i-1].bestcell > 0)
      {
        column->colvalues[i].delcell
          = max2 (column->colvalues[i-1].delcell + scorevalues->gapextend,
                  column->colvalues[i-1].bestcell + scorevalues->gapstart +
                                                    scorevalues->gapextend);
      } else
      {
        column->colvalues[i].delcell = column->colvalues[i-1].delcell +
                                       scorevalues->gapextend;
      }
    } else
    {
      if (column->colvalues[i-1].bestcell > 0)
      {
        column->colvalues[i].delcell = column->colvalues[i-1].bestcell +
                                       scorevalues->gapstart +
                                       scorevalues->gapextend;
      } else
      {
        column->colvalues[i].delcell = MINUSINFTY;
      }
    }
    nw = west;
    column->colvalues[i].bestcell = max3 (column->colvalues[i].repcell,
                                          column->colvalues[i].inscell,
                                          column->colvalues[i].delcell);
    if (column->colvalues[i].bestcell > 0 &&
        column->colvalues[i].bestcell > (long) column->maxvalue)
    {
      column->maxvalue = (unsigned long) column->colvalues[i].bestcell;
      column->pprefixlen = i;
    }
  }
}

static void *locali_allocatedfsconstinfo (GT_UNUSED unsigned int alphasize)
{
  Limdfsconstinfo *lci = gt_malloc (sizeof (Limdfsconstinfo));

  return lci;
}

static void locali_initdfsconstinfo (void *dfsconstinfo,
                                     unsigned int alphasize,...)
{
  va_list ap;
  Limdfsconstinfo *lci = (Limdfsconstinfo *) dfsconstinfo;

  va_start (ap, alphasize);
  lci->scorevalues.matchscore = va_arg (ap, long);
  lci->scorevalues.mismatchscore = va_arg (ap, long);
  lci->scorevalues.gapstart = va_arg (ap, long);
  lci->scorevalues.gapextend = va_arg (ap, long);
  lci->query = va_arg (ap, const Uchar *);
  lci->querylength = va_arg (ap, unsigned long);
  lci->threshold = va_arg (ap, unsigned long);
  va_end(ap);
}

static void locali_freedfsconstinfo (void **dfsconstinfo)
{
  Limdfsconstinfo *lci = (Limdfsconstinfo *) *dfsconstinfo;

  gt_free (lci);
  *dfsconstinfo = NULL;
}

static void locali_initLimdfsstate (DECLAREPTRDFSSTATE (aliasstate),
                                    GT_UNUSED void *dfsconstinfo)
{
  Column *column = (Column *) aliasstate;

  column->colvalues = NULL;
}

static void locali_initLimdfsstackelem (DECLAREPTRDFSSTATE (aliasstate))
{
  ((Column *) aliasstate)->colvalues = NULL;
}

static void locali_freeLimdfsstackelem (DECLAREPTRDFSSTATE (aliasstate))
{
  gt_free(((Column *) aliasstate)->colvalues);
}

static void locali_copyLimdfsstate (DECLAREPTRDFSSTATE(deststate),
                                    const DECLAREPTRDFSSTATE(srcstate),
                                    void *dfsconstinfo)
{
  unsigned long idx;
  Limdfsconstinfo *lci = (Limdfsconstinfo *) dfsconstinfo;
  Column *destcol = (Column *) deststate;
  const Column *srccol = (const Column *) srcstate;

  if (srccol->colvalues != NULL)
  {
    if (destcol->colvalues == NULL)
    {
      destcol->colvalues = gt_malloc (sizeof (Matrixvalue) *
                                      (lci->querylength + 1));
    }
    for (idx = 0; idx<=lci->querylength; idx++)
    {
      destcol->colvalues[idx] = srccol->colvalues[idx];
    }
  }
}

static void locali_fullmatchLimdfsstate (Limdfsresult *limdfsresult,
                                         DECLAREPTRDFSSTATE(aliasstate),
                                         GT_UNUSED Seqpos leftbound,
                                         GT_UNUSED Seqpos rightbound,
                                         GT_UNUSED Seqpos width,
                                         GT_UNUSED unsigned long currentdepth,
                                         void *dfsconstinfo)
{
  Column *column = (Column *) aliasstate;
  const Limdfsconstinfo *lci = (Limdfsconstinfo *) dfsconstinfo;

  if (column->colvalues != NULL)
  {
    if (column->maxvalue >= lci->threshold)
    {
      limdfsresult->status = Limdfssuccess;
      limdfsresult->distance = column->maxvalue;
      limdfsresult->pprefixlen = column->pprefixlen;
    } else
    {
      if (column->maxvalue > 0)
      {
        limdfsresult->status = Limdfscontinue;
      } else
      {
        limdfsresult->status = Limdfsstop;
      }
    }
  } else
  {
    limdfsresult->status = Limdfscontinue;
  }
}

static void locali_nextLimdfsstate (const void *dfsconstinfo,
                                    DECLAREPTRDFSSTATE (aliasoutcol),
                                    GT_UNUSED unsigned long currentdepth,
                                    Uchar currentchar,
                                    const DECLAREPTRDFSSTATE (aliasincol))
{
  const Limdfsconstinfo *lci = (const Limdfsconstinfo *) dfsconstinfo;
  Column *outcol = (Column *) aliasoutcol;
  const Column *incol = (const Column *) aliasincol;

  if (incol->colvalues == NULL)
  {
    gt_assert(currentdepth == 1UL);
    secondcolumn (outcol, &lci->scorevalues, lci->query, lci->querylength,
                  currentchar);
  } else
  {
    nextcolumn (outcol,&lci->scorevalues,currentchar,
                lci->query,lci->querylength,incol);
  }
}

static void locali_inplacenextLimdfsstate (const void *dfsconstinfo,
                                           DECLAREPTRDFSSTATE (aliasstate),
                                           GT_UNUSED unsigned long currentdepth,
                                           Uchar currentchar)
{
  Column *column = (Column *) aliasstate;
  const Limdfsconstinfo *lci = (const Limdfsconstinfo *) dfsconstinfo;

  if (column->colvalues == NULL)
  {
    gt_assert(currentdepth == 1UL);
    secondcolumn (column, &lci->scorevalues, lci->query, lci->querylength,
                  currentchar);
  } else
  {
    inplacenextcolumn (&lci->scorevalues,currentchar,
                       lci->query,lci->querylength,column);
  }
}

const AbstractDfstransformer *locali_AbstractDfstransformer (void)
{
  static const AbstractDfstransformer locali_adfst =
  {
    sizeof (Column),
    locali_allocatedfsconstinfo,
    locali_initdfsconstinfo,
    NULL,
    locali_freedfsconstinfo,
    locali_initLimdfsstate,
    locali_initLimdfsstackelem,
    locali_freeLimdfsstackelem,
    locali_copyLimdfsstate,
    locali_fullmatchLimdfsstate,
    locali_nextLimdfsstate,
    locali_inplacenextLimdfsstate,
#ifdef SKDEBUG
    locali_showLimdfsstate
#endif /*  */
  };
  return &locali_adfst;
}
