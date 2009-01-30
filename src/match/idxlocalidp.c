#include <stdarg.h>
#include "core/ma_api.h"
#include "core/symboldef.h"
#include "core/unused_api.h"
#include "absdfstrans-imp.h"

#define INFTY -1L
#define REPLACEMENTSCORE(A,B) ((A) == (B) ? scorevalues->matchscore\
                                          : scorevalues->mismatchscore)

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

static void firstcolumn (Column *column,
                         const Scorevalues *scorevalues,
                         unsigned long querylength)
{
  unsigned long i;

  if (column->colvalues == NULL)
  {
    column->colvalues = gt_malloc (sizeof (Matrixvalue) * (querylength + 1));
  }
  column->colvalues[0].repcell = 0;
  column->colvalues[0].inscell = scorevalues->gapstart;
  column->colvalues[0].delcell = scorevalues->gapstart;
  column->colvalues[0].bestcell = 0;
  column->maxvalue = 0;
  column->pprefixlen = 0;
  for (i = 1UL; i <= querylength; i++)
  {
    column->colvalues[i].repcell = INFTY;
    column->colvalues[i].inscell = INFTY;
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
        column->colvalues[i].delcell = INFTY;
      }
    }
    column->colvalues[i].bestcell = 0;
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

  if (outcol->colvalues == NULL)
  {
    outcol->colvalues = gt_malloc (sizeof (Matrixvalue) * (querylength + 1));
  }
  outcol->colvalues[0].repcell = outcol->colvalues[0].delcell = INFTY;
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
      outcol->colvalues[0].inscell = INFTY;
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
                                     REPLACEMENTSCORE(dbchar,query[i]);
    } else
    {
      outcol->colvalues[i].repcell = INFTY;
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
        outcol->colvalues[i].inscell = INFTY;
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
        outcol->colvalues[i].delcell = INFTY;
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

  column->colvalues[0].repcell = column->colvalues[0].delcell = INFTY;
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
      column->colvalues[0].inscell = INFTY;
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
                                     REPLACEMENTSCORE(dbchar,query[i]);
    } else
    {
      column->colvalues[i].repcell = INFTY;
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
        column->colvalues[i].inscell = INFTY;
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
        column->colvalues[i].delcell = INFTY;
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
  printf("querylength = %lu\n",lci->querylength);
  printf("threshold = %lu\n",lci->threshold);
}

static void locali_freedfsconstinfo (void **dfsconstinfo)
{
  Limdfsconstinfo *lci = (Limdfsconstinfo *) *dfsconstinfo;

  gt_free (lci);
  *dfsconstinfo = NULL;
}

static void locali_initLimdfsstate (DECLAREPTRDFSSTATE (aliascolumn),
                                    void *dfsconstinfo)
{
  Column *column = (Column *) aliascolumn;
  const Limdfsconstinfo *lci = (Limdfsconstinfo *) dfsconstinfo;

  firstcolumn (column, &lci->scorevalues, lci->querylength);
}

static void locali_initLimdfsstackelem (DECLAREPTRDFSSTATE (aliascolumn))
{
  ((Column *) aliascolumn)->colvalues = NULL;
}

static void locali_freeLimdfsstackelem (DECLAREPTRDFSSTATE (aliascolumn))
{
  gt_free(((Column *) aliascolumn)->colvalues);
}

static void locali_fullmatchLimdfsstate (Limdfsresult *limdfsresult,
                                         DECLAREPTRDFSSTATE(aliascolumn),
                                         GT_UNUSED Seqpos leftbound,
                                         GT_UNUSED Seqpos rightbound,
                                         GT_UNUSED Seqpos width,
                                         GT_UNUSED unsigned long currentdepth,
                                         void *dfsconstinfo)
{
  Column *column = (Column *) aliascolumn;
  const Limdfsconstinfo *lci = (Limdfsconstinfo *) dfsconstinfo;

  if (column->maxvalue >= lci->threshold)
  {
    printf("maxvalue = %lu\n",column->maxvalue);
    limdfsresult->status = Limdfssuccess;
    limdfsresult->distance = column->maxvalue;
    limdfsresult->pprefixlen = column->pprefixlen;
  } else
  {
    if (column->maxvalue > 0)
    {
      printf("maxvalue = %lu\n",column->maxvalue);
      limdfsresult->status = Limdfscontinue;
    } else
    {
      limdfsresult->status = Limdfsstop;
    }
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

  nextcolumn (outcol,&lci->scorevalues,currentchar,
              lci->query,lci->querylength,incol);
}

static void locali_inplacenextLimdfsstate (const void *dfsconstinfo,
                                           DECLAREPTRDFSSTATE (aliascolumn),
                                           GT_UNUSED unsigned long currentdepth,
                                           Uchar currentchar)
{
  Column *column = (Column *) aliascolumn;
  const Limdfsconstinfo *lci = (const Limdfsconstinfo *) dfsconstinfo;

  inplacenextcolumn (&lci->scorevalues,currentchar,
                     lci->query,lci->querylength,column);
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
    locali_fullmatchLimdfsstate,
    locali_nextLimdfsstate,
    locali_inplacenextLimdfsstate,
#ifdef SKDEBUG
    NULL
#endif /*  */
  };
  return &locali_adfst;
}
