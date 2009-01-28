#include <stdarg.h>
#include "core/ma_api.h"
#include "core/symboldef.h"
#include "core/unused_api.h"
#include "absdfstrans-imp.h"

#define INFTY -1L
#define REPLACEMENTSCORE(A,B) ((A) == (B) ? cost->matchscore\
                                          : cost->mismatchscore)

typedef struct
{
  long matchscore,   /* must be positive */
       mismatchscore,/* must be negative */
       gapstart,     /* must be negative */
       gapextend;    /* must be negative */
} Cost;

typedef struct
{
  const Cost *cost;
  const Uchar *qseq;
  unsigned long lengthofqseq,
                threshold;
} Limdfsconstinfo;

typedef struct
{
  long repvalue,
       insvalue,
       delvalue,
       bestvalue;
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
                         const Cost *cost,
                         unsigned long lengthofqseq)
{
  unsigned long i;

  if (column->colvalues == NULL)
  {
    column->colvalues == gt_malloc (sizeof (Matrixvalue) * (lengthofqseq + 1));
  }
  column->colvalues = gt_malloc (sizeof (Matrixvalue) * (lengthofqseq + 1));
  column->colvalues[0].repvalue = 0;
  column->colvalues[0].insvalue = cost->gapstart;
  column->colvalues[0].delvalue = cost->gapstart;
  column->colvalues[0].bestvalue = 0;
  column->maxvalue = 0;
  column->pprefixlen = 0;
  for (i = 1UL; i <= lengthofqseq; i++)
  {
    column->colvalues[i].repvalue = INFTY;
    column->colvalues[i].insvalue = INFTY;
    if (column->colvalues[i-1].delvalue > 0)
    {
      if (column->colvalues[i-1].bestvalue > 0)
      {
        column->colvalues[i].delvalue
          = max2 (column->colvalues[i-1].delvalue + cost->gapextend,
                  column->colvalues[i-1].bestvalue + cost->gapstart +
                                                     cost->gapextend);
      } else
      {
        column->colvalues[i].delvalue = column->colvalues[i-1].delvalue +
                                        cost->gapextend;
      }
    } else
    {
      if (column->colvalues[i-1].bestvalue > 0)
      {
        column->colvalues[i].delvalue = column->colvalues[i-1].bestvalue +
                                        cost->gapstart + cost->gapextend;
      } else
      {
        column->colvalues[i].delvalue = INFTY;
      }
    }
    column->colvalues[i].bestvalue = max3 (column->colvalues[i].repvalue,
                                           column->colvalues[i].insvalue,
                                           column->colvalues[i].delvalue);
    if (column->colvalues[i].bestvalue > (long) column->maxvalue)
    {
      column->maxvalue = (unsigned long) column->colvalues[i].bestvalue;
      column->pprefixlen = i;
    }
  }
}

static void nextcolumn (Column *outcol,
                        const Cost *cost,
                        const Uchar dbchar,
                        const Uchar *qseq,
                        unsigned long lengthofqseq,
                        const Column *incol)
{
  unsigned long i;

  if (outcol->colvalues == NULL)
  {
    outcol->colvalues = gt_malloc (sizeof (Matrixvalue) * (lengthofqseq + 1));
  }
  outcol->colvalues[0].repvalue = outcol->colvalues[0].delvalue = INFTY;
  if (incol->colvalues[0].insvalue > 0)
  {
    if (incol->colvalues[0].bestvalue > 0)
    {
      outcol->colvalues[0].insvalue
        = max2 (incol->colvalues[0].insvalue + cost->gapextend,
                incol->colvalues[0].bestvalue + cost->gapstart +
                                                cost->gapextend);
    } else
    {
      outcol->colvalues[0].insvalue = incol->colvalues[0].insvalue +
                                      cost->gapextend;
    }
  } else
  {
    if (incol->colvalues[0].bestvalue > 0)
    {
      outcol->colvalues[0].insvalue = incol->colvalues[0].bestvalue +
                                      cost->gapstart + cost->gapextend;
    } else
    {
      outcol->colvalues[0].insvalue = INFTY;
    }
  }
  outcol->colvalues[0].bestvalue = max3 (outcol->colvalues[0].repvalue,
                                         outcol->colvalues[0].insvalue,
                                         outcol->colvalues[0].delvalue);
  outcol->maxvalue = (unsigned long) outcol->colvalues[0].bestvalue;
  outcol->pprefixlen = 0;
  for (i = 1UL; i <= lengthofqseq; i++)
  {
    if (incol->colvalues[i-1].bestvalue > 0)
    {
      outcol->colvalues[i].repvalue = incol->colvalues[i-1].bestvalue +
                                      REPLACEMENTSCORE(dbchar,qseq[i]);
    } else
    {
      outcol->colvalues[i].repvalue = INFTY;
    }
    if (incol->colvalues[i].insvalue > 0)
    {
      if (incol->colvalues[i].bestvalue > 0)
      {
        outcol->colvalues[i].insvalue
          = max2 (incol->colvalues[i].insvalue + cost->gapextend,
                  incol->colvalues[i].bestvalue + cost->gapstart +
                                                  cost->gapextend);
      } else
      {
        outcol->colvalues[i].insvalue = incol->colvalues[i].insvalue +
                                        cost->gapextend;
      }
    } else
    {
      if (incol->colvalues[i].bestvalue > 0)
      {
        outcol->colvalues[i].insvalue = incol->colvalues[i].bestvalue +
                                        cost->gapstart + cost->gapextend;
      } else
      {
        outcol->colvalues[i].insvalue = INFTY;
      }
    }
    if (outcol->colvalues[i-1].delvalue > 0)
    {
      if (outcol->colvalues[i-1].bestvalue > 0)
      {
        outcol->colvalues[i].delvalue
          = max2 (outcol->colvalues[i-1].delvalue + cost->gapextend,
                  outcol->colvalues[i-1].bestvalue + cost->gapstart +
                                                     cost->gapextend);
      } else
      {
        outcol->colvalues[i].delvalue = outcol->colvalues[i-1].delvalue +
                                        cost->gapextend;
      }
    } else
    {
      if (outcol->colvalues[i-1].bestvalue > 0)
      {
        outcol->colvalues[i].delvalue = outcol->colvalues[i-1].bestvalue +
                                        cost->gapstart + cost->gapextend;
      } else
      {
        outcol->colvalues[i].delvalue = INFTY;
      }
    }
    outcol->colvalues[i].bestvalue = max3 (outcol->colvalues[i].repvalue,
                                           outcol->colvalues[i].insvalue,
                                           outcol->colvalues[i].delvalue);
    if (outcol->colvalues[i].bestvalue > (long) outcol->maxvalue)
    {
      outcol->maxvalue = (unsigned long) outcol->colvalues[i].bestvalue;
      outcol->pprefixlen = i;
    }
  }
}

static void inplacenextcolumn (const Cost *cost,
                               const Uchar dbchar,
                               const Uchar *qseq,
                               unsigned long lengthofqseq,
                               Column *column)
{
  unsigned long i;
  Matrixvalue nw, west;

  column->colvalues[0].repvalue = column->colvalues[0].delvalue = INFTY;
  if (column->colvalues[0].insvalue > 0)
  {
    if (column->colvalues[0].bestvalue > 0)
    {
      column->colvalues[0].insvalue
        = max2 (column->colvalues[0].insvalue + cost->gapextend,
                column->colvalues[0].bestvalue + cost->gapstart +
                                                 cost->gapextend);
    } else
    {
      column->colvalues[0].insvalue
        = column->colvalues[0].insvalue + cost->gapextend;
    }
  } else
  {
    if (column->colvalues[0].bestvalue > 0)
    {
      column->colvalues[0].insvalue = column->colvalues[0].bestvalue +
                                      cost->gapstart + cost->gapextend;
    } else
    {
      column->colvalues[0].insvalue = INFTY;
    }
  }
  column->colvalues[0].bestvalue = max3 (column->colvalues[0].repvalue,
                                         column->colvalues[0].insvalue,
                                         column->colvalues[0].delvalue);
  column->maxvalue = (unsigned long) column->colvalues[0].bestvalue;
  column->pprefixlen = 0;
  nw = column->colvalues[0];
  for (i = 1UL; i <= lengthofqseq; i++)

  {
    west = column->colvalues[i];
    if (nw.bestvalue > 0)
    {
      column->colvalues[i].repvalue = nw.bestvalue +
                                      REPLACEMENTSCORE(dbchar,qseq[i]);
    } else
    {
      column->colvalues[i].repvalue = INFTY;
    }
    if (west.insvalue > 0)
    {
      if (west.bestvalue > 0)
      {
        column->colvalues[i].insvalue
          = max2 (west.insvalue + cost->gapextend,
                  west.bestvalue + cost->gapstart + cost->gapextend);
      } else
      {
        column->colvalues[i].insvalue = west.insvalue + cost->gapextend;
      }
    } else
    {
      if (west.bestvalue > 0)
      {
        column->colvalues[i].insvalue = west.bestvalue + cost->gapstart +
                                                         cost->gapextend;
      } else
      {
        column->colvalues[i].insvalue = INFTY;
      }
    }
    if (column->colvalues[i-1].delvalue > 0)
    {
      if (column->colvalues[i-1].bestvalue > 0)
      {
        column->colvalues[i].delvalue
          = max2 (column->colvalues[i-1].delvalue + cost->gapextend,
                  column->colvalues[i-1].bestvalue + cost->gapstart +
                                                     cost->gapextend);
      } else
      {
        column->colvalues[i].delvalue = column->colvalues[i-1].delvalue +
                                        cost->gapextend;
      }
    } else
    {
      if (column->colvalues[i-1].bestvalue > 0)
      {
        column->colvalues[i].delvalue = column->colvalues[i-1].bestvalue +
                                        cost->gapstart + cost->gapextend;
      } else
      {
        column->colvalues[i].delvalue = INFTY;
      }
    }
    nw = west;
    column->colvalues[i].bestvalue = max3 (column->colvalues[i].repvalue,
                                           column->colvalues[i].insvalue,
                                           column->colvalues[i].delvalue);
    if (column->colvalues[i].bestvalue > (long) column->maxvalue)
    {
      column->maxvalue = (unsigned long) column->colvalues[i].bestvalue;
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
                                     unsigned int alphasize,
                                     ...)
{
  va_list ap;
  Limdfsconstinfo *lci = (Limdfsconstinfo *) dfsconstinfo;

  va_start (ap, alphasize);
  lci->cost = va_arg (ap, const Cost *);
  lci->qseq = va_arg (ap, const Uchar *);
  lci->lengthofqseq = va_arg (ap, unsigned long);
  lci->threshold = va_arg (ap, unsigned long);
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

  firstcolumn (column, lci->cost, lci->lengthofqseq);
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

  nextcolumn (outcol,lci->cost,currentchar,lci->qseq,lci->lengthofqseq,incol);
}

static void locali_inplacenextLimdfsstate (const void *dfsconstinfo,
                                           DECLAREPTRDFSSTATE (aliascolumn),
                                           GT_UNUSED unsigned long currentdepth,
                                           Uchar currentchar)
{
  Column *column = (Column *) aliascolumn;
  const Limdfsconstinfo *lci = (const Limdfsconstinfo *) dfsconstinfo;

  inplacenextcolumn (lci->cost,currentchar,lci->qseq,lci->lengthofqseq,column);
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
