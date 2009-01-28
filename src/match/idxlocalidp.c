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
  long Rvalue,
       Dvalue,
       Ivalue,
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
  unsigned long j;

  if (column->colvalues == NULL)
  {
    column->colvalues == gt_malloc (sizeof (Matrixvalue) * (lengthofqseq + 1));
  }
  column->colvalues = gt_malloc (sizeof (Matrixvalue) * (lengthofqseq + 1));
  column->colvalues[0].Rvalue = 0;
  column->colvalues[0].Dvalue = cost->gapstart;
  column->colvalues[0].Ivalue = cost->gapstart;
  column->colvalues[0].bestvalue = 0;
  column->maxvalue = 0;
  for (j = 1UL; j <= lengthofqseq; j++)
  {
    column->colvalues[j].Rvalue = INFTY;
    column->colvalues[j].Dvalue = INFTY;
    if (column->colvalues[j - 1].Ivalue > 0)
    {
      if (column->colvalues[j - 1].bestvalue > 0)
      {
        column->colvalues[j].Ivalue
          = max2 (column->colvalues[j - 1].Ivalue + cost->gapextend,
                  column->colvalues[j - 1].bestvalue + cost->gapstart +
                                                       cost->gapextend);
      } else
      {
        column->colvalues[j].Ivalue = column->colvalues[j - 1].Ivalue +
                                      cost->gapextend;
      }
    } else
    {
      if (column->colvalues[j - 1].bestvalue > 0)
      {
        column->colvalues[j].Ivalue = column->colvalues[j - 1].bestvalue +
                                      cost->gapstart + cost->gapextend;
      } else
      {
        column->colvalues[j].Ivalue = INFTY;
      }
    }
    column->colvalues[j].bestvalue = max3 (column->colvalues[j].Rvalue,
                                           column->colvalues[j].Dvalue,
                                           column->colvalues[j].Ivalue);
    if (column->colvalues[j].bestvalue > (long) column->maxvalue)
    {
      column->maxvalue = (unsigned long) column->colvalues[j].bestvalue;
      column->pprefixlen = j;
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
  unsigned long j;

  if (outcol->colvalues == NULL)
  {
    outcol->colvalues = gt_malloc (sizeof (Matrixvalue) * (lengthofqseq + 1));
  }
  outcol->colvalues[0].Rvalue = outcol->colvalues[0].Ivalue = INFTY;
  if (incol->colvalues[0].Dvalue > 0)
  {
    if (incol->colvalues[0].bestvalue > 0)
    {
      outcol->colvalues[0].Dvalue
        = max2 (incol->colvalues[0].Dvalue + cost->gapextend,
                incol->colvalues[0].bestvalue + cost->gapstart +
                                                cost->gapextend);
    } else
    {
      outcol->colvalues[0].Dvalue = incol->colvalues[0].Dvalue +
                                    cost->gapextend;
    }
  } else
  {
    if (incol->colvalues[0].bestvalue > 0)
    {
      outcol->colvalues[0].Dvalue = incol->colvalues[0].bestvalue +
                                    cost->gapstart + cost->gapextend;
    } else
    {
      outcol->colvalues[0].Dvalue = INFTY;
    }
  }
  outcol->colvalues[0].bestvalue = max3 (outcol->colvalues[0].Rvalue,
                                         outcol->colvalues[0].Dvalue,
                                         outcol->colvalues[0].Ivalue);
  outcol->maxvalue = 0;
  outcol->pprefixlen = 0;
  for (j = 1UL; j <= lengthofqseq; j++)
  {
    if (incol->colvalues[j - 1].bestvalue > 0)
    {
      outcol->colvalues[j].Rvalue = incol->colvalues[j - 1].bestvalue +
                                    REPLACEMENTSCORE(dbchar,qseq[j]);
    } else
    {
      outcol->colvalues[j].Rvalue = INFTY;
    }
    if (incol->colvalues[j].Dvalue > 0)
    {
      if (incol->colvalues[j].bestvalue > 0)
      {
        outcol->colvalues[j].Dvalue
          = max2 (incol->colvalues[j].Dvalue + cost->gapextend,
                  incol->colvalues[j].bestvalue + cost->gapstart +
                                                  cost->gapextend);
      } else
      {
        outcol->colvalues[j].Dvalue = incol->colvalues[j].Dvalue +
                                      cost->gapextend;
      }
    } else
    {
      if (incol->colvalues[j].bestvalue > 0)
      {
        outcol->colvalues[j].Dvalue = incol->colvalues[j].bestvalue +
                                      cost->gapstart + cost->gapextend;
      } else
      {
        outcol->colvalues[j].Dvalue = INFTY;
      }
    }
    if (outcol->colvalues[j - 1].Ivalue > 0)
    {
      if (outcol->colvalues[j - 1].bestvalue > 0)
      {
        outcol->colvalues[j].Ivalue
          = max2 (outcol->colvalues[j - 1].Ivalue + cost->gapextend,
                  outcol->colvalues[j - 1].bestvalue + cost->gapstart +
                                                       cost->gapextend);
      } else
      {
        outcol->colvalues[j].Ivalue = outcol->colvalues[j - 1].Ivalue +
                                      cost->gapextend;
      }
    } else
    {
      if (outcol->colvalues[j - 1].bestvalue > 0)
      {
        outcol->colvalues[j].Ivalue = outcol->colvalues[j - 1].bestvalue +
                                      cost->gapstart + cost->gapextend;
      } else
      {
        outcol->colvalues[j].Ivalue = INFTY;
      }
    }
    outcol->colvalues[j].bestvalue = max3 (outcol->colvalues[j].Rvalue,
                                           outcol->colvalues[j].Dvalue,
                                           outcol->colvalues[j].Ivalue);
    if (outcol->colvalues[j].bestvalue > (long) outcol->maxvalue)
    {
      outcol->maxvalue = (unsigned long) outcol->colvalues[j].bestvalue;
      outcol->pprefixlen = j;
    }
  }
}

static void inplacenextcolumn (const Cost *cost,
                               const Uchar dbchar,
                               const Uchar *qseq,
                               unsigned long lengthofqseq,
                               Column *column)
{
  unsigned long j;
  Matrixvalue nw, west;

  column->colvalues[0].Rvalue = column->colvalues[0].Ivalue = INFTY;
  if (column->colvalues[0].Dvalue > 0)
  {
    if (column->colvalues[0].bestvalue > 0)
    {
      column->colvalues[0].Dvalue
        = max2 (column->colvalues[0].Dvalue + cost->gapextend,
                column->colvalues[0].bestvalue + cost->gapstart +
                                                 cost->gapextend);
    } else
    {
      column->colvalues[0].Dvalue
        = column->colvalues[0].Dvalue + cost->gapextend;
    }
  } else
  {
    if (column->colvalues[0].bestvalue > 0)
    {
      column->colvalues[0].Dvalue = column->colvalues[0].bestvalue +
                                    cost->gapstart + cost->gapextend;
    } else
    {
      column->colvalues[0].Dvalue = INFTY;
    }
  }
  column->colvalues[0].bestvalue = max3 (column->colvalues[0].Rvalue,
                                         column->colvalues[0].Dvalue,
                                         column->colvalues[0].Ivalue);
  column->maxvalue = 0;
  nw = column->colvalues[0];
  for (j = 1UL; j <= lengthofqseq; j++)

  {
    west = column->colvalues[j];
    if (nw.bestvalue > 0)
    {
      column->colvalues[j].Rvalue = nw.bestvalue +
                                    REPLACEMENTSCORE(dbchar,qseq[j]);
    } else
    {
      column->colvalues[j].Rvalue = INFTY;
    }

    if (west.Dvalue > 0)
    {
      if (west.bestvalue > 0)
      {
        column->colvalues[j].Dvalue
          = max2 (west.Dvalue + cost->gapextend,
                  west.bestvalue + cost->gapstart + cost->gapextend);
      } else
      {
        column->colvalues[j].Dvalue = west.Dvalue + cost->gapextend;
      }
    } else
    {
      if (west.bestvalue > 0)
      {
        column->colvalues[j].Dvalue = west.bestvalue + cost->gapstart +
                                                       cost->gapextend;
      } else
      {
        column->colvalues[j].Dvalue = INFTY;
      }
    }
    if (column->colvalues[j - 1].Ivalue > 0)
    {
      if (column->colvalues[j - 1].bestvalue > 0)
      {
        column->colvalues[j].Ivalue
          = max2 (column->colvalues[j - 1].Ivalue + cost->gapextend,
                  column->colvalues[j - 1].bestvalue + cost->gapstart +
                  cost->gapextend);
      } else
      {
        column->colvalues[j].Ivalue = column->colvalues[j - 1].Ivalue +
                                      cost->gapextend;
      }
    } else
    {
      if (column->colvalues[j - 1].bestvalue > 0)
      {
        column->colvalues[j].Ivalue = column->colvalues[j - 1].bestvalue +
                                      cost->gapstart + cost->gapextend;
      } else
      {
        column->colvalues[j].Ivalue = INFTY;
      }
    }
    nw = west;
    column->colvalues[j].bestvalue = max3 (column->colvalues[j].Rvalue,
                                           column->colvalues[j].Dvalue,
                                           column->colvalues[j].Ivalue);
    if (column->colvalues[j].bestvalue > (long) column->maxvalue)
    {
      column->maxvalue = (unsigned long) column->colvalues[j].bestvalue;
      column->pprefixlen = j;
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
    limdfsresult->pprefixlen = column->pprefixlen;
    limdfsresult->distance = column->maxvalue;
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
