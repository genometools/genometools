#include <string.h>
#include "core/ma.h"
#include "core/minmax.h"
#include "core/assert_api.h"
#include "core/unused_api.h"
#include "core/divmodmul.h"
#include "match/squarededist.h"
#include "extended/alignment.h"
#include "extended/linearspace.h"
#include "extended/reconstructalignment.h"

static void firstEDtabRtabcolumn(GtUword *EDtabcolumn,
                                 GtUword *Rtabcolumn,
                                 GtUword ulen,
                                 const GtWord gapcost)
{
  GtUword rowindex;
  EDtabcolumn[0] = 0;
  Rtabcolumn[0]  = 0;

  for (rowindex=1; rowindex <= ulen; rowindex++)
  {
    EDtabcolumn[rowindex] =EDtabcolumn[rowindex-1]+ gapcost;
    Rtabcolumn[rowindex]  = rowindex;
  }
}

static void nextEDtabRtabcolumn(GtUword *EDtabcolumn,
                                GtUword *Rtabcolumn,
                                GtUword colindex,
                                GtUword midcolumn, GtUchar b,
                                const GtUchar *useq, GtUword ulen,
                                const GtWord matchcost,
                                const GtWord mismatchcost,
                                const GtWord gapcost)
{
  GtUword rowindex, val,
          northwestEDtabentry,
          westEDtabentry,
          northwestRtabentry,
          westRtabentry=0;
  bool updateRtabcolumn = false;

  gt_assert(EDtabcolumn != NULL);
  westEDtabentry = EDtabcolumn[0]; /* saves the first entry of EDtabcolumn */
  EDtabcolumn[0]+=gapcost;
  if (colindex > midcolumn)
  {
    updateRtabcolumn = true;
    Rtabcolumn[0] = 0;
  }

  for (rowindex = 1UL; rowindex <= ulen; rowindex++) {
    northwestEDtabentry = westEDtabentry;
    northwestRtabentry  = westRtabentry;
    westEDtabentry = EDtabcolumn[rowindex];
    westRtabentry  = Rtabcolumn[rowindex];
    EDtabcolumn[rowindex]+=gapcost; /* 1. recurrence */
    /* Rtabcolumn[rowindex] is unchanged */
    /* 2. recurrence: */
    if ((val = northwestEDtabentry +
        (useq[rowindex-1] == b ? matchcost : mismatchcost))
        <= EDtabcolumn[rowindex])
    {
      EDtabcolumn[rowindex] = val;
      if (updateRtabcolumn)
      {
        Rtabcolumn[rowindex] = northwestRtabentry;
      }
    }
    /* 3. recurrence: */
    if ((val = EDtabcolumn[rowindex-1]+gapcost) < EDtabcolumn[rowindex])
    {
      EDtabcolumn[rowindex] = val;
      if (updateRtabcolumn)
      {
        Rtabcolumn[rowindex] = Rtabcolumn[rowindex-1];
      }
    }

  }
}

static GtUword evaluateallcolumns(GtUword *EDtabcolumn,
                                  GtUword *Rtabcolumn,
                                  GtUword midcol,
                                  const GtUchar *useq,
                                  const GtUchar *vseq,
                                  GtUword ulen, GtUword vlen,
                                  const GtWord matchcost,
                                  const GtWord mismatchcost,
                                  const GtWord gapcost)
{
  GtUword colindex;

  firstEDtabRtabcolumn(EDtabcolumn, Rtabcolumn, ulen, gapcost);
  for (colindex = 1UL; colindex <= vlen; colindex++)
  {
    nextEDtabRtabcolumn(EDtabcolumn, Rtabcolumn, colindex, midcol,
                        vseq[colindex-1], useq, ulen, matchcost,
                        mismatchcost,gapcost);
  }
  return EDtabcolumn[ulen];
}

static GtUword evaluatecrosspoints(const GtUchar *useq,
                                   const GtUchar *vseq,
                                   GtUword ulen,
                                   GtUword vlen,
                                   GtUword *EDtabcolumn,
                                   GtUword *Rtabcolumn,
                                   GtUword *Ctab,
                                   GtUword rowoffset,
                                   const GtWord matchcost,
                                   const GtWord mismatchcost,
                                   const GtWord gapcost)
{
  GtUword midrow, midcol, distance;

  if (vlen >= 2UL)
  {
    midcol = GT_DIV2(vlen);
    distance = evaluateallcolumns(EDtabcolumn, Rtabcolumn, midcol, useq, vseq,
                                  ulen, vlen, matchcost, mismatchcost, gapcost);
    midrow = Rtabcolumn[ulen];
    Ctab[midcol] = rowoffset + midrow;
    (void) evaluatecrosspoints(useq,
                               vseq,
                               midrow,
                               midcol,
                               EDtabcolumn,
                               Rtabcolumn,
                               Ctab,
                               rowoffset,
                               matchcost,
                               mismatchcost,
                               gapcost);
    (void) evaluatecrosspoints(useq+midrow,
                               vseq+midcol,
                               ulen-midrow,
                               vlen-midcol,
                               EDtabcolumn,
                               Rtabcolumn,
                               Ctab+midcol,
                               rowoffset+midrow,
                               matchcost,
                               mismatchcost,
                               gapcost);
    return distance;
  }
  return 0;
}

static GtUword computealignment(const GtUchar *useq,
                                const GtUchar *vseq,
                                GtUword ulen,
                                GtUword vlen,
                                GtAlignment *align,
                                GtUword *Ctab,
                                const GtWord matchcost,
                                const GtWord mismatchcost,
                                const GtWord gapcost)
{
  GtUword distance,
          *EDtabcolumn,
          *Rtabcolumn;
  char *pos, *useq_part;
  useq_part= gt_malloc(sizeof *useq_part * (ulen+1));

  EDtabcolumn = gt_malloc(sizeof *EDtabcolumn * (ulen+1));
  Rtabcolumn = gt_malloc(sizeof *Rtabcolumn * (ulen+1));
  Ctab[vlen] = ulen;
  if (vlen == 1UL) {
    strncpy(useq_part,(const char*)useq,ulen);
    useq_part[ulen] = '\0';
    pos = strrchr(useq_part, vseq[0]);
    if(pos == NULL){
      Ctab[0] = Ctab[1]-1;
      distance = mismatchcost + (ulen-1)*gapcost;
    }
    else{
      Ctab[0] = (pos-useq_part);
      distance = matchcost + (ulen-1)*gapcost;
    }
  }
  else{
    distance = evaluatecrosspoints(useq, vseq, ulen, vlen,
                                   EDtabcolumn, Rtabcolumn,
                                   Ctab, 0, matchcost,
                                   mismatchcost, gapcost);
    Ctab[0] = 0;
  }
  reconstructalignment(align, Ctab, vlen);
  gt_free(EDtabcolumn);
  gt_free(Rtabcolumn);
  return distance;
}

GtUword gt_calc_linearalign_with_costs(const GtUchar *useq, GtUword ulen,
                            const GtUchar *vseq, GtUword vlen,
                            GtAlignment *align,
                            const GtWord matchcost,
                            const GtWord mismatchcost,
                            const GtWord gapcost)
{
  GtUword *Ctab, edist;

  Ctab = gt_malloc(sizeof *Ctab * (vlen+1));
  edist = computealignment(useq, vseq, ulen, vlen, align,
                          Ctab, matchcost, mismatchcost, gapcost);

  gt_free(Ctab);
  return edist;
}

void gt_computelinearspace2(bool showevalue,
                           const GtUchar *useq,
                           GtUword ulen,
                           const GtUchar *vseq,
                           GtUword vlen,
                           const GtWord matchcost,
                           const GtWord mismatchcost,
                           const GtWord gapcost,
                           FILE *fp)
{
  GtAlignment *align;
  GtUword distance;

  gt_assert(useq && ulen && vseq && vlen);
  if (matchcost < 0 || mismatchcost < 0 || gapcost < 0)
  {
    fprintf(stderr,"invalid cost value");
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  align = gt_alignment_new_with_seqs(useq, ulen, vseq, vlen);
  distance =gt_calc_linearalign_with_costs(useq, ulen, vseq, vlen, align,
                                 matchcost, mismatchcost, gapcost);

  gt_alignment_show(align, fp, 80);
  if(showevalue)
    fprintf(fp, "linear costs: "GT_WU"\n", distance);
  gt_alignment_delete(align);
}
