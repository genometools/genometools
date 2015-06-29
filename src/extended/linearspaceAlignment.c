#include "core/ma.h"
#include "core/minmax.h"
#include "core/assert_api.h"
#include "core/unused_api.h"
#include "core/divmodmul.h"
#include "match/squarededist.h"
#include "extended/alignment.h"
#include "extended/linearspaceAlignment.h"

static void firstEDtabRtabcolumn(GtUword *EDtabcolumn,
                                 GtUword *Rtabcolumn,
                                 GtUword ulen,
                                 const GtWord gapcost)
{
  GtUword rowindex;

  for (rowindex=0; rowindex <= ulen; rowindex++)
  {
    EDtabcolumn[rowindex] = gapcost;
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
          westRtabentry = 0;
  bool updateRtabcolumn = false;

  gt_assert(EDtabcolumn != NULL);
  westEDtabentry = EDtabcolumn[0]; /* saves the first entry of EDtabcolumn */
  EDtabcolumn[0]++;
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
    if ((val = northwestEDtabentry + (useq[rowindex-1] == b ? matchcost : mismatchcost)) <
        EDtabcolumn[rowindex])
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
                        vseq[colindex-1], useq, ulen, matchcost, mismatchcost,gapcost);
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

static GtUword determineCtab0(GtUword *Ctab,
                              GtUchar vseq0,
                              const GtUchar *useq)
{
  GtUword rowindex;

  for (rowindex=0; rowindex < Ctab[1]; rowindex++)
  {
    if (vseq0 == useq[rowindex])
    {
      Ctab[0] = rowindex;
      return Ctab[1] - 1;
    }
  }

  Ctab[0] = (Ctab[1] > 0) ?  Ctab[1]-1 : 0;
  return Ctab[1];
}

static void reconstructalignment(GtAlignment *align,
                                 GtUword *Ctab,
                                 GtUword vlen)
{
  GtUword i,j;

  gt_assert(align != NULL && Ctab != NULL);
  for (i = vlen; i > 0; i--) {
    if (Ctab[i] == Ctab[i-1] + 1)
      gt_alignment_add_replacement(align);
    else if (Ctab[i] == Ctab[i-1])
      gt_alignment_add_insertion(align);
    else if (Ctab[i] > Ctab[i-1]) {
      for (j = 0; j < (Ctab[i]-Ctab[i-1])-1; j++)
        gt_alignment_add_deletion(align);
      gt_alignment_add_replacement(align);
    }
  }
  for (j = Ctab[0]; j > 0; j--)
    gt_alignment_add_deletion(align);

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

  EDtabcolumn = gt_malloc(sizeof *EDtabcolumn * (ulen+1));
  Rtabcolumn = gt_malloc(sizeof *Rtabcolumn * (ulen+1));
  Ctab[vlen] = ulen;
  if (vlen == 1UL) {
    distance = determineCtab0(Ctab, vseq[0], useq);
  }
  else{
    distance = evaluatecrosspoints(useq, vseq, ulen, vlen, EDtabcolumn,
                                   Rtabcolumn, Ctab, 0, matchcost, mismatchcost, gapcost);
    (void) determineCtab0(Ctab, vseq[0], useq);
  }
  reconstructalignment(align,Ctab, vlen);
  gt_free(EDtabcolumn);
  gt_free(Rtabcolumn);
  return distance;
}

GtUword gt_computelinearspace_with_costs(const GtUchar *u, GtUword ulen,
                            const GtUchar *v, GtUword vlen,
                            GtAlignment *align,
                            const GtWord matchcost,
                            const GtWord mismatchcost,
                            const GtWord gapcost)
{
  GtUword *Ctab, edist;

  Ctab = gt_malloc(sizeof *Ctab * (vlen+1));
  edist = computealignment(u, v, ulen, vlen, align, Ctab,matchcost, mismatchcost, gapcost);
  gt_free(Ctab);
  return edist;
}

static void alignment_to_stdout(const GtAlignment *align)
{
    gt_assert(align != NULL);
    gt_alignment_show(align, stdout, 80);
}

void gt_computelinearspace_with_output(const GtUchar *useq,
                           GtUword ulen,
                           const GtUchar *vseq,
                           GtUword vlen,
                           const GtWord matchcost,
                           const GtWord mismatchcost,
                           const GtWord gapcost)
{
  GtAlignment *align;

  align = gt_alignment_new_with_seqs(useq, ulen, vseq, vlen);
  gt_computelinearspace_with_costs(useq, ulen, vseq, vlen, align, matchcost, mismatchcost, gapcost);

  alignment_to_stdout(align);
  gt_alignment_delete(align);
}
