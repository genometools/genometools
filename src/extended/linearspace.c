/*
  Copyright (c) 2015 Annika <annika.seidel@studium.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#define LINEAR_EDIST_GAP          ((GtUchar) UCHAR_MAX)

static GtUword alignment_in_square_space(GtAlignment *align,
                                         const GtUchar *useq,
                                         const GtUword ustart,
                                         const GtUword ulen,
                                         const GtUchar *vseq,
                                         const GtUword vstart,
                                         const GtUword vlen,
                                         const GtWord matchcost,
                                         const GtWord mismatchcost,
                                         const GtWord gapcost)
{
  GtUword **E, distance=0;
  GtUword i,j, val;

  E = gt_malloc((sizeof **E)*(ulen+1));
  *E = gt_malloc((sizeof *E)*((vlen+1)*(ulen+1)));
  for (j = 1; j <= ulen; j++)
  {
    E[j] = E[j-1]+vlen+1;
  }

  E[0][0] = 0;
  for (i = 1; i <= ulen; i++)
  {
      E[i][0] = E[i-1][0] + gapcost;
  }

  for (j = 1; j <= vlen; j++)
  {
      E[0][j] = E[0][j-1] + gapcost;
      for (i = 1; i <= ulen; i++)
      {
        E[i][j] = E[i][j-1];

        if ((val = E[i-1][j-1] + (useq[ustart+i-1] == vseq[vstart+j-1] ?
                                  matchcost : mismatchcost))
            <= E[i][j])
        {
          E[i][j] = val;
        }

        if ((val = E[i-1][j] + gapcost) < E[i][j])
        {
          E[i][j] = val;
        }
     }
  }

  i = ulen;
  j = vlen;
  distance = E[i][j];
  while ( i != 0 || j != 0)
  {
    if (i != 0 && j != 0 && E[i][j] == E[i-1][j-1] +
            (useq[ustart+i-1] == vseq[vstart+j-1] ?
                         matchcost : mismatchcost))
    {
      gt_alignment_add_replacement(align);
      i--; j--;
    }
    else if (j != 0 && E[i][j] == E[i][j-1] + gapcost)
    {
      gt_alignment_add_insertion(align);
      j--;
    }
    else if (i != 0 && E[i][j] == E[i-1][j] + gapcost)
    {
      gt_alignment_add_deletion(align);
      i--;
    }
    else
    {
      /*never reach this line*/
      fprintf(stderr,"the impossible happend\n");
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
  gt_free(E[0]);
  gt_free(E);
  return distance;
}
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
    EDtabcolumn[rowindex] = EDtabcolumn[rowindex-1]+ gapcost;
    Rtabcolumn[rowindex]  = rowindex;
  }
}

static void nextEDtabRtabcolumn(GtUword *EDtabcolumn,
                                GtUword *Rtabcolumn,
                                GtUword colindex,
                                GtUword midcolumn, GtUchar b,
                                const GtUchar *useq,
                                const GtUword ustart, const GtUword ulen,
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
    EDtabcolumn[rowindex] += gapcost; /* 1. recurrence */

    /* 2. recurrence: */
    if ((val = northwestEDtabentry +
        (useq[ustart+rowindex-1] == b ? matchcost : mismatchcost))
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
                                  const GtUword ustart, const GtUword ulen,
                                  const GtUchar *vseq,
                                  const GtUword vstart, const GtUword vlen,
                                  const GtWord matchcost,
                                  const GtWord mismatchcost,
                                  const GtWord gapcost)
{
  GtUword colindex;

  firstEDtabRtabcolumn(EDtabcolumn, Rtabcolumn, ulen, gapcost);
  for (colindex = 1UL; colindex <= vlen; colindex++)
  {
    nextEDtabRtabcolumn(EDtabcolumn, Rtabcolumn, colindex, midcol,
                        vseq[vstart+colindex-1], useq, ustart,
                        ulen, matchcost, mismatchcost,gapcost);
  }
  return EDtabcolumn[ulen];
}

static GtUword evaluatecrosspoints(const GtUchar *useq,
                                   const GtUword ustart, GtUword ulen,
                                   const GtUchar *vseq,
                                   const GtUword vstart, GtUword vlen,
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
    distance = evaluateallcolumns(EDtabcolumn, Rtabcolumn, midcol,
                                  useq, ustart, ulen, vseq, vstart, vlen,
                                  matchcost, mismatchcost, gapcost);
    midrow = Rtabcolumn[ulen];
    Ctab[midcol] = rowoffset + midrow;
    (void) evaluatecrosspoints(useq, ustart, midrow,
                               vseq, vstart, midcol,
                               EDtabcolumn,
                               Rtabcolumn,
                               Ctab,
                               rowoffset,
                               matchcost,
                               mismatchcost,
                               gapcost);
    (void) evaluatecrosspoints(useq, ustart + midrow,
                               ulen-midrow,
                               vseq, vstart + midcol,
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

static void determineCtab0(GtUword *Ctab, GtUchar vseq0,
                           const GtUchar *useq, const GtUword ustart)
{
  GtUword rowindex;

  for (rowindex = 0; rowindex < Ctab[1]; rowindex++)
  {
    if (vseq0 == useq[ustart+rowindex])
    {
      Ctab[0] = rowindex;
      return;
    }
  }

  Ctab[0] = (Ctab[1] > 0) ?  Ctab[1]-1 : 0;
}

GtUword gt_calc_linearalign2(const GtUchar *useq, const GtUword ustart,
                             const GtUword ulen,
                             const GtUchar *vseq, const GtUword vstart,
                             const GtUword vlen,
                             GtAlignment *align,
                             const GtWord matchcost,
                             const GtWord mismatchcost,
                             const GtWord gapcost)
{
  GtUword distance, *Ctab, *EDtabcolumn,*Rtabcolumn;

  if (ulen == 0UL)
  {
      distance = construct_trivial_insertion_alignment(align,vlen,
                                                       gapcost);
  }
  else if (vlen == 0UL)
  {
      distance = construct_trivial_deletion_alignment(align,vlen,
                                                       gapcost);
  }
  else if (ulen == 1UL || vlen == 1UL ) {
    distance = alignment_in_square_space(align, useq, ustart, ulen,
                                         vseq, vstart, vlen, matchcost,
                                         mismatchcost, gapcost);
  }
  else
  {
    Ctab = gt_malloc(sizeof *Ctab * (vlen+1));
    EDtabcolumn = gt_malloc(sizeof *EDtabcolumn * (ulen+1));
    Rtabcolumn = gt_malloc(sizeof *Rtabcolumn * (ulen+1));

    Ctab[vlen] = ulen;
    distance = evaluatecrosspoints(useq, ustart, ulen,
                                   vseq, vstart, vlen,
                                   EDtabcolumn, Rtabcolumn,
                                   Ctab, 0, matchcost,
                                   mismatchcost, gapcost);
    determineCtab0(Ctab,vseq[vstart],useq, ustart);

    reconstructalignment_from_Ctab(align, Ctab, vlen);
    gt_free(Ctab);
    gt_free(EDtabcolumn);
    gt_free(Rtabcolumn);
  }
  return distance;
}

GtAlignment *gt_computelinearspace2(const GtUchar *useq,
                            GtUword ustart, GtUword ulen,
                            const GtUchar *vseq,
                            GtUword vstart,GtUword vlen,
                            const GtWord matchcost,
                            const GtWord mismatchcost,
                            const GtWord gapcost)
{
  GtAlignment *align;

  gt_assert(useq && ulen && vseq && vlen);
  if (matchcost < 0 || mismatchcost < 0 || gapcost < 0)
  {
    fprintf(stderr,"invalid cost value");
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  align = gt_alignment_new_with_seqs(useq+ustart, ulen, vseq+vstart, vlen);
  (void)gt_calc_linearalign2(useq, ustart, ulen, vseq, vstart, vlen, align,
                                 matchcost, mismatchcost, gapcost);

  return align;
}
