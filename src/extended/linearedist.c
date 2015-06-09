/*
  Copyright (C) 2015 Annika Seidel, annika.seidel@studium.uni-hamburg.de
  Copyright (C) 2015 Stefan Kurtz, kurtz@zbh.uni-hamburg.de
  Copyright (C) 2015 Joerg Winkler, joerg.winkler@studium.uni-hamburg.de
  Copyright (C) 2014 Dirk Willrodt, willrodt@zbh.uni-hamburg.de
  Copyright (C) 2010 Sascha Steinbiss, steinbiss@zbh.uni-hamburg.de
  Copyright (c) 2006-2009 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "core/minmax.h"
#include "core/assert_api.h"
#include "core/unused_api.h"
#include "core/divmodmul.h"
#include "match/squarededist.h"
#include "extended/linearedist.h"

#define LINEAR_EDIST_GAP          ((GtUchar) UCHAR_MAX)

/*
   The following function computes the first column of the E and R table as
   described in the handout of the lecture on ``A Linear Space Alignmen
   Algorithm''.  EDtabcolumn takes a column of the E table and should
   contain the first column after completion of the function. Analog,
   Rtabcolumn takes a column of the R table and should contain the first column
   after completion of the function (this is the initialization!).
   ulen is the length of the first sequence.
   */
static void firstEDtabRtabcolumn(GtUword *EDtabcolumn,
                                 GtUword *Rtabcolumn,
                                 GtUword ulen)
{
  GtUword rowindex;

  for (rowindex=0; rowindex <= ulen; rowindex++)
  {
    EDtabcolumn[rowindex] = rowindex;
    Rtabcolumn[rowindex]  = rowindex;
  }
}

/*
   The following function computes a column of the E and R table under the
   assumption that EDtabcolumn and Rtabcolumn contain the previous columns
   of the corresponding tables. The function overwrites EDtabcolumn and
   Rtabcolumn with the new columns. b is the character of the
   sequence vseq corresponding to to column which needs to be calculated. Useq
   is the first sequence and ulen its length.

   colindex is the number of the column to be calculated and midcolumn
   is the column, where the tables are divided.
   */
static void nextEDtabRtabcolumn(GtUword *EDtabcolumn,
                                GtUword *Rtabcolumn,
                                GtUword colindex,
                                GtUword midcolumn, GtUchar b,
                                const GtUchar *useq, GtUword ulen)
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
    EDtabcolumn[rowindex]++; /* 1. recurrence */
    /* Rtabcolumn[rowindex] is unchanged */
    /* 2. recurrence: */
    if ((val = northwestEDtabentry + (useq[rowindex-1] == b ? 0 : 1)) <
        EDtabcolumn[rowindex])
    {
      EDtabcolumn[rowindex] = val;
      if (updateRtabcolumn)
      {
        Rtabcolumn[rowindex] = northwestRtabentry;
      }
    }
    /* 3. recurrence: */
    if ((val = EDtabcolumn[rowindex-1]+1) < EDtabcolumn[rowindex])
    {
      EDtabcolumn[rowindex] = val;
      if (updateRtabcolumn)
      {
        Rtabcolumn[rowindex] = Rtabcolumn[rowindex-1];
      }
    }
  }
}

/*
   The following function computes the rightmost EDtabcolumn and Rtabcolumn.
   midcolumn is the actual column, where the tables are divided.
   useq and vseq are the sequences with the corresponding length ulen and vlen.
*/
static GtUword evaluateallcolumns(GtUword *EDtabcolumn,
                                  GtUword *Rtabcolumn,
                                  GtUword midcol,
                                  const GtUchar *useq,
                                  const GtUchar *vseq,
                                  GtUword ulen, GtUword vlen)
{
  GtUword colindex;

  firstEDtabRtabcolumn(EDtabcolumn, Rtabcolumn, ulen);
  for (colindex = 1UL; colindex <= vlen; colindex++)
  {
    nextEDtabRtabcolumn(EDtabcolumn, Rtabcolumn, colindex, midcol,
                        vseq[colindex-1], useq, ulen);
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
                                   GtUword rowoffset)
{
  GtUword midrow, midcol, distance;

  if (vlen >= 2UL)
  {
    midcol = GT_DIV2(vlen);
    distance = evaluateallcolumns(EDtabcolumn, Rtabcolumn, midcol, useq, vseq,
                                  ulen, vlen);
    midrow = Rtabcolumn[ulen];
    Ctab[midcol] = rowoffset + midrow;
    (void) evaluatecrosspoints(useq,
                               vseq,
                               midrow,
                               midcol,
                               EDtabcolumn,
                               Rtabcolumn,
                               Ctab,
                               rowoffset);
    (void) evaluatecrosspoints(useq+midrow,
                               vseq+midcol,
                               ulen-midrow,
                               vlen-midcol,
                               EDtabcolumn,
                               Rtabcolumn,
                               Ctab+midcol,
                               rowoffset+midrow);
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

static GtUword reconstructalignment(GtAlignment *align,
                                       GtUword *Ctab,
                                       GtUword vlen)
{
  GtUword alilen, i,j;

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

  alilen = gt_alignment_get_length(align);
  return alilen;
}

static GtUword computealignment(const GtUchar *useq,
                                const GtUchar *vseq,
                                GtUword ulen,
                                GtUword vlen,
                                GtAlignment *align,
                                GtUword *alilen,
                                GtUword *Ctab)
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
                                   Rtabcolumn, Ctab, 0);
    (void) determineCtab0(Ctab, vseq[0], useq);
  }
  *alilen = reconstructalignment(align, Ctab, vlen);
  gt_free(EDtabcolumn);
  gt_free(Rtabcolumn);
  return distance;
}

GtUword gt_calc_linearalign(const GtUchar *u, GtUword ulen,
                            const GtUchar *v, GtUword vlen,
                            GtAlignment *align, GtUword *alilen)
{
  GtUword *Ctab, edist;

  Ctab = gt_malloc(sizeof *Ctab * (vlen+1));
  edist = computealignment(u, v, ulen, vlen, align, alilen, Ctab);
  gt_free(Ctab);
  return edist;
}

/* just calculate distance, no alignment */
static void fillDPtable(GtUword *dpcolumn,
                        const GtUchar *u, GtUword ulen,
                        const GtUchar *v, GtUword vlen)
{
  GtUword i, j , nw, we;
  for (i = 0; i <= ulen; i++)
    dpcolumn[i] = i;
  for (j = 1UL; j <= vlen; j++) {
    nw = dpcolumn[0];
    dpcolumn[0] = j;
    for (i = 1UL; i <= ulen; i++) {
      we = dpcolumn[i];
      dpcolumn[i] = nw + (u[i-1] == v[j-1] ? 0 : 1); /* replacement */
      if (dpcolumn[i-1] + 1 < dpcolumn[i]) /* deletion */
        dpcolumn[i] = dpcolumn[i-1] + 1;
      if (we + 1 < dpcolumn[i]) /* insertion */
        dpcolumn[i] = we + 1;
      nw = we;
    }
  }
}

GtUword gt_calc_linearedist(const GtUchar *u, GtUword ulen,
                            const GtUchar *v, GtUword vlen)
{
  GtUword *dpcolumn, edist;

  dpcolumn = gt_malloc(sizeof *dpcolumn * (MIN(ulen,vlen) + 1));
  fillDPtable(dpcolumn, ulen <= vlen ? u : v, MIN(ulen,vlen),
                       ulen <= vlen ? v : u, MAX(ulen,vlen));
  edist = dpcolumn[MIN(ulen,vlen)];
  gt_free(dpcolumn);
  return edist;
}

static GtUword evaluate_alcost(const GtAlignment *align)
{
  GtUword alcost = 0;

  gt_assert(align != NULL);
  alcost = gt_alignment_eval(align);

  return alcost;
}

static bool gap_symbol_in_sequence(const GtUchar *seq, GtUword len)
{
  const GtUchar *sptr;

  for (sptr = seq; sptr < seq + len; sptr++)
  {
    if (*sptr == LINEAR_EDIST_GAP)
    {
      return true;
    }
  }
  return false;
}

void gt_checklinearspace(GT_UNUSED bool forward,
                         const GtUchar *useq,
                         GtUword ulen,
                         const GtUchar *vseq,
                         GtUword vlen)
{
  GtAlignment *align;
  GtUword  alcost, alilen, edist1, edist2, edist3;

  if (gap_symbol_in_sequence(useq,ulen))
  {
    fprintf(stderr,"%s: sequence u contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  if (gap_symbol_in_sequence(vseq,vlen))
  {
    fprintf(stderr,"%s: sequence v contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  edist1 = gt_calc_linearedist(useq,ulen,vseq,vlen);
  edist2 = gt_squarededistunit (useq,ulen,vseq,vlen);
  if (edist1 != edist2)
  {
    fprintf(stderr,"gt_calc_linearedist = "GT_WU" != "GT_WU
            " = gt_squarededistunit\n", edist1,edist2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  align = gt_alignment_new_with_seqs(useq, ulen, vseq, vlen);
  edist3 = gt_calc_linearalign(useq, ulen, vseq, vlen, align, &alilen);
  if (edist2 != edist3)
  {
    fprintf(stderr,"gt_calc_linearalign = "GT_WU" != "GT_WU
            " = gt_squarededistunit\n", edist3,edist2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  alcost = evaluate_alcost(align);
  if (edist2 != alcost)
  {
    fprintf(stderr,"evaluate_alcost= "GT_WU" != "GT_WU
            " = gt_squarededistunit\n", alcost, edist2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  gt_alignment_delete(align);
}
