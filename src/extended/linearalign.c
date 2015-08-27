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

#include <string.h>
#include "core/ma.h"
#include "core/minmax.h"
#include "core/assert_api.h"
#include "core/unused_api.h"
#include "core/divmodmul.h"
#include "match/squarededist.h"
#include "extended/alignment.h"
#include "extended/linearalign.h"
#include "extended/reconstructalignment.h"
#include "extended/maxcoordvalue.h"

#define LINEAR_EDIST_GAP          ((GtUchar) UCHAR_MAX)

static void change_score_to_cost_function(const GtWord matchscore,
                                          const GtWord mismatchscore,
                                          const GtWord gapscore,
                                          GtWord *matchcost,
                                          GtWord *mismatchcost,
                                          GtWord *gapcost )
{
  GtWord max;

  max = (MAX3(GT_DIV2(matchscore), GT_DIV2(mismatchscore), 1 + gapscore));

  if (max < 0)
    max = 0;

  *matchcost = 2 * max-matchscore;
  *mismatchcost = 2 * max-mismatchscore;
  *gapcost = max-gapscore;

}

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

  E = gt_malloc((sizeof *E)*(ulen+1));
  *E = gt_malloc((sizeof **E)*((vlen+1)*(ulen+1)));
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
        E[i][j] = E[i][j-1] + gapcost;

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

/*------------------------------global--------------------------------*/
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
  westEDtabentry = EDtabcolumn[0];
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
                                  const GtUword ustart,
                                  const GtUword ulen,
                                  const GtUchar *vseq,
                                  const GtUword vstart,
                                  const GtUword vlen,
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

GtUword gt_calc_linearalign(const GtUchar *useq,
                            const GtUword ustart,
                            const GtUword ulen,
                            const GtUchar *vseq,
                            const GtUword vstart,
                            const GtUword vlen,
                            GtAlignment *align,
                            const GtWord matchcost,
                            const GtWord mismatchcost,
                            const GtWord gapcost)
{
  GtUword distance, *Ctab, *EDtabcolumn,*Rtabcolumn;

  if (ulen == 0UL)
  {
    distance = construct_trivial_insertion_alignment(align,vlen,gapcost);
  }
  else if (vlen == 0UL)
  {
    distance = construct_trivial_deletion_alignment(align,vlen,gapcost);
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

    /* reconstructalignment_from_Ctab(align, Ctab, vlen); */

    reconstructalignment_from_Ctab(align,Ctab,useq,ustart,vseq,
                                   vstart,vlen,matchcost,mismatchcost,
                                   0,gapcost);
    gt_free(Ctab);
    gt_free(EDtabcolumn);
    gt_free(Rtabcolumn);
  }
  return distance;
}

GtAlignment *gt_computelinearspace(const GtUchar *useq,
                            GtUword ustart, GtUword ulen,
                            const GtUchar *vseq,
                            GtUword vstart, GtUword vlen,
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
  (void)gt_calc_linearalign(useq, ustart, ulen, vseq, vstart, vlen, align,
                                 matchcost, mismatchcost, gapcost);

  return align;
}

void gt_print_edist_alignment(const GtUchar *useq, GtUword ustart,
                              GtUword ulen,
                              const GtUchar *vseq,GtUword vstart,
                               GtUword vlen)
{
  GtAlignment *align;

  align = gt_alignment_new_with_seqs(useq+ustart, ulen, vseq+vstart, vlen);
  (void)gt_calc_linearalign(useq, ustart, ulen, vseq, vstart, vlen, align,
                             0,1,1);
  gt_alignment_show(align, stdout, 80);
  gt_alignment_delete(align);
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

/*------------------------------local---------------------------------*/
static void firstLStabcolumn(const GtUword ulen,
                             GtWord *Ltabcolumn,
                             GtUwordPair *Starttabcolumn)
{
  GtUword rowindex;

  gt_assert(Ltabcolumn != NULL && Starttabcolumn != NULL);
  for (rowindex = 0; rowindex <= ulen; rowindex++)
  {
    Ltabcolumn[rowindex] = 0;
    Starttabcolumn[rowindex].a = rowindex;
    Starttabcolumn[rowindex].b = 0;
  }
}

static void nextLStabcolumn(const GtUchar *useq,
                            const GtUword ustart, const GtUword ulen,
                            const GtUchar b, const GtUword colindex,
                            GtWord *Ltabcolumn,
                            GtUwordPair *Starttabcolumn,
                            Gtmaxcoordvalue *max,
                            const GtWord matchscore,
                            const GtWord mismatchscore,
                            const GtWord gapscore)
{
  GtUword rowindex;
  GtUwordPair northwestStarttabentry, westStarttabentry;
  GtWord val, northwestLtabentry, westLtabentry;

  gt_assert(max != NULL);
  westLtabentry = Ltabcolumn[0];
  westStarttabentry = Starttabcolumn[0];

  Ltabcolumn[0] = 0;
  Starttabcolumn[0].a = 0;
  Starttabcolumn[0].b = colindex;
  for (rowindex = 1UL; rowindex <= ulen; rowindex++)
  {
    northwestLtabentry = westLtabentry;
    northwestStarttabentry = westStarttabentry;
    westLtabentry = Ltabcolumn[rowindex];
    westStarttabentry = Starttabcolumn[rowindex];
    Ltabcolumn[rowindex] += gapscore;

    if ((val = northwestLtabentry + (useq[ustart + rowindex-1] ==
         b ? matchscore : mismatchscore)) > Ltabcolumn[rowindex])
    {
      Ltabcolumn[rowindex] = val;
      Starttabcolumn[rowindex] = northwestStarttabentry;
    }
    if ((val = Ltabcolumn[rowindex-1]+gapscore) > Ltabcolumn[rowindex])
    {
      Ltabcolumn[rowindex] = val;
      Starttabcolumn[rowindex]=Starttabcolumn[rowindex-1];
    }
    if (0 > Ltabcolumn[rowindex])
    {
      Ltabcolumn[rowindex] = 0;
      Starttabcolumn[rowindex].a = rowindex;
      Starttabcolumn[rowindex].b = colindex;
    }
    if (Ltabcolumn[rowindex] > gt_max_get_value(max))
    {
      gt_max_coord_update(max, Ltabcolumn[rowindex],
                               Starttabcolumn[rowindex],
                               rowindex, colindex);
    }
  }
}

static Gtmaxcoordvalue *evaluateallLScolumns(GtWord *Ltabcolumn,
                                             GtUwordPair *Starttabcolumn,
                                             const GtUchar *useq,
                                             const GtUword ustart,
                                             const GtUword ulen,
                                             const GtUchar *vseq,
                                             const GtUword vstart,
                                             const GtUword vlen,
                                             const GtWord matchscore,
                                             const GtWord mismatchscore,
                                             const GtWord gapscore)
{
  GtUword colindex;
  Gtmaxcoordvalue *max;

  firstLStabcolumn(ulen, Ltabcolumn, Starttabcolumn);
  max = gt_max_new();
  for (colindex = 1UL; colindex <= vlen; colindex++)
  {
    nextLStabcolumn(useq, ustart, ulen, vseq[vstart+colindex-1], colindex,
                    Ltabcolumn, Starttabcolumn, max,
                    matchscore, mismatchscore, gapscore);
  }
  return max;
}

GtWord gt_calc_linearalign_local(const GtUchar *useq,
                                 const GtUword ustart, const GtUword ulen,
                                 const GtUchar *vseq,
                                 const GtUword vstart, const GtUword vlen,
                                 GtAlignment *align,
                                 const GtWord matchscore,
                                 const GtWord mismatchscore,
                                 const GtWord gapscore)
{
  GtWord *Ltabcolumn;
  GtUwordPair *Starttabcolumn;
  GtUword ulen_part, ustart_part, vlen_part, vstart_part,score;

  Gtmaxcoordvalue *max;
  GtWord matchcost, mismatchcost, gapcost;

  Ltabcolumn = gt_malloc(sizeof *Ltabcolumn * (ulen+1));
  Starttabcolumn = gt_malloc(sizeof *Starttabcolumn * (ulen+1));

  max = evaluateallLScolumns(Ltabcolumn,
                             Starttabcolumn,
                             useq, ustart, ulen,
                             vseq, vstart, vlen,
                             matchscore,
                             mismatchscore,
                             gapscore);
  gt_free(Ltabcolumn);
  gt_free(Starttabcolumn);
  change_score_to_cost_function(matchscore,
                                mismatchscore,
                                gapscore,
                                &matchcost,
                                &mismatchcost,
                                &gapcost );

  if (gt_max_get_length_safe(max))
  {
    ustart_part = ustart+(gt_max_get_start(max)).a;
    vstart_part = vstart+(gt_max_get_start(max)).b;
    ulen_part = gt_max_get_row_length(max);
    vlen_part = gt_max_get_col_length(max);

    gt_alignment_set_seqs(align, &useq[ustart_part], ulen_part,
                                 &vseq[vstart_part], vlen_part);
    gt_calc_linearalign(useq, ustart_part, ulen_part,
                         vseq, vstart_part, vlen_part,
                         align, matchcost, mismatchcost, gapcost);
    score = gt_alignment_eval_with_score(align, matchscore,
                                         mismatchscore, gapscore);
  }else
  {
    gt_alignment_set_seqs(align,(const GtUchar*)"",0,(const GtUchar*)"",0);
    score = 0;
  }

  gt_max_delete(max);

  return score;
}

GtAlignment *gt_computelinearspace_local(const GtUchar *useq,
                                         const GtUword ustart,
                                         const GtUword ulen,
                                         const GtUchar *vseq,
                                         const GtUword vstart,
                                         const GtUword vlen,
                                         const GtWord matchscore,
                                         const GtWord mismatchscore,
                                         const GtWord gapscore)
{
  GtAlignment *align;

  align = gt_alignment_new();
  (void) gt_calc_linearalign_local(useq, ustart, ulen, vseq, vstart, vlen,
                                   align, matchscore, mismatchscore, gapscore);

  return align;
}

/*-------------------------checkfunctions-----------------------------*/
void gt_checklinearspace(GT_UNUSED bool forward,
                         const GtUchar *useq,
                         GtUword ulen,
                         const GtUchar *vseq,
                         GtUword vlen)
{
  GtAlignment *align;
  GtUword edist1, edist2, edist3, edist4;

  /*gt_assert(useq && ulen && vseq && vlen);*/
  if (strchr((const char*)useq, LINEAR_EDIST_GAP))
  {
    fprintf(stderr,"%s: sequence u contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  if (strchr((const char*)vseq, LINEAR_EDIST_GAP))
  {
    fprintf(stderr,"%s: sequence v contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  align = gt_alignment_new_with_seqs(useq, ulen, vseq, vlen);
  edist1 = gt_calc_linearalign(useq, 0, ulen, vseq, 0, vlen,
                               align, 0, 1, 1);
  edist2 = gt_squarededistunit(useq,ulen,vseq,vlen);

  if (edist1 != edist2)
  {
    fprintf(stderr,"gt_calc_linearalign2 = "GT_WU" != "GT_WU
            " = gt_squarededistunit\n", edist1,edist2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  edist3 = gt_alignment_eval_with_score(align,0,1,1);
  if (edist2 != edist3)
  {
    fprintf(stderr,"gt_squarededistunit = "GT_WU" != "GT_WU
            " = gt_alignment_eval_with_score\n", edist2,edist3);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  edist4 = gt_calc_linearedist(useq, ulen, vseq, vlen);
  if (edist3 != edist4)
  {
    fprintf(stderr,"gt_alignment_eval_with_score = "GT_WU" != "GT_WU
            " = gt_calc_linearedist\n", edist3, edist4);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  gt_alignment_delete(align);
}

void gt_checklinearspace_local(GT_UNUSED bool forward,
                               const GtUchar *useq, GtUword ulen,
                               const GtUchar *vseq, GtUword vlen)
{
  GtAlignment *align;
  GtWord score1, score2;

  /*gt_assert(useq && ulen && vseq && vlen);*/
  if (strchr((const char*)useq, LINEAR_EDIST_GAP))
  {
    fprintf(stderr,"%s: sequence u contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  if (strchr((const char*)vseq, LINEAR_EDIST_GAP))
  {
    fprintf(stderr,"%s: sequence v contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  align = gt_alignment_new();
  score1 = gt_calc_linearalign_local(useq, 0, ulen,
                                     vseq, 0, vlen, align, 2, -2, -1);

  score2 = gt_alignment_eval_with_score(align,2, -2, -1);

  if (score1 != score2)
  {
    fprintf(stderr,"gt_calc_linearalign_local = "GT_WD" != "GT_WD
              " = gt_alignment_eval_with_score\n", score1, score2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  gt_alignment_delete(align);
}
