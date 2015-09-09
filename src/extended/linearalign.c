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

#include <ctype.h>
#include <string.h>
#include "core/ma.h"
#include "core/minmax.h"
#include "core/array2dim_api.h"
#include "core/assert_api.h"
#include "core/unused_api.h"
#include "core/divmodmul.h"
#include "match/squarededist.h"
#include "extended/alignment.h"
#include "extended/reconstructalignment.h"
#include "extended/maxcoordvalue.h"

#include "extended/linearalign.h"
#define LINEAR_EDIST_GAP          ((GtUchar) UCHAR_MAX)

static void change_score_to_cost_function(GtWord matchscore,
                                          GtWord mismatchscore,
                                          GtWord gapscore,
                                          GtUword *matchcost,
                                          GtUword *mismatchcost,
                                          GtUword *gapcost )
{
  GtWord max;

  max = (MAX3(GT_DIV2(matchscore), GT_DIV2(mismatchscore), 1 + gapscore));
  if (max < 0)
    max = 0;
  *matchcost = 2 * max-matchscore;
  *mismatchcost = 2 * max-mismatchscore;
  *gapcost = max-gapscore;
}

static void fillDPtab_in_square_space(GtUword **E,
                                      const GtUchar *useq,
                                      GtUword ustart,
                                      GtUword ulen,
                                      const GtUchar *vseq,
                                      GtUword vstart,
                                      GtUword vlen,
                                      GtUword matchcost,
                                      GtUword mismatchcost,
                                      GtUword gapcost)
{
  GtUword i, j;

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
        GtUword val;

        E[i][j] = E[i][j-1] + gapcost;
        if ((val = E[i-1][j-1] + (tolower((int)useq[ustart+i-1]) ==
                                  tolower((int)vseq[vstart+j-1]) ?
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
}
static GtUword alignment_in_square_space(GtAlignment *align,
                                         const GtUchar *useq,
                                         GtUword ustart,
                                         GtUword ulen,
                                         const GtUchar *vseq,
                                         GtUword vstart,
                                         GtUword vlen,
                                         GtUword matchcost,
                                         GtUword mismatchcost,
                                         GtUword gapcost)
{
  GtUword **E, distance, i, j;

  gt_assert(align != NULL);
  gt_array2dim_malloc(E, (ulen+1), (vlen+1));

  fillDPtab_in_square_space(E, useq, ustart, ulen, vseq, vstart, vlen,
                           matchcost, mismatchcost, gapcost);

  i = ulen;
  j = vlen;
  distance = E[i][j];

  /* reconstruct alignment from 2dimarray E */
  while (i > 0 || j > 0)
  {
    if (i > 0 && j > 0 && E[i][j] == E[i-1][j-1] +
            (tolower((int)useq[ustart+i-1]) == tolower((int) vseq[vstart+j-1]) ?
                         matchcost : mismatchcost))
    {
      gt_alignment_add_replacement(align);
      i--; j--;
    }
    else if (j > 0 && E[i][j] == E[i][j-1] + gapcost)
    {
      gt_alignment_add_insertion(align);
      j--;
    }
    else if (i > 0 && E[i][j] == E[i-1][j] + gapcost)
    {
      gt_alignment_add_deletion(align);
      i--;
    }
    else
    {
      gt_assert(false);
    }
  }

  gt_array2dim_delete(E);
  return distance;
}

static void ctab_in_square_space(GtUword *Ctab,
                                 const GtUchar *useq,
                                 GtUword ustart,
                                 GtUword ulen,
                                 const GtUchar *vseq,
                                 GtUword vstart,
                                 GtUword vlen,
                                 GtUword matchcost,
                                 GtUword mismatchcost,
                                 GtUword gapcost,
                                 GtUword rowoffset)
{
  GtUword **E;
  gt_assert(Ctab != NULL);

  gt_array2dim_malloc(E, (ulen+1), (vlen+1));
  fillDPtab_in_square_space(E, useq, ustart, ulen, vseq, vstart, vlen,
                            matchcost, mismatchcost, gapcost);

  evaluate_crosspoints_from_2dimtab(E, Ctab, useq, ustart, ulen,
                                    vseq, vstart, vlen, matchcost,
                                    mismatchcost, gapcost, rowoffset);
  gt_array2dim_delete(E);
}

/*------------------------------global--------------------------------*/
static void firstEDtabRtabcolumn(GtUword *EDtabcolumn,
                                 GtUword *Rtabcolumn,
                                 GtUword ulen,
                                 GtUword gapcost)
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
                                GtUword midcolumn,
                                GtUchar b,
                                const GtUchar *useq,
                                GtUword ustart,
                                GtUword ulen,
                                GtUword matchcost,
                                GtUword mismatchcost,
                                GtUword gapcost)
{
  GtUword rowindex, val,
          northwestEDtabentry,
          westEDtabentry,
          northwestRtabentry,
          westRtabentry = 0;

  gt_assert(EDtabcolumn != NULL);
  westEDtabentry = EDtabcolumn[0];
  EDtabcolumn[0] += gapcost;
  if (colindex > midcolumn)
  {
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
        (tolower((int)useq[ustart+rowindex-1]) == tolower((int)b) ?
         matchcost : mismatchcost)) <= EDtabcolumn[rowindex])
    {
      EDtabcolumn[rowindex] = val;
      if (colindex > midcolumn)
      {
        Rtabcolumn[rowindex] = northwestRtabentry;
      }
    }
    /* 3. recurrence: */
    if ((val = EDtabcolumn[rowindex-1]+gapcost) < EDtabcolumn[rowindex])
    {
      EDtabcolumn[rowindex] = val;
      if (colindex > midcolumn)
      {
        Rtabcolumn[rowindex] = Rtabcolumn[rowindex-1];
      }
    }
  }
}

static GtUword evaluateallEDtabRtabcolumns(GtUword *EDtabcolumn,
                                           GtUword *Rtabcolumn,
                                           GtUword midcol,
                                           const GtUchar *useq,
                                           GtUword ustart,
                                           GtUword ulen,
                                           const GtUchar *vseq,
                                           GtUword vstart,
                                           GtUword vlen,
                                           GtUword matchcost,
                                           GtUword mismatchcost,
                                           GtUword gapcost)
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

static GtUword evaluatelinearcrosspoints(const GtUchar *useq,
                                         GtUword ustart, GtUword ulen,
                                         const GtUchar *vseq,
                                         GtUword vstart, GtUword vlen,
                                         GtUword *EDtabcolumn,
                                         GtUword *Rtabcolumn,
                                         GtUword *Ctab,
                                         GtUword rowoffset,
                                         GtUword matchcost,
                                         GtUword mismatchcost,
                                         GtUword gapcost,
                                         GtUword original_ulen,
                                         GtUword original_vlen)
{
  GtUword midrow, midcol, distance;

  if (vlen >= 2UL)
  {
    if ((ulen+1)*(vlen+1)>(original_ulen+1))
    {
      midcol = GT_DIV2(vlen);
      distance = evaluateallEDtabRtabcolumns(EDtabcolumn, Rtabcolumn, midcol,
                                             useq, ustart, ulen, vseq, vstart,
                                             vlen, matchcost, mismatchcost,
                                             gapcost);
      midrow = Rtabcolumn[ulen];
      Ctab[midcol] = rowoffset + midrow;

       /* upper left corner */
      (void) evaluatelinearcrosspoints(useq, ustart, midrow,
                                       vseq, vstart, midcol,
                                       EDtabcolumn,
                                       Rtabcolumn,
                                       Ctab,
                                       rowoffset,
                                       matchcost,
                                       mismatchcost,
                                       gapcost,
                                       original_ulen,
                                       original_vlen);

      /* bottom right corner */
      (void) evaluatelinearcrosspoints(useq, ustart + midrow,
                                       ulen-midrow,
                                       vseq, vstart + midcol,
                                       vlen-midcol,
                                       EDtabcolumn,
                                       Rtabcolumn,
                                       Ctab+midcol,
                                       rowoffset+midrow,
                                       matchcost,
                                       mismatchcost,
                                       gapcost,
                                       original_ulen,
                                       original_vlen);

      return distance;
    }
    else /* product of subsquences is in O(n) */
    {
      (void) ctab_in_square_space(Ctab, useq, ustart, ulen, vseq, vstart, vlen,
                                  matchcost, mismatchcost, gapcost,rowoffset);
    }
  }
  return 0;
}

static void determineCtab0(GtUword *Ctab, GtUchar vseq0,
                           const GtUchar *useq, GtUword ustart)
{
  GtUword rowindex;

  for (rowindex = 0; rowindex < Ctab[1]; rowindex++)
  {
    if (tolower((int)vseq0) == tolower((int)useq[ustart+rowindex]))
    {
      Ctab[0] = rowindex;
      return;
    }
  }

  Ctab[0] = (Ctab[1] > 0) ?  Ctab[1]-1 : 0;
}

GtUword gt_calc_linearalign(const GtUchar *useq,
                            GtUword ustart,
                            GtUword ulen,
                            const GtUchar *vseq,
                            GtUword vstart,
                            GtUword vlen,
                            GtAlignment *align,
                            GtUword matchcost,
                            GtUword mismatchcost,
                            GtUword gapcost)
{
  GtUword distance, *Ctab, *EDtabcolumn, *Rtabcolumn;

  if (ulen == 0UL)
  {
    return construct_trivial_insertion_alignment(align,vlen,gapcost);
  }
  else if (vlen == 0UL)
  {
    return construct_trivial_deletion_alignment(align,vlen,gapcost);
  }

  if (ulen == 1UL || vlen == 1UL ) {
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
    distance = evaluatelinearcrosspoints(useq, ustart, ulen,
                                         vseq, vstart, vlen,
                                         EDtabcolumn, Rtabcolumn,
                                         Ctab, 0, matchcost,
                                         mismatchcost, gapcost,
                                         ulen, vlen);

    determineCtab0(Ctab,vseq[vstart],useq, ustart);
    reconstructalignment_from_Ctab(align, Ctab, useq, ustart, vseq, vstart,
                                   vlen, matchcost, mismatchcost, 0,gapcost);

    gt_free(Ctab);
    gt_free(EDtabcolumn);
    gt_free(Rtabcolumn);
  }
  return distance;
}

GtUword gt_computelinearspace(GtAlignment *align,
                              const GtUchar *useq,
                              GtUword ustart,
                              GtUword ulen,
                              const GtUchar *vseq,
                              GtUword vstart,
                              GtUword vlen,
                              GtUword matchcost,
                              GtUword mismatchcost,
                              GtUword gapcost)
{
  gt_assert(useq != NULL  && ulen > 0 && vseq != NULL  && vlen > 0);

  gt_alignment_set_seqs(align, useq+ustart, ulen, vseq+vstart, vlen);
  return gt_calc_linearalign(useq, ustart, ulen, vseq, vstart, vlen, align,
                             matchcost, mismatchcost, gapcost);
}

void gt_print_edist_alignment(const GtUchar *useq, GtUword ustart, GtUword ulen,
                              const GtUchar *vseq, GtUword vstart, GtUword vlen)
{
  GtAlignment *align;

  align = gt_alignment_new_with_seqs(useq+ustart, ulen, vseq+vstart, vlen);
  (void)gt_calc_linearalign(useq,ustart,ulen,vseq,vstart,vlen,align,0,1,1);
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
       /* replacement */
      dpcolumn[i] = nw + (tolower((int)u[i-1]) == tolower((int)v[j-1]) ? 0 : 1);
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
static void firstLStabcolumn(GtUword ulen,
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
                            GtUword ustart, GtUword ulen,
                            const GtUchar b, GtUword colindex,
                            GtWord *Ltabcolumn,
                            GtUwordPair *Starttabcolumn,
                            Gtmaxcoordvalue *max,
                            GtWord matchscore,
                            GtWord mismatchscore,
                            GtWord gapscore)
{
  GtUword rowindex;
  GtUwordPair northwestStarttabentry, westStarttabentry;
  GtWord northwestLtabentry, westLtabentry;

  gt_assert(max != NULL);
  westLtabentry = Ltabcolumn[0];
  westStarttabentry = Starttabcolumn[0];

  Ltabcolumn[0] = 0;
  Starttabcolumn[0].a = 0;
  Starttabcolumn[0].b = colindex;
  for (rowindex = 1UL; rowindex <= ulen; rowindex++)
  {
    GtWord val;

    northwestLtabentry = westLtabentry;
    northwestStarttabentry = westStarttabentry;
    westLtabentry = Ltabcolumn[rowindex];
    westStarttabentry = Starttabcolumn[rowindex];
    Ltabcolumn[rowindex] += gapscore;

    if ((val = northwestLtabentry + (tolower((int)useq[ustart + rowindex-1]) ==
                                     tolower((int) b) ?
          matchscore : mismatchscore)) > Ltabcolumn[rowindex])
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
                                             GtUword ustart,
                                             GtUword ulen,
                                             const GtUchar *vseq,
                                             GtUword vstart,
                                             GtUword vlen,
                                             GtWord matchscore,
                                             GtWord mismatchscore,
                                             GtWord gapscore)
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

GtUword gt_computelinearspace_local(GtAlignment *align,
                                    const GtUchar *useq,
                                    GtUword ustart,
                                    GtUword ulen,
                                    const GtUchar *vseq,
                                    GtUword vstart,
                                    GtUword vlen,
                                    GtWord matchscore,
                                    GtWord mismatchscore,
                                    GtWord gapscore)
{
  GtWord *Ltabcolumn;
  GtUwordPair *Starttabcolumn;
  GtUword ulen_part, ustart_part, vlen_part, vstart_part, score;
  Gtmaxcoordvalue *max;
  GtUword matchcost, mismatchcost, gapcost;

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
                                &gapcost);
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
    score = gt_alignment_eval_generic_with_score(false, align, matchscore,
                                                 mismatchscore, gapscore);
  } else
  {
    gt_alignment_set_seqs(align,(const GtUchar*) "",0,(const GtUchar*) "",0);
    score = 0;
  }
  gt_max_delete(max);
  gt_assert(score == gt_alignment_eval_generic_with_score(false, align,
                                                          matchscore,
                                                          mismatchscore,
                                                          gapscore));
  return score;
}

/*-------------------------checkfunctions-----------------------------*/
void gt_checklinearspace(GT_UNUSED bool forward,
                         const GtUchar *useq,
                         GtUword ulen,
                         const GtUchar *vseq,
                         GtUword vlen)
{
  GtAlignment *align;
  GtUword i, edist1, edist2, edist3, edist4,
          matchcost = 0, mismatchcost = 1, gapcost = 1;
  /*immediate result, because squareedistunit cannot handle lower/upper cases*/
  GtUchar *low_useq = malloc(sizeof(*low_useq)*ulen),
          *low_vseq = malloc(sizeof(*low_vseq)*vlen);

  if (memchr(useq, LINEAR_EDIST_GAP,ulen) != NULL)
  {
    fprintf(stderr,"%s: sequence u contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  if (memchr(vseq, LINEAR_EDIST_GAP,vlen) != NULL)
  {
    fprintf(stderr,"%s: sequence v contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  for (i = 0; i < ulen; i++)
    low_useq[i] = tolower((int)useq[i]);
  low_useq[i] = '\0';
  for (i = 0; i < vlen; i++)
    low_vseq[i] = tolower((int)vseq[i]);
  low_vseq[i] = '\0';

  align = gt_alignment_new_with_seqs(low_useq, ulen, low_vseq, vlen);
  edist1 = gt_calc_linearalign(low_useq, 0, ulen, low_vseq, 0, vlen, align,
                               matchcost, mismatchcost, gapcost);

  edist2 = gt_squarededistunit(low_useq, ulen, low_vseq, vlen);

  if (edist1 != edist2)
  {
    fprintf(stderr,"gt_calc_linearalign = "GT_WU" != "GT_WU
            " = gt_squarededistunit\n", edist1,edist2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  edist3 = gt_alignment_eval_generic_with_score(false, align, matchcost,
                                                mismatchcost, gapcost);

  if (edist2 != edist3)
  {printf("useq %s, ulen:"GT_WU", vseq %s, vlen:"GT_WU"\n",useq,ulen,vseq,vlen);
    fprintf(stderr,"gt_squarededistunit = "GT_WU" != "GT_WU
            " = gt_alignment_eval_with_score\n", edist2,edist3);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  edist4 = gt_calc_linearedist(low_useq, ulen, low_vseq, vlen);
  if (edist3 != edist4)
  {
    fprintf(stderr,"gt_alignment_eval_with_score = "GT_WU" != "GT_WU
            " = gt_calc_linearedist\n", edist3, edist4);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_free(low_useq);
  gt_free(low_vseq);
  gt_alignment_delete(align);
}

void gt_checklinearspace_local(GT_UNUSED bool forward,
                               const GtUchar *useq, GtUword ulen,
                               const GtUchar *vseq, GtUword vlen)
{
  GtAlignment *align;
  GtWord i, score1, score2, matchscore = 2, mismatchscore = -2, gapscore = -1;
  GtUchar *low_useq = malloc(sizeof(*low_useq)*(ulen+1)),
          *low_vseq = malloc(sizeof(*low_vseq)*(vlen+1));

  if (memchr(useq, LINEAR_EDIST_GAP,ulen) != NULL)
  {
    fprintf(stderr,"%s: sequence u contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  if (memchr(vseq, LINEAR_EDIST_GAP,vlen) != NULL)
  {
    fprintf(stderr,"%s: sequence v contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  for (i = 0; i < ulen; i++)
    low_useq[i] = tolower((int)useq[i]);
  low_useq[i] = '\0';
  for (i = 0; i < vlen; i++)
    low_vseq[i] = tolower((int)vseq[i]);
  low_vseq[i] = '\0';

  align = gt_alignment_new();
  score1 = gt_computelinearspace_local(align,useq, 0, ulen, vseq, 0, vlen,
                                       matchscore, mismatchscore, gapscore);
  score2 = gt_alignment_eval_generic_with_score(false, align, matchscore,
                                                mismatchscore, gapscore);
  gt_alignment_delete(align);
  if (score1 != score2)
  {
    fprintf(stderr,"gt_computelinearspace_local = "GT_WU" != "GT_WU
            " = gt_alignment_eval_generic_with_score\n", score1, score2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_free(low_useq);
  gt_free(low_vseq);
}
