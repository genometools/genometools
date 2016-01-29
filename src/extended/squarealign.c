/*
  Copyright (c) 2015 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
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
#include <ctype.h>
#include "core/array2dim_api.h"
#include "extended/maxcoordvalue.h"
#include "core/minmax.h"
#include "extended/reconstructalignment.h"

#include "extended/squarealign.h"

/*----------------------------global alignment--------------------------------*/
static void fillDPtab_in_square_space(GtUword **E,
                                      const GtUchar *useq,
                                      GtUword ustart,
                                      GtUword ulen,
                                      const GtUchar *vseq,
                                      GtUword vstart,
                                      GtUword vlen,
                                      const GtScoreHandler *scorehandler)
{
  GtUword i, j, gapcost;
  gapcost =  gt_scorehandler_get_gapscore(scorehandler);

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
      val = E[i-1][j-1] + gt_scorehandler_get_replacement(scorehandler,
                                                          useq[ustart+i-1],
                                                          vseq[vstart+j-1]);
      if (val <= E[i][j])
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

/* create an global alignment in square space, to use it in linear context you
 * have to generate an spacemanager before, in any other case it can be NULL */
GtUword gt_squarealign_calculate_generic (GtLinspaceManagement *spacemanager,
                                          GtAlignment *align,
                                          const GtUchar *useq,
                                          GtUword ustart,
                                          GtUword ulen,
                                          const GtUchar *vseq,
                                          GtUword vstart,
                                          GtUword vlen,
                                          const GtScoreHandler *scorehandler)
{
  GtUword **E, distance;
  gt_assert(align && scorehandler);

  if (spacemanager == NULL)
  {
    /*use it in normally case*/
    gt_array2dim_malloc(E, (ulen+1), (vlen+1));
  }
  else
  {
    /*use it in lineraspace context*/
    E = gt_linspace_management_change_to_square(spacemanager,ulen,vlen);
  }
  fillDPtab_in_square_space(E, useq, ustart, ulen,
                            vseq, vstart, vlen, scorehandler);

  distance = E[ulen][vlen];
  /* reconstruct alignment from 2dimarray E */
  gt_reconstructalignment_from_EDtab(align, E, useq, ustart, ulen, vseq, vstart,
                                     vlen, scorehandler);
  if (spacemanager == NULL)
  {
    gt_array2dim_delete(E);
  }
  return distance;
}

GtUword gt_squarealign_global_distance_only(const GtUchar *useq,
                                            GtUword ustart,
                                            GtUword ulen,
                                            const GtUchar *vseq,
                                            GtUword vstart,
                                            GtUword vlen,
                                            const GtScoreHandler *scorehandler)
{
  GtUword **E, distance;

  gt_assert(scorehandler);

  gt_array2dim_malloc(E, (ulen+1), (vlen+1));
  fillDPtab_in_square_space(E, useq, ustart, ulen,
                            vseq, vstart, vlen, scorehandler);
  distance = E[ulen][vlen];
  gt_array2dim_delete(E);
  return distance;
}

/* create a global alignment in square space with constant cost values,
 * to use it in linear context you have to generate an spacemanager before,
 * in any other case it can be NULL */
GtUword gt_squarealign_calculate(GtLinspaceManagement *spacemanager,
                                 GtAlignment *align,
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
  GtUword distance;
  GtScoreHandler *scorehandler;

  gt_assert(align);
  scorehandler = gt_scorehandler_new(matchcost, mismatchcost, 0, gapcost);
  distance = gt_squarealign_calculate_generic (spacemanager, align,
                                               useq, ustart,  ulen,
                                               vseq,  vstart, vlen,
                                               scorehandler);
  gt_scorehandler_delete(scorehandler);
  return distance;
}

void gt_squarealign_print_edit_alignment(const GtUchar *useq, GtUword ustart,
                                         GtUword ulen, const GtUchar *vseq,
                                         GtUword vstart, GtUword vlen)
{
  GtAlignment *align;
  align = gt_alignment_new_with_seqs(useq+ustart, ulen, vseq+vstart, vlen);

  gt_squarealign_calculate(NULL, align, useq, ustart, ulen,
                           vseq, vstart, vlen, 0, 1, 1);
  gt_alignment_show(align, true, stdout, 80);
  gt_alignment_delete(align);
}

static void evaluate_crosspoints_from_2dimtab(
                                       GtUword **E,
                                       GtUword *Ctab,
                                       const GtScoreHandler *scorehandler,
                                       const GtUchar *useq,
                                       GtUword ustart,
                                       GtUword ulen,
                                       const GtUchar *vseq,
                                       GtUword vstart,
                                       GtUword vlen,
                                       GtUword rowoffset)

{
  GtUword gapcost, idx, jdx;

  gt_assert(E && Ctab);
  gapcost = gt_scorehandler_get_gapscore(scorehandler);
  idx = ulen;
  jdx = vlen;
  while (jdx > 1 || idx > 0)
  {
    if (idx > 0 && jdx > 0 && E[idx][jdx] == E[idx-1][jdx-1] +
                                   gt_scorehandler_get_replacement(scorehandler,
                                                             useq[ustart+idx-1],
                                                            vseq[vstart+jdx-1]))
    {
      idx--;
      jdx--;
      Ctab[jdx] = idx + rowoffset;
    }
    else if (idx > 0 && E[idx][jdx] == E[idx-1][jdx] + gapcost)
    {
      idx--;
      /*continue*/
    }
    else if (jdx > 0 && E[idx][jdx] == E[idx][jdx-1] + gapcost)
    {
      jdx--;
      Ctab[jdx] = idx + rowoffset;
    }
    else
      gt_assert(false);
  }
}

/* fill crosspointtable ctab for part of sequences useq and vseq in square
 * space, use it to combine square calculating with linear calculating */
GtUword gt_squarealign_ctab(GtLinspaceManagement *spacemanager,
                            const GtScoreHandler *scorehandler,
                            GtUword *Ctab,
                            const GtUchar *useq,
                            GtUword ustart,
                            GtUword ulen,
                            const GtUchar *vseq,
                            GtUword vstart,
                            GtUword vlen,
                            GtUword rowoffset)
{
  GtUword **E, distance;

  gt_assert(Ctab && spacemanager && scorehandler);

  E = gt_linspace_management_change_to_square(spacemanager,ulen,vlen);
  E = gt_linspace_management_change_to_square(spacemanager,ulen, vlen);

  fillDPtab_in_square_space(E, useq, ustart, ulen,
                            vseq, vstart, vlen, scorehandler);
  distance = E[ulen][vlen];

  evaluate_crosspoints_from_2dimtab(E, Ctab, scorehandler, useq, ustart, ulen,
                                    vseq, vstart, vlen,  rowoffset);
  return distance;
}

/*----------------------------local alignment---------------------------------*/
static GtWord fillDPtab_in_square_space_local(GtWord **Ltabcolumn,
                                              GtMaxcoordvalue *max,
                                              const GtUchar *useq,
                                              GtUword ustart,
                                              GtUword ulen,
                                              const GtUchar *vseq,
                                              GtUword vstart,
                                              GtUword vlen,
                                             const GtScoreHandler *scorehandler)
{
  GtUword i, j;
  GtWord gapscore, maxscore, overall_maxscore = 0;

  gapscore = gt_scorehandler_get_gapscore(scorehandler);
  Ltabcolumn[0][0] = 0;
  for (i = 1; i <= ulen; i++)
  {
      Ltabcolumn[i][0] = 0;
  }
  for (j = 1; j <= vlen; j++)
  {
      Ltabcolumn[0][j] = 0;
      for (i = 1; i <= ulen; i++)
      {
        GtWord val;

        Ltabcolumn[i][j] = Ltabcolumn[i][j-1] + gapscore;

        val = Ltabcolumn[i-1][j-1] +
                                  gt_scorehandler_get_replacement(scorehandler,
                                           useq[ustart+i-1], vseq[vstart+j-1]);

        if (val >= Ltabcolumn[i][j])
        {
          Ltabcolumn[i][j] = val;
        }
        if ((val = Ltabcolumn[i-1][j] + gapscore) > Ltabcolumn[i][j])
        {
          Ltabcolumn[i][j] = val;
        }
        Ltabcolumn[i][j] = MAX(Ltabcolumn[i][j],0);
        maxscore = Ltabcolumn[i][j];

        if (maxscore > (GtWord) overall_maxscore)
        {
          overall_maxscore = maxscore;
          gt_maxcoordvalue_coord_update_without_start (max, maxscore, i, j);
        }
     }
  }

  return overall_maxscore;
}

/* create an local alignment in square space, to use it in linear context you
 * have to generate an spacemanager before, in any other case it can be NULL */
GtWord gt_squarealign_calculate_local_generic(GtLinspaceManagement
                                              *spacemanager,
                                              GtAlignment *align,
                                              const GtUchar *useq,
                                              GtUword ustart,
                                              GtUword ulen,
                                              const GtUchar *vseq,
                                              GtUword vstart,
                                              GtUword vlen,
                                              const GtScoreHandler
                                              *scorehandler)
{
  GtWord score = 0, **Ltabcolumn;
  GtMaxcoordvalue *max;

  gt_assert(align != NULL);
  if (spacemanager == NULL)
  {
    /*use it in normally case*/
    gt_array2dim_malloc(Ltabcolumn, (ulen+1), (vlen+1));
    max = gt_maxcoordvalue_new();
  }
  else
  {
    /*use it in lineraspace context*/
    Ltabcolumn = (GtWord **)
                 gt_linspace_management_change_to_square(spacemanager,
                                                         ulen, vlen);
    max = gt_linspace_management_get_maxspace(spacemanager);
  }

  score = fillDPtab_in_square_space_local(Ltabcolumn, max, useq, ustart, ulen,
                                          vseq, vstart, vlen, scorehandler);

  /* reconstruct local alignment from 2dimarray Ltabcolumn */
  gt_reconstructalignment_from_Ltab(align, Ltabcolumn, max,
                                    useq, ustart, ulen,
                                    vseq, vstart, vlen,
                                    scorehandler);

  if (gt_maxcoordvalue_get_length_safe(max))
  {
    ustart = ustart + (gt_maxcoordvalue_get_start(max)).a;
    vstart = vstart + (gt_maxcoordvalue_get_start(max)).b;
    ulen = gt_maxcoordvalue_get_row_length(max);
    vlen = gt_maxcoordvalue_get_col_length(max);

    gt_alignment_set_seqs(align, &useq[ustart], ulen,
                                 &vseq[vstart], vlen);
  }

  if (spacemanager == NULL)
  {
    gt_array2dim_delete(Ltabcolumn);
    gt_maxcoordvalue_delete(max);
  }
  return score;
}

/* create an local alignment in square space with constant score values,
 * to use it in linear context you have to generate an spacemanager before,
 * in any other case it can be NULL */
GtWord gt_squarealign_calculate_local(GtLinspaceManagement *spacemanager,
                                      GtAlignment *align,
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
  GtWord score;
  gt_assert(align);
  GtScoreHandler *scorehandler = gt_scorehandler_new(matchscore,
                                                     mismatchscore, 0,
                                                     gapscore);
  score = gt_squarealign_calculate_local_generic(spacemanager, align,
                                                 useq, ustart, ulen,
                                                 vseq, vstart, vlen,
                                                 scorehandler);
  gt_scorehandler_delete(scorehandler);
  return score;
}
