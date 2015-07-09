/*
  Copyright (C) 2015 Annika Seidel, annika.seidel@studium.uni-hamburg.de

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

#include "core/assert_api.h"
#include "core/unused_api.h"
#include "core/error.h"
#include "core/array.h"
#include "core/ma.h"
#include "core/minmax.h"

#include "extended/checksequence.h"
#include "extended/linearspace_local.h"
#include "extended/linearspace.h"
#include "extended/reconstructalignment.h"
#include "extended/alignment.h"
#include "extended/max.h"

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

static void nextLStabcolumn(const GtUchar *useq, const GtUword ulen,
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

    if ((val = northwestLtabentry + (useq[rowindex-1] ==
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
    if ((val = 0) > Ltabcolumn[rowindex])
    {
      Ltabcolumn[rowindex] = val;
      Starttabcolumn[rowindex].a = rowindex;
      Starttabcolumn[rowindex].b = colindex;
    }
    if (Ltabcolumn[rowindex] > gt_max_get_value(max))
    {
      gt_max_set_value(max, Ltabcolumn[rowindex]);
      gt_max_set_start(max, Starttabcolumn[rowindex]);
      gt_max_set_end(max, rowindex, colindex);
    }
  }
}

static Gtmaxcoordvalue *evaluateallLScolumns(const GtUchar *useq,
                                             const GtUword ulen,
                                             const GtUchar *vseq,
                                             const GtUword vlen,
                                             GtWord *Ltabcolumn,
                                             GtUwordPair *Starttabcolumn,
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
    nextLStabcolumn(useq, ulen, vseq[colindex-1], colindex,
                    Ltabcolumn, Starttabcolumn, max,
                    matchscore, mismatchscore, gapscore);
  }
  return max;
}

static GtWord gt_calc_linearscore(const GtUchar *useq, GtUword ulen,
                                  const GtUchar *vseq, GtUword vlen,
                                  const GtWord matchscore,
                                  const GtWord mismatchscore,
                                  const GtWord gapscore)
{
    Gtmaxcoordvalue *max;
    GtWord score, *Ltabcolumn;
    GtUwordPair *Starttabcolumn;

    Ltabcolumn = gt_malloc(sizeof *Ltabcolumn * (ulen+1));
    Starttabcolumn = gt_malloc(sizeof *Starttabcolumn * (ulen+1));

    max = evaluateallLScolumns(useq, ulen, vseq,vlen,
                               Ltabcolumn, Starttabcolumn,
                               matchscore, mismatchscore, gapscore);
    score = gt_max_get_value(max);
    gt_max_delete(max);

    return(score);
}

static void change_score_to_cost_function(const GtWord matchscore,
                                          const GtWord mismatchscore,
                                          const GtWord gapscore,
                                          GtWord *matchcost,
                                          GtWord *mismatchcost,
                                          GtWord *gapcost )
{
  GtWord temp;
//TODO:hier noch max und div nutzen und aufrunden
  temp=0;
  if (matchscore/2 > temp)
    temp=matchscore/2;
  else if (mismatchscore/2 > temp)
    temp=mismatchscore/2;
  else if (gapscore > temp)
    temp=gapscore;

  *matchcost = 2*temp-matchscore;
  *mismatchcost = 2*temp-mismatchscore;
  *gapcost = temp-gapscore;

}

static GtAlignment *gt_calc_linearalign_local(const GtUchar *useq,GtUword ulen,
                                              const GtUchar *vseq,GtUword vlen,
                                              const GtWord matchscore,
                                              const GtWord mismatchscore,
                                              const GtWord gapscore)
{
  GtWord *Ltabcolumn;
  GtUwordPair *Starttabcolumn;
  GtUword ulen_part, vlen_part;
  const GtUchar *useq_part, *vseq_part;
  GtAlignment *align;
  Gtmaxcoordvalue *max;
  GtWord matchcost, mismatchcost, gapcost;

  /*gt_assert(ulen > 0 && vlen > 0);*/

  Ltabcolumn = gt_malloc(sizeof *Ltabcolumn * (ulen+1));
  Starttabcolumn = gt_malloc(sizeof *Starttabcolumn * (ulen+1));

  max = evaluateallLScolumns(useq, ulen,
                             vseq, vlen,
                             Ltabcolumn,
                             Starttabcolumn,
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
    useq_part = &useq[(gt_max_get_start(max)).a];
    vseq_part = &vseq[(gt_max_get_start(max)).b];
    ulen_part = gt_max_get_row_length(max);
    vlen_part = gt_max_get_col_length(max);
    align = gt_alignment_new_with_seqs(useq_part,ulen_part,vseq_part,vlen_part);
    gt_calc_linearalign_with_costs(useq_part, ulen_part, vseq_part, vlen_part,
                                   align, matchcost, mismatchcost, gapcost);
  }else
  {
    align = gt_alignment_new_with_seqs(useq,ulen,vseq,vlen);
    gt_calc_linearalign_with_costs(useq, ulen, vseq, vlen,
                                   align, matchcost, mismatchcost, gapcost);
  }

  gt_max_delete(max);

  return align;
}

static GtWord fillLtable(GtWord *lcolumn,
                         const GtUchar *u, GtUword ulen,
                         const GtUchar *v, GtUword vlen,
                         const GtWord matchscore,
                         const GtWord mismatchscore,
                         const GtWord gapscore)
{
  GtUword i, j;
  GtWord nw, we, max = 0;
  for (i = 0; i <= ulen; i++)
    lcolumn[i] = 0;
  for (j = 1UL; j <= vlen; j++)
  {
    nw = lcolumn[0];
    lcolumn[0] = 0;
    for (i = 1UL; i <= ulen; i++)
    {
      we = lcolumn[i];
      lcolumn[i] = nw + (u[i-1] == v[j-1] ?
                        matchscore : mismatchscore); /* replacement */
      if (lcolumn[i-1] + gapscore > lcolumn[i]) /* deletion */
        lcolumn[i] = lcolumn[i-1] + gapscore;
      if (we + gapscore >lcolumn[i]) /* insertion */
        lcolumn[i] = we + gapscore;
      if (0 > lcolumn[i])
        lcolumn[i]=0;
      nw = we;
      if (lcolumn[i] > max)
        max = lcolumn[i];
    }
  }
  return max;
}

/* use this function to calculate score if no char matches*/
static GtUword gt_calc_linearscore_safe(const GtUword ulen,
                                        const GtUword vlen,
                                        const GtWord mismatchscore,
                                        const GtWord gapscore)
{
  GtUword idx, length, score=0;

  length = MAX(ulen,vlen);
  score += mismatchscore;
  for (idx = 1; idx < length ; idx++)
  {
    score += gapscore;
  }
  return score;
}

static GtUword gt_calc_linearscore_with_table(const GtUchar *useq, GtUword ulen,
                                              const GtUchar *vseq, GtUword vlen,
                                              const GtWord matchscore,
                                              const GtWord mismatchscore,
                                              const GtWord gapscore)
{
  GtWord *lcolumn, score;

  lcolumn = gt_malloc(sizeof *lcolumn * (MIN(ulen,vlen) + 1));
  score = fillLtable(lcolumn, ulen <= vlen ? useq : vseq, MIN(ulen,vlen),
                     ulen <= vlen ? vseq : useq, MAX(ulen,vlen),
                     matchscore, mismatchscore, gapscore);

  gt_free(lcolumn);
  return score;
}

void gt_computelinearspace_local(bool showevalue,
                                 const GtUchar *useq, GtUword ulen,
                                 const GtUchar *vseq, GtUword vlen,
                                 const GtWord matchscore,
                                 const GtWord mismatchscore,
                                 const GtWord gapscore,
                                 FILE *fp)
{
  GtAlignment *align;
  GtUword score;

  /*gt_assert(useq && ulen && vseq && vlen);*/
  align = gt_calc_linearalign_local(useq, ulen, vseq, vlen,
                                    matchscore, mismatchscore, gapscore);
 /* fprintf(fp, ">> local alignment: u = %s, v = %s\n",useq,vseq);*/
  gt_alignment_show(align, fp, 80);
  if(showevalue)
  {
    score = gt_alignment_eval_with_score(align, matchscore,
                                         mismatchscore, gapscore);
    fprintf(fp, "linear score: "GT_WD"\n", score);
  }
  gt_alignment_delete(align);
}

void gt_checklinearspace_local(GT_UNUSED bool forward,
                               const GtUchar *useq, GtUword ulen,
                               const GtUchar *vseq, GtUword vlen)
{
  /*printf("useq %s, vseq %s, ulen:"GT_WU", vlen:"GT_WU"\n",
   * useq,vseq,ulen,vlen);*/
  GtWord score1, score2, score3, score4;
  GtAlignment *align;
  bool no_align = false;

  gt_assert(useq && ulen && vseq && vlen);

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

  score1 = gt_calc_linearscore(useq, ulen, vseq, vlen, 2, -2, -1);
  score2 = gt_calc_linearscore_with_table(useq, ulen, vseq, vlen, 2, -2, -1);

  if (score1 != score2)
  {
    fprintf(stderr,"gt_calc_linearscore = "GT_WD" != "GT_WD
            " = gt_calc_linearscore_with_table\n", score1, score2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  if (score1 == 0)
     no_align =true;

  align = gt_calc_linearalign_local(useq, ulen, vseq, vlen, 2, -2, -1);
  score3 = gt_alignment_eval_with_score(align, 2,-2,-1);

  if (no_align)
  {
    score4 = gt_calc_linearscore_safe(ulen,vlen, -2, -1);
    if (score3 != score4)
    {
      fprintf(stderr,"gt_calc_linearalign_local = "GT_WD" != "GT_WD
              " = gt_calc_linearscore_safe\n", score3, score4);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
  else{
    if (score1 != score3)
    {
      fprintf(stderr,"gt_calc_linearscore = "GT_WD" != "GT_WD
              " = gt_calc_linearalign_local\n", score1, score3);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }

  gt_alignment_delete(align);
}
