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
#include "core/array.h"
#include "core/assert_api.h"
#include "core/divmodmul.h"
#include "core/error.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/unused_api.h"

#include "extended/alignment.h"
#include "extended/maxcoordvalue.h"
#include "extended/linearspace_local.h"
#include "extended/linearspace.h"
#include "extended/reconstructalignment.h"

#define LINEAR_EDIST_GAP          ((GtUchar) UCHAR_MAX)

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

static Gtmaxcoordvalue *evaluateallLScolumns(const GtUchar *useq,
                                             const GtUword ustart,
                                             const GtUword ulen,
                                             const GtUchar *vseq,
                                             const GtUword vstart,
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
    nextLStabcolumn(useq, ustart, ulen, vseq[vstart+colindex-1], colindex,
                    Ltabcolumn, Starttabcolumn, max,
                    matchscore, mismatchscore, gapscore);
  }
  return max;
}

static GtWord gt_calc_linearscore(const GtUchar *useq,
                                  const GtUword ustart, const GtUword ulen,
                                  const GtUchar *vseq,
                                  const GtUword vstart, const GtUword vlen,
                                  const GtWord matchscore,
                                  const GtWord mismatchscore,
                                  const GtWord gapscore)
{
    Gtmaxcoordvalue *max;
    GtWord score, *Ltabcolumn;
    GtUwordPair *Starttabcolumn;

    Ltabcolumn = gt_malloc(sizeof *Ltabcolumn * (ulen + 1));
    Starttabcolumn = gt_malloc(sizeof *Starttabcolumn * (ulen + 1));

    max = evaluateallLScolumns(useq, ustart, ulen, vseq, vstart, vlen,
                               Ltabcolumn, Starttabcolumn,
                               matchscore, mismatchscore, gapscore);
    score = gt_max_get_value(max);
    gt_max_delete(max);
    gt_free(Ltabcolumn);
    gt_free(Starttabcolumn);

    return(score);
}

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

  /*gt_assert(ulen > 0 && vlen > 0);*/

  Ltabcolumn = gt_malloc(sizeof *Ltabcolumn * (ulen+1));
  Starttabcolumn = gt_malloc(sizeof *Starttabcolumn * (ulen+1));

  max = evaluateallLScolumns(useq, ustart, ulen,
                             vseq, vstart, vlen,
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
    ustart_part = ustart+(gt_max_get_start(max)).a;
    vstart_part = vstart+(gt_max_get_start(max)).b;
    ulen_part = gt_max_get_row_length(max);
    vlen_part = gt_max_get_col_length(max);

    gt_alignment_set_seqs(align, &useq[ustart_part], ulen_part,
                                 &vseq[vstart_part], vlen_part);
    gt_calc_linearalign2(useq, ustart_part, ulen_part,
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

  /*gt_assert(useq && ulen && vseq && vlen);*/
  align = gt_alignment_new();
  (void) gt_calc_linearalign_local(useq, ustart, ulen, vseq, vstart,vlen, align,
                                    matchscore, mismatchscore, gapscore);

  return align;
}

static GtWord fillLtable(GtWord *lcolumn,
                         const GtUchar *u, GtUword ustart, GtUword ulen,
                         const GtUchar *v, GtUword vstart, GtUword vlen,
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
      lcolumn[i] = nw + (u[ustart+i-1] == v[vstart+j-1] ?
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
static GtUword gt_calc_linearscore_with_table(const GtUchar *useq,
                                              const GtUword ustart,
                                              const GtUword ulen,
                                              const GtUchar *vseq,
                                              const GtUword vstart,
                                              const GtUword vlen,
                                              const GtWord matchscore,
                                              const GtWord mismatchscore,
                                              const GtWord gapscore)
{
  GtWord *lcolumn, score;

  lcolumn = gt_malloc(sizeof *lcolumn * (MIN(ulen,vlen) + 1));
  score = fillLtable(lcolumn, ulen <= vlen ? useq : vseq,
                     ulen <= vlen ? ustart:vstart,MIN(ulen,vlen),
                     ulen <= vlen ? vseq : useq,
                     ulen <= vlen ? vstart : ustart,MAX(ulen,vlen),
                     matchscore, mismatchscore, gapscore);

  gt_free(lcolumn);
  return score;
}

void gt_checklinearspace_local(GT_UNUSED bool forward,
                               const GtUchar *useq, GtUword ulen,
                               const GtUchar *vseq, GtUword vlen)
{
  GtWord score1, score2, score3;
  GtAlignment *align;

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

  score1 = gt_calc_linearscore(useq, 0, ulen, vseq, 0, vlen, 2, -2, -1);
  score2 = gt_calc_linearscore_with_table(useq, 0, ulen,
                                          vseq, 0, vlen, 2, -2, -1);

  if (score1 != score2)
  {
    fprintf(stderr,"gt_calc_linearscore = "GT_WD" != "GT_WD
            " = gt_calc_linearscore_with_table\n", score1, score2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  align = gt_alignment_new();
  score3 = gt_calc_linearalign_local(useq, 0, ulen,
                                     vseq, 0, vlen, align, 2, -2, -1);

    if (score1 != score3)
    {
      fprintf(stderr,"gt_calc_linearscore = "GT_WD" != "GT_WD
              " = gt_calc_linearalign_local\n", score1, score3);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }

  gt_alignment_delete(align);
}
