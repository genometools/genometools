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
#include "core/assert_api.h"
#include "core/divmodmul.h"
#include "core/ma_api.h"
#include "core/minmax.h"

#include "extended/alignment.h"
#include "extended/affinealign_linear.h"
#include "extended/affinealign_linear_local.h"
#include "extended/alignment.h"
#include "extended/maxcoordvalue.h"

#define LINEAR_EDIST_GAP          ((GtUchar) UCHAR_MAX)
typedef enum {
  R,
  D,
  I,
  X /*unknown*/
} Edge;

typedef struct {
  GtWord Rvalue, Dvalue, Ivalue, totalvalue;
}Atabentry;

typedef struct {
  GtUwordPair Rstart, Dstart, Istart;
}Starttabentry;

static GtUword add_safe(const GtWord val1, const GtWord val2)
{
  if (val1 != GT_WORD_MIN && val2 != GT_WORD_MIN)
     return val1+val2;

  return GT_WORD_MIN;
}

static void firstAStabcolumn(const GtUword ulen,
                             Atabentry *Atabcolumn,
                             Starttabentry *Starttabcolumn,
                             const GtWord gap_opening,
                             const GtWord gap_extension)
{
  GtUword rowindex;

  Atabcolumn[0].Rvalue = GT_WORD_MIN;
  Atabcolumn[0].Dvalue = GT_WORD_MIN;
  Atabcolumn[0].Ivalue = GT_WORD_MIN;
  Atabcolumn[0].totalvalue = 0;

  Starttabcolumn[0].Rstart.a = 0;
  Starttabcolumn[0].Rstart.b = 0;
  Starttabcolumn[0].Dstart.a = 0;
  Starttabcolumn[0].Dstart.b = 0;
  Starttabcolumn[0].Istart.a = 0;
  Starttabcolumn[0].Istart.b = 0;

  for (rowindex = 1; rowindex <= ulen; rowindex++)
  {
    Atabcolumn[rowindex].Rvalue = GT_WORD_MIN;
    Atabcolumn[rowindex].Dvalue = (gap_opening + gap_extension);
    Atabcolumn[rowindex].Ivalue = GT_WORD_MIN;
    Atabcolumn[rowindex].totalvalue = 0;

    Starttabcolumn[rowindex].Rstart.a = rowindex;
    Starttabcolumn[rowindex].Rstart.b = 0;
    Starttabcolumn[rowindex].Dstart.a = rowindex;
    Starttabcolumn[rowindex].Dstart.b = 0;
    Starttabcolumn[rowindex].Istart.a = rowindex;
    Starttabcolumn[rowindex].Istart.b = 0;
  }
}

static GtUwordPair setStarttabentry(GtWord entry, Atabentry aTab,
                                    Starttabentry sTab,
                                    const GtWord replacement,
                                    const GtWord gap_opening,
                                    const GtWord gap_extension,
                                    const Edge edge)
{
  GtUwordPair start;
  switch (edge) {
  case R:
    if (entry == aTab.Rvalue + replacement)
       start = sTab.Rstart;
    else if (entry == aTab.Dvalue + replacement)
       start = sTab.Dstart;
    else if (entry == aTab.Ivalue + replacement)
       start = sTab.Istart;
    else
      start = sTab.Rstart;
    break;
  case D:
    if (entry == aTab.Rvalue + gap_opening + gap_extension)
       start = sTab.Rstart;
    else if (entry == aTab.Dvalue + gap_extension)
       start = sTab.Dstart;
    else if (entry == aTab.Ivalue + gap_opening + gap_extension)
       start = sTab.Istart;
    else
      start = sTab.Rstart;
    break;
  case I:
    if (entry == aTab.Rvalue + gap_opening + gap_extension)
       start = sTab.Rstart;
    else if (entry == aTab.Dvalue + gap_opening + gap_extension)
       start = sTab.Dstart;
    else if (entry == aTab.Ivalue + gap_extension)
       start = sTab.Istart;
    else
      start = sTab.Rstart;
    break;
  default:
    start.a = 0;
    start.b = 0;
  }
  return start;
}

static void nextAStabcolumn(const GtUchar *useq, GtUword ustart,
                            const GtUword ulen,
                            const GtUchar b,
                            Atabentry *Atabcolumn,
                            Starttabentry *Starttabcolumn,
                            const GtWord matchscore,
                            const GtWord mismatchscore,
                            const GtWord gap_opening,
                            const GtWord gap_extension,
                            const GtUword colindex,
                            Gtmaxcoordvalue *max)
{
  Atabentry Anw, Awe;
  Starttabentry Snw, Swe;
  GtUword rowindex;
  GtWord replacement, temp, val1, val2;
  GtUwordPair start;

  Anw = Atabcolumn[0];
  Snw = Starttabcolumn[0];
  Atabcolumn[0].Rvalue = GT_WORD_MIN;
  Atabcolumn[0].Dvalue = GT_WORD_MIN;
  Atabcolumn[0].Ivalue = (gap_opening + gap_extension);
  temp = MAX3(Atabcolumn[0].Rvalue,
              Atabcolumn[0].Dvalue,
              Atabcolumn[0].Ivalue);
  Atabcolumn[0].totalvalue = ((temp > 0)? temp : 0);

  if (Atabcolumn[0].totalvalue == 0) {
    Starttabcolumn[0].Rstart.a = 0;
    Starttabcolumn[0].Rstart.b = colindex;
    Starttabcolumn[0].Dstart.a = 0;
    Starttabcolumn[0].Dstart.b = colindex;
    Starttabcolumn[0].Istart.a = 0;
    Starttabcolumn[0].Istart.b = colindex;
  }

  if (Atabcolumn[0].totalvalue > gt_max_get_value(max))
    {
      if (Atabcolumn[0].totalvalue == Atabcolumn[0].Rvalue)
         start = Starttabcolumn[0].Rstart;
      else if (Atabcolumn[0].totalvalue == Atabcolumn[0].Dvalue)
         start = Starttabcolumn[0].Dstart;
      else if (Atabcolumn[0].totalvalue == Atabcolumn[0].Ivalue)
         start = Starttabcolumn[0].Istart;

      gt_max_coord_update(max, Atabcolumn[0].totalvalue,
                          start, 0, colindex);
    }
  for (rowindex = 1; rowindex <= ulen; rowindex++)
  {
    Awe = Atabcolumn[rowindex];
    Swe = Starttabcolumn[rowindex];

    /*calculate Rvalue*/
    replacement = (useq[ustart+rowindex-1]==b? matchscore:mismatchscore);
    Atabcolumn[rowindex].Rvalue = add_safe(Anw.totalvalue,replacement);
    Starttabcolumn[rowindex].Rstart =
    setStarttabentry(Atabcolumn[rowindex].Rvalue, Anw, Snw,
                     replacement,gap_opening,gap_extension,R);

    /*calculate Dvalue*/
    val1 = add_safe(Atabcolumn[rowindex-1].Dvalue,gap_extension);
    val2 = add_safe(Atabcolumn[rowindex-1].totalvalue,
                   (gap_opening+gap_extension));
    Atabcolumn[rowindex].Dvalue = MAX(val1,val2);
    Starttabcolumn[rowindex].Dstart =
    setStarttabentry(Atabcolumn[rowindex].Dvalue, Atabcolumn[rowindex-1],
                     Starttabcolumn[rowindex-1], replacement,gap_opening,
                     gap_extension,D);

    /*calculate Ivalue*/
    val1=(add_safe(Awe.Ivalue,gap_extension));
    val2=(add_safe(Awe.totalvalue,(gap_opening+gap_extension)));
    Atabcolumn[rowindex].Ivalue = MAX(val1,val2);
    Starttabcolumn[rowindex].Istart =
    setStarttabentry(Atabcolumn[rowindex].Ivalue, Awe, Swe, replacement,
                     gap_opening, gap_extension,I);

    /*calculate totalvalue*/
    temp = MAX3(Atabcolumn[rowindex].Rvalue,
                Atabcolumn[rowindex].Dvalue,
                Atabcolumn[rowindex].Ivalue);
    Atabcolumn[rowindex].totalvalue = ((temp > 0)? temp : 0);

    /* set start indices for Atab-values*/
    if (Atabcolumn[rowindex].totalvalue == 0)
    {
      Starttabcolumn[rowindex].Rstart.a = rowindex;
      Starttabcolumn[rowindex].Rstart.b = colindex;
      Starttabcolumn[rowindex].Dstart.a = rowindex;
      Starttabcolumn[rowindex].Dstart.b = colindex;
      Starttabcolumn[rowindex].Istart.a = rowindex;
      Starttabcolumn[rowindex].Istart.b = colindex;
    }

    /*set new max*/
    if (Atabcolumn[rowindex].totalvalue > gt_max_get_value(max))
    {
      if (Atabcolumn[rowindex].totalvalue == Atabcolumn[rowindex].Rvalue)
         start = Starttabcolumn[rowindex].Rstart;
      else if (Atabcolumn[rowindex].totalvalue == Atabcolumn[rowindex].Dvalue)
         start = Starttabcolumn[rowindex].Dstart;
      else if (Atabcolumn[rowindex].totalvalue == Atabcolumn[rowindex].Ivalue)
         start = Starttabcolumn[rowindex].Istart;

      gt_max_coord_update(max, Atabcolumn[rowindex].totalvalue,
                          start, rowindex, colindex);
    }
    Anw=Awe;
    Snw=Swe;
  }
}

static Gtmaxcoordvalue *evaluateallAStabcolumns(const GtUchar *useq,
                                                const GtUword ustart,
                                                const GtUword ulen,
                                                const GtUchar *vseq,
                                                const GtUword vstart,
                                                const GtUword vlen,
                                                Atabentry *Atabcolumn,
                                                Starttabentry *Starttabcolumn,
                                                const GtWord matchscore,
                                                const GtWord mismatchscore,
                                                const GtWord gap_opening,
                                                const GtWord gap_extension)
{
  GtUword colindex;
  Gtmaxcoordvalue *max;
  firstAStabcolumn(ulen, Atabcolumn, Starttabcolumn,
                      gap_opening, gap_extension);

  max = gt_max_new();
  for (colindex = 1UL; colindex <= vlen; colindex++)
  {
    nextAStabcolumn(useq, ustart, ulen,
                       vseq[vstart+colindex-1],
                       Atabcolumn,
                       Starttabcolumn,
                       matchscore,mismatchscore,
                       gap_opening,
                       gap_extension,
                       colindex,max);
  }
  return max;
}

static void change_score_to_cost_affine_function(const GtWord matchscore,
                                                 const GtWord mismatchscore,
                                                 const GtWord gap_opening,
                                                 const GtWord gap_extension,
                                                 GtWord *replacement_cost,
                                                 GtWord *gap_opening_cost,
                                                 GtWord *gap_extension_cost)
{
  GtWord temp1, temp2, max, matchcost;

  temp1 = MAX(GT_DIV2(matchscore), GT_DIV2(mismatchscore));

  temp2 = MAX(1 + gap_opening, 1 + gap_extension);

  max = MAX(temp1, temp2);
  if (max < 0)
    max = 0;

  matchcost =  2 * max-matchscore; /*set matchcost to zero*/
  *replacement_cost = 2 * max-mismatchscore-matchcost;
  *gap_opening_cost = max-gap_opening;
  *gap_extension_cost = max-gap_extension-matchcost;
}

static GtUword gt_calc_affinealign_linear_local(const GtUchar *useq,
                                                const GtUword ustart,
                                                const GtUword ulen,
                                                const GtUchar *vseq,
                                                const GtUword vstart,
                                                const GtUword vlen,
                                                GtAlignment *align,
                                                const GtWord matchscore,
                                                const GtWord mismatchscore,
                                                const GtWord gap_opening,
                                                const GtWord gap_extension)
{
  GtUword score, ulen_part, ustart_part, vlen_part, vstart_part;
  GtWord replacement_cost,
         gap_opening_cost,
         gap_extension_cost;
  Atabentry *Atabcolumn;
  Starttabentry *Starttabcolumn;
  Gtmaxcoordvalue *max;

  Atabcolumn = gt_malloc(sizeof *Atabcolumn * (ulen+1));
  Starttabcolumn = gt_malloc(sizeof *Starttabcolumn * (ulen+1));

  max = evaluateallAStabcolumns(useq, ustart, ulen, vseq, vstart, vlen,
                                Atabcolumn, Starttabcolumn,
                                matchscore, mismatchscore,
                                gap_opening, gap_extension);

  score = gt_max_get_value(max);

  if (gt_max_get_length_safe(max))
  {

    ustart_part = ustart+(gt_max_get_start(max)).a;
    vstart_part = vstart+(gt_max_get_start(max)).b;
    ulen_part = gt_max_get_row_length(max);
    vlen_part = gt_max_get_col_length(max);

    gt_alignment_set_seqs(align,&useq[ustart_part],ulen_part,
                                &vseq[vstart_part],vlen_part);

    change_score_to_cost_affine_function(matchscore,mismatchscore,
                                         gap_opening,gap_extension,
                                         &replacement_cost,
                                         &gap_opening_cost,
                                         &gap_extension_cost);

    gt_calc_affinealign_linear(useq, ustart_part, ulen_part,
                               vseq, vstart_part, vlen_part,
                               align, replacement_cost,
                               gap_opening_cost,gap_extension_cost);
  }else
  {
     gt_alignment_set_seqs(align,(const GtUchar*)"",0,
                                        (const GtUchar*)"",0);
     score = 0;
  }

  gt_max_delete(max);
  gt_free(Atabcolumn);
  gt_free(Starttabcolumn);

  return(score);
}

GtAlignment *gt_computeaffinelinearspace_local(const GtUchar *useq,
                                       const GtUword ustart,
                                       const GtUword ulen,
                                       const GtUchar *vseq,
                                       const GtUword vstart,
                                       const GtUword vlen,
                                       const GtWord matchscore,
                                       const GtWord mismatchscore,
                                       const GtWord gap_opening,
                                       const GtWord gap_extension)
{
  GtAlignment *align;

  align = gt_alignment_new();
  (void) gt_calc_affinealign_linear_local(useq, ustart, ulen,
                                           vseq, vstart, vlen,
                                           align, matchscore,mismatchscore,
                                           gap_opening, gap_extension);
  return align;
}

void gt_checkaffinelinearspace_local(GT_UNUSED bool forward,
                               const GtUchar *useq,
                               GtUword ulen,
                               const GtUchar *vseq,
                               GtUword vlen)
{
  GtAlignment *align;
  GtUword affine_score1, affine_score2;

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
  affine_score1 = gt_calc_affinealign_linear_local(useq,0, ulen,
                                                vseq,0,vlen,align,
                                                6,-3,-2,-1);
  affine_score2 = gt_alignment_eval_with_affine_score(align, 6,-3,-2,-1);

  if (affine_score1 != affine_score2)
  {
    fprintf(stderr,"gt_calc_affinealign_linear_local = "GT_WU" != "GT_WU
            " = gt_alignment_eval_with_affine_score\n", affine_score1,
                                                        affine_score2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  gt_alignment_delete(align);
}
