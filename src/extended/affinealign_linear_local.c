#include "core/assert_api.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#include "extended/affinealign_linear.h"
#include "extended/affinealign_linear_local.h"
#include "extended/alignment.h"
#include "extended/checksequence.h"
#include "extended/max.h"

typedef struct {
  GtWord Rvalue, Dvalue, Ivalue, totalvalue;
}Atabentry;

static GtUword add_safe(const GtWord val1, const GtWord val2)
{
  if (val1 != GT_WORD_MIN && val2 != GT_WORD_MIN)
     return val1+val2;

  return GT_WORD_MIN;
}

static void firstAStabcolumn(const GtUword ulen,
                             Atabentry *Atabcolumn,
                             GtUwordPair *Starttabcolumn,
                             const GtWord gap_opening,
                             const GtWord gap_extension)
{
  GtUword rowindex;

  Atabcolumn[0].Rvalue = GT_WORD_MIN;
  Atabcolumn[0].Dvalue = GT_WORD_MIN;
  Atabcolumn[0].Ivalue = GT_WORD_MIN;
  Atabcolumn[0].totalvalue = 0;

  Starttabcolumn[0].a = 0;
  Starttabcolumn[0].b = 0;

  for (rowindex = 1; rowindex <= ulen; rowindex++)
  {
    Atabcolumn[rowindex].Rvalue = GT_WORD_MIN;
    Atabcolumn[rowindex].Dvalue = (gap_opening + gap_extension);
    Atabcolumn[rowindex].Ivalue = GT_WORD_MIN;
    Atabcolumn[rowindex].totalvalue = 0;

    Starttabcolumn[rowindex].a = rowindex;
    Starttabcolumn[rowindex].b = 0;
  }
}

static void nextAStabcolumn(const GtUchar *useq,
                               const GtUword ulen,
                               const GtUchar b,
                               Atabentry *Atabcolumn,
                               GtUwordPair *Starttabcolumn,
                               const GtWord replacement_score,
                               const GtWord gap_opening,
                               const GtWord gap_extension,
                               const GtUword colindex,
                               Gtmaxcoordvalue *max)
{
  Atabentry Anw, Awe;
  GtUwordPair Snw, Swe;
  GtUword rowindex;
  GtWord temp, val1, val2;

  Anw = Atabcolumn[0];
  Snw = Starttabcolumn[0];
  Atabcolumn[0].Rvalue = GT_WORD_MIN;
  Atabcolumn[0].Dvalue = GT_WORD_MIN;
  Atabcolumn[0].Ivalue = (gap_opening + gap_extension);
  temp = MAX3(Atabcolumn[0].Rvalue,
              Atabcolumn[0].Dvalue,
              Atabcolumn[0].Ivalue);
  Atabcolumn[0].totalvalue = ((temp > 0)? temp : 0);
  Starttabcolumn[0].a = 0;
  Starttabcolumn[0].b = colindex;

  if (Atabcolumn[0].totalvalue > gt_max_get_value(max))
    {
      gt_max_set_value(max, Atabcolumn[0].totalvalue);
      gt_max_set_start(max, Starttabcolumn[0]);
      gt_max_set_end(max, 0, colindex);
    }
  for (rowindex = 1; rowindex <= ulen; rowindex++)
  {
    Awe = Atabcolumn[rowindex];
    Swe = Starttabcolumn[rowindex];

    /*calculate Rvalue*/
    Atabcolumn[rowindex].Rvalue = add_safe(Anw.totalvalue,
                             (useq[rowindex-1]==b? 3:replacement_score));
    /*calculate Dvalue*/
    val1 = add_safe(Atabcolumn[rowindex-1].Dvalue,gap_extension);
    val2 = add_safe(Atabcolumn[rowindex-1].totalvalue, 
                   (gap_opening+gap_extension));
    Atabcolumn[rowindex].Dvalue = MAX(val1,val2);
    
    /*calculate Ivalue*/
    val1=(add_safe(Awe.Ivalue,gap_extension));
    val2=(add_safe(Awe.totalvalue,(gap_opening+gap_extension)));
    Atabcolumn[rowindex].Ivalue = MAX(val1,val2);
    
    /*calculate totalvalue*/
    temp = MAX3(Atabcolumn[rowindex].Rvalue,
                Atabcolumn[rowindex].Dvalue,
                Atabcolumn[rowindex].Ivalue);
    Atabcolumn[rowindex].totalvalue = ((temp > 0)? temp : 0);

    /* set start indices for Atab-values*/
    if (Atabcolumn[rowindex].totalvalue == Atabcolumn[rowindex].Rvalue)
    {
      Starttabcolumn[rowindex].a = Snw.a;
      Starttabcolumn[rowindex].b = Snw.b;
    }
    else if (Atabcolumn[rowindex].totalvalue == Atabcolumn[rowindex].Dvalue)
    {
      Starttabcolumn[rowindex].a = Starttabcolumn[rowindex-1].a;
      Starttabcolumn[rowindex].b = Starttabcolumn[rowindex-1].b;
    }
    else if (Atabcolumn[rowindex].totalvalue == Atabcolumn[rowindex].Ivalue)
    {//ignore
      Starttabcolumn[rowindex].a = Swe.a;
      Starttabcolumn[rowindex].b = Swe.b;
    }
    else
    {
      Starttabcolumn[rowindex].a = rowindex;
      Starttabcolumn[rowindex].b = colindex;
    }
    /*set new max*/
    if (Atabcolumn[rowindex].totalvalue > gt_max_get_value(max))
    {
      gt_max_set_value(max, Atabcolumn[rowindex].totalvalue);
      gt_max_set_start(max, Starttabcolumn[rowindex]);
      gt_max_set_end(max, rowindex, colindex);
    }
    Anw=Awe;
    Snw=Swe;
  }
}

static Gtmaxcoordvalue *evaluateallAStabcolumns(const GtUchar *useq,
                                                GtUword ulen,
                                                const GtUchar *vseq,
                                                GtUword vlen,
                                                Atabentry *Atabcolumn,
                                                GtUwordPair *Starttabcolumn,
                                                const GtWord replacement_score,
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
    nextAStabcolumn(useq, ulen,
                       vseq[colindex-1],
                       Atabcolumn,
                       Starttabcolumn,
                       replacement_score,
                       gap_opening,
                       gap_extension,
                       colindex,max);
  }
  return max;
}

static GtUword gt_calc_affinealign_linear_local(const GtUchar *useq,
                                          const GtUword ulen,
                                          const GtUchar *vseq,
                                          const GtUword vlen,
                                          GtAlignment *align,
                                          const GtWord replacement_score,
                                          const GtWord gap_opening,
                                          const GtWord gap_extension)
{
  GtUword score;
  GtUword ulen_part, vlen_part;
  const GtUchar *useq_part, *vseq_part;
  Atabentry *Atabcolumn;
  GtUwordPair *Starttabcolumn;
  Gtmaxcoordvalue *max;

  Atabcolumn = gt_malloc(sizeof *Atabcolumn * (ulen+1));
  Starttabcolumn = gt_malloc(sizeof *Starttabcolumn * (ulen+1));

  max = evaluateallAStabcolumns(useq, ulen, vseq,vlen,
                                Atabcolumn, Starttabcolumn,
                                replacement_score, gap_opening,
                                gap_extension);
  score = gt_max_get_value(max);
  
  //change_affine_score_to_cost_function TODO

  if (gt_max_get_length_safe(max))
  {
    useq_part = &useq[(gt_max_get_start(max)).a];
    vseq_part = &vseq[(gt_max_get_start(max)).b];
    ulen_part = gt_max_get_row_length(max);
    vlen_part = gt_max_get_col_length(max);
    gt_alignment_set_seqs(align,useq_part,ulen_part,vseq_part,vlen_part);
    score = gt_calc_affinealign_linear(useq_part, ulen_part,
                                       vseq_part, vlen_part,
                                       align, 0,1,3);
  }else
  {
    gt_alignment_set_seqs(align,useq,ulen,vseq,vlen);
    score = gt_calc_affinealign_linear(useq, ulen, vseq, vlen,
                                       align, 6,1,3);
  }
  gt_max_delete(max);
  gt_free(Atabcolumn);
  gt_free(Starttabcolumn);
  return(score);
}

void gt_computeaffinelinearspace_local(bool showevalue,
                                 const GtUchar *useq, GtUword ulen,
                                 const GtUchar *vseq, GtUword vlen,
                                 const GtWord replacement_score,
                                 const GtWord gap_opening,
                                 const GtWord gap_extension,
                                 FILE *fp)
{
  GtAlignment *align;
  GtUword score;

  /*gt_assert(useq && ulen && vseq && vlen);*/
  align = gt_alignment_new();
  score = gt_calc_affinealign_linear_local(useq, ulen, vseq, vlen,
                                           align, replacement_score,
                                           gap_opening, gap_extension);
  gt_alignment_show(align, fp, 80);
  if (showevalue)
    fprintf(fp, "local affine score: "GT_WU"\n", score);
  gt_alignment_delete(align);
}
