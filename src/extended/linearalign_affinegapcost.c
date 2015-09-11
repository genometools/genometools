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

#include <ctype.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/minmax.h"
#include "core/error.h"
#include "core/types_api.h"
#include "core/divmodmul.h"
#include "core/ma_api.h"
#include "extended/affinealign.h"
#include "extended/linearalign_affinegapcost.h"
#include "extended/linearalign_utilities.h"
#include "extended/maxcoordvalue.h"
#include "extended/reconstructalignment.h"

#define LINEAR_EDIST_GAP          ((GtUchar) UCHAR_MAX)

typedef struct {
  GtUwordPair Rstart, Dstart, Istart;
} Starttabentry;

static void change_score_to_cost_affine_function(GtWord matchscore,
                                                 GtWord mismatchscore,
                                                 GtWord gap_opening,
                                                 GtWord gap_extension,
                                                 GtUword *match_cost,
                                                 GtUword *mismatch_cost,
                                                 GtUword *gap_opening_cost,
                                                 GtUword *gap_extension_cost)
{
  GtWord temp1, temp2, max;

  temp1 = MAX(GT_DIV2(matchscore), GT_DIV2(mismatchscore));
  temp2 = MAX(0, 1 + gap_extension);
  max = MAX(temp1, temp2);
  *match_cost = 2 * max-matchscore;
  *mismatch_cost = 2 * max-mismatchscore;
  *gap_opening_cost = -gap_opening;
  *gap_extension_cost = max-gap_extension;
}

/*------------------------------global--------------------------------*/
inline AffineAlignEdge set_edge(GtWord Rdist,
                                GtWord Ddist,
                                GtWord Idist)
{
  GtUword minvalue;
  minvalue = MIN3(Rdist, Ddist, Idist);

  if (Rdist == minvalue)
    return Affine_R;
  else if (Ddist == minvalue)
    return Affine_D;
  else if (Idist == minvalue)
    return Affine_I;

  return Affine_X;
}

static Rnode get_Rtabentry(const Rtabentry *rtab, const AffineAlignEdge edge)
{
  switch (edge) {
  case Affine_R:
    return rtab->val_R;
  case Affine_D:
    return rtab->val_D;
  case Affine_I:
    return rtab->val_I;
  default:
    gt_assert(false);
  }
}

static inline void firstAtabRtabentry(AffinealignDPentry *Atabcolumn,
                               GtUword gap_opening,
                               AffineAlignEdge edge)
{
  switch (edge) {
  case Affine_R:
    Atabcolumn[0].Rvalue = 0;
    Atabcolumn[0].Dvalue = GT_WORD_MAX;
    Atabcolumn[0].Ivalue = GT_WORD_MAX;
    break;
  case Affine_D:
    Atabcolumn[0].Rvalue = GT_WORD_MAX;
    Atabcolumn[0].Dvalue = 0;
    Atabcolumn[0].Ivalue = GT_WORD_MAX;
    break;
  case Affine_I:
    Atabcolumn[0].Rvalue = GT_WORD_MAX;
    Atabcolumn[0].Dvalue = GT_WORD_MAX;
    Atabcolumn[0].Ivalue = 0;
    break;
  default:
    Atabcolumn[0].Rvalue = 0;
    Atabcolumn[0].Dvalue = gap_opening;
    Atabcolumn[0].Ivalue = gap_opening;
  }
}

static void firstAtabRtabcolumn(GtUword ulen,
                                AffinealignDPentry *Atabcolumn,
                                Rtabentry *Rtabcolumn,
                                GtUword gap_opening,
                                GtUword gap_extension,
                                AffineAlignEdge edge)
{
  GtUword rowindex;
  firstAtabRtabentry(Atabcolumn, gap_opening, edge);

  Atabcolumn[0].Redge = Affine_X;
  Atabcolumn[0].Dedge = Affine_X;
  Atabcolumn[0].Iedge = Affine_X;

  Rtabcolumn[0].val_R.idx = 0;
  Rtabcolumn[0].val_D.idx = 0;
  Rtabcolumn[0].val_I.idx = 0;

  Rtabcolumn[0].val_R.edge = Affine_R;
  Rtabcolumn[0].val_D.edge = Affine_D;
  Rtabcolumn[0].val_I.edge = Affine_I;

  for (rowindex = 1; rowindex <= ulen; rowindex++)
  {
    Atabcolumn[rowindex].Rvalue = GT_WORD_MAX;
    Atabcolumn[rowindex].Dvalue = add_safe_max(Atabcolumn[rowindex-1].Dvalue,
                                               gap_extension);
    Atabcolumn[rowindex].Ivalue = GT_WORD_MAX;

    Atabcolumn[rowindex].Redge = Affine_X;
    Atabcolumn[rowindex].Dedge = Affine_D;
    Atabcolumn[rowindex].Iedge = Affine_X;

    Rtabcolumn[rowindex].val_R.idx = rowindex;
    Rtabcolumn[rowindex].val_D.idx = rowindex;
    Rtabcolumn[rowindex].val_I.idx = rowindex;

    Rtabcolumn[rowindex].val_R.edge = Affine_R;
    Rtabcolumn[rowindex].val_D.edge = Affine_D;
    Rtabcolumn[rowindex].val_I.edge = Affine_I;
  }
}

static void nextAtabRtabcolumn(const GtUchar *useq,
                               GtUword ustart,
                               GtUword ulen,
                               const GtUchar b,
                               AffinealignDPentry *Atabcolumn,
                               Rtabentry *Rtabcolumn,
                               GtUword matchcost,
                               GtUword mismatchcost,
                               GtUword gap_opening,
                               GtUword gap_extension,
                               GtUword midcolumn,
                               GtUword colindex)
{
  AffinealignDPentry northwestAffinealignDPentry, westAffinealignDPentry;
  Rtabentry northwestRtabentry, westRtabentry;
  GtWord rowindex, rcost, rdist, ddist, idist, minvalue;

  northwestAffinealignDPentry = Atabcolumn[0];
  northwestRtabentry = Rtabcolumn[0];

  rdist = add_safe_max(Atabcolumn[0].Rvalue, gap_extension + gap_opening);
  ddist = add_safe_max(Atabcolumn[0].Dvalue, gap_extension + gap_opening);
  idist = add_safe_max(Atabcolumn[0].Ivalue, gap_extension);

  minvalue = MIN3(rdist, ddist, idist);
  Atabcolumn[0].Ivalue = minvalue;
  Atabcolumn[0].Rvalue = GT_WORD_MAX;
  Atabcolumn[0].Dvalue = GT_WORD_MAX;

  Atabcolumn[0].Redge = Affine_X;
  Atabcolumn[0].Dedge = Affine_X;
  Atabcolumn[0].Iedge = Affine_I;

  if (colindex > midcolumn)
  {
    northwestRtabentry = Rtabcolumn[0];
    Rtabcolumn[0].val_R.idx = Rtabcolumn[0].val_I.idx;
    Rtabcolumn[0].val_D.idx = Rtabcolumn[0].val_I.idx;
    Rtabcolumn[0].val_I.idx = Rtabcolumn[0].val_I.idx;

    Rtabcolumn[0].val_R.edge = Affine_X;
    Rtabcolumn[0].val_D.edge = Affine_X;
    Rtabcolumn[0].val_I.edge = Rtabcolumn[0].val_I.edge;
  }

  for (rowindex = 1; rowindex <= ulen; rowindex++)
  {
    westAffinealignDPentry = Atabcolumn[rowindex];
    westRtabentry = Rtabcolumn[rowindex];

    rcost = tolower((int)useq[ustart+rowindex-1]) == tolower((int)b) ?
            matchcost:mismatchcost;
    rdist = add_safe_max(northwestAffinealignDPentry.Rvalue, rcost);
    ddist = add_safe_max(northwestAffinealignDPentry.Dvalue, rcost);
    idist = add_safe_max(northwestAffinealignDPentry.Ivalue, rcost);

    minvalue = MIN3(rdist, ddist, idist);
    Atabcolumn[rowindex].Rvalue = minvalue;
    Atabcolumn[rowindex].Redge = set_edge(rdist, ddist, idist);

    rdist = add_safe_max(Atabcolumn[rowindex-1].Rvalue,
                         gap_extension + gap_opening);
    ddist = add_safe_max(Atabcolumn[rowindex-1].Dvalue,gap_extension);
    idist = add_safe_max(Atabcolumn[rowindex-1].Ivalue,
                         gap_extension + gap_opening);

    minvalue = MIN3(rdist, ddist, idist);
    Atabcolumn[rowindex].Dvalue = minvalue;
    Atabcolumn[rowindex].Dedge = set_edge(rdist, ddist, idist);

    rdist = add_safe_max(westAffinealignDPentry.Rvalue,
                         gap_extension + gap_opening);
    ddist = add_safe_max(westAffinealignDPentry.Dvalue,
                         gap_extension + gap_opening);
    idist = add_safe_max(westAffinealignDPentry.Ivalue, gap_extension);

    minvalue = MIN3(rdist, ddist, idist);
    Atabcolumn[rowindex].Ivalue = minvalue;
    Atabcolumn[rowindex].Iedge = set_edge(rdist, ddist, idist);

    if (colindex > midcolumn)
    {
      Rtabcolumn[rowindex].val_R = get_Rtabentry(&northwestRtabentry,
                                                  Atabcolumn[rowindex].Redge);
      Rtabcolumn[rowindex].val_D = get_Rtabentry(&Rtabcolumn[rowindex-1],
                                                  Atabcolumn[rowindex].Dedge);
      Rtabcolumn[rowindex].val_I = get_Rtabentry(&westRtabentry,
                                                  Atabcolumn[rowindex].Iedge);
    }
    northwestAffinealignDPentry = westAffinealignDPentry;
    northwestRtabentry = westRtabentry;
  }
}

static GtUword evaluateallAtabRtabcolumns(const GtUchar *useq,
                                          GtUword ustart,
                                          GtUword ulen,
                                          const GtUchar *vseq,
                                          GtUword vstart,
                                          GtUword vlen,
                                          AffinealignDPentry *Atabcolumn,
                                          Rtabentry *Rtabcolumn,
                                          GtUword matchcost,
                                          GtUword mismatchcost,
                                          GtUword gap_opening,
                                          GtUword gap_extension,
                                          GtUword midcolumn,
                                          AffineAlignEdge edge)
{
  GtUword colindex;
  firstAtabRtabcolumn(ulen, Atabcolumn, Rtabcolumn,
                      gap_opening, gap_extension, edge);

  for (colindex = 1UL; colindex <= vlen; colindex++)
  {
    nextAtabRtabcolumn(useq, ustart,ulen,
                       vseq[vstart+colindex-1],
                       Atabcolumn,
                       Rtabcolumn,
                       matchcost,
                       mismatchcost,
                       gap_opening,
                       gap_extension,
                       midcolumn,
                       colindex);
  }

  return MIN3(Atabcolumn[ulen].Rvalue,
              Atabcolumn[ulen].Dvalue,
              Atabcolumn[ulen].Ivalue);
}

AffineAlignEdge minAdditionalCosts(const AffinealignDPentry *entry,
                                   const AffineAlignEdge edge,
                                   GtUword gap_opening)
{
  GtUword rdist, ddist, idist;

  switch (edge) {
    case Affine_D:
      rdist = add_safe_max(entry->Rvalue, gap_opening);
      ddist = entry->Dvalue;
      idist = add_safe_max(entry->Ivalue, gap_opening);
     break;
    case Affine_I:
      rdist = add_safe_max(entry->Rvalue, gap_opening);
      ddist = add_safe_max(entry->Dvalue, gap_opening);
      idist = entry->Ivalue;
      break;
    default:
      rdist = entry->Rvalue;
      ddist = entry->Dvalue;
      idist = entry->Ivalue;
  }

  return set_edge(rdist, ddist, idist);
}

static GtUword evaluateaffinecrosspoints(const GtUchar *useq,
                                         GtUword ustart,
                                         GtUword ulen,
                                         GtUword original_ulen,
                                         const GtUchar *vseq,
                                         GtUword vstart,
                                         GtUword vlen,
                                         GtUword original_vlen,
                                         AffinealignDPentry *Atabcolumn,
                                         Rtabentry *Rtabcolumn,
                                         GtUword *Ctab,
                                         GtUword rowoffset,
                                         GtUword matchcost,
                                         GtUword mismatchcost,
                                         GtUword gap_opening,
                                         GtUword gap_extension,
                                         AffineAlignEdge from_edge,
                                         AffineAlignEdge to_edge)
{
  if (vlen >= 2UL)
  {
    if ((ulen+1)*(vlen+1)>(original_ulen+1)) {
    GtUword  midrow = 0, midcol = GT_DIV2(vlen), distance, colindex;
    AffineAlignEdge bottomtype, midtype = Affine_X;

    distance = evaluateallAtabRtabcolumns(useq, ustart, ulen,
                                          vseq, vstart, vlen,
                                          Atabcolumn, Rtabcolumn,
                                          matchcost, mismatchcost,
                                          gap_opening,
                                          gap_extension,
                                          midcol, from_edge);

    bottomtype = minAdditionalCosts(&Atabcolumn[ulen], to_edge, gap_opening);
    switch (bottomtype) {
      case Affine_R:
        midrow = Rtabcolumn[ulen].val_R.idx;
        midtype = Rtabcolumn[ulen].val_R.edge;
        break;
      case Affine_D:
        midrow = Rtabcolumn[ulen].val_D.idx;
        midtype = Rtabcolumn[ulen].val_D.edge;
        break;
      case Affine_I:
        midrow = Rtabcolumn[ulen].val_I.idx;
        midtype = Rtabcolumn[ulen].val_I.edge;
        break;
      case Affine_X: /*never reach this line*/
        gt_assert(false);
    }
    Ctab[midcol] = rowoffset + midrow;
    gt_assert(midcol > 0);
    if (midrow == 0) {
      for (colindex = midcol-1; colindex > 0; colindex--)
        Ctab[colindex] = Ctab[midcol];
    }
    else{/* upper left corner */
      switch (midtype) {
        case Affine_R:
          if (midcol > 1)
            Ctab[midcol-1] = Ctab[midcol] == 0 ? 0: Ctab[midcol] - 1;

            (void) evaluateaffinecrosspoints(useq, ustart, midrow-1,
                                             original_ulen,
                                             vseq, vstart, midcol-1,
                                             original_vlen,
                                             Atabcolumn,Rtabcolumn,
                                             Ctab,rowoffset,
                                             matchcost, mismatchcost,
                                             gap_opening,
                                             gap_extension,
                                             from_edge,midtype);
          break;
        case Affine_D:
          (void) evaluateaffinecrosspoints(useq,ustart,midrow-1,
                                           original_ulen,
                                           vseq,vstart,midcol,
                                           original_vlen,
                                           Atabcolumn,Rtabcolumn,
                                           Ctab,rowoffset,
                                           matchcost, mismatchcost,
                                           gap_opening,
                                           gap_extension,
                                           from_edge,midtype);
          break;
        case Affine_I:
          if (midcol>1)
            Ctab[midcol-1] = (Ctab[midcol]);
          (void) evaluateaffinecrosspoints(useq,ustart,midrow,
                                           original_ulen,
                                           vseq,vstart,midcol-1,
                                           original_vlen,
                                           Atabcolumn,Rtabcolumn,
                                           Ctab,rowoffset,
                                           matchcost, mismatchcost,
                                           gap_opening,
                                           gap_extension,
                                           from_edge,midtype);
          break;
        case Affine_X: /*never reach this line*/
                gt_assert(false);
      }
    }
   /*bottom right corner */
   evaluateaffinecrosspoints(useq, ustart+midrow, ulen-midrow, original_ulen,
                             vseq, vstart+midcol,vlen-midcol, original_vlen,
                             Atabcolumn,Rtabcolumn,
                             Ctab+midcol,rowoffset+midrow,
                             matchcost, mismatchcost,
                             gap_opening,
                             gap_extension,
                             midtype, to_edge);
    return distance;}
    else /* product of subsquences is in O(n) */
    {
      affine_ctab_in_square_space(Ctab, useq, ustart, ulen, vseq, vstart, vlen,
                                  matchcost, mismatchcost, gap_opening,
                                  gap_extension, rowoffset, from_edge, to_edge);
    }
  }
  return 0;
}

static void affine_determineCtab0(GtUword *Ctab, GtUchar vseq0,
                                  const GtUchar *useq,
                                  GtUword ustart,
                                  const GtWord matchcost,
                                  const GtWord mismatchcost,
                                  const GtWord gap_opening)
{
  GtUword rowindex;

  if (Ctab[1] == 1 || Ctab[1] == 0)
  {
    Ctab[0] = 0;
    return;
  }
  else
  {
    if (Ctab[2]-Ctab[1] > 1)
    {
      if (gap_opening > (mismatchcost-matchcost))
      {
        Ctab[0] = 0;
        return;
      }
      else
      {
        if (tolower((int)vseq0) == tolower((int)useq[ustart]))
        {
          Ctab[0] = 0;
          return;
        }
        for (rowindex = Ctab[1]-1; rowindex >= 0; rowindex--)
        {
          if (tolower((int)vseq0) == tolower((int)useq[ustart+rowindex]))
          {
            Ctab[0] = rowindex;
            return;
          }
        }
        Ctab[0] = 0;
        return;
      }
    }
    else
    {
      if (tolower((int)vseq0) == tolower((int)useq[ustart+Ctab[1]-1]))
      {
          Ctab[0] = Ctab[1]-1;
          return;
      }
      else if (tolower((int)vseq0) == tolower((int)useq[ustart]))
      {
          Ctab[0] = 0;
          return;
      }
      if (gap_opening > (mismatchcost-matchcost))
      {
        Ctab[0] = Ctab[1]-1;
        return;
      }
      else
      {
        for (rowindex = 0; rowindex < Ctab[1]; rowindex++)
        {
          if (tolower((int)vseq0) == tolower((int)useq[ustart+rowindex]))
          {
            Ctab[0] = rowindex;
            return;
          }
        }
         Ctab[0] = Ctab[1]-1;
         return;
      }
    }
  }

  Ctab[0] = (Ctab[1] > 0) ?  Ctab[1]-1 : 0;

}

GtUword gt_calc_affinealign_linear(const GtUchar *useq, GtUword ustart,
                                   GtUword ulen,
                                   const GtUchar *vseq, GtUword vstart,
                                   GtUword vlen,
                                   GtAlignment *align,
                                   GtUword matchcost,
                                   GtUword mismatchcost,
                                   GtUword gap_opening,
                                   GtUword gap_extension)
{
  GtUword distance, *Ctab;
  AffinealignDPentry *Atabcolumn;
  Rtabentry *Rtabcolumn;
  GtAlignment *square_align;
  GtUchar *low_useq, *low_vseq;

  if (ulen == 0UL)
  {
      distance = construct_trivial_insertion_alignment(align, vlen,
                                                      gap_extension);
      distance += gap_opening;
      return distance;
  }
  else if (vlen == 0UL)
  {
      distance = construct_trivial_deletion_alignment(align, ulen,
                                                      gap_extension);
      distance += gap_opening;
      return distance;
  }

  low_useq = sequence_to_lower_case(useq, ulen);
  low_vseq = sequence_to_lower_case(vseq, vlen);
  if (ulen == 1UL || vlen == 1UL )
  {
    square_align = gt_affinealign(useq+ustart, ulen,
                                  vseq+vstart, vlen,
                                  (int)matchcost, (int)mismatchcost,
                                  (int)gap_opening,
                                  (int)gap_extension);
    gt_alignment_clone(square_align, align);

    distance = gt_alignment_eval_generic_with_affine_score(false, align,
                                                           matchcost,
                                                           mismatchcost,
                                                           gap_opening,
                                                           gap_extension);
    gt_alignment_delete(square_align);
  }
  else
  {
    Ctab = gt_malloc(sizeof *Ctab * (vlen+1));
    Atabcolumn = gt_malloc(sizeof *Atabcolumn * (ulen+1));
    Rtabcolumn = gt_malloc(sizeof *Rtabcolumn * (ulen+1));

    Ctab[vlen] = ulen;
    distance = evaluateaffinecrosspoints(useq, ustart, ulen, ulen,
                                         vseq, vstart, vlen, vlen,
                                         Atabcolumn, Rtabcolumn,
                                         Ctab, 0, matchcost, mismatchcost,
                                         gap_opening,gap_extension,
                                         Affine_X,Affine_X);

    affine_determineCtab0(Ctab, vseq[vstart],useq, ustart,
                          matchcost, mismatchcost, gap_opening);

    reconstructalignment_from_Ctab(align,Ctab,useq,ustart,vseq,
                                   vstart,vlen,matchcost,mismatchcost,
                                   gap_opening,gap_extension);

    gt_free(Ctab);
    gt_free(Atabcolumn);
    gt_free(Rtabcolumn);
  }
  gt_free(low_useq);
  gt_free(low_vseq);
  return distance;
}

void gt_computeaffinelinearspace(GtAlignment *align,
                                 const GtUchar *useq,
                                 GtUword ustart,
                                 GtUword ulen,
                                 const GtUchar *vseq,
                                 GtUword vstart,
                                 GtUword vlen,
                                 GtUword matchcost,
                                 GtUword mismatchcost,
                                 GtUword gap_opening,
                                 GtUword gap_extension)
{
  gt_assert(useq != NULL  && ulen > 0 && vseq != NULL  && vlen > 0);

  gt_alignment_set_seqs(align,useq+ustart, ulen, vseq+vstart, vlen);
  (void) gt_calc_affinealign_linear(useq, ustart, ulen,
                                vseq, vstart, vlen,
                                align, matchcost, mismatchcost,
                                gap_opening,gap_extension);

}

/*------------------------------local--------------------------------*/
static void firstAStabcolumn(GtUword ulen,
                             AffinealignDPentry *Atabcolumn,
                             Starttabentry *Starttabcolumn,
                             GtWord gap_opening,
                             GtWord gap_extension)
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

static GtUwordPair setStarttabentry(GtWord entry, AffinealignDPentry *Atab,
                                    Starttabentry *Stab,
                                    GtWord replacement,
                                    GtWord gap_opening,
                                    GtWord gap_extension,
                                    const AffineAlignEdge edge)
{
  GtUwordPair start;
  switch (edge) {
    case Affine_R:
      if (entry == Atab->Rvalue + replacement)
         start = Stab->Rstart;
      else if (entry == Atab->Dvalue + replacement)
         start = Stab->Dstart;
      else if (entry == Atab->Ivalue + replacement)
         start = Stab->Istart;
      else
        start = Stab->Rstart;
      break;
    case Affine_D:
      if (entry == Atab->Rvalue + gap_opening + gap_extension)
         start = Stab->Rstart;
      else if (entry == Atab->Dvalue + gap_extension)
         start = Stab->Dstart;
      else if (entry == Atab->Ivalue + gap_opening + gap_extension)
         start = Stab->Istart;
      else
        start = Stab->Rstart;
      break;
    case Affine_I:
      if (entry == Atab->Rvalue + gap_opening + gap_extension)
         start = Stab->Rstart;
      else if (entry == Atab->Dvalue + gap_opening + gap_extension)
         start = Stab->Dstart;
      else if (entry == Atab->Ivalue + gap_extension)
         start = Stab->Istart;
      else
        start = Stab->Rstart;
      break;
    default:
      start.a = 0;
      start.b = 0;
  }
  return start;
}

static void nextAStabcolumn(const GtUchar *useq, GtUword ustart,
                            GtUword ulen,
                            const GtUchar b,
                            AffinealignDPentry *Atabcolumn,
                            Starttabentry *Starttabcolumn,
                            GtWord matchscore,
                            GtWord mismatchscore,
                            GtWord gap_opening,
                            GtWord gap_extension,
                            GtUword colindex,
                            Gtmaxcoordvalue *max)
{
  AffinealignDPentry northwestAffinealignDPentry, westAffinealignDPentry;
  Starttabentry Snw, Swe;
  GtUword rowindex;
  GtWord replacement, temp, val1, val2;
  GtUwordPair start;

  northwestAffinealignDPentry = Atabcolumn[0];
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
    westAffinealignDPentry = Atabcolumn[rowindex];
    Swe = Starttabcolumn[rowindex];

    /*calculate Rvalue*/
    replacement = (tolower((int)useq[ustart+rowindex-1]) == tolower((int)b) ?
                   matchscore : mismatchscore);
    Atabcolumn[rowindex].Rvalue =
              add_safe_min(northwestAffinealignDPentry.totalvalue, replacement);
    Starttabcolumn[rowindex].Rstart =
    setStarttabentry(Atabcolumn[rowindex].Rvalue, &northwestAffinealignDPentry,
                      &Snw, replacement,gap_opening, gap_extension, Affine_R);

    /*calculate Dvalue*/
    val1 = add_safe_min(Atabcolumn[rowindex-1].Dvalue, gap_extension);
    val2 = add_safe_min(Atabcolumn[rowindex-1].totalvalue,
                       (gap_opening+gap_extension));
    Atabcolumn[rowindex].Dvalue = MAX(val1,val2);
    Starttabcolumn[rowindex].Dstart =
    setStarttabentry(Atabcolumn[rowindex].Dvalue, &Atabcolumn[rowindex-1],
                     &Starttabcolumn[rowindex-1], replacement, gap_opening,
                     gap_extension,Affine_D);

    /*calculate Ivalue*/
    val1=(add_safe_min(westAffinealignDPentry.Ivalue,gap_extension));
    val2=(add_safe_min(westAffinealignDPentry.totalvalue,
                              gap_opening+gap_extension));
    Atabcolumn[rowindex].Ivalue = MAX(val1,val2);
    Starttabcolumn[rowindex].Istart =
    setStarttabentry(Atabcolumn[rowindex].Ivalue, &westAffinealignDPentry, &Swe,
                     replacement, gap_opening, gap_extension, Affine_I);

    /*calculate totalvalue*/
    temp = MAX3(Atabcolumn[rowindex].Rvalue,
                Atabcolumn[rowindex].Dvalue,
                Atabcolumn[rowindex].Ivalue);
    Atabcolumn[rowindex].totalvalue = temp > 0 ? temp : 0;

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
    northwestAffinealignDPentry=westAffinealignDPentry;
    Snw=Swe;
  }
}

static Gtmaxcoordvalue *evaluateallAStabcolumns(const GtUchar *useq,
                                                GtUword ustart,
                                                GtUword ulen,
                                                const GtUchar *vseq,
                                                GtUword vstart,
                                                GtUword vlen,
                                                AffinealignDPentry *Atabcolumn,
                                                Starttabentry *Starttabcolumn,
                                                GtWord matchscore,
                                                GtWord mismatchscore,
                                                GtWord gap_opening,
                                                GtWord gap_extension)
{
  GtUword colindex;
  Gtmaxcoordvalue *max;

  firstAStabcolumn(ulen, Atabcolumn, Starttabcolumn,
                   gap_opening, gap_extension);

  max = gt_max_new();
  for (colindex = 1UL; colindex <= vlen; colindex++)
  {
    nextAStabcolumn(useq, ustart, ulen, vseq[vstart+colindex-1],
                    Atabcolumn, Starttabcolumn, matchscore,mismatchscore,
                    gap_opening, gap_extension, colindex, max);
  }
  return max;
}

static GtUword gt_calc_affinealign_linear_local(const GtUchar *useq,
                                                GtUword ustart,
                                                GtUword ulen,
                                                const GtUchar *vseq,
                                                GtUword vstart,
                                                GtUword vlen,
                                                GtAlignment *align,
                                                GtWord matchscore,
                                                GtWord mismatchscore,
                                                GtWord gap_opening,
                                                GtWord gap_extension)
{
  GtUword score, ulen_part, ustart_part, vlen_part, vstart_part,
          match_cost, mismatch_cost, gap_opening_cost, gap_extension_cost;
  AffinealignDPentry *Atabcolumn;
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
                                         &match_cost,
                                         &mismatch_cost,
                                         &gap_opening_cost,
                                         &gap_extension_cost);

    gt_calc_affinealign_linear(useq, ustart_part, ulen_part,
                               vseq, vstart_part, vlen_part,
                               align, match_cost, mismatch_cost,
                               gap_opening_cost,gap_extension_cost);
  }else
  {
     gt_alignment_set_seqs(align,(const GtUchar*)"", 0, (const GtUchar*)"", 0);
     score = 0;
  }

  gt_max_delete(max);
  gt_free(Atabcolumn);
  gt_free(Starttabcolumn);

  return(score);
}

void gt_computeaffinelinearspace_local(GtAlignment *align,
                                       const GtUchar *useq,
                                       GtUword ustart,
                                       GtUword ulen,
                                       const GtUchar *vseq,
                                       GtUword vstart,
                                       GtUword vlen,
                                       GtWord matchscore,
                                       GtWord mismatchscore,
                                       GtWord gap_opening,
                                       GtWord gap_extension)
{
  (void) gt_calc_affinealign_linear_local(useq, ustart, ulen,
                                          vseq, vstart, vlen,
                                          align, matchscore,mismatchscore,
                                          gap_opening, gap_extension);
}

/*----------------------------checkfunctions--------------------------*/
void gt_checkaffinelinearspace(GT_UNUSED bool forward,
                               const GtUchar *useq,
                               GtUword ulen,
                               const GtUchar *vseq,
                               GtUword vlen)
{
  GtAlignment *align_linear, *align_square;
  GtUword affine_score1, affine_score2, affine_score3,
          matchcost = 0, mismatchcost = 4, gap_opening = 4, gap_extension = 1;
  /* immediate result, because affinealign (square) cannot
   * handle lower/upper cases*/
  GtUchar *low_useq, *low_vseq;

  gt_assert(useq && vseq);
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

  low_useq = sequence_to_lower_case(useq, ulen);
  low_vseq = sequence_to_lower_case(vseq, vlen);

  align_linear = gt_alignment_new_with_seqs(low_useq, ulen, low_vseq, vlen);

  affine_score1 = gt_calc_affinealign_linear(low_useq, 0, ulen,
                                             low_vseq, 0, vlen,
                                             align_linear, matchcost,
                                             mismatchcost, gap_opening,
                                             gap_extension);
  affine_score2 = gt_alignment_eval_generic_with_affine_score(false,
                                                      align_linear, matchcost,
                                                      mismatchcost, gap_opening,
                                                      gap_extension);

  if (affine_score1 != affine_score2)
  {
    fprintf(stderr,"gt_calc_affinealign_linear = "GT_WU" != "GT_WU
            " = gt_alignment_eval_with_affine_score\n", affine_score1,
                                                        affine_score2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  align_square = gt_affinealign(low_useq, ulen, low_vseq, vlen, matchcost,
                                mismatchcost, gap_opening, gap_extension);
  affine_score3 = gt_alignment_eval_generic_with_affine_score(false,
                                                      align_square, matchcost,
                                                      mismatchcost, gap_opening,
                                                      gap_extension);

  if (affine_score1 != affine_score3)
  {
    fprintf(stderr,"gt_calc_affinealign_linear = "GT_WU" != "GT_WU
            " = gt_affinealign\n", affine_score1, affine_score3);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_free(low_useq);
  gt_free(low_vseq);
  gt_alignment_delete(align_linear);
  gt_alignment_delete(align_square);
}

void gt_checkaffinelinearspace_local(GT_UNUSED bool forward,
                                     const GtUchar *useq,
                                     GtUword ulen,
                                     const GtUchar *vseq,
                                     GtUword vlen)
{
  GtAlignment *align;
  GtUword affine_score1, affine_score2;
  GtWord matchscore = 6, mismatchscore = -3,
         gap_opening = -2, gap_extension = -1;
  GtUchar *low_useq, *low_vseq;

  gt_assert(useq && vseq);
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

  low_useq = sequence_to_lower_case(useq, ulen);
  low_vseq = sequence_to_lower_case(vseq, vlen);

  align = gt_alignment_new();
  affine_score1 = gt_calc_affinealign_linear_local(low_useq, 0, ulen, low_vseq,
                                                   0, vlen, align, matchscore,
                                                   mismatchscore, gap_opening,
                                                   gap_extension);

  affine_score2 = gt_alignment_eval_generic_with_affine_score(false, align,
                                                 matchscore, mismatchscore,
                                                 gap_opening, gap_extension);

  if (affine_score1 != affine_score2)
  {
    fprintf(stderr,"gt_calc_affinealign_linear_local = "GT_WU" != "GT_WU
            " = gt_alignment_eval_with_affine_score\n", affine_score1,
                                                        affine_score2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  gt_free(low_useq);
  gt_free(low_vseq);
  gt_alignment_delete(align);
}
