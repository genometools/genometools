/*
  Copyright (C) 2015 Annika Seidel, annika.seidel@studium.uni-hamburg.de
  Copyright (c) 2007-2009 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include <limits.h>
#include "core/assert_api.h"
#include "core/array2dim_api.h"
#include "core/minmax.h"
#include "extended/linearalign_affinegapcost.h"
#include "extended/linearalign_utilities.h"

#include "extended/affinealign.h"

static void affinealign_fill_table(AffinealignDPentry **dptable,
                                   const GtUchar *u, GtUword ulen,
                                   const GtUchar *v, GtUword vlen,
                                   GtUword matchcost, GtUword mismatchcost,
                                   GtUword gap_opening, GtUword gap_extension,
                                   AffineAlignEdge edge)
{
  GtUword i, j, Rvalue, Dvalue, Ivalue, minvalue;
  int rcost;
  gt_assert(dptable && u && v);
  /*gt_assert(ulen && vlen);*/
  for (i = 0; i <= ulen; i++) {
    for (j = 0; j <= vlen; j++) {
      if (!i && !j) {
        switch (edge) {
          case Affine_R:
            dptable[0][0].Rvalue = 0;
            dptable[0][0].Dvalue = GT_WORD_MAX;
            dptable[0][0].Ivalue = GT_WORD_MAX;
            break;
          case Affine_D:
            dptable[0][0].Rvalue = GT_WORD_MAX;
            dptable[0][0].Dvalue = 0;
            dptable[0][0].Ivalue = GT_WORD_MAX;
            break;
          case Affine_I:
            dptable[0][0].Rvalue = GT_WORD_MAX;
            dptable[0][0].Dvalue = GT_WORD_MAX;
            dptable[0][0].Ivalue = 0;
            break;
          default:
            dptable[0][0].Rvalue = 0;
            dptable[0][0].Dvalue = gap_opening;
            dptable[0][0].Ivalue = gap_opening;
        }
      }
      else {
        /* compute A_affine(i,j,R) */
        if (!i || !j)
          dptable[i][j].Rvalue = LONG_MAX;
        else {
          rcost  = (u[i-1] == v[j-1]) ? matchcost : mismatchcost;
          Rvalue = add_safe_max(dptable[i-1][j-1].Rvalue, rcost);
          Dvalue = add_safe_max(dptable[i-1][j-1].Dvalue, rcost);
          Ivalue = add_safe_max(dptable[i-1][j-1].Ivalue, rcost);
          minvalue = MIN3(Rvalue, Dvalue, Ivalue);
          gt_assert(minvalue != LONG_MAX);
          dptable[i][j].Rvalue = minvalue;
          /* set backtracing edge */
          if (Rvalue == minvalue)
            dptable[i][j].Redge = Affine_R;
          else if (Dvalue == minvalue)
            dptable[i][j].Redge = Affine_D;
          else /* Ivalue == minvalue */
            dptable[i][j].Redge = Affine_I;
        }
        /* compute A_affine(i,j,D) */
        if (!i)
          dptable[i][j].Dvalue = LONG_MAX;
        else {
          Rvalue = add_safe_max(dptable[i-1][j].Rvalue,
                                gap_opening + gap_extension);
          Dvalue = add_safe_max(dptable[i-1][j].Dvalue, gap_extension);
          Ivalue = add_safe_max(dptable[i-1][j].Ivalue,
                                gap_opening + gap_extension);
          minvalue = MIN3(Rvalue, Dvalue, Ivalue);
          gt_assert(minvalue != ULONG_MAX);
          dptable[i][j].Dvalue = minvalue;
          /* set backtracing edge */
          if (Rvalue == minvalue)
            dptable[i][j].Dedge = Affine_R;
          else if (Dvalue == minvalue)
            dptable[i][j].Dedge = Affine_D;
          else /* Ivalue == minvalue */
            dptable[i][j].Dedge = Affine_I;
        }
        /* compute A_affine(i,j,I) */
        if (!j)
          dptable[i][j].Ivalue = LONG_MAX;
        else {
          Rvalue = add_safe_max(dptable[i][j-1].Rvalue,
                                gap_opening + gap_extension);
          Dvalue = add_safe_max(dptable[i][j-1].Dvalue,
                                gap_opening + gap_extension);
          Ivalue = add_safe_max(dptable[i][j-1].Ivalue, gap_extension);
          minvalue = MIN3(Rvalue, Dvalue, Ivalue);
          gt_assert(minvalue != LONG_MAX);
          dptable[i][j].Ivalue = minvalue;
          /* set backtracing edge */
          if (Rvalue == minvalue)
            dptable[i][j].Iedge = Affine_R;
          else if (Dvalue == minvalue)
            dptable[i][j].Iedge = Affine_D;
          else /* Ivalue == minvalue */
            dptable[i][j].Iedge = Affine_I;
        }
      }
    }
  }
}

void affinealign_traceback(GtAlignment *a,
                           AffinealignDPentry * const *dptable,
                           GtUword i, GtUword j)
{
  GtWord minvalue;
  AffineAlignEdge edge;
  gt_assert(a && dptable);
  /* determine min{A_affine(m,n,x) | x in {R,D,I}} */
  minvalue = MIN3(dptable[i][j].Rvalue, dptable[i][j].Dvalue,
                  dptable[i][j].Ivalue);
  if (dptable[i][j].Rvalue == minvalue)
    edge = Affine_R;
  else if (dptable[i][j].Dvalue == minvalue)
    edge = Affine_D;
  else /* dptable[i][j].Ivalue == minvalue */
    edge = Affine_I;
  /* backtracing */
  while (i > 0 || j > 0) {
    switch (edge) {
      case Affine_R:
        gt_assert(dptable[i][j].Rvalue != LONG_MAX);
        gt_alignment_add_replacement(a);
        edge = dptable[i][j].Redge;
        gt_assert(i > 0 && j > 0);
        i--;
        j--;
        break;
      case Affine_D:
        gt_alignment_add_deletion(a);
        edge = dptable[i][j].Dedge;
        gt_assert(i);
        i--;
        break;
      case Affine_I:
        gt_alignment_add_insertion(a);
        edge = dptable[i][j].Iedge;
        gt_assert(j);
        j--;
        break;
      default:
        gt_assert(false);
    }
  }
}

GtAlignment* gt_affinealign(const GtUchar *u, GtUword ulen,
                            const GtUchar *v, GtUword vlen,
                            GtUword matchcost, GtUword mismatchcost,
                            GtUword gap_opening,
                            GtUword gap_extension)
{
  AffinealignDPentry **dptable;
  GtAlignment *align;

  gt_assert(u && v);
  gt_array2dim_malloc(dptable, ulen+1, vlen+1);
  affinealign_fill_table(dptable, u, ulen, v, vlen, matchcost, mismatchcost,
                         gap_opening, gap_extension, Affine_X);
  align = gt_alignment_new_with_seqs(u, ulen,  v, vlen);
  affinealign_traceback(align, dptable, ulen, vlen);
  gt_array2dim_delete(dptable);
  return align;
}

static void evaluate_affinecrosspoints_from_2dimtab(GtUword *Ctab,
                                                AffinealignDPentry **Atabcolumn,
                                                GtUword ulen, GtUword vlen,
                                                GtUword gap_opening,
                                                GtUword rowoffset,
                                                AffineAlignEdge edge)
{
  GtUword i, j;
  gt_assert(Atabcolumn != NULL);

  i = ulen;
  j = vlen;
  edge = minAdditionalCosts(&Atabcolumn[i][j], edge, gap_opening);

  while (i > 0 || j > 1) {
    switch (edge) {
      case Affine_R:
        gt_assert(Atabcolumn[i][j].Rvalue != LONG_MAX);
        Ctab[j-1] = i-1 + rowoffset;
        edge = Atabcolumn[i][j].Redge;
        gt_assert(i > 0 && j > 0);
        i--;
        j--;
        break;
      case Affine_D:
        edge = Atabcolumn[i][j].Dedge;
        gt_assert(i);
        i--;
        break;
      case Affine_I:
        Ctab[j-1] = i + rowoffset;
        edge = Atabcolumn[i][j].Iedge;
        gt_assert(j);
        j--;
        break;
      default:
        gt_assert(false);
    }
  }
}

void affine_ctab_in_square_space(GtUword *Ctab,
                                 const GtUchar *useq,
                                 GtUword ustart,
                                 GtUword ulen,
                                 const GtUchar *vseq,
                                 GtUword vstart,
                                 GtUword vlen,
                                 GtUword matchcost,
                                 GtUword mismatchcost,
                                 GtUword gap_opening,
                                 GtUword gap_extension,
                                 GtUword rowoffset,
                                 AffineAlignEdge from_edge,
                                 AffineAlignEdge to_edge)
{
  AffinealignDPentry **Atabcolumn;
  gt_assert(Ctab != NULL);

  gt_array2dim_malloc(Atabcolumn, (ulen+1), (vlen+1));
  affinealign_fill_table(Atabcolumn, &useq[ustart], ulen, &vseq[vstart], vlen,
                         matchcost, mismatchcost, gap_opening,
                         gap_extension, from_edge);

  evaluate_affinecrosspoints_from_2dimtab(Ctab, Atabcolumn, ulen, vlen,
                                          gap_opening, rowoffset, to_edge);
  gt_array2dim_delete(Atabcolumn);
}
