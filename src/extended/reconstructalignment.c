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
#include "extended/alignment.h"
#include "extended/linspace_management.h"
#include "extended/reconstructalignment.h"

GtUword gt_reconstructalignment_trivial_deletion(GtAlignment *align,
                                                 GtUword len,
                                                 GtUword gapcost)
{
  GtUword idx;

  for (idx = 0; idx < len; idx ++)
  {
    gt_alignment_add_deletion(align);
  }
  return len * gapcost;
}

GtUword gt_reconstructalignment_trivial_insertion(GtAlignment *align,
                                                  GtUword len,
                                                  GtUword gapcost)
{
  GtUword idx;

  for (idx = 0; idx < len; idx ++)
  {
    gt_alignment_add_insertion(align);
  }

  return len * gapcost;
}

/* reconstruct global alignment from square space table E */
void gt_reconstructalignment_from_EDtab(GtAlignment *align,
                                        GtUword * const *E,
                                        const GtUchar *useq,
                                        GtUword ustart,
                                        GtUword ulen,
                                        const GtUchar *vseq,
                                        GtUword vstart,
                                        GtUword vlen,
                                        const GtScoreHandler *scorehandler)
{
  GtUword i = ulen, j = vlen, gapcost;
  const GtUchar *uptr = useq + ustart - 1, *vptr = vseq + vstart - 1;

  gt_assert(align && E && scorehandler);
  gapcost = gt_scorehandler_get_gapscore(scorehandler);
  while (i > 0 || j > 0)
  {
    if (i > 0 && j > 0 && E[i][j] == E[i-1][j-1] +
                                     gt_scorehandler_get_replacement(
                                                          scorehandler,
                                                          uptr[i],vptr[j]))
    {
      gt_alignment_add_replacement(align);
      i--;
      j--;
      continue;
    }
    if (j > 0 && E[i][j] == E[i][j-1] + gapcost)
    {
      gt_alignment_add_insertion(align);
      j--;
      continue;
    }
    if (i > 0 && E[i][j] == E[i-1][j] + gapcost)
    {
      gt_alignment_add_deletion(align);
      i--;
      continue;
    }
    gt_assert(false);
  }
}

/* reconstruct local alignment from square space table Ltab */
void gt_reconstructalignment_from_Ltab(GtAlignment *align,
                                       GtWord **Ltabcolumn,
                                       GtMaxcoordvalue *max,
                                       const GtUchar *useq,
                                       GtUword ustart,
                                       GT_UNUSED GtUword ulen,
                                       const GtUchar *vseq,
                                       GtUword vstart,
                                       GT_UNUSED GtUword vlen,
                                       const GtScoreHandler *scorehandler)
{
  GtUword i, j;
  GtWord gapscore;
  GtUwordPair max_end;

  gt_assert(align && Ltabcolumn && scorehandler);
  max_end = gt_maxcoordvalue_get_end(max);
  i = max_end.a;
  j = max_end.b;
  gt_assert(i <= ulen && j <= vlen);
  gapscore = gt_scorehandler_get_gapscore(scorehandler);
  while ((i > 0 || j > 0) && Ltabcolumn[i][j] != 0)
  {
    if (i > 0 && j > 0 && Ltabcolumn[i][j] == Ltabcolumn[i-1][j-1] +
                                              gt_scorehandler_get_replacement(
                                                    scorehandler,
                                                    useq[ustart+i-1],
                                                    vseq[vstart+j-1]))
    {
      gt_alignment_add_replacement(align);
      i--;
      j--;
      continue;
    }
    if (j > 0 && Ltabcolumn[i][j] == Ltabcolumn[i][j-1] + gapscore)
    {
      gt_alignment_add_insertion(align);
      j--;
      continue;
    }
    if (i > 0 && Ltabcolumn[i][j] == Ltabcolumn[i-1][j] + gapscore)
    {
      gt_alignment_add_deletion(align);
      i--;
      continue;
    }
    gt_assert(false);
  }
  gt_maxcoordvalue_set_start(max,i,j);
}

/* reconstruct alignment from crosspoint table
   crosspoints relating to midolumn */
void gt_reconstructalignment_from_Ctab(GtAlignment *align,
                                       const GtUword *Ctab,
                                       const GtUchar *useq,
                                       GtUword ustart,
                                       const GtUchar *vseq,
                                       GtUword vstart,
                                       GtUword vlen,
                                       const GtScoreHandler *scorehandler)
{
  GtUword i,j, repl, indel, gap_opening, gap_extension;
  gt_assert(align && Ctab && scorehandler);

  gap_opening = gt_scorehandler_get_gap_opening(scorehandler);
  gap_extension = gt_scorehandler_get_gapscore(scorehandler);

  for (i = vlen; i > 0; i--) {
    if (Ctab[i] == Ctab[i-1] + 1)
    {
      if (i > 1 && Ctab[i-2] == Ctab[i-1])
        indel = 2*gap_extension + gap_opening;
      else
        indel = (2*gap_extension + 2*gap_opening);

      repl = gt_scorehandler_get_replacement(scorehandler,
                                             vseq[vstart+i-1],
                                             useq[ustart+Ctab[i]-1]);

      if (indel>repl)
        gt_alignment_add_replacement(align);
      else
      {
        gt_alignment_add_deletion(align);
        gt_alignment_add_insertion(align);
      }
    }
    else if (Ctab[i] == Ctab[i-1])
      gt_alignment_add_insertion(align);
    else if (Ctab[i] > Ctab[i-1]) {
      indel = 0; repl = 0;
      for (j = 0; j < (Ctab[i]-Ctab[i-1])-1; j++)
       gt_alignment_add_deletion(align);
      /*replacmente or insertion+deletion*/
      if (i > 1 && Ctab[i-2] == Ctab[i-1])
        indel = 2*gap_extension;
      else
        indel = (2*gap_extension+gap_opening);

      repl = gt_scorehandler_get_replacement (scorehandler,
                                              vseq[vstart+i-1],
                                              useq[ustart+Ctab[i]-j-1]);
      if (indel > repl)
        gt_alignment_add_replacement(align);
      else
      {
        gt_alignment_add_deletion(align);
        gt_alignment_add_insertion(align);
      }
    }
  }
  for (j = Ctab[0]; j > 0; j--)
    gt_alignment_add_deletion(align);
}

/*reconstruct alignment from crosspoints, crosspoints relating to diagonalband*/
void gt_reconstructalignment_from_Dtab(GtAlignment *align,
                                       const GtDiagAlignentry *Dtab,
                                       GtUword ulen, GtUword vlen)
{
  GtUword i,j;

  gt_assert(align && Dtab);

  for (j = ulen; j > Dtab[vlen].currentrowindex; j--)
  {
    gt_alignment_add_deletion(align);
  }
  for (i = vlen; i > 0; i--) {
    gt_assert(Dtab[i].currentrowindex != GT_UWORD_MAX);
    if (Dtab[i].currentrowindex == Dtab[i-1].currentrowindex + 1)
    {
      if (Dtab[i].last_type == Linear_R)
       gt_alignment_add_replacement(align);

      else if (Dtab[i].last_type == Linear_D)
      {
         gt_alignment_add_deletion(align);
         gt_alignment_add_insertion(align);
      }
      else if (Dtab[i].last_type == Linear_I)
      {
         gt_alignment_add_insertion(align);
         gt_alignment_add_deletion(align);
      }
    }
    else if (Dtab[i].currentrowindex == Dtab[i-1].currentrowindex)
      gt_alignment_add_insertion(align);

    else if (Dtab[i].currentrowindex > Dtab[i-1].currentrowindex)
    {
      if (Dtab[i].last_type == Linear_R)
      {
        gt_alignment_add_replacement(align);
        for (j = 0; j < (Dtab[i].currentrowindex -
                         Dtab[i-1].currentrowindex)-1; j++)
        {
          gt_alignment_add_deletion(align);
        }
      }
      else if (Dtab[i].last_type == Linear_I)
      {
        gt_alignment_add_insertion(align);
        for (j = 0; j < (Dtab[i].currentrowindex -
                         Dtab[i-1].currentrowindex); j++)
        {
          gt_alignment_add_deletion(align);
        }
      }
      else
      {
        gt_assert(false);
        /*for (j = 0; j < (Dtab[i].currentrowindex -
                         Dtab[i-1].currentrowindex)-1; j++)
        {
          gt_alignment_add_deletion(align);
        }
          gt_alignment_add_replacement(align);*/

      }
    }
  }
  for (j = Dtab[0].currentrowindex; j > 0; j--)
  {
    gt_alignment_add_deletion(align);
  }
}

/* reconstruct alignment from crosspoints (affine gapcosts),
 * crosspoints relating to diagonalband */
void gt_reconstructalignment_from_affineDtab(GtAlignment *align,
                                             const GtAffineDiagAlignentry *Dtab,
                                             GtAffineAlignEdge edge,
                                             GT_UNUSED const GtUchar *useq,
                                             GtUword ulen,
                                             GT_UNUSED const GtUchar *vseq,
                                             GtUword vlen)
{
  GtUword i,j;
  GtDiagAlignentry node, prevnode;
  GtAffineAlignEdge prevedge;
  gt_assert(align != NULL && Dtab != NULL);

  switch (edge) {
    case Affine_R:
      node = Dtab[vlen].val_R;
      break;
    case Affine_D:
      node = Dtab[vlen].val_D;
      break;
    case Affine_I:
      node = Dtab[vlen].val_I;
      break;
    default:
      gt_assert(false);
#ifdef NDEBUG
      exit(GT_EXIT_PROGRAMMING_ERROR);
#endif
  }

  for (j = ulen; j > node.currentrowindex; j--)
  {
    gt_alignment_add_deletion(align);
  }
  prevedge = edge;
  /*consider pairwise distance between prevnode and node */
  /*prevedge->prevnode, prevnode.edge->node, node.edge->nextnode*/
  for (i = vlen; i > 0; i--)
  {
    prevnode = node;
    switch (prevnode.last_type) {
      case Affine_R:
        node = Dtab[i-1].val_R;
        break;
      case Affine_D:
        node = Dtab[i-1].val_D;
        break;
      case Affine_I:
        node = Dtab[i-1].val_I;
        break;
      default:
        gt_assert(false);
#ifdef NDEBUG
        exit(GT_EXIT_PROGRAMMING_ERROR);
#endif
    }

    gt_assert(prevnode.currentrowindex != GT_UWORD_MAX);
    if (prevnode.currentrowindex == node.currentrowindex + 1)
    {
      if (prevedge == Affine_R)
      {
        gt_alignment_add_replacement(align);
      }
      else if (prevedge == Affine_D)
      {
         gt_alignment_add_deletion(align);
         gt_alignment_add_insertion(align);
      }
      else if (prevedge == Affine_I)
      {
         gt_alignment_add_insertion(align);
         gt_alignment_add_deletion(align);
      }
    }
    else if (prevnode.currentrowindex == node.currentrowindex)
    {
      gt_alignment_add_insertion(align);
    }
    else if (prevnode.currentrowindex > node.currentrowindex)
    {
      if (prevedge == Affine_R)
      {
        gt_alignment_add_replacement(align);

        for (j = 0; j < prevnode.currentrowindex - node.currentrowindex - 1;
             j++)
        {
          gt_alignment_add_deletion(align);
        }
      }
      else if (prevedge == Affine_I)
      {
        gt_alignment_add_insertion(align);
        for (j = 0; j < prevnode.currentrowindex - node.currentrowindex; j++)
        {
          gt_alignment_add_deletion(align);
        }
      }
      else
      {
        gt_assert(false);
      }
    }
    prevedge = prevnode.last_type;
  }
  for (j = node.currentrowindex; j > 0; j--)
  {
    gt_alignment_add_deletion(align);
  }
}
