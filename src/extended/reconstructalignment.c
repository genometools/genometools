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
#include "extended/alignment.h"
#include "extended/linearalign_utilities.h"

#include "extended/reconstructalignment.h"

GtUword construct_trivial_deletion_alignment(GtAlignment *align,
                                             GtUword len,
                                             GtUword gapcost)
{
  GtUword idx;

  for (idx = 0; idx < len; idx ++)
  {
    gt_alignment_add_deletion(align);
  }

  return (len*gapcost);
}

GtUword construct_trivial_insertion_alignment(GtAlignment *align,
                                              GtUword len,
                                              GtUword gapcost)
{
  GtUword idx;

  for (idx = 0; idx < len; idx ++)
  {
    gt_alignment_add_insertion(align);
  }

  return (len*gapcost);
}

/* reconstruct alignment from square space table ED */
void reconstructalignment_from_EDtab(GtAlignment *align, GtUword **E,
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
  gt_assert(align && E);
  i = ulen;
  j = vlen;

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
}

/* reconstruct alignment from crosspoint table
 * crosspoints relating to midolumn
 */
void reconstructalignment_from_Ctab(GtAlignment *align,
                                    const GtUword *Ctab,
                                    const GtUchar *useq,
                                    GtUword ustart,
                                    const GtUchar *vseq,
                                    GtUword vstart,
                                    GtUword vlen,
                                    GtUword matchcost,
                                    GtUword mismatchcost,
                                    GtUword gap_opening,
                                    GtUword gap_extension)
{
  GtUword i,j, indel, repl;
  gt_assert(align != NULL && Ctab != NULL);
  for (i = vlen; i > 0; i--) {
    if (Ctab[i] == Ctab[i-1] + 1)
    {
      if (i > 1 && Ctab[i-2] == Ctab[i-1])
        indel = 2*gap_extension + gap_opening;
      else
        indel = (2*gap_extension + 2*gap_opening);
      if (tolower((int)vseq[vstart+i-1])==tolower((int)useq[ustart+Ctab[i]-1]))
        repl = matchcost;
      else
        repl = mismatchcost;
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
      if (tolower((int)vseq[vstart+i-1]) ==
                                         tolower((int)useq[ustart+Ctab[i]-j-1]))
      {
        repl = matchcost;
      }
      else
        repl = mismatchcost;
      if (indel>repl)
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
void reconstructalignment_from_Dtab(GtAlignment *align,
                                    const Diagentry *Dtab,GtUword ulen,
                                    GtUword vlen)
{
  GtUword i,j;

  gt_assert(align != NULL && Dtab != NULL);

  for (j = ulen; j > Dtab[vlen].currentrowindex; j--)
  {
    gt_alignment_add_deletion(align);
  }
  for (i = vlen; i > 0; i--) {
    gt_assert(Dtab[i].currentrowindex != GT_UWORD_MAX);
    if (Dtab[i].currentrowindex == Dtab[i-1].currentrowindex + 1)
    {
      if (Dtab[i].edge == Linear_R)
       gt_alignment_add_replacement(align);

      else if (Dtab[i].edge == Linear_D)
      {
         gt_alignment_add_deletion(align);
         gt_alignment_add_insertion(align);
      }
      else if (Dtab[i].edge == Linear_I)
      {
         gt_alignment_add_insertion(align);
         gt_alignment_add_deletion(align);
      }
    }
    else if (Dtab[i].currentrowindex == Dtab[i-1].currentrowindex)
      gt_alignment_add_insertion(align);

    else if (Dtab[i].currentrowindex > Dtab[i-1].currentrowindex)
    {
      if (Dtab[i].edge == Linear_R)
      {
        gt_alignment_add_replacement(align);
        for (j = 0; j < (Dtab[i].currentrowindex -
                         Dtab[i-1].currentrowindex)-1; j++)
        {
          gt_alignment_add_deletion(align);
        }
      }
      else if (Dtab[i].edge == Linear_I)
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
        for (j = 0; j < (Dtab[i].currentrowindex -
                         Dtab[i-1].currentrowindex)-1; j++)
        {
          gt_alignment_add_deletion(align);
        }
          gt_alignment_add_replacement(align);

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
void reconstructalignment_from_affineDtab(GtAlignment *align,
                                          const AffineDiagentry *Dtab,
                                          AffineAlignEdge edge,
                                          const GtUchar *useq, GtUword ulen,
                                          const GtUchar *vseq, GtUword vlen)
{
  GtUword i,j;
  Diagentry node, prevnode;
  AffineAlignEdge prevedge;
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
    switch (prevnode.edge) {
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
    }

    gt_assert(prevnode.currentrowindex != GT_UWORD_MAX);
    if (prevnode.currentrowindex == node.currentrowindex + 1)
    {
      if (prevedge == Affine_R)
        {gt_alignment_add_replacement(align);}
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
      {gt_alignment_add_insertion(align);}

    else if (prevnode.currentrowindex > node.currentrowindex)
    {
      if (prevedge == Affine_R)
      {
        gt_alignment_add_replacement(align);

        for (j = 0; j < (prevnode.currentrowindex -
                         node.currentrowindex)-1; j++)
        {
          gt_alignment_add_deletion(align);
        }
      }
      else if (prevedge == Affine_I)
      {
        gt_alignment_add_insertion(align);
        for (j = 0; j < (prevnode.currentrowindex -
                         node.currentrowindex); j++)
        {
          gt_alignment_add_deletion(align);
        }
      }
      else
      {
        for (j = 0; j < (prevnode.currentrowindex -
                         node.currentrowindex)-1; j++)
        {
          gt_alignment_add_deletion(align);
        }
        if (prevnode.edge == Affine_I)
        {
          if (tolower((int)vseq[i-1]) ==
                                      tolower((int)useq[node.currentrowindex]))
          {
            gt_alignment_add_replacement(align);
          }
          else
          {
            gt_alignment_add_deletion(align);
            gt_alignment_add_insertion(align);
          }
        }
        else
          gt_alignment_add_replacement(align);
      }
    }
    prevedge = prevnode.edge;
  }
  for (j = node.currentrowindex; j > 0; j--)
  {
    gt_alignment_add_deletion(align);
  }
}
