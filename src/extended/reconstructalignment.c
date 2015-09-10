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
#include "extended/reconstructalignment.h"
#include "extended/alignment.h"

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

/* reconstruct alignment from crosspoint table */
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

void evaluate_crosspoints_from_2dimtab(GtUword **E,
                                       GtUword *Ctab,
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
  GtUword idx, jdx;

  idx = ulen;
  jdx = vlen;
  while (jdx > 1 || idx > 0)
  {
    if (idx > 0 && jdx > 0 && E[idx][jdx] == E[idx-1][jdx-1] +
       (tolower((int)useq[ustart+idx-1]) == tolower((int) vseq[vstart+jdx-1]) ?
                                                     matchcost : mismatchcost))
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
