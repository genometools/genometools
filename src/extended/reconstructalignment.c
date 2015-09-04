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

#include "extended/reconstructalignment.h"

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
      if (vseq[vstart+i-1]==useq[ustart+Ctab[i]-1])
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
      if (vseq[vstart+i-1]==useq[ustart+Ctab[i]-j-1])
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
