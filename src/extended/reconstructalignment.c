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
                                    const GtUword vlen)
{
  GtUword i,j;

  gt_assert(align != NULL && Ctab != NULL);
  for (i = vlen; i > 0; i--) {
    if (Ctab[i] == Ctab[i-1] + 1)
      gt_alignment_add_replacement(align);
    else if (Ctab[i] == Ctab[i-1])
      gt_alignment_add_insertion(align);
    else if (Ctab[i] > Ctab[i-1]) {
      for (j = 0; j < (Ctab[i]-Ctab[i-1])-1; j++)
        gt_alignment_add_deletion(align);
      gt_alignment_add_replacement(align);

    }
  }
  for (j = Ctab[0]; j > 0; j--)
    gt_alignment_add_deletion(align);
}

GtUword construct_trivial_deletion_alignment(GtAlignment *align,
                                             const GtUword len,
                                             const GtWord gapcost)
{
  GtUword idx;

  for (idx = 0; idx < len; idx ++)
  {
    gt_alignment_add_deletion(align);
  }

  return (len*gapcost);
}

GtUword construct_trivial_insertion_alignment(GtAlignment *align,
                                              const GtUword len,
                                              const GtWord gapcost)
{
  GtUword idx;

  for (idx = 0; idx < len; idx ++)
  {
    gt_alignment_add_insertion(align);
  }

  return (len*gapcost);
}
