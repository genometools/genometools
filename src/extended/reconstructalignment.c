#include "extended/reconstructalignment.h"

void reconstructalignment(GtAlignment *align,
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
