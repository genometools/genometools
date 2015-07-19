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
   /* for (i = vlen; i > 1; i--) {
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
  if (Ctab[1]==0)
     gt_alignment_add_insertion(align);
  if (Ctab[1]>0)
    gt_alignment_add_replacement(align); 
  for (j = Ctab[1]; j > 1; j--)
    gt_alignment_add_deletion(align);*/

}


GtUword construct_trivial_alignment(GtAlignment *align, GtUword len, 
                                    const GtWord gapcost,
                                    void (*indel)(GtAlignment*))
{
  GtUword idx, distance=0;
  
  for (idx = 0; idx < len; idx ++)
  {
    indel(align);
    distance += gapcost;
  }
  return distance;
}
