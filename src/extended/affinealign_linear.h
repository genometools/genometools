#ifndef AFFINEALIGN_LINEAR_H
#define AFFINEALIGN_LINEAR_H

#include "core/unused_api.h"
#include "core/types_api.h"
#include "extended/alignment.h"

void gt_checkaffinelinearspace(GT_UNUSED bool forward,
                               const GtUchar *useq,
                               GtUword ulen,
                               const GtUchar *vseq,
                               GtUword vlen);

void gt_computeaffinelinearspace(bool showevalue, 
                                 const GtUchar *useq,
                                 GtUword ulen,
                                 const GtUchar *vseq,
                                 GtUword vlen,
                                 const GtWord replacement_cost,
                                 const GtWord gap_opening,
                                 const GtWord gap_extension,
                                 FILE *fp);
GtUword gt_calc_affinealign_linear(const GtUchar *useq,
                                          const GtUword ulen,
                                          const GtUchar *vseq,
                                          const GtUword vlen,
                                          GtAlignment *align,
                                          const GtWord replacement_cost,
                                          const GtWord gap_opening,
                                          const GtWord gap_extension);
#endif
