#ifndef AFFINEALIGN_LINEAR_H
#define AFFINEALIGN_LINEAR_H

#include "core/types_api.h"
void gt_computeaffinelinearspace(bool showevalue, 
                                 const GtUchar *useq,
                                 GtUword ulen,
                                 const GtUchar *vseq,
                                 GtUword vlen,
                                 const GtWord replacement_cost,
                                 const GtWord gap_opening,
                                 const GtWord gap_extension,
                                 FILE *fp);
#endif
