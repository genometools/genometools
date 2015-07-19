#ifndef AFFINEALIGN_LINEAR_LOCAL_H
#define AFFINEALIGN_LINEAR_LOCAL_H

#include "core/types_api.h"

void gt_computeaffinelinearspace_local(
                                 const GtUchar *useq, GtUword ulen,
                                 const GtUchar *vseq, GtUword vlen,
                                 const GtWord replacement_score,
                                 const GtWord gap_opening,
                                 const GtWord gap_extension,
                                 FILE *fp);
#endif
