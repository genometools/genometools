#include "core/types_api.h"

GtUword diagonalband_affine_distance_only(const GtUchar *useq,
                                         const GtUword ustart,
                                         const GtUword ulen,
                                         const GtUchar *vseq,
                                         const GtUword vstart,
                                         const  GtUword vlen,
                                         const GtWord left_dist,
                                         const GtWord right_dist,
                                         const GtWord matchcost,
                                         const GtWord mismatchcost,
                                         const GtWord gap_opening,
                                         const GtWord gap_extension);
