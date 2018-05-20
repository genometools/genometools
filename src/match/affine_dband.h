#ifndef AFFINE_DBAND_H
#define AFFINE_DBAND_H
#include <stdbool.h>
#include <inttypes.h>
#include "ft-eoplist.h"

typedef struct GtAffineDPreservoir GtAffineDPreservoir;

GtAffineDPreservoir *gt_affine_diagonalband_new(bool opt_memory,
                                                bool keepcolumns,
                                                GtUword max_ulen,
                                                GtUword max_vlen);

void gt_affine_diagonalband_delete(GtAffineDPreservoir *adpr);

GtUword gt_affine_iter_diagonalband_align(GtEoplist *eoplist,
                                       GtAffineDPreservoir *adpr,
                                       int8_t gap_opening, /* > 0 */
                                       int8_t gap_extension, /* > 0 */
                                       const int8_t * const *scorematrix2D,
                                       int8_t smallest_score,
                                       const GtUchar *useq,
                                       GtUword ulen,
                                       const GtUchar *vseq,
                                       GtUword vlen,
                                       bool no_score_run,
                                       GtUword expected_score);

#endif
