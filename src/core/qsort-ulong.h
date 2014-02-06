#ifndef QSORT_ULONG_H
#define QSORT_ULONG_H

#include <stdbool.h>
#include "core/types_api.h"

void gt_direct_qsort_ulong (GtUword insertionsortthreshold,
                            bool handlenotswapped,
                            GtUword *arr,
                            GtUword len);

#endif
