#ifndef QSORT_ULONG_H
#define QSORT_ULONG_H

#include <stdbool.h>

void gt_direct_qsort_ulong (unsigned long insertionsortthreshold,
                            bool handlenotswapped,
                            unsigned long *arr,
                            unsigned long len);

#endif
