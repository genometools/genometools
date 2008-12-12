#ifndef MSORT_H
#define MSORT_H

#include <stdlib.h>
#include "core/fptr_api.h"

/* merge sort; the interface equals qsort(3) */
void gt_msort(void *base, size_t nmemb, size_t size,
              GtCompare compar);

void gt_msort_r(void *base, size_t nmemb, size_t size, void *comparinfo,
                GtCompareWithData compar);

#endif
