#ifndef MSORT_H
#define MSORT_H

#include <stdlib.h>

/* merge sort; the interface equals qsort(3) */
void gt_msort(void *base, size_t nmemb, size_t size,
              int (*compar)(const void *, const void *));

void gt_msort_r(void *base, size_t nmemb, size_t size, void *comparinfo,
                int(*compar)(void *, const void *, const void *));

#endif
