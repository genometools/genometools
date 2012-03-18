#ifdef  QSORTNAME
#undef  QSORTNAME
#endif

#define QSORTNAME(NAME)                   NAME##_ulong
#define ARRAY_GET_ulong(ARR,RELIDX)       ARR[RELIDX]
#define ARRAY_SET_ulong(ARR,RELIDX,VALUE) ARR[RELIDX] = VALUE

typedef unsigned long QSORTNAME(Sorttype);

#include "core/qsort-direct.gen"
