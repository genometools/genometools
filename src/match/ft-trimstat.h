#ifndef FT_TRIMSTAT_H
#define FT_TRIMSTAT_H
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#ifndef OUTSIDE_OF_GT
#else
#include "gt-defs.h"
#endif

typedef struct GtFtTrimstat GtFtTrimstat;

GtFtTrimstat *gt_ft_trimstat_new(void);

#ifndef NDEBUG
void gt_ft_trimstat_add(GtFtTrimstat *trimstat,
                        bool diedout,
                        GtUword sumvalid,
                        GtUword maxvalid,
                        GtUword d,
                        size_t spaceforfront,
                        GtUword cache_size);
void gt_ft_trimstat_add_matchlength(GtFtTrimstat *trimstat,
                                    uint32_t matchlength);
#else
void gt_ft_trimstat_add(GT_UNUSED GtFtTrimstat *trimstat,
                        GT_UNUSED bool diedout,
                        GT_UNUSED GtUword sumvalid,
                        GT_UNUSED GtUword maxvalid,
                        GT_UNUSED GtUword d,
                        GT_UNUSED size_t spaceforfront,
                        GT_UNUSED GtUword cache_size);
void gt_ft_trimstat_add_matchlength(GT_UNUSED GtFtTrimstat *trimstat,
                                    GT_UNUSED uint32_t matchlength);
#endif

void gt_ft_trimstat_delete(GtFtTrimstat *trimstat,bool verbose);

#endif
