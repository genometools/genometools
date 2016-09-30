#ifndef FT_TRIMSTAT_H
#define FT_TRIMSTAT_H
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "core/types_api.h"

typedef struct GtFtTrimstat GtFtTrimstat;

GtFtTrimstat *gt_ft_trimstat_new(void);

void gt_ft_trimstat_add(GtFtTrimstat *trimstat,
                        bool diedout,
                        GtUword sumvalid,
                        GtUword maxvalid,
                        GtUword d,
                        size_t spaceforfront,
                        GtUword cache_size);

void gt_ft_trimstat_delete(GtFtTrimstat *trimstat,bool verbose);

#endif
