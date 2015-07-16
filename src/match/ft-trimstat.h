#ifndef FT_TRIMSTAT_H
#define FT_TRIMSTAT_H
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "core/types_api.h"

typedef struct Trimstat Trimstat;

Trimstat *trimstat_new(double errorpercentage,
                       GtUword minmatchpercentage,
                       GtUword maxalignedlendifference);

void trimstat_add(Trimstat *trimstat,bool diedout,
                  GtUword sumvalid,
                  GtUword maxvalid,
                  GtUword d,
                  size_t spaceforfront,
                  GtUword cache_size);

void trimstat_delete(Trimstat *trimstat,double total_time,bool verbose);

#endif
