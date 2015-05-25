#ifndef FT_POLISH_H
#define FT_POLISH_H
#include <stdint.h>
#include "core/types_api.h"

typedef struct Polishing_info Polishing_info;

Polishing_info *polishing_info_new(GtUword cut_depth,
                                   double errorpercentage);

void polishing_info_delete(Polishing_info *pol_info);

bool history_is_polished(const Polishing_info *pol_info,uint64_t matchhistory);

bool history_is_polished_brute_force(const Polishing_info *pol_info,
                                     uint64_t matchhistory);

uint64_t polishing_info_maxvalue(const Polishing_info *pol_info);

#endif
