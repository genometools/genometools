#ifndef FT_POLISH_H
#define FT_POLISH_H
#include <stdint.h>
#include "core/types_api.h"

#define DEFAULT_MATCHSCORE_BIAS 1.0  /* has no effect */

typedef struct
{
  int16_t score_sum, diff_from_max;
} Polishing_value;

typedef struct
{
  GtUword entries, cut_depth, mask;
  GtWord difference_score, match_score;
  Polishing_value *values;
} Polishing_info;

Polishing_info *polishing_info_new(GtUword cut_depth,
                                   double errorpercentage,
                                   double matchscore_bias);

void polishing_info_delete(Polishing_info *pol_info);

bool history_is_polished_brute_force(const Polishing_info *pol_info,
                                     uint64_t matchhistory);

uint64_t polishing_info_maxvalue(const Polishing_info *pol_info);

#define HISTORY_IS_POLISHED(POL_INFO,MATCHHISTORY,LSB)\
        ((POL_INFO)->values[LSB].diff_from_max >= 0 &&\
         (POL_INFO)->values[LSB].score_sum +\
         (POL_INFO)->values[((MATCHHISTORY) >> (POL_INFO)->cut_depth) &\
                             (POL_INFO)->mask].diff_from_max >= 0)

#endif
