#ifndef FT_POLISH_H
#define FT_POLISH_H
#include <stdint.h>
#include "core/types_api.h"
/* The following values are used when no pruning is performed */
#define GT_MIN_PERC_MAT_HISTORY 1
#define GT_MAX_ALI_LEN_DIFF     UINT32_MAX

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

Polishing_info *polishing_info_new(double errorpercentage,
                                   GtUword history_size);

Polishing_info *polishing_info_new_with_bias(double errorpercentage,
                                             double matchscore_bias,
                                             GtUword history_size);

void polishing_info_delete(Polishing_info *pol_info);

bool history_is_polished_brute_force(const Polishing_info *pol_info,
                                     uint64_t matchhistory,
                                     bool withoutput);

uint64_t polishing_info_maxvalue(const Polishing_info *pol_info);

#define GT_HISTORY_IS_POLISHED(POL_INFO,MATCHHISTORY)\
        ((POL_INFO)->values[(MATCHHISTORY) & (POL_INFO)->mask].diff_from_max \
         >= 0 &&\
         (POL_INFO)->values[(MATCHHISTORY) & (POL_INFO)->mask].score_sum +\
         (POL_INFO)->values[((MATCHHISTORY) >> (POL_INFO)->cut_depth) &\
                             (POL_INFO)->mask].diff_from_max >= 0)

#endif
