#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "core/assert_api.h"
#include "core/ma_api.h"
#include "core/types_api.h"
#include "ft-polish.h"

static void fill_polishing_info(Polishing_info *pol_info,
                                GtUword currentdepth,
                                GtUword prefix, GtWord score, GtWord maxscore)
{
  gt_assert(currentdepth <= pol_info->cut_depth);
  if (currentdepth == pol_info->cut_depth)
  {
    gt_assert(prefix < pol_info->entries && score >= INT16_MIN + maxscore);
    pol_info->values[prefix].diff_from_max = (int16_t) (score - maxscore);
    pol_info->values[prefix].score_sum = (int16_t) score;
  } else
  {
    if (score > maxscore)
    {
      maxscore = score;
    }
    gt_assert(score >= INT16_MIN - pol_info->difference_score);
    fill_polishing_info(pol_info,currentdepth+1, prefix << 1,
                        score - pol_info->difference_score,maxscore);
    gt_assert(score <= INT16_MAX - pol_info->match_score);
    fill_polishing_info(pol_info,currentdepth+1,(prefix << 1) | 1UL,
                        score + pol_info->match_score,maxscore);
  }
}

Polishing_info *polishing_info_new(GtUword cut_depth,
                                   double errorpercentage)
{
  Polishing_info *pol_info = gt_malloc(sizeof *pol_info);

  gt_assert(pol_info != NULL);
  pol_info->entries = 1UL << cut_depth;
  pol_info->mask = pol_info->entries - 1;
  pol_info->values = gt_malloc(sizeof *pol_info->values * pol_info->entries);
  gt_assert(pol_info->values != NULL);
  pol_info->cut_depth = cut_depth;
  pol_info->match_score = 20.0 * errorpercentage;
  gt_assert(pol_info->match_score <= 1000.0);
  pol_info->difference_score = 1000.0 - pol_info->match_score;
  fill_polishing_info(pol_info,0,0,0,0);
  return pol_info;
}

uint64_t polishing_info_maxvalue(const Polishing_info *pol_info)
{
  return (((uint64_t) 1) << (2 * pol_info->cut_depth)) - 1;
}

static char *intbits2string(GtUword bits,GtUword bs)
{
  char *csptr;
  GtUword mask;
  static char cs[64+1];

  gt_assert(bits > 0);
  for (csptr = cs, mask = 1 << (bits-1); mask > 0; mask >>= 1, csptr++)
  {
     *csptr = (bs & mask) ? '1' : '0';
  }
  *csptr = '\0';
  return cs;
}

void polishing_info_show(const Polishing_info *pol_info)
{
  GtUword idx;

  printf("pi->cut_depth=" GT_WU "\n",pol_info->cut_depth);
  printf("pi->entries=" GT_WU "\n",pol_info->entries);
  printf("pi->match_score=" GT_WD "\n",pol_info->match_score);
  printf("pi->difference_score=" GT_WD "\n",pol_info->difference_score);
  printf("pi->mask=%s\n",intbits2string(pol_info->cut_depth,pol_info->mask));
  for (idx = 0; idx < pol_info->entries; idx++)
  {
    printf("[%s]=%+hd/%+hd\n",intbits2string(pol_info->cut_depth,idx),
                              pol_info->values[idx].score_sum,
                              pol_info->values[idx].diff_from_max);
  }
}

void polishing_info_delete(Polishing_info *pol_info)
{
  if (pol_info != NULL)
  {
    gt_free(pol_info->values);
    gt_free(pol_info);
  }
}

bool history_is_polished_brute_force(const Polishing_info *pol_info,
                                     uint64_t matchhistory)
{
  GtUword idx;
  uint64_t mask;
  GtWord sum_score = 0;

  for (mask = (uint64_t) 1, idx = 0; idx < 2 * pol_info->cut_depth;
       idx++, mask <<= 1)
  {
    if (matchhistory & mask)
    {
      sum_score += pol_info->match_score;
    } else
    {
      sum_score -= pol_info->difference_score;
    }
    if (sum_score < 0)
    {
      return false;
    }
  }
  return true;
}
