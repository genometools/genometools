#include "core/unused_api.h"
#include "core/assert_api.h"
#include "core/arraydef_api.h"
#include "core/minmax_api.h"
#include "match/ft-trimstat.h"

struct GtFtTrimstat
{
  GtUword diedout, *trim_dist, *matchlength_dist;
  GtArrayGtUword distance_dist, maxvalid_dist;
  size_t spaceforfront_total;
  double sum_meanvalid;
};

GtFtTrimstat *gt_ft_trimstat_new(void)
{
  GtFtTrimstat *trimstat = gt_malloc(sizeof *trimstat);

  gt_assert(trimstat != 0);
  trimstat->trim_dist = gt_calloc(101,sizeof *trimstat->trim_dist);
  gt_assert(trimstat->trim_dist != NULL);
  trimstat->matchlength_dist = gt_calloc(101,
                                         sizeof *trimstat->matchlength_dist);
  gt_assert(trimstat->matchlength_dist != NULL);
  trimstat->diedout = 0;
  GT_INITARRAY(&trimstat->distance_dist,GtUword);
  GT_INITARRAY(&trimstat->maxvalid_dist,GtUword);
  trimstat->spaceforfront_total = 0;
  trimstat->sum_meanvalid = 0.0;
  return trimstat;
}

#ifndef NDEBUG
void gt_ft_trimstat_add(GtFtTrimstat *trimstat,
                        bool diedout,
                        GtUword sumvalid,
                        GtUword maxvalid,
                        GtUword d,
                        size_t spaceforfront)
{
  if (trimstat == NULL)
  {
    return;
  }
  while (maxvalid >= trimstat->maxvalid_dist.allocatedGtUword)
  {
    GtUword idx;
    const GtUword allocated = trimstat->maxvalid_dist.allocatedGtUword;

    trimstat->maxvalid_dist.allocatedGtUword
      = trimstat->maxvalid_dist.allocatedGtUword * 1.2 + 128UL;
    trimstat->maxvalid_dist.spaceGtUword
      = gt_realloc(trimstat->maxvalid_dist.spaceGtUword,
                   sizeof *trimstat->maxvalid_dist.spaceGtUword *
                   trimstat->maxvalid_dist.allocatedGtUword);
    for (idx = allocated; idx < trimstat->maxvalid_dist.allocatedGtUword; idx++)
    {
      trimstat->maxvalid_dist.spaceGtUword[idx] = 0;
    }
  }
  gt_assert(maxvalid < trimstat->maxvalid_dist.allocatedGtUword);
  trimstat->maxvalid_dist.spaceGtUword[maxvalid]++;
  if (diedout)
  {
    trimstat->diedout++;
  } else
  {
    const GtUword fullfronts = (d+1) * (d+1);
    GtUword percentage;

    percentage
      = (GtUword) (100.0 * (double) (fullfronts - sumvalid)/fullfronts);
    gt_assert(percentage <= 100UL);
    trimstat->sum_meanvalid += (double) sumvalid/(d+1);
    trimstat->trim_dist[percentage]++;
    GT_CHECKARRAYSPACE(&trimstat->distance_dist,GtUword,32);
    trimstat->distance_dist.spaceGtUword[trimstat->
                                         distance_dist.nextfreeGtUword++] = d;
    trimstat->spaceforfront_total += spaceforfront;
  }
}

void gt_ft_trimstat_add_matchlength(GtFtTrimstat *trimstat,
                                    uint32_t matchlength)
{
  gt_assert(trimstat != NULL && trimstat->matchlength_dist != NULL);
  trimstat->matchlength_dist[GT_MIN(100,matchlength)]++;
}
#else
void gt_ft_trimstat_add(GT_UNUSED GtFtTrimstat *trimstat,
                        GT_UNUSED bool diedout,
                        GT_UNUSED GtUword sumvalid,
                        GT_UNUSED GtUword maxvalid,
                        GT_UNUSED GtUword d,
                        GT_UNUSED size_t spaceforfront)
{
  return;
}

void gt_ft_trimstat_add_matchlength(GT_UNUSED GtFtTrimstat *trimstat,
                                    GT_UNUSED uint32_t matchlength)
{
  return;
}
#endif

static int gt_ft_trimstat_compare_GtUword(const void *va, const void *vb)
{
  GtUword a = *((const GtUword *) va);
  GtUword b = *((const GtUword *) vb);

  if (a < b)
  {
    return -1;
  }
  if (a > b)
  {
    return 1;
  }
  return 0;
}

#define MEGABYTES(X) ((double) (X)/(1UL << 20))

void gt_ft_trimstat_delete(GtFtTrimstat *trimstat)
{
  if (trimstat != NULL)
  {
    gt_free(trimstat->trim_dist);
    GT_FREEARRAY(&trimstat->distance_dist,GtUword);
    GT_FREEARRAY(&trimstat->maxvalid_dist,GtUword);
    gt_free(trimstat->matchlength_dist);
    gt_free(trimstat);
  }
}

void gt_ft_trimstat_out(const GtFtTrimstat *trimstat,bool verbose)
{
  if (trimstat != NULL)
  {
    printf("died_out=" GT_WU "\t",trimstat->diedout);
    if (trimstat->distance_dist.nextfreeGtUword > 0)
    {
      printf("mean_valid=%.2f\t",
             trimstat->sum_meanvalid/trimstat->distance_dist.nextfreeGtUword);
      printf("frontspace=%.2f\n",
             MEGABYTES((double) trimstat->spaceforfront_total/
                       trimstat->distance_dist.nextfreeGtUword));
    } else
    {
      printf("mean_valid=undef\t");
      printf("frontspace=undef\n");
    }
    if (verbose)
    {
      GtUword idx, count = 1UL, matchlength_sum = 0, matchlength_cum = 0;

      for (idx = 0; idx <= 100UL; idx++)
      {
        matchlength_sum += trimstat->matchlength_dist[idx];
      }
      for (idx = 0; idx <= 100UL; idx++)
      {
        if (trimstat->matchlength_dist[idx] > 0)
        {
          matchlength_cum += trimstat->matchlength_dist[idx];
          printf("# matchlength%s" GT_WU ": " GT_WU " times, "
                 "total=" GT_WU " (%.2f), "
                 "cum=%.2f%%\n",
                 idx < 100UL ? "=" : ">=",
                 idx,trimstat->matchlength_dist[idx],
                 idx * trimstat->matchlength_dist[idx],
                 (double) trimstat->matchlength_dist[idx]/matchlength_sum,
                 100.0 * (double) matchlength_cum/matchlength_sum);
        }
      }
      for (idx = 0; idx <= 100UL; idx++)
      {
        if (trimstat->trim_dist[idx] > 0)
        {
          printf("# trim by " GT_WU "%%: " GT_WU " times\n",
                 idx,trimstat->trim_dist[idx]);
        }
      }
      qsort(trimstat->distance_dist.spaceGtUword,
            trimstat->distance_dist.nextfreeGtUword,
            sizeof *trimstat->distance_dist.spaceGtUword,
            gt_ft_trimstat_compare_GtUword);
      if (trimstat->distance_dist.nextfreeGtUword > 0)
      {
        GtUword previous = trimstat->distance_dist.spaceGtUword[0];
        for (idx = 1UL; idx < trimstat->distance_dist.nextfreeGtUword; idx++)
        {
          if (previous == trimstat->distance_dist.spaceGtUword[idx])
          {
            count++;
          } else
          {
            printf("distance " GT_WU ": " GT_WU " times\n",previous,count);
            count = 1UL;
            previous = trimstat->distance_dist.spaceGtUword[idx];
          }
        }
        printf("distance " GT_WU ": " GT_WU " times\n",previous,count);
      }
      for (idx = 0; idx < trimstat->maxvalid_dist.allocatedGtUword; idx++)
      {
        if (trimstat->maxvalid_dist.spaceGtUword[idx] > 0)
        {
          printf("maxvalid=" GT_WU ": " GT_WU " times\n",idx,
                 trimstat->maxvalid_dist.spaceGtUword[idx]);
        }
      }
    }
  }
}
