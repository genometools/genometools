#include "core/ma_api.h"
#include "core/types_api.h"
#include "core/assert_api.h"
#include "ft-trimstat.h"

struct Trimstat
{
  GtUword diedout, dist_nextfree, dist_allocated, allocated_maxvalid,
          *dist_maxvalid, *trimdist, *distarray;
  size_t spaceforfront_total;
  double sum_meanvalid;
  /* the following are to create the output at successive statements. */
  double errorpercentage;
  GtUword minmatchpercentage;
  GtUword maxalignedlendifference;
  GtUword max_cache_size;
};

Trimstat *trimstat_new(double errorpercentage,
                       GtUword minmatchpercentage,
                       GtUword maxalignedlendifference)
{
  Trimstat *trimstat = gt_malloc(sizeof *trimstat);

  gt_assert(trimstat != 0);
  trimstat->trimdist = gt_calloc(101,sizeof *trimstat->trimdist);
  gt_assert(trimstat->trimdist != NULL);
  trimstat->diedout = 0;
  trimstat->distarray = NULL;
  trimstat->dist_nextfree = 0;
  trimstat->dist_allocated = 0;
  trimstat->allocated_maxvalid = 0;
  trimstat->dist_maxvalid = NULL;
  trimstat->spaceforfront_total = 0;
  trimstat->sum_meanvalid = 0.0;
  trimstat->max_cache_size = 0;
  trimstat->errorpercentage = errorpercentage;
  trimstat->minmatchpercentage = minmatchpercentage;
  trimstat->maxalignedlendifference = maxalignedlendifference;
  return trimstat;
}

void trimstat_add(Trimstat *trimstat,bool diedout,
                  GtUword sumvalid,
                  GtUword maxvalid,
                  GtUword d,
                  size_t spaceforfront,
                  GtUword cache_size)
{
  if (trimstat == NULL)
  {
    return;
  }
  while (maxvalid >= trimstat->allocated_maxvalid)
  {
    GtUword idx;
    const GtUword allocated = trimstat->allocated_maxvalid;

    trimstat->allocated_maxvalid = trimstat->allocated_maxvalid * 1.2 + 128UL;
    trimstat->dist_maxvalid = gt_realloc(trimstat->dist_maxvalid,
                                         sizeof *trimstat->dist_maxvalid *
                                         trimstat->allocated_maxvalid);
    for (idx = allocated; idx < trimstat->allocated_maxvalid; idx++)
    {
      trimstat->dist_maxvalid[idx] = 0;
    }
  }
  gt_assert(maxvalid < trimstat->allocated_maxvalid);
  trimstat->dist_maxvalid[maxvalid]++;
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
    trimstat->trimdist[percentage]++;
    if (trimstat->dist_nextfree >= trimstat->dist_allocated)
    {
      trimstat->dist_allocated = trimstat->dist_allocated +
                                 trimstat->dist_allocated/5 + 1024;
      trimstat->distarray
        = gt_realloc(trimstat->distarray,
                     sizeof *trimstat->distarray * trimstat->dist_allocated);
      gt_assert(trimstat->distarray != NULL);
    }
    trimstat->distarray[trimstat->dist_nextfree++] = d;
    trimstat->spaceforfront_total += spaceforfront;
  }
  if (trimstat->max_cache_size < cache_size)
  {
    trimstat->max_cache_size = cache_size;
  }
}

static int compare_ulong(const void *va, const void *vb)
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

void trimstat_delete(Trimstat *trimstat,double total_time,bool verbose)
{
  if (trimstat != NULL)
  {
    printf("erp=%.1f\t",trimstat->errorpercentage);
    printf("mmp=" GT_WU "\t",trimstat->minmatchpercentage);
    printf("mad=" GT_WU "\t",trimstat->maxalignedlendifference);
    printf("died_out=" GT_WU "\t",trimstat->diedout);
    if (trimstat->dist_nextfree > 0)
    {
      printf("mean_valid=%.2f\t",
             trimstat->sum_meanvalid/trimstat->dist_nextfree);
      printf("frontspace=%.2f\t",
             MEGABYTES((double) trimstat->spaceforfront_total/
                       trimstat->dist_nextfree));
    } else
    {
      printf("mean_valid=undef\t");
      printf("frontspace=undef\t");
    }
    printf("time=%.2f\n",total_time);
    if (verbose)
    {
      GtUword idx, count = 1UL;

      printf("max_cache_size = " GT_WU " bytes\n",trimstat->max_cache_size);
      for (idx = 0; idx <= 100UL; idx++)
      {
        if (trimstat->trimdist[idx] > 0)
        {
          printf("# trim by " GT_WU "%%: " GT_WU " times\n",
                 idx,trimstat->trimdist[idx]);
        }
      }
      qsort(trimstat->distarray,trimstat->dist_nextfree,
            sizeof *trimstat->distarray,
            compare_ulong);
      if (trimstat->dist_nextfree > 0)
      {
        GtUword previous = trimstat->distarray[0];
        for (idx = 1UL; idx < trimstat->dist_nextfree; idx++)
        {
          if (previous == trimstat->distarray[idx])
          {
            count++;
          } else
          {
            printf("distance " GT_WU ": " GT_WU " times\n",previous,count);
            count = 1UL;
            previous = trimstat->distarray[idx];
          }
        }
        printf("distance " GT_WU ": " GT_WU " times\n",previous,count);
      }
      for (idx = 0; idx < trimstat->allocated_maxvalid; idx++)
      {
        if (trimstat->dist_maxvalid[idx] > 0)
        {
          printf("maxvalid=" GT_WU ": " GT_WU " times\n",idx,
                 trimstat->dist_maxvalid[idx]);
        }
      }
    }
    gt_free(trimstat->trimdist);
    gt_free(trimstat->distarray);
    gt_free(trimstat->dist_maxvalid);
    gt_free(trimstat);
  }
}
