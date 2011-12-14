/*
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include "core/disc_distri.h"
#include "core/format64.h"
#include "core/ma.h"
#include "extended/assembly_stats_calculator.h"

struct GtAssemblyStatsCalculator
{
  uint64_t numofseq;
  uint64_t sumlength;
  unsigned long minlength;
  unsigned long maxlength;
  GtDiscDistri *lengths;
};

GtAssemblyStatsCalculator *gt_assembly_stats_calculator_new(void)
{
  GtAssemblyStatsCalculator *asc;

  asc = gt_malloc(sizeof (GtAssemblyStatsCalculator));
  asc->lengths = gt_disc_distri_new();
  asc->numofseq = 0;
  asc->sumlength = 0;
  asc->minlength = 0;
  asc->maxlength = 0;
  return asc;
}

void gt_assembly_stats_calculator_delete(GtAssemblyStatsCalculator *asc)
{
  if (asc != NULL)
  {
    gt_disc_distri_delete(asc->lengths);
    gt_free(asc);
  }
}

void gt_assembly_stats_calculator_add(GtAssemblyStatsCalculator *asc,
    unsigned long length)
{
  gt_assert(asc != NULL);
  gt_assert(length != 0);
  gt_disc_distri_add(asc->lengths, length);
  (asc->numofseq)++;
  asc->sumlength += length;
  if (asc->minlength == 0 || length < asc->minlength)
    asc->minlength = length;
  if (length > asc->maxlength)
    asc->maxlength = length;
}

#define NOF_NSTATS 2 /* N50, N80 */
typedef struct
{
  unsigned long       nvalue[NOF_NSTATS];
  unsigned long long  min[NOF_NSTATS];
  bool                done[NOF_NSTATS];
  char                *name[NOF_NSTATS];
  unsigned long long  current;
} Nstats;

static void calcNstats(unsigned long key, unsigned long long value,
                        void *data)
{
  Nstats *nstats = data;
  int i;
  nstats->current += (key*value);
  for (i = 0; i < NOF_NSTATS; i++)
  {
    if (!nstats->done[i])
    {
      if (nstats->current >= nstats->min[i])
      {
        nstats->done[i] = true;
        nstats->nvalue[i] = key;
      }
    }
  }
}

#define initNstat(INDEX, NAME, PORTION)\
  nstats.name[INDEX] = (NAME);\
  nstats.min[INDEX] = (unsigned long long) (asc->sumlength * (PORTION) / 100);\
  nstats.nvalue[INDEX] = 0;\
  nstats.done[INDEX] = false

void gt_assembly_stats_calculator_show(GtAssemblyStatsCalculator *asc,
    GtLogger *logger)
{
  Nstats nstats;
  int i;

  gt_assert(asc != NULL);
  initNstat(0, "N50", 50);
  initNstat(1, "N80", 80);
  nstats.current = 0;

  gt_disc_distri_foreach_in_reverse_order(asc->lengths, calcNstats, &nstats);

  gt_logger_log(logger, "number of contigs: "Formatuint64_t"",
      PRINTuint64_tcast(asc->numofseq));
  gt_logger_log(logger, "total length:      "Formatuint64_t"",
      PRINTuint64_tcast(asc->sumlength));
  gt_logger_log(logger, "average size:      %.2f",
      (double) asc->sumlength/asc->numofseq);
  gt_logger_log(logger, "longest contig:    %lu", asc->maxlength);
  for (i = 0; i < NOF_NSTATS ; i++)
    gt_logger_log(logger, "%s:               %lu", nstats.name[i],
        nstats.nvalue[i]);
  gt_logger_log(logger, "smallest contig:   %lu", asc->minlength);
}
