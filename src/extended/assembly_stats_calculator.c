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

#include "core/disc_distri_api.h"
#include "core/format64.h"
#include "core/ma.h"
#include "extended/assembly_stats_calculator.h"

struct GtAssemblyStatsCalculator
{
  uint64_t numofseq;
  uint64_t sumlength;
  unsigned long minlength;
  unsigned long maxlength;
  unsigned long genome_length;
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
  asc->genome_length = 0;
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

void gt_assembly_stats_calculator_set_genome_length(
    GtAssemblyStatsCalculator *asc, unsigned long genome_length)
{
  gt_assert(asc != NULL);
  asc->genome_length = genome_length;
}

#define NOF_LIMITS 5

#define MAX_NOF_NSTATS 4
typedef struct
{
  unsigned long       nvalue[MAX_NOF_NSTATS];
  unsigned long       lvalue[MAX_NOF_NSTATS];
  unsigned long long  min[MAX_NOF_NSTATS];
  bool                done[MAX_NOF_NSTATS];
  char                *name[MAX_NOF_NSTATS];
  unsigned long long  limit[NOF_LIMITS];
  unsigned long long  larger_than_limit[NOF_LIMITS];
  unsigned long long  current_len;
  unsigned long long  current_num;
  unsigned int        nofstats;
  unsigned long long  median;
  unsigned long long  half_num;
} Nstats;

static void calcNstats(unsigned long key, unsigned long long value,
                        void *data)
{
  Nstats *nstats = data;
  unsigned int i;
  nstats->current_len += (key * value);
  nstats->current_num += value;
  for (i = 0; i < (unsigned int) NOF_LIMITS; i++)
  {
    if ((unsigned long long) key > nstats->limit[i])
      nstats->larger_than_limit[i] = nstats->current_num;
  }
  if (nstats->median == 0 && nstats->current_num >= nstats->half_num)
  {
    nstats->median = (unsigned long long) key;
  }
  for (i = 0; i < nstats->nofstats; i++)
  {
    if (!nstats->done[i] && (nstats->current_len >= nstats->min[i]))
    {
      nstats->done[i] = true;
      nstats->nvalue[i] = key;
      nstats->lvalue[i] = (unsigned long) nstats->current_num;
    }
  }
}

#define initNstat(INDEX, NAME, LENGTH)\
  nstats.name[INDEX] = (NAME);\
  nstats.min[INDEX] = (unsigned long long) (LENGTH);\
  nstats.nvalue[INDEX] = 0;\
  nstats.lvalue[INDEX] = 0;\
  nstats.done[INDEX] = false;\
  nstats.nofstats++

#define initLimit(INDEX, LENGTH)\
  nstats.limit[INDEX] = (unsigned long long) (LENGTH);\
  nstats.larger_than_limit[INDEX] = 0

void gt_assembly_stats_calculator_show(GtAssemblyStatsCalculator *asc,
    GtLogger *logger)
{
  Nstats nstats;
  unsigned int i;

  gt_assert(asc != NULL);
  nstats.nofstats = 0;
  initNstat(0, "50: ", asc->sumlength * 0.5);
  initNstat(1, "80: ", asc->sumlength * 0.8);
  if (asc->genome_length > 0)
  {
    initNstat(2, "G50:", asc->genome_length * 0.5);
    initNstat(3, "G80:", asc->genome_length * 0.8);
  }

  initLimit(0, 500ULL);
  initLimit(1, 1000ULL);
  initLimit(2, 10000ULL);
  initLimit(3, 100000ULL);
  initLimit(4, 1000000ULL);
  nstats.current_len = 0;
  nstats.current_num = 0;
  nstats.half_num = (unsigned long long) (asc->numofseq >> 1);
  nstats.median = 0;
  gt_disc_distri_foreach_in_reverse_order(asc->lengths, calcNstats, &nstats);

  gt_logger_log(logger, "number of contigs:     "Formatuint64_t"",
      PRINTuint64_tcast(asc->numofseq));
  if (asc->genome_length > 0)
  {
    gt_logger_log(logger, "genome length:         %lu",
        PRINTuint64_tcast(asc->genome_length));
  }
  gt_logger_log(logger, "total contigs length:  "Formatuint64_t"",
      PRINTuint64_tcast(asc->sumlength));
  if (asc->genome_length > 0)
  {
    gt_logger_log(logger, "   as %% of genome:     %.2f %%",
        (double) asc->sumlength * 100 / asc->genome_length);
  }
  gt_logger_log(logger, "mean contig size:      %.2f",
      (double) asc->sumlength / asc->numofseq);
  gt_logger_log(logger, "median contig size:    %llu", nstats.median);
  gt_logger_log(logger, "longest contig:        %lu", asc->maxlength);
  gt_logger_log(logger, "shortest contig:       %lu", asc->minlength);
  gt_logger_log(logger, "contigs > 500 nt:      %llu (%.2f %%)",
      nstats.larger_than_limit[0], (double) nstats.larger_than_limit[0] * 100
      / asc->numofseq);
  gt_logger_log(logger, "contigs > 1K nt:       %llu (%.2f %%)",
      nstats.larger_than_limit[1], (double) nstats.larger_than_limit[1] * 100
      / asc->numofseq);
  gt_logger_log(logger, "contigs > 10K nt:      %llu (%.2f %%)",
      nstats.larger_than_limit[2], (double) nstats.larger_than_limit[2] * 100
      / asc->numofseq);
  gt_logger_log(logger, "contigs > 100K nt:     %llu (%.2f %%)",
      nstats.larger_than_limit[3], (double) nstats.larger_than_limit[3] * 100
      / asc->numofseq);
  gt_logger_log(logger, "contigs > 1M nt:       %llu (%.2f %%)",
      nstats.larger_than_limit[4], (double) nstats.larger_than_limit[4] * 100
      / asc->numofseq);
  for (i = 0; i < nstats.nofstats ; i++)
  {
    if (nstats.nvalue[i] > 0)
    {
      gt_logger_log(logger, "N%s                  %lu", nstats.name[i],
          nstats.nvalue[i]);
      gt_logger_log(logger, "L%s                  %lu", nstats.name[i],
          nstats.lvalue[i]);
    }
    else
    {
      gt_logger_log(logger, "N%s                  n.a.",
          nstats.name[i]);
      gt_logger_log(logger, "L%s                  n.a.",
          nstats.name[i]);
    }
  }
}
