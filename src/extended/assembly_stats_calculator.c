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

#include "core/compat.h"
#include "core/disc_distri_api.h"
#include "core/format64.h"
#include "core/ma.h"
#include "extended/assembly_stats_calculator.h"

struct GtAssemblyStatsCalculator
{
  uint64_t numofseq;
  uint64_t sumlength;
  GtUword minlength;
  GtUword maxlength;
  GtUword genome_length;
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
    GtUword length)
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
    GtAssemblyStatsCalculator *asc, GtUword genome_length)
{
  gt_assert(asc != NULL);
  asc->genome_length = genome_length;
}

#define NOF_LIMITS 5

#define MAX_NOF_NSTATS 4
typedef struct
{
  char         *name[MAX_NOF_NSTATS];
  bool         done[MAX_NOF_NSTATS];
  GtUint64     current_len,
               current_num,
               first_quartile,
               fourth_num,
               half_num,
               larger_than_limit[NOF_LIMITS],
               limit[NOF_LIMITS],
               median,
               min[MAX_NOF_NSTATS],
               third_quartile,
               three_fourth_num;
  GtUword      nvalue[MAX_NOF_NSTATS],
               lvalue[MAX_NOF_NSTATS];
  unsigned int nofstats;
} Nstats;

static void calcNstats(GtUword key, GtUint64 value,
                        void *data)
{
  Nstats *nstats = data;
  unsigned int i;
  nstats->current_len += (key * value);
  nstats->current_num += value;
  for (i = 0; i < (unsigned int) NOF_LIMITS; i++)
  {
    if ((GtUint64) key > nstats->limit[i])
      nstats->larger_than_limit[i] = nstats->current_num;
  }
  if (nstats->third_quartile == 0 && nstats->current_num >= nstats->fourth_num)
  {
    nstats->third_quartile = (GtUint64) key;
  }
  if (nstats->median == 0 && nstats->current_num >= nstats->half_num)
  {
    nstats->median = (GtUint64) key;
  }
  if (nstats->first_quartile == 0 && nstats->current_num >=
      nstats->three_fourth_num)
  {
    nstats->first_quartile = (GtUint64) key;
  }
  for (i = 0; i < nstats->nofstats; i++)
  {
    if (!nstats->done[i] && (nstats->current_len >= nstats->min[i]))
    {
      nstats->done[i] = true;
      nstats->nvalue[i] = key;
      nstats->lvalue[i] = (GtUword) nstats->current_num;
    }
  }
}

#define initNstat(INDEX, NAME, LENGTH)\
  nstats.name[INDEX] = (NAME);\
  nstats.min[INDEX] = (GtUint64) (LENGTH);\
  nstats.nvalue[INDEX] = 0;\
  nstats.lvalue[INDEX] = 0;\
  nstats.done[INDEX] = false;\
  nstats.nofstats++

#define initLimit(INDEX, LENGTH)\
  nstats.limit[INDEX] = (GtUint64) (LENGTH);\
  nstats.larger_than_limit[INDEX] = 0

GtUword gt_assembly_stats_calculator_nstat(GtAssemblyStatsCalculator *asc,
    GtUword n)
{
  Nstats nstats;
  gt_assert(n > 0);
  gt_assert(n < (GtUword)100UL);
  nstats.min[0] = (GtUint64) (asc->sumlength * ((float)n / 100U));
  nstats.nvalue[0] = 0;
  nstats.lvalue[0] = 0;
  nstats.done[0] = false;
  nstats.nofstats = 1U;
  nstats.current_len = 0;
  nstats.current_num = 0;
  nstats.half_num = (GtUint64) (asc->numofseq >> 1);
  nstats.fourth_num = nstats.half_num >> 1;
  nstats.three_fourth_num = nstats.fourth_num + nstats.half_num;
  nstats.median = 0;
  nstats.first_quartile = 0;
  nstats.third_quartile = 0;
  gt_disc_distri_foreach_in_reverse_order(asc->lengths, calcNstats, &nstats);
  return nstats.nvalue[0];
}

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
  nstats.half_num = (GtUint64) (asc->numofseq >> 1);
  nstats.fourth_num = nstats.half_num >> 1;
  nstats.three_fourth_num = nstats.fourth_num + nstats.half_num;
  nstats.median = 0;
  nstats.first_quartile = 0;
  nstats.third_quartile = 0;
  gt_disc_distri_foreach_in_reverse_order(asc->lengths, calcNstats, &nstats);

  gt_logger_log(logger, "number of contigs:     "Formatuint64_t"",
                PRINTuint64_tcast(asc->numofseq));
  if (asc->genome_length > 0)
  {
    gt_logger_log(logger, "genome length:         "GT_WU"",
                  PRINTuint64_tcast(asc->genome_length));
  }
  gt_logger_log(logger, "total contigs length:  "Formatuint64_t"",
                PRINTuint64_tcast(asc->sumlength));
  if (asc->genome_length > 0)
  {
    gt_logger_log(logger, "   as %% of genome:     %.2f %%", (double)
                  asc->sumlength * 100 / asc->genome_length);
  }
  gt_logger_log(logger, "mean contig size:      %.2f",
                (double) asc->sumlength / asc->numofseq);
  gt_logger_log(logger, "contig size first quartile: " GT_LLU,
                nstats.first_quartile);
  gt_logger_log(logger, "median contig size:         " GT_LLU, nstats.median);
  gt_logger_log(logger, "contig size third quartile: " GT_LLU,
                nstats.third_quartile);
  gt_logger_log(logger, "longest contig:             " GT_WU, asc->maxlength);
  gt_logger_log(logger, "shortest contig:            " GT_WU, asc->minlength);
  gt_logger_log(logger, "contigs > 500 nt:           " GT_LLU " (%.2f %%)",
                nstats.larger_than_limit[0], (double)
                nstats.larger_than_limit[0] * 100 / asc->numofseq);
  gt_logger_log(logger, "contigs > 1K nt:            " GT_LLU " (%.2f %%)",
                nstats.larger_than_limit[1], (double)
                nstats.larger_than_limit[1] * 100 / asc->numofseq);
  gt_logger_log(logger, "contigs > 10K nt:           " GT_LLU " (%.2f %%)",
                nstats.larger_than_limit[2], (double)
                nstats.larger_than_limit[2] * 100 / asc->numofseq);
  gt_logger_log(logger, "contigs > 100K nt:          " GT_LLU " (%.2f %%)",
                nstats.larger_than_limit[3], (double)
                nstats.larger_than_limit[3] * 100 / asc->numofseq);
  gt_logger_log(logger, "contigs > 1M nt:            " GT_LLU " (%.2f %%)",
                nstats.larger_than_limit[4], (double)
                nstats.larger_than_limit[4] * 100 / asc->numofseq);
  for (i = 0; i < nstats.nofstats ; i++)
  {
    if (nstats.nvalue[i] > 0)
    {
      gt_logger_log(logger, "N%s                  "GT_WU"", nstats.name[i],
          nstats.nvalue[i]);
      gt_logger_log(logger, "L%s                  "GT_WU"", nstats.name[i],
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
