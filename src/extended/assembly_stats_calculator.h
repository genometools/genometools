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

#ifndef ASSEMBLY_STATS_CALCULATOR_H
#define ASSEMBLY_STATS_CALCULATOR_H

#include "core/logger.h"

typedef struct GtAssemblyStatsCalculator GtAssemblyStatsCalculator;

GtAssemblyStatsCalculator *gt_assembly_stats_calculator_new(void);
void gt_assembly_stats_calculator_delete(GtAssemblyStatsCalculator *asc);

void gt_assembly_stats_calculator_add(GtAssemblyStatsCalculator *asc,
    unsigned long length);

/* if the following is set to a value > 0 (which is the default value)
 * then the NG50 and NG80 are also calculated */
void gt_assembly_stats_calculator_set_genome_length(
    GtAssemblyStatsCalculator *asc, unsigned long genome_length);

void gt_assembly_stats_calculator_show(GtAssemblyStatsCalculator *asc,
    GtLogger *logger);

#endif
