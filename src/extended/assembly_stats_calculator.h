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
#include "core/types_api.h"

/* GtAssemblyStatsCalculator computes statistics for sequence sets,
   which are in particular useful for assemblies (contigs sets) */
typedef struct GtAssemblyStatsCalculator GtAssemblyStatsCalculator;

/* Create a new GtAssemblyStatsCalculator object */
GtAssemblyStatsCalculator *gt_assembly_stats_calculator_new(void);

/* Delete the GtAssemblyStatsCalculator <asc> */
void                       gt_assembly_stats_calculator_delete(
                                                GtAssemblyStatsCalculator *asc);

/* Add to the GtAssemblyStatsCalculator <asc> the sequence length <length>;
   the calculator (only) needs to know about the length of each sequence in the
   set in order to compute the statistics*/
void                       gt_assembly_stats_calculator_add(
                                                GtAssemblyStatsCalculator *asc,
                                                GtUword length);

/* Compute the N statistics <n> for the GtAssemblyStatsCalculator <asc>;
   <n> is an integer between 0 and 100 (extremes excluded);
   e.g. for N50 use n = 50 */
GtUword                    gt_assembly_stats_calculator_nstat(
                                                 GtAssemblyStatsCalculator *asc,
                                                 GtUword n);

/* Set the genome length for the GtAssemblyStatsCalculator <asc>;
   this is not necessary, but if the genome length is set;
   then the NG/LG statistics can also be computed;
   to disable NG statistics computation set the genome lenght to 0
   (this is the default) */
void                       gt_assembly_stats_calculator_set_genome_length(
                                                 GtAssemblyStatsCalculator *asc,
                                                 GtUword genome_length);

/* Print the statistics for GtAssemblyStatsCalculator <asc>
   to the GtLogger <logger> in a pre-formatted table */
void                       gt_assembly_stats_calculator_show(
                                                 GtAssemblyStatsCalculator *asc,
                                                 GtLogger *logger);

#endif
