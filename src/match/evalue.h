/*
  Copyright (c) 2016 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
  Copyright (c) 2016 Center for Bioinformatics, University of Hamburg

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

#ifndef EVALUE_H
#define EVALUE_H

#include "core/types_api.h"
#include "extended/scorehandler.h"
#include "match/karlin_altschul_stat.h"
#include "querymatch.h"

/*
 * calculates evalue for an alignment
 * <ma> = number of matches
 * <mm> = number of mismatches
 * <id> = number of indels
 */
double gt_evalue_calculate(const GtKarlinAltschulStat *ka,
                           const GtScoreHandler *scorehandler,
                           GtUword ma,
                           GtUword mm,
                           GtUword id,
                           GtUword searchspace);

double gt_evalue_calculate_on_bitscore(const GtKarlinAltschulStat *ka,
                                       double bit_score,
                                       GtUword searchspace);

GtUword gt_evalue_calculate_searchspace(const GtKarlinAltschulStat *ka,
                                        const GtEncseq *dbencseq,
                                        GtUword query_idx_length);

/* use it only if the total number of sequences contained in <queryencseq>
 * is one, better use  gt_evalue_calculate_searchspace instead */
GtUword gt_evalue_calculate_searchspace_by_encseqs(
                                                 const GtKarlinAltschulStat *ka,
                                                 const GtEncseq *dbencseq,
                                                 const GtEncseq *queryencseq);

#endif
