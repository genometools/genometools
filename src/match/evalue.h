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

#include "core/unused_api.h"
#include "core/types_api.h"
#include "extended/scorehandler.h"
#include "match/karlin_altschul_stat.h"

/*
  calculates effective searchspace for query sequence of length
  <query_idx_length> and a set of database sequences.
  <total_length_of_db> is the total number of characters in all db sequences
  including separators and wildcards
 */
GtUword gt_evalue_searchspace(const GtKarlinAltschulStat *ka,
                              GtUword query_idx_length);

GtWord gt_evalue_raw_score(const GtKarlinAltschulStat *ka,
                           GtUword matches,
                           GtUword mismatches,
                           GtUword indels);

double gt_evalue_raw_score2bit_score(const GtKarlinAltschulStat *ka,
                                     GtWord raw_score);

GtWord gt_evalue_bit_score2raw_score(const GtKarlinAltschulStat *ka,
                                     double bit_score);

double gt_evalue_from_raw_score(const GtKarlinAltschulStat *ka,
                                GtWord raw_score,
                                GtUword searchspace);

double gt_evalue_from_bitscore(const GtKarlinAltschulStat *ka,
                               double bit_score,
                               GtUword searchspace);

double gt_evalue_from_eop_count(const GtKarlinAltschulStat *ka,
                                GtUword matches,
                                GtUword mismatches,
                                GtUword indels,
                                GtUword searchspace);

int gt_evalue_unit_test(GT_UNUSED GtError *err);

#endif
