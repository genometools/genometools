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

#ifndef KARLIN_ALTSCHUL_STAT_H
#define KARLIN_ALTSCHUL_STAT_H

#include <stdbool.h>
#include "extended/scorehandler.h"
#include "core/error_api.h"
typedef struct GtKarlinAltschulStat GtKarlinAltschulStat;

/*
  determine karlin altschul parameters lambda, H, K, alpha and beta for
  scorefunction given as <scorehandler>

  If <numofchars> is 0, then the precomputation is done for
  gapped alignments. Otherwise, the precomputation is
  based on an alphabet of size <numchars> and calculated by using scoring
  frequency statistics for ungapped alignments,

  returns a GtKarlinAltschulStat object if no error occured, otherwise returns
  NULL and <err> is set
 */
GtKarlinAltschulStat *gt_karlin_altschul_stat_new(unsigned int numchars,
                                           const GtScoreHandler *scorehandler);

GtKarlinAltschulStat *gt_karlin_altschul_stat_new_gapped(void);

void gt_karlin_altschul_stat_add_keyvalues(
                             GtKarlinAltschulStat *karlin_altschul_stat,
                             GtUword db_totallength,
                             GtUword db_numofsequences);

void gt_karlin_altschul_stat_delete(GtKarlinAltschulStat *ka);
double gt_karlin_altschul_stat_get_lambda(const GtKarlinAltschulStat *ka);
double gt_karlin_altschul_stat_get_K(const GtKarlinAltschulStat *ka);
double gt_karlin_altschul_stat_get_logK(const GtKarlinAltschulStat *ka);
double gt_karlin_altschul_stat_get_alphadlambda(const GtKarlinAltschulStat *ka);
double gt_karlin_altschul_stat_get_beta(const GtKarlinAltschulStat *ka);
GtWord gt_karlin_altschul_get_num_of_db_seqs(const GtKarlinAltschulStat *ka);
GtWord gt_karlin_altschul_get_total_length_db(const GtKarlinAltschulStat *ka);
GtWord gt_karlin_altschul_stat_mismatchscore(const GtKarlinAltschulStat *ka);
GtWord gt_karlin_altschul_stat_matchscore(const GtKarlinAltschulStat *ka);
GtWord gt_karlin_altschul_stat_gapscore(const GtKarlinAltschulStat *ka);
int gt_karlin_altschul_stat_unit_test(GtError *err);

#endif
