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

#include "core/alphabet.h"
#include "extended/scorehandler.h"
typedef struct GtKarlinAltschulStat GtKarlinAltschulStat;

GtKarlinAltschulStat *gt_karlin_altschul_stat_new(void);
void gt_karlin_altschul_stat_delete(GtKarlinAltschulStat *ka);
double gt_karlin_altschul_stat_get_lambda(const GtKarlinAltschulStat *ka);
double gt_karlin_altschul_stat_get_K(const GtKarlinAltschulStat *ka);
double gt_karlin_altschul_stat_get_logK(const GtKarlinAltschulStat *ka);
double gt_karlin_altschul_stat_get_alpha_div_lambda(const GtKarlinAltschulStat *ka);
double gt_karlin_altschul_stat_get_beta(const GtKarlinAltschulStat *ka);

/* determine karlin altschul parameters lambda, H, K, alpha and beta for 
 * a given <alphabet> and scorefunction given as <scorehandler>
 * 
 * If <gapped_alignments> is not set, paraemters are calculated by using scoring 
 * frequency statistics for ungapped alignments,
 * else precomputed values are used for gapped alignment
 * 
 * returns 0 if no error occured, otherwise returns 1
 * analog to BLAST
 */
int gt_karlin_altschul_stat_calculate_params(GtKarlinAltschulStat *ka,
                                             bool gapped_alignment,
                                             GtAlphabet *alphabet,
                                             GtScoreHandler *scorehandler,
                                             GtError *err);

#endif
