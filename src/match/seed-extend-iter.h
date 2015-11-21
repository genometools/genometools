/*
  Copyright (c) 2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#ifndef SEED_EXTEND_ITER_H
#define SEED_EXTEND_ITER_H

#include "core/str_api.h"
#include "core/types_api.h"
#include "core/error_api.h"

typedef struct GtSeedextendMatchIterator GtSeedextendMatchIterator;

void gt_seedextend_match_iterator_delete(GtSeedextendMatchIterator *semi);

GtSeedextendMatchIterator *gt_seedextend_match_iterator_new(
                                            const GtStr *matchfilename,
                                            GtError *err);

int gt_seedextend_match_iterator_next(GtSeedextendMatchIterator *semi,
                                      GtError *err);

const GtEncseq *gt_seedextend_match_iterator_aencseq(
                        const GtSeedextendMatchIterator *semi);

const GtEncseq *gt_seedextend_match_iterator_bencseq(
                        const GtSeedextendMatchIterator *semi);

GtUword gt_seedextend_match_iterator_history_size(
                        const GtSeedextendMatchIterator *semi);

GtUword gt_seedextend_match_iterator_errorpercentage(
                        const GtSeedextendMatchIterator *semi);

bool gt_seedextend_match_iterator_bias_parameters(
                        const GtSeedextendMatchIterator *semi);

bool gt_seedextend_match_iterator_has_seedline(
                        const GtSeedextendMatchIterator *semi);

GtQuerymatch *gt_seedextend_match_iterator_querymatch_ptr(
                        GtSeedextendMatchIterator *semi);

GtUword gt_seedextend_match_iterator_seedlen(
                             const GtSeedextendMatchIterator *semi);

GtUword gt_seedextend_match_iterator_seedpos1(
                             const GtSeedextendMatchIterator *semi);

GtUword gt_seedextend_match_iterator_seedpos2(
                            const GtSeedextendMatchIterator *semi);

#endif
