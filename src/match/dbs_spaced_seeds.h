/*
  Copyright (c) 2017 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2017 Center for Bioinformatics, University of Hamburg

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

#ifndef DBS_SPACED_SEEDS_H
#define DBS_SPACED_SEEDS_H
#include "core/codetype.h"

#define GT_SPACED_SEED_FIRST_SPAN 15

int gt_spaced_seed_span(GtCodetype spaced_seed);

int gt_spaced_seed_weight(GtCodetype spaced_seed);

void gt_spaced_seed_weight_range(int *min_weight,int *max_weight, int span);

typedef struct GtSpacedSeedSpec GtSpacedSeedSpec;

GtSpacedSeedSpec *gt_spaced_seed_spec_new(GtCodetype spacedseed);

GtSpacedSeedSpec *gt_spaced_seed_spec_new_from_ws(int weight,int span);

void gt_spaced_seed_spec_delete(GtSpacedSeedSpec *seed_spec);

GtCodetype gt_spaced_seed_extract_generic(const GtSpacedSeedSpec *seed_spec,
                                          GtCodetype kmer);

#endif
