/*
  Copyright (c) 2017 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
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

#ifndef WEIGHTED_LIS_FILTER_H
#define WEIGHTED_LIS_FILTER_H

#include <stdbool.h>
#include "core/types_api.h"
#include "core/arraydef.h"

/*
  This module implements a filter algorithm to leave only those alignments,
  which forms the longest mutually consistent set. The algorithm computes
  the identity^2 * len weighted longest increasing subset (LIS) analogous to
  the delta filter of MUMMER Software.

  @article{kurtz2004versatile,
  title={Versatile and open software for comparing large genomes},
  author={Kurtz, Stefan and Phillippy, Adam and Delcher, Arthur L and
         Smoot, Michael and Shumway, Martin and Antonescu,Corina and
         Salzberg, Steven L},
  journal={Genome biology},
  volume={5},
  number={2},
  pages={R12},
  year={2004},
  publisher={BioMed Central}
}*/

/* the following type is used for the table of matches */
typedef struct GtWLisFilterMatches GtWLisFilterMatches;

/* the constructor for the table of matches; */
GtWLisFilterMatches *gt_wlis_filter_matches_new(void);

/* function to reset the table of matches, to reuse it. */
void          gt_wlis_filter_matches_reset(GtWLisFilterMatches *allmatches);

/* the destructor for the table of matches */
void          gt_wlis_filter_matches_delete(GtWLisFilterMatches *allmatches);

/* function to add a match, described by its coordinates <s_start>, <s_end>,
   <q_start>,<q_end> and <distance> value, to the table of matches */
void           gt_wlis_filter_matches_add(GtWLisFilterMatches *allmatches,
                                          GtUword s_start, GtUword s_end,
                                          GtUword q_start, GtUword q_end,
                                          GtUword distance,
                                          bool store_querymatch);

/* function that applies the filter-algorithm to a given table of matches
   <allmatches> and stores the result in <result>. If <forward> is set, the
   algorithm interprets the q-coordinates regarding to the forward strand,
   else regarding to the reverse strand */
void gt_wlis_filter_evaluate(GtArrayGtUword *result,
                             GtUword *sum_distance,
                             GtUword *sum_aligned_len_chain,
                             GtUword *chain_weighted_score,
                             GtWLisFilterMatches *allmatches);

#endif
