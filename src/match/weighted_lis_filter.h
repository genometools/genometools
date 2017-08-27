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

/* the following type is used to store a chain of indices, whose corresponding
   alignments have passed the filter */
typedef struct
{
  GtUword *chain,
           nextfree,
           allocated;
}GtFilter;

/* the constructor for the filter */
GtFilter      *gt_filter_new();

/* function to reset the filter, to reuse it */
void          gt_filter_reset(GtFilter *filter);

/* the destructor for the filter */
void          gt_filter_delete(GtFilter *filter);

/* the following type is used for the table of matches */
typedef struct GtAllMatches GtAllMatches;

/* the constructor for the table of matches, setting the orientation to
   reverse if <reverse> is set */
GtAllMatches *gt_filter_all_matches_new(bool reverse);

/* function to reset the table of matches, to reuse it, setting the orientation
   for the next use to be reverse if <reverse> is set */
void          gt_filter_all_matches_reset(GtAllMatches *allmatches,
                                          bool reverse);
/* the destructor for the filter */
void          gt_filter_all_matches_delete(GtAllMatches *allmatches);

/* function to add matches to the table of matches */
void           gt_filter_all_matches_add(GtAllMatches *allmatches,
                                         GtUword s_start, GtUword s_end,
                                         GtUword q_start, GtUword q_end,
                                         float weight);

/* function that applies the filter-algorithm to a given table of matches
   <allmatches> and stores the result in <filter> */
void          gt_filter_call(GtAllMatches *allmatches,
                             GtFilter *filter);

#endif
