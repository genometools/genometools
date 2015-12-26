/*
  Copyright (c) 2015 Fabian Sobanski <0sobansk@informatik.uni-hamburg.de>
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

#ifndef ONE_DIM_CHAINER_MATCH_H
#define ONE_DIM_CHAINER_MATCH_H

#include <inttypes.h>

#include "core/ma.h"
#include "match/querymatch.h"

struct GtOneDimChainerMatch;

/* A struct that defines a match object that mightbecome part of the chain.
It contains its succesor <suc> in the maximal chain. <suc> is NULL if this
match is the last in the chain. the sequence ID <queryseqnum>, a counter for
references <refcount>, an index <querystart> for the beginning of this match,
an index <queryend> for its end, as well as a variable <chainweight>. */
typedef struct GtOneDimChainerMatch {
  struct GtOneDimChainerMatch *suc;
  uint64_t queryseqnum;
  unsigned refcount;

  GtUword querystart;
  GtUword queryend;
  GtUword chainweight; /* weight of max chain up to, excluding, this match */
  GtUword dist;
} GtOneDimChainerMatch;

/* Allocates memory for a GtOneDimChainerMatch <match> and returns a pointer to
it. It assigns important variables from a GtQueryMatch reference provided by the
match iterator to the new <match> as well as the currently calculated maximum
chain start <maxchainstart> and the current maximum chain weight <chainweight>.
*/
GtOneDimChainerMatch *gt_1d_chainer_match_new(
    GtQuerymatch *querymatchptr, GtOneDimChainerMatch *maxchainstart,
    GtUword chainweight);

void gt_1d_chainer_match_delete(GtOneDimChainerMatch *one_dim_chainer_match);

/*
 * in: matchfilename, overlap
 * out: chainstart, err
 * returns: had_err
 */
int gt_1d_chainer_calc_chain(const GtStr *matchfilename, GtUword overlap,
                             GtOneDimChainerMatch *chainstart,
                             GtError *err);

#endif
