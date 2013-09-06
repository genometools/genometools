/*
  Copyright (c) 2013 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#ifndef RANDOMCODES_HPCORRECT_H
#define RANDOMCODES_HPCORRECT_H

#include <stdint.h>
#include "core/encseq.h"
#include "match/hplstore.h"
#include "match/seqnumrelpos.h"

int gt_randomcodes_hpcorrect_process_bucket(void *data,
    const GtUword *bucketofsuffixes, const GtSeqnumrelpos *snrp,
    const uint16_t *lcptab_bucket, GtUword numberofsuffixes,
    unsigned int sortingdepth, GtError *err);

typedef struct GtRandomcodesHpcorrectData GtRandomcodesHpcorrectData;

GtRandomcodesHpcorrectData *gt_randomcodes_hpcorrect_data_new(
    GtEncseq *encseq, GtHplstore *hplstore, unsigned int k,
    bool firstpass, unsigned int r /* ignored if !firstpass */,
    unsigned int mintrustp /* use 0 to disable */,
    unsigned int maxuntrustp /* use 100 to disable */,
    bool greedy_clustering, bool skip_read_ends, bool skip_hmer_ends,
    bool skip_rc, bool non_redundant, bool best_score_clustering,
    bool manhattan, /* true: score = - manhattan distance,
                       false: use hard-coded distance matrix  */
    GtWord clustering_param, /* best_score_clustering: minimal scores
                                percentile;
                                otherwise: mimimal score */
    GtUword maxwidth, int rext_cl_minscore,
    int rext_I_minscore, int rext_J_minscore, int rext_R_minscore,
    int rext_D_minscore, int rext_J_lminscore, GtUword rext_J_lwidth,
    GtStr *indexname, unsigned int threadnum, GtError *err);

void gt_randomcodes_hpcorrect_data_delete(GtRandomcodesHpcorrectData *sdata);

#endif
