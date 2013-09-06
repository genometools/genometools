/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef RANDOMCODES_FIND_SELDOM_H
#define RANDOMCODES_FIND_SELDOM_H

#include <stdint.h>
#include "core/encseq.h"
#include "match/seqnumrelpos.h"

int gt_randomcodes_find_seldom_process_bucket(void *data,
    const GtUword *bucketofsuffixes, const GtSeqnumrelpos *snrp,
    const uint16_t *lcptab_bucket, GtUword numberofsuffixes,
    unsigned int sortingdepth, GtError *err);

typedef struct GtRandomcodesFindSeldomData GtRandomcodesFindSeldomData;

GtRandomcodesFindSeldomData *gt_randomcodes_find_seldom_data_new(
    GtEncseq *encseq, unsigned int k, unsigned int c,
    GtBitsequence *seldom_reads);

void gt_randomcodes_find_seldom_data_collect_stats(
    GtRandomcodesFindSeldomData *sdata, unsigned int threadnum,
    GtUword *nofseldomkmers);

void gt_randomcodes_find_seldom_data_delete(GtRandomcodesFindSeldomData *sdata);

#endif
