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

#ifndef RANDOMSAMPLES_H
#define RANDOMSAMPLES_H

#include "core/encseq_api.h"

typedef struct GtRandomSamples GtRandomSamples;
GtRandomSamples* gt_randomsamples_new(GtEncseq *encseq, GtTimer *timer,
    unsigned int correction_kmersize, GtLogger *default_logger,
    GtLogger *verbose_logger, GtError *err);
int gt_randomsamples_run(GtRandomSamples *rs, unsigned long sampling_factor,
    unsigned int numofparts, unsigned long maximumspace);
void gt_randomsamples_delete(GtRandomSamples* rs);

#endif
