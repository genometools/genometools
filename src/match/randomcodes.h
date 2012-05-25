/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#ifndef RANDOMCODES_H
#define RANDOMCODES_H

#include <inttypes.h>
#include "core/log_api.h"
#include "core/error_api.h"
#include "core/encseq_api.h"
#include "seqnumrelpos.h"

typedef int (*GtRandomcodesintervalprocess)(void *,
                                           const unsigned long *,
                                           const GtSeqnumrelpos *,
                                           const uint16_t *,
                                           unsigned long,
                                           unsigned int,
                                           GtError *);

typedef void (*GtRandomcodesintervalprocess_end)(void *);

void gt_rungetencseqkmers_rc(const GtEncseq *encseq,unsigned int kmersize);

int storerandomcodes_getencseqkmers_twobitencoding(
                    const GtEncseq *encseq,
                    unsigned int kmersize,
                    unsigned int numofparts,
                    unsigned long maximumspace,
                    unsigned int correction_kmersize,
                    bool usefirstcodes,
                    unsigned int sampling_factor,
                    bool usemaxdepth,
                    bool withsuftabcheck, /* set to false, only for tests */
                    bool onlyaccumulation, /* set to false, only for tests */
                    bool onlyallrandomcodes, /* set to false, only for tests */
                    GT_UNUSED unsigned int addbscache_depth, /* set to 5U */
                    unsigned long phase2extra, /* extra space needed in proc.
                                                  intervals */
                    bool radixsmall,      /* set to true */
                    unsigned int radixparts, /* set to 2U */
                    GtRandomcodesintervalprocess itvprocess,
                    GtRandomcodesintervalprocess_end itvprocess_end,
                    void *itvprocessdatatab,
                    GtLogger *logger,
                    GtTimer *timer,
                    GtError *err);

#endif
