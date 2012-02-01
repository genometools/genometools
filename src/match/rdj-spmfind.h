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

#ifndef RDJ_SPMFIND_H
#define RDJ_SPMFIND_H

#include <stdint.h>
#include "core/error_api.h"
#include "core/encseq_api.h"
#include "match/seqnumrelpos.h"

/*
 * elimtrans: if false, the blindtrie are never actually used and
 *            all SPM are output, not only the irreducible ones
 * showspm: if true, SPM are shown on stdout in text format; otherwise
 *          they are saved to file in binary format
 */

typedef struct GtBUstate_spm GtBUstate_spmeq;

GtBUstate_spmeq *gt_spmfind_eqlen_state_new(const GtEncseq *encseq,
    unsigned long minmatchlength, unsigned long w_maxsize, bool elimtrans,
    bool showspm, const char *indexname, unsigned int threadnum,
    GtLogger *default_logger, GtLogger *verbose_logger, GtError *err);

void gt_spmfind_eqlen_state_delete(GtBUstate_spmeq *state);

int gt_spmfind_eqlen_process(void *data,
    const unsigned long *seqnum_relpos_bucket, const GtSeqnumrelpos *snrp,
    const uint16_t *lcptab_bucket, unsigned long nonspecials,
    unsigned long spaceforbucketprocessing, GtError *err);

void gt_spmfind_eqlen_process_end(void *data);

typedef struct GtBUstate_spm GtBUstate_spmvar;

GtBUstate_spmvar *gt_spmfind_varlen_state_new(const GtEncseq *encseq,
    unsigned long minmatchlength, unsigned long w_maxsize, bool elimtrans,
    bool showspm, const char *indexname, unsigned int threadnum,
    GtLogger *default_logger, GtLogger *verbose_logger, GtError *err);

void gt_spmfind_varlen_state_delete(GtBUstate_spmvar *state);

int gt_spmfind_varlen_process(void *data,
    const unsigned long *seqnum_relpos_bucket, const GtSeqnumrelpos *snrp,
    const uint16_t *lcptab_bucket, unsigned long nonspecials,
    unsigned long spaceforbucketprocessing, GtError *err);

void gt_spmfind_varlen_process_end(void *data);

unsigned long gt_spmfind_varlen_nof_trans_spm(GtBUstate_spmvar *state);

unsigned long gt_spmfind_varlen_nof_irr_spm(GtBUstate_spmvar *state);

unsigned long gt_spmfind_eqlen_nof_trans_spm(GtBUstate_spmeq *state);

unsigned long gt_spmfind_eqlen_nof_irr_spm(GtBUstate_spmeq *state);

#endif
