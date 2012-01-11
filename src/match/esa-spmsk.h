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

#ifndef ESA_SPMSK_H
#define ESA_SPMSK_H

#include <stdint.h>
#include "core/unused_api.h"
#include "core/error_api.h"
#include "core/encseq_api.h"
#include "seqnumrelpos.h"

typedef struct GtBUstate_spmsk GtBUstate_spmsk;

GtBUstate_spmsk *gt_spmsk_inl_new(const GtEncseq *encseq,
                                  GtReadmode readmode,
                                  unsigned long minmatchlength,
                                  bool countspms,
                                  bool outputspms,
                                  GT_UNUSED const char *indexname);

unsigned long gt_spmsk_inl_delete(GtBUstate_spmsk *state);

int gt_spmsk_inl_process(void *data,
                         const unsigned long *seqnum_relpos_bucket,
                         const GtSeqnumrelpos *snrp,
                         const uint16_t *lcptab_bucket,
                         unsigned long nonspecials,
                         unsigned long spaceforbucketprocessing,
                         /* this parameter can be 0 in case where
                            user defined memlimit or derived
                            memlimit is not available. If it is available
                            then, value is > 0 and gives the amount of space
                            in bytes which can be used in the bucket processing
                            without exceeding the given memlimit. Note
                            that after each part, the bucket processing
                            machinery must delete its space to a minimum,
                            as all space required adds up to the maximum.
                         */
                         GtError *err);

void gt_spmsk_inl_process_end(GT_UNUSED void *data);

#endif
