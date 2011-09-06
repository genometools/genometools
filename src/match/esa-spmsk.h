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
#include "core/error_api.h"
#include "core/encseq_api.h"

typedef struct GtSpmsk_state GtSpmsk_state;

GtSpmsk_state *gt_spmsk_new(const GtEncseq *encseq,
                            GtReadmode readmode,
                            unsigned long minmatchlength);

void gt_spmsk_delete(GtSpmsk_state *state);

int gt_spmsk_process(GtSpmsk_state *state,
                     const unsigned long *suftab_bucket,
                     const uint16_t *lcptab_bucket,
                     unsigned long nonspecials,
                     GtError *err);

#endif
