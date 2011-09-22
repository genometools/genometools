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

void gt_spmsk_inl_delete(GtBUstate_spmsk *state);

int gt_spmsk_inl_process(void *data,
                         const unsigned long *suftab_bucket,
                         const GtSeqnumrelpostab *snrp,
                         const uint16_t *lcptab_bucket,
                         unsigned long nonspecials,
                         GtError *err);

#endif
