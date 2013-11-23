/*
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#ifndef SFX_SAIN_H
#define SFX_SAIN_H

#include "core/timer_api.h"
#include "core/logger_api.h"
#include "core/error_api.h"
#include "core/encseq.h"
#include "bare-encseq.h"

typedef unsigned int GtUsainindextype;

GtUsainindextype *gt_sain_encseq_sortsuffixes(const GtEncseq *encseq,
                                              GtReadmode readmode,
                                              bool intermediatecheck,
                                              bool finalcheck,
                                              GtLogger *logger,
                                              GtTimer *timer);

GtUsainindextype *gt_sain_bare_encseq_sortsuffixes(
                                      const GtBareEncseq *bare_encseq,
                                      GtReadmode readmode,
                                      bool intermediatecheck,
                                      bool finalcheck,
                                      GtLogger *logger,
                                      GtTimer *timer);

GtUsainindextype *gt_sain_plain_sortsuffixes(const GtUchar *plainseq,
                                             GtUword len,
                                             bool intermediatecheck,
                                             bool finalcheck,
                                             GtLogger *logger,
                                             GtTimer *timer);

int gt_sain_checkmaxsequencelength(GtUword len,bool forencseq,GtError *err);

typedef struct GtSainSufLcpIterator GtSainSufLcpIterator;

GtSainSufLcpIterator *gt_sain_suf_lcp_iterator_new(bool withlcp,
                                                   GtUchar *sequence,
                                                   GtUword len,
                                                   GtReadmode readmode,
                                                   GtUword numofchars,
                                                   GtError *err);

GtUword gt_sain_suf_lcp_iterator_nonspecials(const GtSainSufLcpIterator
                                                   *suflcpiterator);

void gt_sain_suf_lcp_iterator_delete(GtSainSufLcpIterator *suflcpiterator);

GtUword gt_sain_suf_lcp_iterator_nonspecials(const GtSainSufLcpIterator
                                                   *suflcpiterator);

GtReadmode gt_sain_suf_lcp_iterator_readmode(const GtSainSufLcpIterator
                                             *suflcpiterator);

GtUword gt_sain_suf_lcp_iterator_next(GtUword *lcpvalue,
                                      GtSainSufLcpIterator *suflcpiterator);

const GtBareEncseq *gt_sain_suf_lcp_iterator_bare_encseq(
         const GtSainSufLcpIterator *suflcpiterator);

#endif
