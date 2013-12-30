/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef SFX_BLTRIE_H
#define SFX_BLTRIE_H

#include "core/encseq.h"
#include "core/readmode.h"
#include "sfx-suffixgetset.h"
#include "sfx-lcpvalues.h"

typedef struct GtBlindtrie GtBlindtrie;

GtBlindtrie *gt_blindtrie_new(GtSuffixsortspace *suffixsortspace,
                              GtUword maxnumofsuffixes,
                              unsigned int nodenumberincrement,
                              const GtEncseq *encseq,
                              bool cmpcharbychar,
                              GtEncseqReader *esr1,
                              GtEncseqReader *esr2,
                              GtReadmode readmode);

void gt_blindtrie_resize(GtBlindtrie *blindtrie,unsigned int maxnumofnodes);

size_t gt_blindtrie_size(GtUword maxnumofsuffixes);

size_t gt_blindtrie_current_size(const GtBlindtrie *blindtrie);

void gt_blindtrie_reset(GtBlindtrie *blindtrie);

void gt_blindtrie_suffixsort(GtBlindtrie *blindtrie,
                             GtUword subbucketleft,
                             GtLcpvalues *tableoflcpvalues,
                             GtUword numberofsuffixes,
                             GtUword offset,
                             GtUword sortmaxdepth,
                             GtProcessunsortedsuffixrange
                               processunsortedsuffixrange,
                             void *processunsortedsuffixrangeinfo);

void gt_blindtrie_delete(GtBlindtrie *blindtrie);

bool gt_blindtrie_retrieve(GtBlindtrie *blindtrie,
                           GtUword currentstartpos,
                           GtUword currenttwobitencodingstoppos);

#endif
