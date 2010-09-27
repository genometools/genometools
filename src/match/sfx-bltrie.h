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

typedef struct Blindtrie Blindtrie;

Blindtrie *gt_blindtrie_new(GtSuffixsortspace *suffixsortspace,
                            unsigned long numofsuffixes,
                            const GtEncseq *encseq,
                            bool cmpcharbychar,
                            GtEncseqReader *esr1,
                            GtEncseqReader *esr2,
                            GtReadmode readmode);

unsigned long gt_blindtrie_suffixsort(
                            Blindtrie *blindtrie,
                            unsigned long subbucketleft,
                            unsigned long *lcpsubtab,
                            unsigned long numberofsuffixes,
                            unsigned long offset,
                            unsigned long maxdepth,
                            void *voiddcov,
                            Dc_processunsortedrange dc_processunsortedrange);

void gt_blindtrie_delete(Blindtrie *blindtrie);

#endif
