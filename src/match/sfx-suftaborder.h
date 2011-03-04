/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef SFX_SUFTABORDER_H
#define SFX_SUFTABORDER_H

#include "core/error_api.h"
#include "core/encseq.h"
#include "core/readmode.h"
#include "sfx-suffixgetset.h"

void gt_checksortedsuffixes(const char *filename,
                            int line,
                            const GtEncseq *encseq,
                            GtReadmode readmode,
                            const GtSuffixsortspace *suffixsortspace,
                            unsigned long subbucketleft,
                            unsigned long numberofsuffixes,
                            bool specialsareequal,
                            bool specialsareequalatdepth0,
                            unsigned long depth);

void gt_checkifprefixesareidentical(const char *filename,
                                    int line,
                                    const GtEncseq *encseq,
                                    GtReadmode readmode,
                                    const GtSuffixsortspace *suffixsortspace,
                                    unsigned long subbucketleft,
                                    unsigned long width,
                                    unsigned long depth);
#endif
