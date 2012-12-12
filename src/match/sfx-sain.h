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
#include "core/encseq.h"

void gt_sain_encseq_sortsuffixes(const GtEncseq *encseq,
                                 GtReadmode readmode,
                                 bool intermediatecheck,
                                 bool finalcheck,
                                 GtLogger *logger,
                                 GtTimer *timer);

void gt_sain_plain_sortsuffixes(const GtUchar *plainseq,
                                unsigned long len,
                                bool intermediatecheck,
                                GtLogger *logger,
                                GtTimer *timer);

#endif
