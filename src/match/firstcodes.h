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

#ifndef FIRSTCODES_H
#define FIRSTCODES_H

#include "core/log_api.h"
#include "core/error_api.h"
#include "core/encseq_api.h"

int storefirstcodes_getencseqkmers_twobitencoding(const GtEncseq *encseq,
                                                  unsigned int kmersize,
                                                  unsigned int numofparts,
                                                  unsigned long maximumspace,
                                                  unsigned int minmatchlength,
                                                  bool withsuftabcheck,
                                                  bool countspms,
                                                  bool outputspms,
                                                  GtLogger *logger,
                                                  GtError *err);

#endif
