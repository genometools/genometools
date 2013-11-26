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

#ifndef FMI_SUFBWTSTREAM_H
#define FMI_SUFBWTSTREAM_H

#include <stdbool.h>
#include "core/chardef.h"
#include "core/str_array.h"
#include "core/logger_api.h"
#include "core/error_api.h"
#include "fmindex.h"

int gt_sufbwt2fmindex(Fmindex *fmindex,
                      GtSpecialcharinfo *specialcharinfo,
                      unsigned int log2bsize,
                      unsigned int log2markdist,
                      const char *outfmindex,
                      const GtStrArray *indexnametab,
                      bool storeindexpos,
                      GtLogger *logger,
                      GtError *err);

#endif
