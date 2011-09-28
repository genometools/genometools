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

#ifndef TYR_OCCRATIO_H
#define TYR_OCCRATIO_H

#include <stdbool.h>
#include "core/str.h"
#include "core/error_api.h"
#include "core/arraydef.h"
#include "core/logger.h"

int gt_tyr_occratio_func(const char *inputindex,
                         bool scanfile,
                         unsigned long minmersize,
                         unsigned long maxmersize,
                         GtArrayuint64_t *uniquedistribution,
                         GtArrayuint64_t *nonuniquedistribution,
                         GtArrayuint64_t *nonuniquemultidistribution,
                         GtLogger *logger,
                         GtError *err);

#endif
