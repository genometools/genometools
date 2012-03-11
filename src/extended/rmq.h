/*
  Copyright (c) 2008 Johannes Fischer <johannes.fischer@kit.edu>
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>

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

#ifndef RMQ_H
#define RMQ_H

#include <inttypes.h>
#include "core/error_api.h"

#ifdef GT_LONGLCPVALUES
typedef unsigned long GtRMQvaluetype;
#else
typedef uint32_t GtRMQvaluetype;
#endif

typedef struct GtRMQ GtRMQ;

GtRMQ*        gt_rmq_new(const GtRMQvaluetype *data, unsigned long size);

size_t gt_rmq_size(const GtRMQ *rmq);

unsigned long gt_rmq_find_min_index(const GtRMQ *rmq, unsigned long start,
                                    unsigned long end);
GtRMQvaluetype gt_rmq_find_min_value(const GtRMQ *rmq, unsigned long start,
                                     unsigned long end);
void          gt_rmq_delete(GtRMQ *rmq);

int           gt_rmq_unit_test(GtError *err);
#endif
