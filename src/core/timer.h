/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef TIMER_H
#define TIMER_H

#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>
#include "core/error.h"

/* the timer class */
typedef struct GtTimer GtTimer;

GtTimer* gt_timer_new(void);
void     gt_timer_start(GtTimer*);
void     gt_timer_stop(GtTimer*);
void     gt_timer_show(GtTimer*, FILE*);
/* <fmt> must be a format string for four %ld numbers, which are filled with:
   elapsed seconds, elapsed microseconds, used usertime in seconds,
   systemtime in seconds. */
void     gt_timer_show_formatted(GtTimer*, const char *fmt, FILE*);
void     gt_timer_delete(GtTimer*);

#endif
