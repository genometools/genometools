/*
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#include <stdio.h>
#include "core/spacepeak.h"
#include "core/ma.h"
#include "core/spacecalc.h"
#include "core/thread_api.h"

typedef struct
{
  unsigned long current,
                max;
  GtMutex *mutex;
} GtSpacepeakLogger;

static GtSpacepeakLogger *peaklogger = NULL;

void gt_spacepeak_init(void)
{
  gt_assert(!peaklogger);
  peaklogger = gt_malloc(sizeof (GtSpacepeakLogger));
  peaklogger->current = gt_ma_get_space_current();
  peaklogger->max = 0;
  peaklogger->mutex = gt_mutex_new();
}

void gt_spacepeak_add(unsigned long size)
{
  gt_assert(peaklogger);
  gt_mutex_lock(peaklogger->mutex);
  peaklogger->current += size;
  if (peaklogger->current > peaklogger->max)
    peaklogger->max = peaklogger->current;
  gt_mutex_unlock(peaklogger->mutex);
}

void gt_spacepeak_free(unsigned long size)
{
  gt_assert(peaklogger && size <= peaklogger->current);
  gt_mutex_lock(peaklogger->mutex);
  peaklogger->current -= size;
  gt_mutex_unlock(peaklogger->mutex);
}
unsigned long gt_spacepeak_get_space_peak(void)
{
  gt_assert(peaklogger);
  return peaklogger->max;
}

void gt_spacepeak_show_space_peak(FILE *outfp)
{
  gt_assert(peaklogger);
  fprintf(outfp, "# combined space peak in megabytes: %.2f\n",
          GT_MEGABYTES(peaklogger->max));
}

void gt_spacepeak_clean()
{
  if (!peaklogger) return;
  gt_mutex_delete(peaklogger->mutex);
  gt_free(peaklogger);
}
