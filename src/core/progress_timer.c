/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg

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
#include <time.h>
#include "core/assert_api.h"
#include "core/ma_api.h"
#include "core/progress_timer_api.h"

struct GtProgressTimer
{
  clock_t startclock, overalltime;
  const char *eventdescription;
};

GtProgressTimer* gt_progress_timer_new(const char *desc)
{
  GtProgressTimer *pt;

  pt = gt_malloc(sizeof (GtProgressTimer));
  if (desc == NULL)
  {
    pt->startclock = pt->overalltime = 0;
    pt->eventdescription = NULL;
  } else
  {
    pt->startclock = clock();
    pt->overalltime = 0;
    pt->eventdescription = desc;
  }
  return pt;
}

void gt_progress_timer_start_new_state(GtProgressTimer *pt,
                                       const char *newevent,
                                       FILE *fp)
{
  clock_t stopclock;
  gt_assert(pt);

  stopclock = clock();
  fprintf(fp,"# TIME %s %.2f\n",pt->eventdescription,
           (double) (stopclock-pt->startclock)/
           (double) CLOCKS_PER_SEC);
  (void) fflush(fp);
  pt->overalltime += (stopclock - pt->startclock);
  if (newevent == NULL)
  {
    fprintf(fp,"# TIME overall %.2f\n",
                (double) pt->overalltime/(double) CLOCKS_PER_SEC);
    (void) fflush(fp);
  } else
  {
    pt->startclock = stopclock;
    pt->eventdescription = newevent;
  }
}

void gt_progress_timer_delete(GtProgressTimer *pt)
{
  if (!pt) return;
  gt_free(pt);
}
