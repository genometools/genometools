/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/ma.h"
#include "libgtcore/timer.h"
#include "libgtcore/xansi.h"
#include "libgtcore/xposix.h"

typedef enum {
  TIMER_RUNNING,
  TIMER_STOPPED
} Timerstate;

struct Timer {
  struct timeval start_tv,
                 stop_tv;
  struct rusage  start_ru,
                 stop_ru;
  Timerstate state;
};

Timer* timer_new(void)
{
  Timer *t;
  t = ma_malloc(sizeof (Timer));
  t->state = TIMER_RUNNING;
  return t;
}

void timer_start(Timer *t)
{
  assert(t);
  gettimeofday(&t->start_tv, NULL);
  xgetrusage(RUSAGE_SELF, &t->start_ru);
  t->state = TIMER_RUNNING;
}

void timer_stop(Timer *t)
{
  assert(t);
  if (t->state == TIMER_RUNNING) {
    gettimeofday(&t->stop_tv, NULL);
    xgetrusage(RUSAGE_SELF, &t->stop_ru);
    t->state = TIMER_STOPPED;
  }
}

void timer_show(Timer *t, FILE *fp)
{
  if (t->state == TIMER_RUNNING)
    timer_stop(t);
  assert(t->state == TIMER_STOPPED);
  fprintf(fp, "%lds real %lds user %lds system\n",
          (long)(t->stop_tv.tv_sec - t->start_tv.tv_sec),
          (long)(t->stop_ru.ru_utime.tv_sec - t->stop_ru.ru_utime.tv_sec),
          (long)(t->stop_ru.ru_stime.tv_sec - t->stop_ru.ru_stime.tv_sec));
}

void timer_del(Timer *t)
{
  if (!t) return;
  ma_free(t);
}
