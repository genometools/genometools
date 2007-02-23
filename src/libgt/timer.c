/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "timer.h"
#include "xansi.h"
#include "xposix.h"

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

Timer* timer_new(Env *env)
{
  Timer *t;
  t = env_ma_malloc(env, sizeof (Timer));
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
          t->stop_tv.tv_sec - t->start_tv.tv_sec,
          t->stop_ru.ru_utime.tv_sec - t->stop_ru.ru_utime.tv_sec,
          t->stop_ru.ru_stime.tv_sec - t->stop_ru.ru_stime.tv_sec);
}

void timer_del(Timer *t, Env *env)
{
  if (!t) return;
  env_ma_free(t, env);
}
