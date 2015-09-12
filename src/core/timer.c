/*
  Copyright (c) 2006-2008 Gordon Gremme <gordon@gremme.org>
  Copyright (c)      2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2006-2012 Center for Bioinformatics, University of Hamburg

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
#include <sys/time.h>
#include "core/cstr_api.h"
#include "core/ma.h"
#include "core/str_api.h"
#include "core/timer_api.h"
#include "core/unused_api.h"
#include "core/xposix.h"

typedef enum {
  TIMER_RUNNING,
  TIMER_STOPPED
} Timerstate;

struct GtTimer {
#ifndef _WIN32
  struct timeval gstart_tv,
                 start_tv,
                 stop_tv;
  struct rusage  gstart_ru,
                 start_ru,
                 stop_ru;
#endif
  Timerstate state;
  char *statedesc;
  bool has_desc;
  bool omit_last_stage;
  bool show_cpu_time;
};

GtTimer* gt_timer_new(void)
{
  GtTimer *t = gt_malloc(sizeof *t);
  t->state = TIMER_RUNNING;
  t->statedesc = NULL;
  t->has_desc = false;
  t->omit_last_stage = false;
  t->show_cpu_time = false;
  return t;
}

GtTimer* gt_timer_new_with_progress_description(const char* desc)
{
  GtTimer *t = gt_timer_new();
  t->statedesc = gt_cstr_dup(desc);
  t->has_desc = true;
  return t;
}

void gt_timer_start(GtTimer *t)
{
#ifndef _WIN32
  gt_assert(t);
  gettimeofday(&t->gstart_tv, NULL);
  gettimeofday(&t->start_tv, NULL);
  gt_xgetrusage(RUSAGE_SELF, &t->start_ru);
  gt_xgetrusage(RUSAGE_SELF, &t->gstart_ru);
  t->state = TIMER_RUNNING;
#else
  /* XXX */
  fprintf(stderr, "gt_timer_start() not implemented\n");
  exit(EXIT_FAILURE);
#endif
}

void gt_timer_stop(GtTimer *t)
{
  gt_assert(t);
#ifndef _WIN32
  if (t->state == TIMER_RUNNING) {
    gettimeofday(&t->stop_tv, NULL);
    gt_xgetrusage(RUSAGE_SELF, &t->stop_ru);
    t->state = TIMER_STOPPED;
  }
#else
  /* XXX */
  fprintf(stderr, "gt_timer_stop() not implemented\n");
  exit(EXIT_FAILURE);
#endif
}

#ifndef _WIN32
static int timeval_subtract(struct timeval *result,
                            struct timeval *x,
                            struct timeval *y)
{
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;
  return x->tv_sec < y->tv_sec;
}

GT_UNUSED
static int timeval_add(struct timeval *result,
                       const struct timeval *x,
                       const struct timeval *y)
{
  result->tv_sec = x->tv_sec + y->tv_sec;
  result->tv_usec = x->tv_usec + y->tv_usec;
  while (result->tv_usec > 1000000) {
    result->tv_usec -= 1000000;
    result->tv_sec++;
  }
  return 0;
}
#endif

void gt_timer_show_formatted(GtTimer *t, const char *fmt, FILE *fp)
{
#ifndef _WIN32
  struct timeval elapsed_tv;
  if (t->state == TIMER_RUNNING)
    gt_timer_stop(t);
  gt_assert(t->state == TIMER_STOPPED);
  timeval_subtract(&elapsed_tv, &t->stop_tv, &t->gstart_tv);
  fprintf(fp, fmt,
          (GtWord)(elapsed_tv.tv_sec),
          (GtWord)(elapsed_tv.tv_usec),
          (GtWord)(t->stop_ru.ru_utime.tv_sec - t->start_ru.ru_utime.tv_sec),
          (GtWord)(t->stop_ru.ru_stime.tv_sec - t->start_ru.ru_stime.tv_sec));
#else
  /* XXX */
  fprintf(stderr, "gt_timer_show_formatted() not implemented\n");
  exit(EXIT_FAILURE);
#endif
}

void gt_timer_get_formatted(GtTimer *t, const char *fmt, GtStr *str)
{
#ifndef _WIN32
  struct timeval elapsed_tv;
  char buf[BUFSIZ];
  if (t->state == TIMER_RUNNING)
    gt_timer_stop(t);
  gt_assert(t->state == TIMER_STOPPED);
  timeval_subtract(&elapsed_tv, &t->stop_tv, &t->gstart_tv);
  (void) snprintf(buf, BUFSIZ-1, fmt,
          (GtWord)(elapsed_tv.tv_sec),
          (GtWord)(elapsed_tv.tv_usec),
          (GtWord)(t->stop_ru.ru_utime.tv_sec - t->start_ru.ru_utime.tv_sec),
          (GtWord)(t->stop_ru.ru_stime.tv_sec - t->start_ru.ru_stime.tv_sec));
  gt_str_append_cstr(str, buf);
#else
  /* XXX */
  fprintf(stderr, "gt_timer_get_formatted() not implemented\n");
  exit(EXIT_FAILURE);
#endif
}

void gt_timer_show(GtTimer *t, FILE *fp)
{
  gt_timer_show_formatted(t, GT_WD ".%06lds real " GT_WD "s user " GT_WD
                          "s system\n", fp);
}

#ifndef _WIN32
static void gt_timer_print_progress_report(GtTimer *t,
    struct timeval *elapsed_tv, struct timeval *elapsed_user_tv,
    struct timeval *elapsed_sys_tv, const char *desc, FILE *fp)
{
  fprintf(fp,"# TIME %s "GT_WD".%02ld",
          desc,
          (GtWord)(elapsed_tv->tv_sec),
          (GtWord)(elapsed_tv->tv_usec)/10000);
  if (t->show_cpu_time) {
    fprintf(fp, " (user: "GT_WD".%02ld; sys: "GT_WD".%02ld)\n",
            (GtWord)(elapsed_user_tv->tv_sec),
            (GtWord)(elapsed_user_tv->tv_usec)/10000,
            (GtWord)(elapsed_sys_tv->tv_sec),
            (GtWord)(elapsed_sys_tv->tv_usec)/10000);
  }
  else {
    fprintf(fp, "\n");
  }
}
#endif

void gt_timer_show_progress(GtTimer *t, const char *desc, FILE *fp)
{
#ifndef _WIN32
  gt_timer_show_progress_formatted(t, fp, "%s", desc);
#else
  /* XXX */
  fprintf(stderr, "gt_timer_show_progress() not implemented\n");
  exit(EXIT_FAILURE);
#endif
}

void gt_timer_show_progress_formatted(GtTimer *t, FILE *fp, const char *desc,
                                      ...)
{
#ifndef _WIN32
  va_list ap;
  gt_assert(t && desc);
  va_start(ap, desc);
  gt_timer_show_progress_va(t, fp, desc, ap);
  va_end(ap);
#else
  /* XXX */
  fprintf(stderr, "gt_timer_show_progress_formatted() not implemented\n");
  exit(EXIT_FAILURE);
#endif
}

void gt_timer_show_progress_va(GtTimer *t, FILE *fp, const char *desc,
                               va_list ap)
{
#ifndef _WIN32
  char buf[BUFSIZ];
  struct timeval elapsed_tv, elapsed_user_tv, elapsed_sys_tv;
  gt_assert(t && desc);

  gettimeofday(&t->stop_tv, NULL);
  gt_xgetrusage(RUSAGE_SELF, &t->stop_ru);
  timeval_subtract(&elapsed_tv, &t->stop_tv, &t->start_tv);
  timeval_subtract(&elapsed_user_tv, &t->stop_ru.ru_utime,
    &t->start_ru.ru_utime);
  timeval_subtract(&elapsed_sys_tv, &t->stop_ru.ru_stime,
    &t->start_ru.ru_stime);
  gt_timer_print_progress_report(t, &elapsed_tv, &elapsed_user_tv,
    &elapsed_sys_tv, t->statedesc, fp);
  if (t->statedesc)
    gt_free(t->statedesc);
  (void) vsnprintf(buf, BUFSIZ, desc, ap);
  t->statedesc = gt_cstr_dup(buf);
  gettimeofday(&t->start_tv, NULL);
  gt_xgetrusage(RUSAGE_SELF, &t->start_ru);
#else
  /* XXX */
  fprintf(stderr, "gt_timer_show_progress_va() not implemented\n");
  exit(EXIT_FAILURE);
#endif
}

void gt_timer_show_progress_final(GtTimer *t, FILE *fp)
{
#ifndef _WIN32
  struct timeval elapsed_tv, elapsed_user_tv, elapsed_sys_tv;
  const char overall_desc[] = "overall";

  gt_timer_stop(t);
  if (!t->omit_last_stage) {
    timeval_subtract(&elapsed_tv, &t->stop_tv, &t->start_tv);
    timeval_subtract(&elapsed_user_tv, &t->stop_ru.ru_utime,
      &t->start_ru.ru_utime);
    timeval_subtract(&elapsed_sys_tv, &t->stop_ru.ru_stime,
      &t->start_ru.ru_stime);
    gt_timer_print_progress_report(t, &elapsed_tv, &elapsed_user_tv,
      &elapsed_sys_tv, t->statedesc, fp);
  }
  timeval_subtract(&elapsed_tv, &t->stop_tv, &t->gstart_tv);
  timeval_subtract(&elapsed_user_tv, &t->stop_ru.ru_utime,
    &t->gstart_ru.ru_utime);
  timeval_subtract(&elapsed_sys_tv, &t->stop_ru.ru_stime,
    &t->gstart_ru.ru_stime);
  gt_timer_print_progress_report(t, &elapsed_tv, &elapsed_user_tv,
    &elapsed_sys_tv, overall_desc, fp);
#else
  /* XXX */
  fprintf(stderr, "gt_timer_show_progress_final() not implemented\n");
  exit(EXIT_FAILURE);
#endif
}

void gt_timer_show_cpu_time_by_progress(GtTimer *t)
{
  t->show_cpu_time = true;
}

void gt_timer_omit_last_stage(GtTimer *t)
{
  t->omit_last_stage = true;
}

void gt_timer_delete(GtTimer *t)
{
  if (!t) return;
  if (t->statedesc)
    gt_free(t->statedesc);
  gt_free(t);
}
