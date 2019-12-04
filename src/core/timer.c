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
#ifdef _WIN32
  #include <windows.h>
  #define _SECOND 10000000
#endif

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
#else
  double gstart_tv,
         start_tv,
         stop_tv;
  GtUword gstart_u, /* user time in 100ns */
          start_u,
          stop_u,
          gstart_k, /* kernel time in 100ns */
          start_k,
          stop_k;
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

#ifdef _WIN32
static void gt_timer_win_getvalues(double *time_val, GtUword *kernel_val,
                                   GtUword *user_val)
{
  GT_UNUSED FILETIME _a, _b;
  FILETIME _c, _d;
  LARGE_INTEGER _tmp, freq;
  QueryPerformanceCounter(&_tmp);
  QueryPerformanceFrequency(&freq);
  *time_val = (double)_tmp.QuadPart / (double) freq.QuadPart;
  GetProcessTimes(GetCurrentProcess(), &_a, &_b, &_c, &_d);
  *kernel_val = ((((GtUword) _c.dwHighDateTime) << 32) + _c.dwLowDateTime);
  *user_val = ((((GtUword) _d.dwHighDateTime) << 32) + _d.dwLowDateTime);
}
#endif

void gt_timer_start(GtTimer *t)
{
  gt_assert(t);
  t->state = TIMER_RUNNING;
#ifndef _WIN32
  gettimeofday(&t->gstart_tv, NULL);
  gettimeofday(&t->start_tv, NULL);
  gt_xgetrusage(RUSAGE_SELF, &t->start_ru);
  gt_xgetrusage(RUSAGE_SELF, &t->gstart_ru);
#else
  gt_timer_win_getvalues(&t->gstart_tv, &t->gstart_k, &t->gstart_u);
  t->start_tv = t->gstart_tv;
  t->start_k = t->gstart_k;
  t->start_u = t->gstart_u;
#endif
}

void gt_timer_stop(GtTimer *t)
{
  gt_assert(t);
  if (t->state == TIMER_RUNNING) {
    #ifndef _WIN32
      gettimeofday(&t->stop_tv, NULL);
      gt_xgetrusage(RUSAGE_SELF, &t->stop_ru);
    #else
      gt_timer_win_getvalues(&t->stop_tv, &t->stop_k, &t->stop_u);
    #endif
    t->state = TIMER_STOPPED;
  }
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
  if (t->state == TIMER_RUNNING)
    gt_timer_stop(t);
  gt_assert(t->state == TIMER_STOPPED);
#ifndef _WIN32
  struct timeval elapsed_tv;
  timeval_subtract(&elapsed_tv, &t->stop_tv, &t->gstart_tv);
  fprintf(fp, fmt,
          (GtWord)(elapsed_tv.tv_sec),
          (GtWord)(elapsed_tv.tv_usec),
          (GtWord)(t->stop_ru.ru_utime.tv_sec - t->start_ru.ru_utime.tv_sec),
          (GtWord)(t->stop_ru.ru_stime.tv_sec - t->start_ru.ru_stime.tv_sec));
#else
  double realtime = (t->stop_tv - t->gstart_tv);
  fprintf(fp, fmt,
          (GtWord)(realtime),
          (GtWord)(1000000.0 * (realtime - (GtWord)realtime)),
          (GtWord)((t->stop_u - t->start_u) / _SECOND),
          (GtWord)((t->stop_k - t->start_k) / _SECOND));
#endif
}

GtWord gt_timer_elapsed_usec(GtTimer *t)
{
  if (t->state == TIMER_RUNNING)
    gt_timer_stop(t);
  gt_assert(t->state == TIMER_STOPPED);
#ifndef _WIN32
  struct timeval elapsed_tv;
  timeval_subtract(&elapsed_tv, &t->stop_tv, &t->gstart_tv);
  return (GtWord) (elapsed_tv.tv_usec + elapsed_tv.tv_sec * 1000000L);
#else
  double realtime = (double)(t->stop_tv - t->gstart_tv);
  return (realtime * 1000000L);
#endif
}

void gt_timer_get_formatted(GtTimer *t,
                            const char *fmt,
                            GtStr *str)
{
  if (t->state == TIMER_RUNNING)
    gt_timer_stop(t);
  gt_assert(t->state == TIMER_STOPPED);
#ifndef _WIN32
  struct timeval elapsed_tv;
  char buf[BUFSIZ];
  timeval_subtract(&elapsed_tv, &t->stop_tv, &t->gstart_tv);
  (void) snprintf(buf, BUFSIZ-1, fmt,
          (GtWord)(elapsed_tv.tv_sec),
          (GtWord)(elapsed_tv.tv_usec),
          (GtWord)(t->stop_ru.ru_utime.tv_sec - t->start_ru.ru_utime.tv_sec),
          (GtWord)(t->stop_ru.ru_stime.tv_sec - t->start_ru.ru_stime.tv_sec));
#else
  double realtime;
  char buf[BUFSIZ];
  realtime = (double)(t->stop_tv - t->gstart_tv);
  (void) snprintf(buf, BUFSIZ-1, fmt,
          (GtWord)(realtime),
          (GtWord)(1000000.0 * (realtime - (GtWord)realtime)),
          (GtWord)((t->stop_u - t->start_u) / _SECOND),
          (GtWord)((t->stop_k - t->start_k) / _SECOND));
#endif
  gt_str_append_cstr(str, buf);
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
#else
static void gt_timer_print_progress_report(GtTimer *t, double realtime,
                                           GtUword utime, GtUword ktime,
                                           const char *desc, FILE *fp)
{
  fprintf(fp,"# TIME %s %.02f", desc, realtime);
  if (t->show_cpu_time) {
    fprintf(fp, " (user: %.02f; sys: %.02f)\n",
               (double)utime / (double)_SECOND,
               (double)ktime / (double)_SECOND);
  }
  else {
    fprintf(fp, "\n");
  }
}
#endif

void gt_timer_show_progress(GtTimer *t,
                            const char *desc,
                            FILE *fp)
{
  gt_timer_show_progress_formatted(t, fp, "%s", desc);
}

void gt_timer_show_progress_formatted(GtTimer *t,
                                      FILE *fp,
                                      const char *desc,
                                      ...)
{
  va_list ap;
  gt_assert(t && desc);
  va_start(ap, desc);
  gt_timer_show_progress_va(t, fp, desc, ap);
  va_end(ap);
}

void gt_timer_show_progress_va(GtTimer *t,
                               FILE *fp,
                               const char *desc,
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
  char buf[BUFSIZ];
  gt_assert(t && desc);
  gt_timer_win_getvalues(&t->stop_tv, &t->stop_k, &t->stop_u);
  gt_timer_print_progress_report(t, t->stop_tv - t->start_tv,
                                 t->stop_u - t->start_u,
                                 t->stop_k - t->start_k, t->statedesc, fp);
  if (t->statedesc)
    gt_free(t->statedesc);
  (void) vsnprintf(buf, BUFSIZ, desc, ap);
  t->statedesc = gt_cstr_dup(buf);
  gt_timer_win_getvalues(&t->start_tv, &t->start_k, &t->start_u);
#endif
}

void gt_timer_show_progress_final(GtTimer *t,
                                  FILE *fp)
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
  const char overall_desc[] = "overall";
  gt_timer_stop(t);
  if (!t->omit_last_stage) {
    gt_timer_print_progress_report(t, t->stop_tv - t->start_tv,
                                   t->stop_u - t->start_u,
                                   t->stop_k - t->start_k, t->statedesc, fp);
  }
  gt_timer_print_progress_report(t, t->stop_tv - t->gstart_tv,
                                 t->stop_u - t->gstart_u,
                                 t->stop_k - t->gstart_k, overall_desc, fp);
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
