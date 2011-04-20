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

#ifndef TIMER_API_H
#define TIMER_API_H

#include <stdio.h>

/* The <GtTimer> class encapsulates a timer which can be used for run-time
   measurements. */
typedef struct GtTimer GtTimer;

/* Creates a new <GtTimer> object. */
GtTimer* gt_timer_new(void);
/* Creates a new <GtTimer> object and stores the first description. */
GtTimer* gt_timer_new_with_progress_description(const char* desc);
/* Starts the time measurement on <t>. */
void     gt_timer_start(GtTimer *t);
/* Stops the time measurement on <t>. */
void     gt_timer_stop(GtTimer *t);
/* Outputs the current state of <t> in the format
   "%ld.%06lds real %lds user %lds system" to file
   pointer <fp> (see <gt_timer_show_formatted>).
   The timer is then stopped. */
void     gt_timer_show(GtTimer *t, FILE *fp);
/* Outputs the current state of <t> in a user-defined format given by <fmt>.
   <fmt> must be a format string for four %ld numbers, which are filled with:
   elapsed seconds, elapsed microseconds, used usertime in seconds,
   system time in seconds. The output is written to <fp>. */
void     gt_timer_show_formatted(GtTimer *t, const char *fmt, FILE *fp);
/* Outputs the current state of <t> on <fp> since the last call of
   <gt_timer_show_progress()> or the last start of <t>, along with the current
   description. The timer is not stopped, but updated with <desc> to be the
   next description. */
void     gt_timer_show_progress(GtTimer *t, const char *desc, FILE *fp);
/* Outputs the overall time measured with <t> from start to now on <fp>. */
void     gt_timer_show_progress_final(GtTimer *t, FILE *fp);
/* Show also user and sys time in output of gt_timer_show_progress[_final] */
void     gt_timer_show_cpu_time_by_progress(GtTimer *t);
/* Hide output of last stage time in gt_timer_show_progress_final */
void     gt_timer_omit_last_stage(GtTimer *t);
/* Deletes <t>. */
void     gt_timer_delete(GtTimer *t);

#endif
