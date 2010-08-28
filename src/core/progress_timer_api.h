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

#ifndef PROGRESS_TIMER_API_H
#define PROGRESS_TIMER_API_H
#include <stdio.h>

/* The <GtProgressTimer> class implements a timer which can be used to
   measure the time needed for several consecutive steps.  */
typedef struct GtProgressTimer GtProgressTimer;

/* Creates a new <GtProgressTimer> with the intitial state <desc>.
   Optionally, a flag <with_bar> can be set describing the use of a progress
   bar, so the progress timer will not interfere in that case. */
GtProgressTimer* gt_progress_timer_new(const char *desc);
/* Announce the end of the current state and the beginning of the new state
   to <pt>. The new state will be described by <newevent>. The time needed for
   the now completed state is written to <fp>. 
   To announce the end of timing and print the overall-time use 
   <newevent> = NULL */
void             gt_progress_timer_start_new_state(GtProgressTimer *pt,
                                                   const char *newevent,
                                                   FILE *fp);
/* Deletes <pt> and frees all associated space. */
void             gt_progress_timer_delete(GtProgressTimer *pt);

#endif
