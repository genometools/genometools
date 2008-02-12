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

#include <assert.h>
#include <errno.h>
#include <signal.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/ioctl.h>
#include <sys/termios.h>
#include "libgtcore/unused.h"
#include "libgtcore/xposix.h"

#define DEFAULT_WINDOW_SIZE  80
#define UPDATE_INTERVAL      1U  /* update the progress bar every second */
#define MAXIMUM_WINDOW_SIZE  512 /* the maximum window size
                                    do not change this without adjusting the
                                    snprintf() statment producing the bar! */

static int window_size; /* the window size of the terminal */
static unsigned long long last_computation;
static unsigned long long processed_counter;
static unsigned long long computed_eta;
/* progress counter */
static volatile const unsigned long long *computation_counter;
static volatile sig_atomic_t window_resized;
static time_t computation_start,
              computed_eta_time,
              eta;

static int output_is_possible(void)
{
  return getpgrp() == tcgetpgrp(STDOUT_FILENO);
}

static void set_window_size(void)
{
  struct winsize winsize;
  if (ioctl(STDOUT_FILENO, (unsigned long) TIOCGWINSZ, &winsize) != -1 &&
      winsize.ws_col != 0) {
    if ((int) winsize.ws_col > MAXIMUM_WINDOW_SIZE)
      window_size = MAXIMUM_WINDOW_SIZE;
    else
      window_size = (int) winsize.ws_col;
  }
  else
    window_size = DEFAULT_WINDOW_SIZE;
  window_size += 1; /* for the trailing '\0' */
}

/*
 27% |*************                                              |    00:07 ETA
 27% |*************                                              | 01:00:07 ETA
100% |***********************************************************|    00:10
100% |***********************************************************| 02:00:10
--5--                                                             ------14------
*/

static void refresh_progressbar(void)
{
  time_t seconds, minutes, hours, eta_diff;
  char buf[MAXIMUM_WINDOW_SIZE + 1];
  int bar_length, width;
  double percent;
  processed_counter = *computation_counter;
  /* percent */
  if (last_computation)
    percent = ((double) processed_counter / (double) last_computation);
  else
    percent = 1.0;
  (void) snprintf(buf, window_size, "\r%3d%% ", (int) (percent * 100.0));
  /* the bar */
  bar_length = window_size - 22; /* 5 + 2 + 14 + 1 (for '\0') */
  if (bar_length >= 2) { /* make sure the bar is at least 2 '*' wide */
    width = (int) (percent * bar_length);
    (void) snprintf(buf + strlen(buf), window_size - strlen(buf),
                    "|%.*s%*s|", width,
                    "**************************************************"
                    "**************************************************"
                    "**************************************************"
                    "**************************************************"
                    "**************************************************"
                    "**************************************************"
                    "**************************************************"
                    "**************************************************"
                    "**************************************************"
                    "**************************************************",
                    bar_length - width, "");
  }

#define SPLIT_SECONDS            \
        hours = seconds / 3600;  \
        seconds -= hours * 3600; \
        minutes = seconds / 60;  \
        seconds -= minutes * 60;

  /* ETA */
  if (processed_counter == last_computation) {
    /* computation is finished -> show elapsed time */
    seconds = xtime(NULL) - computation_start;
    SPLIT_SECONDS;
    if (hours) {
      (void) snprintf(buf + strlen(buf), window_size - strlen(buf),
                      " %2u:%02u:%02u", (unsigned int) hours,
                      (unsigned int) minutes, (unsigned int) seconds);
    }
    else {
      (void) snprintf(buf + strlen(buf), window_size - strlen(buf),
                      "    %02u:%02u", (unsigned int) minutes,
                      (unsigned int) seconds);
    }
  }
  else if (processed_counter == 0) {
    (void) snprintf(buf + strlen(buf), window_size - strlen(buf),
                    "    --:-- ETA");
  }
  else {
    if (computed_eta != processed_counter) {
      /* compute new ETA */
      computed_eta = processed_counter;
      computed_eta_time = xtime(NULL);
      eta = (time_t)
            (((double) (computed_eta_time - computation_start) /
              (double) (processed_counter)) *
             (last_computation - processed_counter));
      seconds = eta;
    }
    else {
      /* use previous ETA */
      eta_diff = xtime(NULL) - computed_eta_time;
      if (eta_diff >= eta)
        seconds = 0;
      else
        seconds = eta - eta_diff;
    }
    if (seconds) {
      SPLIT_SECONDS;
      if (hours) {
        (void) snprintf(buf + strlen(buf), window_size - strlen(buf),
                        " %2u:%02u:%02u ETA", (unsigned int) hours,
                        (unsigned int) minutes, (unsigned int) seconds);
      }
      else {
        (void) snprintf(buf + strlen(buf), window_size - strlen(buf),
                        "    %02u:%02u ETA", (unsigned int) minutes,
                        (unsigned int) seconds);
      }
    }
    else {
      (void) snprintf(buf + strlen(buf), window_size - strlen(buf),
                      "    --:-- ETA");
    }
  }

  /* file line with blanks */
  (void) snprintf(buf + strlen(buf), window_size - strlen(buf), "%*s",
                  window_size - (int) strlen(buf), "");
  xwrite(STDOUT_FILENO, buf, strlen(buf));
}

static void update_progressbar(UNUSED int sigraised)
{
  int last_errno = errno;
  assert(sigraised == SIGALRM);
  if (window_resized) {
    set_window_size();
    window_resized = 0;
  }
  if (output_is_possible())
    refresh_progressbar();
  (void) xsignal(SIGALRM, update_progressbar); /* set signal handler again (for
                                                  systems which switch it back
                                                  to SIG_DFL)*/
  (void) alarm(UPDATE_INTERVAL);
  errno = last_errno;
}

static void sig_winch(UNUSED int sigraised)
{
  assert(sigraised == SIGWINCH);
  window_resized = 1;
}

void progressbar_start(const unsigned long long *current_computation,
                       unsigned long long number_of_computations)
{
  computation_counter = current_computation;
  last_computation = number_of_computations;
  computed_eta = 0;
  assert(*current_computation == 0);
  computation_start = xtime(NULL);
  set_window_size();
  if (output_is_possible())
    refresh_progressbar();
  /* register signal handlers */
  (void) xsignal(SIGALRM, update_progressbar); /* the timer */
  (void) xsignal(SIGWINCH, sig_winch);         /* window resizing */
  (void) alarm(UPDATE_INTERVAL);               /* set alarm */
}

void progressbar_stop(void)
{
  (void) alarm(0); /* reset alarm */
  if (!output_is_possible())
    return;
  /* ensure the complete bar has been shown */
  if (processed_counter != last_computation) {
    last_computation = *computation_counter;
    refresh_progressbar();
  }
  xwrite(STDOUT_FILENO, "\n", 1); /* trailing newline */
  /* unregister signal handlers */
  (void) xsignal(SIGALRM, SIG_DFL);  /* the timer */
  (void) xsignal(SIGWINCH, SIG_DFL); /* window resizing */
}
