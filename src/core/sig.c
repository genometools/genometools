/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg

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

#include <signal.h>
#include "core/sig.h"
#include "core/xposix.h"

void gt_sig_register_all(void (*func)(int sigraised))
{
  /* POSIX */
  (void) gt_xsignal(SIGABRT, func);
  (void) gt_xsignal(SIGBUS, func);
  (void) gt_xsignal(SIGFPE, func);
  (void) gt_xsignal(SIGHUP, func);
  (void) gt_xsignal(SIGILL, func);
  (void) gt_xsignal(SIGINT, func);
  (void) gt_xsignal(SIGPIPE, func);
  (void) gt_xsignal(SIGQUIT, func);
  (void) gt_xsignal(SIGSEGV, func);
  (void) gt_xsignal(SIGTERM, func);
  (void) gt_xsignal(SIGTSTP, func);
  (void) gt_xsignal(SIGTTIN, func);
  (void) gt_xsignal(SIGTTOU, func);
  (void) gt_xsignal(SIGUSR1, func);
  (void) gt_xsignal(SIGUSR2, func);
  (void) gt_xsignal(SIGSYS, func);
  (void) gt_xsignal(SIGXCPU, func);
  (void) gt_xsignal(SIGXFSZ, func);

  /* OpenBSD */
#ifdef SIGEMT
  (void) gt_xsignal(SIGEMT, func);
#endif
}

void gt_sig_unregister_all(void)
{
  /* POSIX */
  (void) gt_xsignal(SIGABRT, SIG_DFL);
  (void) gt_xsignal(SIGBUS, SIG_DFL);
  (void) gt_xsignal(SIGFPE, SIG_DFL);
  (void) gt_xsignal(SIGHUP, SIG_DFL);
  (void) gt_xsignal(SIGILL, SIG_DFL);
  (void) gt_xsignal(SIGINT, SIG_DFL);
  (void) gt_xsignal(SIGPIPE, SIG_DFL);
  (void) gt_xsignal(SIGQUIT, SIG_DFL);
  (void) gt_xsignal(SIGSEGV, SIG_DFL);
  (void) gt_xsignal(SIGTERM, SIG_DFL);
  (void) gt_xsignal(SIGTSTP, SIG_DFL);
  (void) gt_xsignal(SIGTTIN, SIG_DFL);
  (void) gt_xsignal(SIGTTOU, SIG_DFL);
  (void) gt_xsignal(SIGUSR1, SIG_DFL);
  (void) gt_xsignal(SIGUSR2, SIG_DFL);
  (void) gt_xsignal(SIGSYS, SIG_DFL);
  (void) gt_xsignal(SIGXCPU, SIG_DFL);
  (void) gt_xsignal(SIGXFSZ, SIG_DFL);

  /* OpenBSD */
#ifdef SIGEMT
  (void) gt_xsignal(SIGEMT, SIG_DFL);
#endif
}
