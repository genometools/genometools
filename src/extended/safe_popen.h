/*
  Copyright (c) 2015 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#ifndef SAFE_POPEN_H
#define SAFE_POPEN_H

#include <stdio.h>
#include <sys/types.h>
#include "core/error_api.h"

typedef struct {
  FILE  *read_fd;
  FILE  *write_fd;
  pid_t child_pid;
} GtSafePipe;

/* open a child process, execute <path> with <argv> arguments and <envp> by
   calling execve(), returning a pipe to read and write to the new processes
   stdin or stdout.
   <path> should hold the complete path to the executable, <argv> is a NULL
   terminated array of strings, <envp> likewise, where the strings should be of
   the form key=value (see execve()).
   returns NULL on error. */
/* Will not work on windows, as it uses fork() */
GtSafePipe *gt_safe_popen(const char *path,
                          char *const argv[],
                          char *const envp[],
                          GtError *err);

/* Closes pipe <p>, waits for child to exit return 0 on success. */
int         gt_safe_pclose(GtSafePipe *p);

#endif
