/*
  Copyright (c) 2005-2010 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef XPOSIX_API_H
#define XPOSIX_API_H

#ifndef S_SPLINT_S
#ifndef _WIN32
#include <glob.h>
#endif
#include <stdio.h>
#include <time.h>
#include <sys/types.h>
#ifndef _WIN32
#include <sys/resource.h>
#endif
#include <sys/stat.h>

/* XPOSIX module */

/*
  This module contains wrappers for the POSIX functions we use.
  These functions always terminate the program if an error occurs.
  That is, one can use this functions without the need to check for errors.
*/

/* Wrapper around <close()>, terminating on error. */
void   gt_xclose(int d);
/* Wrapper around <fdopen()>, terminating on error. */
FILE*  gt_xfdopen(int filedes, const char *mode);
/* Wrapper around <fstat()>, terminating on error. */
void   gt_xfstat(int fd, struct stat *sb);
#ifndef _WIN32
/* Wrapper around <getrusage()>, terminating on error. */
void   gt_xgetrusage(int who, struct rusage *rusage);
/* Wrapper around <glob()>, terminating on error. */
void   gt_xglob(const char *pattern, int flags,
                int (*errfunc)(const char*, int), glob_t *pglob);
#endif
/* Wrapper around <open()>, terminating on error. */
int    gt_xopen(const char *path, int flags, mode_t mode);
/* Wrapper around <mkdir()>, terminating on error. */
void   gt_xmkdir(const char *path);
/* Wrapper around <mkstemp()>, terminating on error. */
int    gt_xmkstemp(char *temp);
#ifndef _WIN32
/* Low-level wrapper for the <mmap()> routine, terminating on error. */
void*  gt_xmmap(void *addr, size_t len, int prot, int flags, int fd,
                off_t offset);
/* Generic unmapping routine, terminating on error. */
void   gt_xmunmap(void *addr, size_t len);
#endif
/* Wrapper around <raise()>, terminating on error. */
void   gt_xraise(int sig);
/* Wrapper around <signal()>, terminating on error. */
void (*gt_xsignal(int sigcatch, void (*func)(int sigraised)))(int);
/* Wrapper around <stat()>, terminating on error. */
void   gt_xstat(const char *path, struct stat *sb);
/* Wrapper around <time()>, terminating on error. */
time_t gt_xtime(time_t *tloc);
/* Wrapper around <unlink()>, terminating on error. */
void   gt_xunlink(const char *path);
/* Wrapper around <write()>, terminating on error. */
void   gt_xwrite(int d, const void *buf, size_t nbytes);

#endif /* S_SPLINT_S */
#endif
