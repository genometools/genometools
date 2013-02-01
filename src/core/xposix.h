/*
  Copyright (c) 2005-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef XPOSIX_H
#define XPOSIX_H

#ifndef S_SPLINT_S
#include <glob.h>
#include <stdio.h>
#include <time.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/stat.h>

/*
  This module contains wrappers for the POSIX functions we use.
  These functions always terminate the program if an error occurs.
  That is, one can use this functions without the need to check for errors.
*/

void   gt_xclose(int d);
FILE*  gt_xfdopen(int filedes, const char *mode);
void   gt_xfstat(int fd, struct stat *sb);
void   gt_xgetrusage(int who, struct rusage *rusage);
void   gt_xglob(const char *pattern, int flags,
                int (*errfunc)(const char*, int), glob_t *pglob);
int    gt_xopen(const char *path, int flags, mode_t mode);
void   gt_xmkdir(const char *path);
int    gt_xmkstemp(char *temp);
/* low-level wrapper for the mmap() routine */
void*  gt_xmmap(void *addr, size_t len, int prot, int flags, int fd,
                off_t offset);
/* generic unmapping routine */
void   gt_xmunmap(void *addr, size_t len);
void   gt_xraise(int sig);
void (*gt_xsignal(int sigcatch, void (*func)(int sigraised)))(int);
void   gt_xstat(const char *path, struct stat *sb);
time_t gt_xtime(time_t *tloc);
void   gt_xunlink(const char *path);
void   gt_xwrite(int d, const void *buf, size_t nbytes);

#endif
#endif
