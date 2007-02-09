/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef XPOSIX_H
#define XPOSIX_H

#include <sys/stat.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <glob.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

/*
  This module contains wrappers for the POSIX functions we use.
  These functions always terminate the program if an error occurs.
  That is, one can use this functions without the need to check for errors.
*/

void   xclose(int d);
void   xfstat(int fd, struct stat *sb);
void   xgetrusage(int who, struct rusage *rusage);
void   xglob(const char *pattern, int flags,
             int (*errfunc)(const char*, int), glob_t *pglob);
int    xopen(const char *path, int flags, mode_t mode);
/* low level wrappers for the mmap() and munmap() routines */
void*  xmmap(void *addr, size_t len, int prot, int flags, int fd, off_t offset);
/* high level routine for memory mapping of files (read only) */
void*  xmap_read(const char *filename, size_t *len);
/* generic unmapping routine */
void   xmunmap(void *addr, size_t len);
void   xraise(int sig);
void (*xsignal(int sigcatch, void (*func)(int sigraised)))(int);
void   xstat(const char *path, struct stat *sb);
time_t xtime(time_t *tloc);
void   xunlink(const char *path);
void   xwrite(int d, const void *buf, size_t nbytes);

#endif
