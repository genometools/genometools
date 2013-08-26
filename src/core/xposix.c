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

#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifndef _WIN32
#include <sys/mman.h>
#else
#include <direct.h>
#endif
#include "core/compat.h"
#include "core/xposix.h"

void gt_xclose(int d)
{
  if (close(d)) {
    perror("cannot close file descriptor");
    exit(EXIT_FAILURE);
  }
}

FILE* gt_xfdopen(int filedes, const char *mode)
{
  FILE *fp;
  if (!(fp = fdopen(filedes, mode))) {
    perror("cannot fdopen");
    exit(EXIT_FAILURE);
  }
  return fp;
}

void gt_xfstat(int fd, struct stat *sb)
{
  if (fstat(fd, sb)) {
    perror("cannot fstat");
    exit(EXIT_FAILURE);
  }
}

#ifndef _WIN32
void gt_xgetrusage(int who, struct rusage *rusage)
{
  if (getrusage(who, rusage)) {
    perror("cannot getrusage");
    exit(EXIT_FAILURE);
  }
}
#endif

#ifndef _WIN32
void gt_xglob(const char *pattern, int flags,
              int (*errfunc)(const char*, int), glob_t *pglob)
{
  int rval;
  errno = 0;
  if ((rval = glob(pattern, flags, errfunc, pglob))) {
    fprintf(stderr, "cannot glob: ");
    switch (rval) {
      case GLOB_NOSPACE:
        fprintf(stderr, "out of memory");
        break;
#ifdef GLOB_NOMATCH
      case GLOB_NOMATCH:
        fprintf(stderr, "pattern not found");
        break;
#endif
      default:
        fprintf(stderr, "reason unknown");
    }
    if (errno)
      fprintf(stderr, " (%s)\n", strerror(errno));
    else
      (void) fputc('\n', stderr);
    exit(EXIT_FAILURE);
  }
}
#endif

int gt_xopen(const char *path, int flags, mode_t mode)
{
  int fd;

  if ((fd = open(path, flags, mode)) == -1) {
    fprintf(stderr, "open(): cannot open file descriptor '%s': %s\n", path,
            strerror(errno));
    exit(EXIT_FAILURE);
  }

  return fd;
}

void gt_xmkdir(const char *path)
{
#ifndef _WIN32
  if (mkdir(path, S_IRWXU | S_IRWXG | S_IRWXO)) {
#else
  if (_mkdir(path)) {
#endif
    fprintf(stderr, "mkdir(): cannot make directory '%s': %s\n", path,
            strerror(errno));
    exit(EXIT_FAILURE);
  }
}

int gt_xmkstemp(char *temp)
{
  int fd;
  if ((fd = gt_mkstemp(temp)) == -1) {
    perror("cannot mkstemp");
    exit(EXIT_FAILURE);
  }
  return fd;
}

#ifndef _WIN32
void* gt_xmmap(void *addr, size_t len, int prot, int flags, int fd,
               off_t offset)
{
  void *map;
  if ((map = mmap(addr, len, prot, flags, fd, offset)) == MAP_FAILED) {
    perror("cannot mmap");
    exit(EXIT_FAILURE);
  }
  return map;
}

void gt_xmunmap(void *addr, size_t len)
{
  if (munmap(addr, len)) {
    perror("cannot munmap");
    exit(EXIT_FAILURE);
  }
}
#endif

void gt_xraise(int sig)
{
  if (raise(sig)) {
    perror("cannot raise");
    exit(EXIT_FAILURE);
  }
}

void (*gt_xsignal(int sigcatch, void (*func)(int sigraised)))(int)
{
  void (*rval)(int);
  if ((rval = signal(sigcatch, func)) == SIG_ERR) {
    perror("cannot register signal");
    exit(EXIT_FAILURE);
  }
  return rval;
}

void gt_xstat(const char *path, struct stat *sb)
{
  if (stat(path, sb)) {
    fprintf(stderr, "cannot stat() file '%s': %s\n", path, strerror(errno));
    exit(EXIT_FAILURE);
  }
}

time_t gt_xtime(time_t *tloc)
{
  time_t t;
  if ((t = time(tloc)) == -1) {
    perror("cannot determine time");
    exit(EXIT_FAILURE);
  }
  return t;
}

void gt_xunlink(const char *path)
{
  if (unlink(path)) {
    fprintf(stderr, "cannot unlink \'%s\': %s\n", path, strerror(errno));
    exit(EXIT_FAILURE);
  }
}

void gt_xwrite(int d, const void *buf, size_t nbytes)
{
  size_t pos = 0;
  ssize_t rval;
  int error = 0;
  while (nbytes > pos) {
    rval = write(d, (char*) buf + pos, nbytes - pos);
    switch (rval) {
      case -1:
        if (errno == EINTR || errno == EAGAIN)
          continue;
        error = 1;
        break;
      case 0:
        errno = EPIPE;
        error = 1;
        break;
      default:
        pos += rval;
    }
    if (error)
      break;
  }
  if (error) {
    perror("cannot write to file descriptor");
    exit(EXIT_FAILURE);
  }
}
