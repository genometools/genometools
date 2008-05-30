/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libgtcore/xbzlib.h"

BZFILE* xbzopen(const char *path, const char *mode)
{
  BZFILE* file;
  if (!(file = BZ2_bzopen(path, mode))) {
    fprintf(stderr, "BZ2_bzopen(): cannot open file '%s': %s\n", path,
            strerror(errno));
    exit(EXIT_FAILURE);
  }
  return file;
}

int xbzfgetc(BZFILE *bzfile)
{
  char c;
  return xbzread(bzfile, &c, 1) ? c : EOF;
}

static int bzputc(int c, BZFILE *bzfile)
{
  char cc = (char) c; /* required for big endian systems */
  return BZ2_bzwrite(bzfile, &cc, 1) == 1 ? cc : -1;
}

void xbzfputc(int c, BZFILE *bzfile)
{
  if (bzputc(c, bzfile) == -1) {
    fprintf(stderr, "cannot put character to compressed file\n");
    exit(EXIT_FAILURE);
  }
}

static int bzputs(const char *str, BZFILE *bzfile)
{
  int len = strlen(str);
  return BZ2_bzwrite(bzfile, (char*) str, len) == len ? len : -1;
}

void xbzfputs(const char *str, BZFILE *bzfile)
{
  if (bzputs(str, bzfile) == -1) {
    fprintf(stderr, "cannot put string to compressed file\n");
    exit(EXIT_FAILURE);
  }
}

int xbzread(BZFILE *file, void *buf, unsigned len)
{
  int rval;
  if ((rval = BZ2_bzread(file, buf, len)) == -1) {
    fprintf(stderr, "cannot read from compressed file\n");
    exit(EXIT_FAILURE);
  }
  return rval;
}

void xbzwrite(BZFILE *file, void *buf, unsigned len)
{
  assert(buf && len);
  if (BZ2_bzwrite(file, buf, len) != len) {
    fprintf(stderr, "cannot write it compressed file\n");
    exit(EXIT_FAILURE);
  }
}

void xbzrewind(BZFILE **file, const char *orig_path, const char *orig_mode)
{
  /* simulate a rewind with close/open */
  BZ2_bzclose(*file);
  *file = xbzopen(orig_path, orig_mode);
}
