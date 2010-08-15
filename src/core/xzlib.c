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

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/xzlib.h"

gzFile gt_xgzopen(const char *path, const char *mode)
{
  gzFile file;
  if (!(file = gzopen(path, mode))) {
    fprintf(stderr, "gzopen(): cannot open file '%s': %s\n", path,
            strerror(errno));
    exit(EXIT_FAILURE);
  }
  return file;
}

int gt_xgzfgetc(gzFile file)
{
  char c;
  return gt_xgzread(file, &c, 1) ? c : EOF;
}

void gt_xgzfputc(int c, gzFile file)
{
  int errnum;
  if (gzputc(file, c) == -1) {
    fprintf(stderr, "cannot put character to compressed file: %s\n",
            gzerror(file, &errnum));
    exit(EXIT_FAILURE);
  }
}

void gt_xgzfputs(const char *str, gzFile file)
{
  int errnum;
  if (gzputs(file, str) == -1) {
    fprintf(stderr, "cannot put string to compressed file: %s\n",
            gzerror(file, &errnum));
    exit(EXIT_FAILURE);
  }
}

int gt_xgzread(gzFile file, void *buf, unsigned len)
{
  int errnum, rval;
  if ((rval = gzread(file, buf, len)) == -1) {
    fprintf(stderr, "cannot read from compressed file: %s\n",
            gzerror(file, &errnum));
    exit(EXIT_FAILURE);
  }
  return rval;
}

void gt_xgzwrite(gzFile file, void *buf, unsigned len)
{
  int errnum;
  gt_assert(buf);
  if (gzwrite(file, buf, len) != len) {
    fprintf(stderr, "cannot write to compressed file: %s\n",
            gzerror(file, &errnum));
    exit(EXIT_FAILURE);
  }
}

void gt_xgzrewind(gzFile file)
{
  if (gzrewind(file) == -1) {
    fprintf(stderr, "cannot rewind compressed file\n");
    exit(EXIT_FAILURE);
  }
}

void gt_xgzclose(gzFile file)
{
  const char *msg;
  int errnum;
  if (gzclose(file)) {
    msg = gzerror(file, &errnum);
    if (errnum == Z_ERRNO) { /* error in file system */
      perror("cannot close compressed file");
      exit(EXIT_FAILURE);
    }
    else { /* error in zlib */
      fprintf(stderr, "cannot close compressed file: %s\n", msg);
      exit(EXIT_FAILURE);
    }
  }
}
