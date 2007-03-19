/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgtcore/xbzlib.h>

BZFILE* xbzopen(const char *path, const char *mode)
{
  BZFILE* file;
  if (!(file = BZ2_bzopen(path, mode))) {
    fprintf(stderr, "cannot open file '%s': %s\n", path, strerror(errno));
    exit(EXIT_FAILURE);
  }
  return file;
}

int xbzread(BZFILE* file, void *buf, unsigned len)
{
  int rval;
  if ((rval = BZ2_bzread(file, buf, len)) == -1) {
    fprintf(stderr, "cannod read from compressed file\n");
    exit(EXIT_FAILURE);
  }
  return rval;
}

void xbzrewind(BZFILE **file, const char *orig_path, const char *orig_mode)
{
  /* simulate a rewind with close/open */
  BZ2_bzclose(*file);
  *file = xbzopen(orig_path, orig_mode);
}
