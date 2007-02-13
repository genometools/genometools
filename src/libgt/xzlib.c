/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "xzlib.h"

gzFile xgzopen(const char *path, const char *mode)
{
  gzFile file;
  if (!(file = gzopen(path, mode))) {
    fprintf(stderr, "cannot open file '%s': %s\n", path, strerror(errno));
    exit(EXIT_FAILURE);
  }
  return file;
}

void xgzclose(gzFile file)
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
