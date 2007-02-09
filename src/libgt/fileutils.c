/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "fileutils.h"
#include "xansi.h"
#include "xposix.h"

unsigned int file_exists(const char *path)
{
  FILE *file;
  if ((file = fopen(path, "r")) == NULL)
    return 0;
  xfclose(file);
  return 1;
}

unsigned int file_is_newer(const char *a, const char *b)
{
  struct stat stat_a, stat_b;
  assert(a && b);
  xstat(a, &stat_a);
  xstat(b, &stat_b);
  if (stat_a.st_mtime > stat_b.st_mtime ||
      (stat_a.st_mtime == stat_b.st_mtime &&
       stat_a.st_mtime > stat_b.st_mtime)) {
    return 1;
  }
  return 0;
}

unsigned long file_number_of_lines(FILE *fp)
{
  unsigned long number_of_lines = 0;
  fpos_t current_pos;
  int cc;
  assert(fp);
  xfgetpos(fp, &current_pos);
  xfseek(fp, SEEK_SET, 0);
  while ((cc = getc(fp)) != EOF)
    if (cc == '\n') number_of_lines++;
  if (ferror(fp)) {
    perror("cannot read char");
    exit(EXIT_FAILURE);
  }
  xfsetpos(fp, &current_pos);
  return number_of_lines;
}
