/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "error.h"
#include "fileutils.h"
#include "splitter.h"
#include "xansi.h"
#include "xposix.h"

bool file_exists(const char *path)
{
  FILE *file;
  if ((file = fopen(path, "r")) == NULL)
    return false;
  xfclose(file);
  return true;
}

bool file_is_newer(const char *a, const char *b)
{
  struct stat stat_a, stat_b;
  assert(a && b);
  xstat(a, &stat_a);
  xstat(b, &stat_b);
  if (stat_a.st_mtime > stat_b.st_mtime ||
      (stat_a.st_mtime == stat_b.st_mtime &&
       stat_a.st_mtime > stat_b.st_mtime)) {
    return true;
  }
  return false;
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

Str* file_dirname(const char *file)
{
  long i;
  Str *path;
  for (i = strlen(file) - 1; i >= 0; i--) {
    if (file[i] == '/')
      break;
  }
  if (i > 0) {
    path = str_new();
    str_append_cstr_nt(path, file, i);
    return path;
  }
  return NULL;
}

Str* file_find_in_path(const char *file)
{
  char *pathvariable, *pathcomponent = NULL;
  Splitter *splitter;
  unsigned long i;
  Str *path;

  assert(file);

  /* check if 'file' has dirname */
  path = file_dirname(file);
  if (path)
    return path;
  /* 'file' has no dirname -> scan $PATH */
  pathvariable = getenv("PATH");
  if (pathvariable)
    pathvariable = xstrdup(pathvariable); /* make writeable copy */
  else
    error("environment variable $PATH is not defined");
  splitter = splitter_new();
  splitter_split(splitter, pathvariable, strlen(pathvariable), ':');
  path = str_new();
  for (i = 0; i < splitter_size(splitter); i++) {
    pathcomponent = splitter_get_token(splitter, i);
    str_set_length(path, 0);
    str_append_cstr(path, pathcomponent);
    str_append_char(path, '/');
    str_append_cstr(path, file);
    if (file_exists(str_get(path)))
      break;
  }

  if (i < splitter_size(splitter)) { /* file found in path -> return path */
    str_set_length(path, 0);
    str_append_cstr(path, pathcomponent);
  }
  else { /* file not found in path -> return NULL */
    str_free(path);
    path = NULL;
  }

  /* free */
  free(pathvariable);
  splitter_free(splitter);

  return path;
}
