/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c)      2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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
#include <stdio.h>
#include <stdlib.h>
#include "libgtcore/array.h"
#include "libgtcore/cstr.h"
#include "libgtcore/fileutils.h"
#include "libgtcore/ma.h"
#include "libgtcore/splitter.h"
#include "libgtcore/xansi.h"
#include "libgtcore/xposix.h"

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

unsigned long file_number_of_lines(const char *path)
{
  unsigned long number_of_lines = 0;
  GenFile *fp;
  int cc;
  assert(path);
  fp = genfile_xopen(path, "r");
  while ((cc = genfile_xfgetc(fp)) != EOF)
    if (cc == '\n') number_of_lines++;
  genfile_close(fp);
  return number_of_lines;
}

const char* file_suffix(const char *path)
{
  const char *suffixptr;
  assert(path);
  suffixptr = path + genfile_basename_length(path) - 1;
  while (suffixptr > path) {
    if (*suffixptr == '/')
      return "";
    else if (*suffixptr == '.')
      break;
    suffixptr--;
  }
  return suffixptr;
}

void file_dirname(Str *path, const char *file)
{
  long i;
  str_reset(path);
  for (i = strlen(file) - 1; i >= 0; i--) {
    if (file[i] == '/')
      break;
  }
  if (i > 0)
    str_append_cstr_nt(path, file, i);
}

int file_find_in_path(Str *path, const char *file, Error *err)
{
  char *pathvariable, *pathcomponent = NULL;
  Splitter *splitter = NULL;
  unsigned long i;
  int had_err = 0;

  error_check(err);
  assert(file);

  /* check if 'file' has dirname */
  file_dirname(path, file);
  if (str_length(path))
    return had_err;
  /* 'file' has no dirname -> scan $PATH */
  pathvariable = getenv("PATH");
  if (pathvariable)
    pathvariable = cstr_dup(pathvariable); /* make writeable copy */
  else {
    error_set(err, "environment variable $PATH is not defined");
    had_err = -1;
  }

  if (!had_err) {
    splitter = splitter_new();
    splitter_split(splitter, pathvariable, strlen(pathvariable), ':');
    for (i = 0; i < splitter_size(splitter); i++) {
      pathcomponent = splitter_get_token(splitter, i);
      str_set_length(path, 0);
      str_append_cstr(path, pathcomponent);
      str_append_char(path, '/');
      str_append_cstr(path, file);
      if (file_exists(str_get(path)))
        break;
    }
    if (i < splitter_size(splitter)) {
      /* file found in path */
      str_set_length(path, 0);
      str_append_cstr(path, pathcomponent);
    }
    else {
      /* file not found in path  */
      str_reset(path);
    }
  }

  /* free */
  ma_free(pathvariable);
  splitter_delete(splitter);

  return had_err;
}

off_t files_estimate_total_size(const StrArray *filenames)
{
  unsigned long filenum;
  off_t totalsize = 0;
  struct stat sb;
  GenFileMode gfm;
  int fd;

  for (filenum = 0; filenum < strarray_size(filenames); filenum++)
  {
    fd = xopen(strarray_get(filenames,filenum), O_RDONLY, 0);
    xfstat(fd, &sb);
    gfm = genfilemode_determine(strarray_get(filenames,filenum));
    if (gfm == GFM_UNCOMPRESSED)
    {
      totalsize += sb.st_size;
    } else
    {
      totalsize += (4*sb.st_size); /* expected compression rate for
                                      sequence is 0.25 */
    }
    xclose(fd);
  }
  return totalsize;
}
