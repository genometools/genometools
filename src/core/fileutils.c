/*
  Copyright (c) 2006-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <stdio.h>
#include <stdlib.h>
#include "core/assert_api.h"
#include "core/array.h"
#include "core/cstr_api.h"
#include "core/fileutils_api.h"
#include "core/file.h"
#include "core/ma.h"
#include "core/sequence_buffer.h"
#include "core/splitter.h"
#include "core/types_api.h"
#include "core/xansi.h"
#include "core/xposix.h"

bool gt_file_exists(const char *path)
{
  FILE *file;
  if ((file = fopen(path, "r")) == NULL)
    return false;
  gt_xfclose(file);
  return true;
}

bool gt_file_is_newer(const char *a, const char *b)
{
  struct stat stat_a, stat_b;
  gt_assert(a && b);
  gt_xstat(a, &stat_a);
  gt_xstat(b, &stat_b);
  if (stat_a.st_mtime > stat_b.st_mtime ||
      (stat_a.st_mtime == stat_b.st_mtime &&
       stat_a.st_mtime > stat_b.st_mtime)) {
    return true;
  }
  return false;
}

unsigned long gt_file_number_of_lines(const char *path)
{
  unsigned long number_of_lines = 0;
  GtFile *fp;
  int cc;
  gt_assert(path);
  fp = gt_file_xopen(path, "r");
  while ((cc = gt_file_xfgetc(fp)) != EOF)
    if (cc == '\n') number_of_lines++;
  gt_file_delete(fp);
  return number_of_lines;
}

const char* gt_file_suffix(const char *path)
{
  const char *suffixptr;
  gt_assert(path);
  suffixptr = path + gt_file_basename_length(path) - 1;
  while (suffixptr > path) {
    if (*suffixptr == '/')
      return "";
    else if (*suffixptr == '.')
      break;
    suffixptr--;
  }
  return suffixptr;
}

void gt_file_dirname(GtStr *path, const char *file)
{
  long i;
  gt_str_reset(path);
  for (i = strlen(file) - 1; i >= 0; i--) {
    if (file[i] == '/')
      break;
  }
  if (i > 0)
    gt_str_append_cstr_nt(path, file, i);
}

int gt_file_find_in_path(GtStr *path, const char *file, GtError *err)
{
  return gt_file_find_in_env(path, file, "PATH", err);
}

int gt_file_find_in_env(GtStr *path, const char *file, const char *env,
                        GtError *err)
{
  char *pathvariable, *pathcomponent = NULL;
  GtSplitter *splitter = NULL;
  unsigned long i;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(file);

  /* check if 'file' has dirname */
  gt_file_dirname(path, file);
  if (gt_str_length(path))
    return had_err;
  /* 'file' has no dirname -> scan $env */
  pathvariable = getenv(env);
  if (pathvariable)
    pathvariable = gt_cstr_dup(pathvariable); /* make writeable copy */
  else {
    gt_error_set(err, "environment variable $%s is not defined", env);
    had_err = -1;
  }

  if (!had_err) {
    splitter = gt_splitter_new();
    gt_splitter_split(splitter, pathvariable, strlen(pathvariable), ':');
    for (i = 0; i < gt_splitter_size(splitter); i++) {
      pathcomponent = gt_splitter_get_token(splitter, i);
      gt_str_reset(path);
      gt_str_append_cstr(path, pathcomponent);
      gt_str_append_char(path, '/');
      gt_str_append_cstr(path, file);
      if (gt_file_exists(gt_str_get(path)))
        break;
    }
    if (i < gt_splitter_size(splitter)) {
      /* file found in path */
      gt_str_reset(path);
      gt_str_append_cstr(path, pathcomponent);
    }
    else {
      /* file not found in path  */
      gt_str_reset(path);
    }
  }

  /* free */
  gt_free(pathvariable);
  gt_splitter_delete(splitter);

  return had_err;
}

off_t gt_file_estimate_size(const char *filename)
{
  off_t size;
  struct stat sb;
  GtFileMode gfm;
  int fd;

  gt_assert(filename);

  fd = gt_xopen(filename, O_RDONLY, 0);
  gt_xfstat(fd, &sb);
  gfm = gt_file_mode_determine(filename);
  if (gfm == GFM_UNCOMPRESSED)
    size = sb.st_size;
  else
    size = sb.st_size  * 4; /* expected compression rate for sequence is 0.25 */
  gt_xclose(fd);

  return size;
}

off_t gt_files_estimate_total_size(const GtStrArray *filenames)
{
  unsigned long filenum;
  off_t totalsize = 0;

  for (filenum = 0; filenum < gt_str_array_size(filenames); filenum++)
    totalsize += gt_file_estimate_size(gt_str_array_get(filenames, filenum));

  return totalsize;
}

int gt_files_guess_if_protein_sequences(const GtStrArray *filenames,
                                        GtError *err)
{
  unsigned int countnonbases = 0,
               currentposition;
  GtUchar currentchar;
  GtSequenceBuffer *fb;
  int retval;

  gt_error_check(err);
  fb = gt_sequence_buffer_new_guess_type(filenames, err);
  if (!fb) return -1;

  for (currentposition = 0; currentposition < 1000U;
       currentposition++)
  {
    retval = gt_sequence_buffer_next(fb,&currentchar,err);
    if (retval < 0)
    {
      gt_sequence_buffer_delete(fb);
      return -1;
    }
    if (retval == 0)
    {
      break;
    }
    switch (currentchar)
    {
      case 'L':
      case 'I':
      case 'F':
      case 'E':
      case 'Q':
      case 'P':
      case 'X':
      case 'Z': countnonbases++;
                break;
      default:  break;
    }
    if (countnonbases > 0)
    {
      break;
    }
  }
  gt_sequence_buffer_delete(fb);
  if (countnonbases > 0)
  {
    return 1; /* guess it is a protein sequence */
  }
  return 0; /* guess it is a dna sequence */
}
