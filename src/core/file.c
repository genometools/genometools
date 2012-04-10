/*
  Copyright (c) 2005-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <stdio.h>
#include <string.h>
#include "core/cstr_api.h"
#include "core/fa.h"
#include "core/ma.h"
#include "core/xansi_api.h"
#include "core/xbzlib.h"
#include "core/xzlib.h"

struct GtFile {
  GtFileMode mode;
  union {
    FILE *file;
    gzFile gzfile;
    BZFILE *bzfile;
  } fileptr;
  char *orig_path,
       *orig_mode,
       unget_char;
  bool is_stdin,
       unget_used;
};

GtFileMode gt_file_mode_determine(const char *path)
{
  size_t path_length;
  if (!path)
    return GT_FILE_MODE_UNCOMPRESSED;
  path_length = strlen(path);
  if (path_length >= 4 && strcmp(".gz", path + path_length - 3) == 0)
    return GT_FILE_MODE_GZIP;
  if (path_length >= 5 && strcmp(".bz2", path + path_length - 4) == 0)
    return GT_FILE_MODE_BZIP2;
  return GT_FILE_MODE_UNCOMPRESSED;
}

const char* gt_file_mode_suffix(GtFileMode mode)
{
  switch (mode) {
    case GT_FILE_MODE_UNCOMPRESSED:
      return "";
    case GT_FILE_MODE_GZIP:
      return ".gz";
    case GT_FILE_MODE_BZIP2:
      return ".bz2";
    default:
      gt_assert(0);
      return "";
  }
  /* due do warning on solaris:
     warning: control reaches end of non-void function */
  return "";
}

size_t gt_file_basename_length(const char *path)
{
  size_t path_length;

  gt_assert(path);
  path_length = strlen(path);

  if (path_length >= 4 && strcmp(".gz", path + path_length - 3) == 0)
    return path_length - 3;
  if (path_length >= 5 && strcmp(".bz2", path + path_length - 4) == 0)
    return path_length - 4;
  return path_length;
}

GtFile* gt_file_new(const char *path, const char *mode, GtError *err)
{
  gt_error_check(err);
  gt_assert(mode);
  return gt_file_open(gt_file_mode_determine(path), path, mode, err);
}

GtFile* gt_file_open(GtFileMode file_mode, const char *path, const char *mode,
                     GtError *err)
{
  GtFile *file;
  gt_error_check(err);
  gt_assert(mode);
  file = gt_calloc(1, sizeof (GtFile));
  file->mode = file_mode;
  if (path) {
    switch (file_mode) {
      case GT_FILE_MODE_UNCOMPRESSED:
        file->fileptr.file = gt_fa_fopen(path, mode, err);
        if (!file->fileptr.file) {
          gt_file_delete_without_handle(file);
          return NULL;
        }
        break;
      case GT_FILE_MODE_GZIP:
        file->fileptr.gzfile = gt_fa_gzopen(path, mode, err);
        if (!file->fileptr.gzfile) {
          gt_file_delete_without_handle(file);
          return NULL;
        }
        break;
      case GT_FILE_MODE_BZIP2:
        file->fileptr.bzfile = gt_fa_bzopen(path, mode, err);
        if (!file->fileptr.bzfile) {
          gt_file_delete_without_handle(file);
          return NULL;
        }
        file->orig_path = gt_cstr_dup(path);
        file->orig_mode = gt_cstr_dup(path);
        break;
      default: gt_assert(0);
    }
  }
  else {
    gt_assert(file_mode == GT_FILE_MODE_UNCOMPRESSED);
    file->fileptr.file = stdin;
    file->is_stdin = true;
  }
  return file;
}

GtFile* gt_file_xopen_file_mode(GtFileMode file_mode, const char *path,
                                const char *mode)
{
  GtFile *file;
  gt_assert(mode);
  file = gt_calloc(1, sizeof (GtFile));
  file->mode = file_mode;
  if (path) {
    switch (file_mode) {
      case GT_FILE_MODE_UNCOMPRESSED:
        file->fileptr.file = gt_fa_xfopen(path, mode);
        break;
      case GT_FILE_MODE_GZIP:
        file->fileptr.gzfile = gt_fa_xgzopen(path, mode);
        break;
      case GT_FILE_MODE_BZIP2:
        file->fileptr.bzfile = gt_fa_xbzopen(path, mode);
        file->orig_path = gt_cstr_dup(path);
        file->orig_mode = gt_cstr_dup(path);
        break;
      default: gt_assert(0);
    }
  }
  else {
    gt_assert(file_mode == GT_FILE_MODE_UNCOMPRESSED);
    file->fileptr.file = stdin;
    file->is_stdin = true;
  }
  return file;
}

GtFile* gt_file_xopen(const char *path, const char *mode)
{
  gt_assert(mode);
  return gt_file_xopen_file_mode(gt_file_mode_determine(path), path, mode);
}

GtFile* gt_file_new_from_fileptr(FILE *fp)
{
  GtFile *file;
  gt_assert(fp);
  file = gt_calloc(1, sizeof (GtFile));
  file->mode = GT_FILE_MODE_UNCOMPRESSED;
  file->fileptr.file = fp;
  return file;
}

GtFileMode gt_file_mode(const GtFile *file)
{
  gt_assert(file);
  return file->mode;
}

int gt_file_xfgetc(GtFile *file)
{
  int c = -1;
  if (file) {
    if (file->unget_used) {
      c = file->unget_char;
      file->unget_used = false;
    }
    else {
      switch (file->mode) {
        case GT_FILE_MODE_UNCOMPRESSED:
          c = gt_xfgetc(file->fileptr.file);
          break;
        case GT_FILE_MODE_GZIP:
          c = gt_xgzfgetc(file->fileptr.gzfile);
          break;
        case GT_FILE_MODE_BZIP2:
          c = gt_xbzfgetc(file->fileptr.bzfile);
          break;
        default: gt_assert(0);
      }
    }
  }
  else
    c = gt_xfgetc(stdin);
  return c;
}

void gt_file_unget_char(GtFile *file, char c)
{
  if (file) {
    gt_assert(!file->unget_used); /* only one char can be unget at a time */
    file->unget_char = c;
    file->unget_used = true;
  }
  else
    gt_xungetc(c, stdin);
}

static int vgzprintf(gzFile file, const char *format, va_list va, int buflen)
{
  int len;
  if (!buflen) {
    char buf[BUFSIZ];
    /* no buffer length given -> try static buffer */
    len = gt_xvsnprintf(buf, sizeof (buf), format, va);
    if (len >= BUFSIZ) {
      return len; /* unsuccessful trial -> return buffer length for next call */
    }
    gt_xgzwrite(file, buf, len);
  }
  else {
    char *dynbuf;
    /* buffer length given -> use dynamic buffer */
    dynbuf = gt_malloc((buflen + 1) * sizeof (char));
    len = gt_xvsnprintf(dynbuf, (buflen + 1) * sizeof (char), format, va);
    gt_assert(len == buflen);
    gt_xgzwrite(file, dynbuf, buflen);
    gt_free(dynbuf);
  }
  return 0; /* success */
}

static int vbzprintf(BZFILE *file, const char *format, va_list va, int buflen)
{
  int len;
  if (!buflen) {
    char buf[BUFSIZ];
    /* no buffer length given -> try static buffer */
    len = gt_xvsnprintf(buf, sizeof (buf), format, va);
    if (len >= BUFSIZ)
      return len; /* unsuccessful trial -> return buffer length for next call */
    gt_xbzwrite(file, buf, len);
  }
  else {
    char *dynbuf;
    /* buffer length given -> use dynamic buffer */
    dynbuf = gt_malloc((buflen + 1) * sizeof (char));
    len = gt_xvsnprintf(dynbuf, (buflen + 1) * sizeof (char), format, va);
    gt_assert(len == buflen);
    gt_xbzwrite(file, dynbuf, buflen);
    gt_free(dynbuf);
  }
  return 0; /* success */
}

static int xvprintf(GtFile *file, const char *format, va_list va, int buflen)
{
  int rval = 0;

  if (!file) /* implies stdout */
    gt_xvfprintf(stdout, format, va);
  else {
    switch (file->mode) {
      case GT_FILE_MODE_UNCOMPRESSED:
        gt_xvfprintf(file->fileptr.file, format, va);
        break;
      case GT_FILE_MODE_GZIP:
        rval = vgzprintf(file->fileptr.gzfile, format, va, buflen);
        break;
      case GT_FILE_MODE_BZIP2:
        rval = vbzprintf(file->fileptr.bzfile, format, va, buflen);
        break;
      default: gt_assert(0);
    }
  }
  return rval;
}

void gt_file_xprintf(GtFile *file, const char *format, ...)
{
  va_list va;
  int rval;
  va_start(va, format);
  if ((rval = xvprintf(file, format, va, 0))) {
    gt_assert(rval > 0); /* negative return is not possible, rval should denote
                            the necessary buffer length -> try again with it */
    /* reset variable arguments */
    va_end(va);
    va_start(va, format);
    /* try again */
    rval = xvprintf(file, format, va, rval);
    gt_assert(!rval); /* xvprintf() should not fail with given buffer length */
  }
  va_end(va);
}

void gt_file_xfputc(int c, GtFile *file)
{
  if (!file)
    return gt_xfputc(c, stdout);
  switch (file->mode) {
    case GT_FILE_MODE_UNCOMPRESSED:
      gt_xfputc(c, file->fileptr.file);
      break;
    case GT_FILE_MODE_GZIP:
      gt_xgzfputc(c, file->fileptr.gzfile);
      break;
    case GT_FILE_MODE_BZIP2:
      gt_xbzfputc(c, file->fileptr.bzfile);
      break;
    default: gt_assert(0);
  }
}

void gt_file_xfputs(const char *cstr, GtFile *file)
{
  if (!file)
    return gt_xfputs(cstr, stdout);
  switch (file->mode) {
    case GT_FILE_MODE_UNCOMPRESSED:
      gt_xfputs(cstr, file->fileptr.file);
      break;
    case GT_FILE_MODE_GZIP:
      gt_xgzfputs(cstr, file->fileptr.gzfile);
      break;
    case GT_FILE_MODE_BZIP2:
      gt_xbzfputs(cstr, file->fileptr.bzfile);
      break;
    default: gt_assert(0);
  }
}

int gt_file_xread(GtFile *file, void *buf, size_t nbytes)
{
  int rval = -1;
  if (file) {
    switch (file->mode) {
      case GT_FILE_MODE_UNCOMPRESSED:
        rval = gt_xfread(buf, 1, nbytes, file->fileptr.file);
        break;
      case GT_FILE_MODE_GZIP:
        rval = gt_xgzread(file->fileptr.gzfile, buf, nbytes);
        break;
      case GT_FILE_MODE_BZIP2:
        rval = gt_xbzread(file->fileptr.bzfile, buf, nbytes);
        break;
      default: gt_assert(0);
    }
  }
  else
    rval = gt_xfread(buf, 1, nbytes, stdin);
  return rval;
}

void gt_file_xwrite(GtFile *file, void *buf, size_t nbytes)
{
  if (!file) {
    gt_xfwrite(buf, 1, nbytes, stdout);
    return;
  }
  switch (file->mode) {
    case GT_FILE_MODE_UNCOMPRESSED:
      gt_xfwrite(buf, 1, nbytes, file->fileptr.file);
      break;
    case GT_FILE_MODE_GZIP:
      gt_xgzwrite(file->fileptr.gzfile, buf, nbytes);
      break;
    case GT_FILE_MODE_BZIP2:
      gt_xbzwrite(file->fileptr.bzfile, buf, nbytes);
      break;
    default: gt_assert(0);
  }
}

void gt_file_xrewind(GtFile *file)
{
  gt_assert(file);
  switch (file->mode) {
    case GT_FILE_MODE_UNCOMPRESSED:
      rewind(file->fileptr.file);
      break;
    case GT_FILE_MODE_GZIP:
      gt_xgzrewind(file->fileptr.gzfile);
      break;
    case GT_FILE_MODE_BZIP2:
      gt_xbzrewind(&file->fileptr.bzfile, file->orig_path, file->orig_mode);
      break;
    default: gt_assert(0);
  }
}

void gt_file_delete_without_handle(GtFile *file)
{
  if (!file) return;
  gt_free(file->orig_path);
  gt_free(file->orig_mode);
  gt_free(file);
}

void gt_file_delete(GtFile *file)
{
  if (!file) return;
  switch (file->mode) {
    case GT_FILE_MODE_UNCOMPRESSED:
        if (!file->is_stdin)
          gt_fa_fclose(file->fileptr.file);
      break;
    case GT_FILE_MODE_GZIP:
        gt_fa_gzclose(file->fileptr.gzfile);
      break;
    case GT_FILE_MODE_BZIP2:
        gt_fa_bzclose(file->fileptr.bzfile);
      break;
    default: gt_assert(0);
  }
  gt_file_delete_without_handle(file);
}
