/*
  Copyright (c) 2005-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "libgtcore/cstr.h"
#include "libgtcore/fa.h"
#include "libgtcore/genfile.h"
#include "libgtcore/ma.h"
#include "libgtcore/xansi.h"
#include "libgtcore/xbzlib.h"
#include "libgtcore/xzlib.h"

struct GenFile {
  GenFileMode mode;
  union {
    FILE *file;
    gzFile gzfile;
    BZFILE *bzfile;
  } fileptr;
  char *orig_path,
       *orig_mode;
};

GenFileMode genfilemode_determine(const char *path)
{
  size_t path_length;

  assert(path);
  path_length = strlen(path);

  if (path_length >= 4 && strcmp(".gz", path + path_length - 3) == 0)
    return GFM_GZIP;
  if (path_length >= 5 && strcmp(".bz2", path + path_length - 4) == 0)
    return GFM_BZIP2;
  return GFM_UNCOMPRESSED;
}

const char* genfilemode_suffix(GenFileMode mode)
{
  switch (mode) {
    case GFM_UNCOMPRESSED:
      return "";
    case GFM_GZIP:
      return ".gz";
    case GFM_BZIP2:
      return ".bz2";
    default:
      assert(0);
      return "";
  }
  /* due do warning on solaris:
     warning: control reaches end of non-void function */
  return "";
}

size_t genfile_basename_length(const char *path)
{
  size_t path_length;

  assert(path);
  path_length = strlen(path);

  if (path_length >= 4 && strcmp(".gz", path + path_length - 3) == 0)
    return path_length - 3;
  if (path_length >= 5 && strcmp(".bz2", path + path_length - 4) == 0)
    return path_length - 4;
  return path_length;
}

GenFile* genfile_open(GenFileMode genfilemode, const char *path,
                      const char *mode)
{
  GenFile *genfile;
  assert(path && mode);
  genfile = ma_calloc(1, sizeof (GenFile));
  genfile->mode = genfilemode;
  switch (genfilemode) {
    case GFM_UNCOMPRESSED:
      genfile->fileptr.file = fa_fopen(path, mode);
      if (!genfile->fileptr.file) {
        genfile_delete(genfile);
        return NULL;
      }
      break;
    case GFM_GZIP:
      genfile->fileptr.gzfile = fa_gzopen(path, mode);
      if (!genfile->fileptr.gzfile) {
        genfile_delete(genfile);
        return NULL;
      }
      break;
    case GFM_BZIP2:
      genfile->fileptr.bzfile = fa_bzopen(path, mode);
      if (!genfile->fileptr.bzfile) {
        genfile_delete(genfile);
        return NULL;
      }
      genfile->orig_path = cstr_dup(path);
      genfile->orig_mode = cstr_dup(path);
      break;
    default: assert(0);
  }
  return genfile;
}

GenFile* genfile_xopen_w_gfmode(GenFileMode genfilemode, const char *path,
                                const char *mode)
{
  GenFile *genfile;
  assert(path && mode);
  genfile = ma_calloc(1, sizeof (GenFile));
  genfile->mode = genfilemode;
  switch (genfilemode) {
    case GFM_UNCOMPRESSED:
      genfile->fileptr.file = fa_xfopen(path, mode);
      break;
    case GFM_GZIP:
      genfile->fileptr.gzfile = fa_xgzopen(path, mode);
      break;
    case GFM_BZIP2:
      genfile->fileptr.bzfile = fa_xbzopen(path, mode);
      genfile->orig_path = cstr_dup(path);
      genfile->orig_mode = cstr_dup(path);
      break;
    default: assert(0);
  }
  return genfile;
}

GenFile* genfile_xopen(const char *path, const char *mode)
{
  assert(path && mode);
  return genfile_xopen_w_gfmode(genfilemode_determine(path), path, mode);
}

GenFile* genfile_new(FILE *fp)
{
  GenFile *genfile;
  assert(fp);
  genfile = ma_calloc(1, sizeof (GenFile));
  genfile->mode = GFM_UNCOMPRESSED;
  genfile->fileptr.file = fp;
  return genfile;
}

GenFileMode genfile_mode(GenFile *genfile)
{
  assert(genfile);
  return genfile->mode;
}

int genfile_xfgetc(GenFile *genfile)
{
  int c = -1;
  if (genfile) {
    switch (genfile->mode) {
      case GFM_UNCOMPRESSED:
        c = xfgetc(genfile->fileptr.file);
        break;
      case GFM_GZIP:
        c = xgzfgetc(genfile->fileptr.gzfile);
        break;
      case GFM_BZIP2:
        c = xbzfgetc(genfile->fileptr.bzfile);
        break;
      default: assert(0);
    }
  }
  else
    c = xfgetc(stdin);
  return c;
}

static int vgzprintf(gzFile file, const char *format, va_list va)
{
  char buf[BUFSIZ];
  int len;
  len = vsnprintf(buf, sizeof (buf), format, va);
  assert(len <= BUFSIZ);
  return gzwrite(file, buf, (unsigned) len);
}

static int vbzprintf(BZFILE *file, const char *format, va_list va)
{
  char buf[BUFSIZ];
  int len;
  len = vsnprintf(buf, sizeof (buf), format, va);
  assert(len <= BUFSIZ);
  return BZ2_bzwrite(file, buf, len);
}

static int xvprintf(GenFile *genfile, const char *format, va_list va)
{
  int rval = -1;

  if (!genfile) /* implies stdout */
    rval = vfprintf(stdout, format, va);
  else {
    switch (genfile->mode) {
      case GFM_UNCOMPRESSED:
        rval = vfprintf(genfile->fileptr.file, format, va);
        break;
      case GFM_GZIP:
        rval = vgzprintf(genfile->fileptr.gzfile, format, va);
        break;
      case GFM_BZIP2:
        rval = vbzprintf(genfile->fileptr.bzfile, format, va);
        break;
      default: assert(0);
    }
  }
  return rval;
}

void genfile_xprintf(GenFile *genfile, const char *format, ...)
{
  va_list va;
  va_start(va, format);
  if (xvprintf(genfile, format, va) < 0) {
    fprintf(stderr, "genfile_xprintf(): xvprintf() returned negative value\n");
    exit(EXIT_FAILURE);
  }
  va_end(va);
}

void genfile_xfputc(int c, GenFile *genfile)
{
  if (!genfile)
    return xfputc(c, stdout);
  switch (genfile->mode) {
    case GFM_UNCOMPRESSED:
      xfputc(c, genfile->fileptr.file);
      break;
    case GFM_GZIP:
      xgzfputc(c, genfile->fileptr.gzfile);
      break;
    case GFM_BZIP2:
      xbzfputc(c, genfile->fileptr.bzfile);
      break;
    default: assert(0);
  }
}

void genfile_xfputs(const char *str, GenFile *genfile)
{
  if (!genfile)
    return xfputs(str, stdout);
  switch (genfile->mode) {
    case GFM_UNCOMPRESSED:
      xfputs(str, genfile->fileptr.file);
      break;
    case GFM_GZIP:
      xgzfputs(str, genfile->fileptr.gzfile);
      break;
    case GFM_BZIP2:
      xbzfputs(str, genfile->fileptr.bzfile);
      break;
    default: assert(0);
  }
}

int genfile_xread(GenFile *genfile, void *buf, size_t nbytes)
{
  int rval = -1;
  if (genfile) {
    switch (genfile->mode) {
      case GFM_UNCOMPRESSED:
        rval = xfread(buf, 1, nbytes, genfile->fileptr.file);
        break;
      case GFM_GZIP:
        rval = xgzread(genfile->fileptr.gzfile, buf, nbytes);
        break;
      case GFM_BZIP2:
        rval = xbzread(genfile->fileptr.bzfile, buf, nbytes);
        break;
      default: assert(0);
    }
  }
  else
    rval = xfread(buf, 1, nbytes, stdin);
  return rval;
}

void genfile_xwrite(GenFile *genfile, void *buf, size_t nbytes)
{
  if (!genfile) {
    xfwrite(buf, 1, nbytes, stdout);
    return;
  }
  switch (genfile->mode) {
    case GFM_UNCOMPRESSED:
      xfwrite(buf, 1, nbytes, genfile->fileptr.file);
      break;
    case GFM_GZIP:
      xgzwrite(genfile->fileptr.gzfile, buf, nbytes);
      break;
    case GFM_BZIP2:
      xbzwrite(genfile->fileptr.bzfile, buf, nbytes);
      break;
    default: assert(0);
  }
}

void genfile_xrewind(GenFile *genfile)
{
  assert(genfile);
  switch (genfile->mode) {
    case GFM_UNCOMPRESSED:
      rewind(genfile->fileptr.file);
      break;
    case GFM_GZIP:
      xgzrewind(genfile->fileptr.gzfile);
      break;
    case GFM_BZIP2:
      xbzrewind(&genfile->fileptr.bzfile, genfile->orig_path,
                genfile->orig_mode);
      break;
    default: assert(0);
  }
}

void genfile_delete(GenFile *genfile)
{
  if (!genfile) return;
  ma_free(genfile->orig_path);
  ma_free(genfile->orig_mode);
  ma_free(genfile);
}

void genfile_close(GenFile *genfile)
{
  if (!genfile) return;
  switch (genfile->mode) {
    case GFM_UNCOMPRESSED:
        fa_fclose(genfile->fileptr.file);
      break;
    case GFM_GZIP:
        fa_gzclose(genfile->fileptr.gzfile);
      break;
    case GFM_BZIP2:
        fa_bzclose(genfile->fileptr.bzfile);
      break;
    default: assert(0);
  }
  genfile_delete(genfile);
}
